#!/usr/bin/env perl
# samtools view sample.bam chr4:55589744-55589911 | perl complexindel.pl
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;

my $refseq = "ucsc.hg19.fasta";
my $samtools = "samtools";

##
# internal logic below
my $rawoutput = 0;
my $debug = 0;
# Phred score offset
# default: 33 (Sanger, Illumina 1.8+)
# common alternative: 64 (Illumina 1.3+ / 1.5+)
my $phredoffset = 33;
my $quality_threshold = 20;
my $min_depth = 50;
my $min_af = 0;
my $depth_bam = undef;
my $skip_lowqual = 0;
my $skip_lowdepth = 0;
my $skip_lowaf = 0;
my $max_samtools_depth = 500000; # overriding samtools depth default cap of 8000x depth
GetOptions ("refseq=s" => \$refseq,
	    "samtools=s" => \$samtools,
	    "rawoutput" => \$rawoutput,
	    "debug" => \$debug,
	    "phredoffset=i" => \$phredoffset,
	    "quality_threshold=i" => \$quality_threshold,
	    "min_depth=i" => \$min_depth,
	    "min_af=f" => \$min_af,
	    "depth_bam=s" => \$depth_bam,
	    "skip_lowqual" => \$skip_lowqual,
	    "skip_lowdepth" => \$skip_lowdepth,
	    "skip_lowaf" => \$skip_lowaf,);

my %detected_mutations;

LINE:while (<>) {
    chomp;
    next LINE if length($_) == 0 || substr($_,0,1) eq "@";
    my @fields = split("\t", $_);
    next LINE if ($fields[2] eq "*" || $fields[3] eq "*" || $fields[5] eq "*");
    my ($cigar_total_length, @cigar) = &parse_cigar($fields[5]);
    my $rawrefseq = faidx(sprintf("%s:%d-%d", $fields[2], $fields[3], $fields[3] + $cigar_total_length - 1));
    my $readseq = $fields[9];
    my $readqual = $fields[10];
    my $direction = $fields[1] & 0x0010 ? "-" : "+";
    print join("\t", $fields[0], scalar @fields, $fields[5], $readseq, $readqual)."\n" if $debug;
    my @candidate_mutations = reconstruct_alignment($cigar_total_length, \@cigar, uc($readseq), $readqual, uc($rawrefseq), $fields[3]);
    foreach my $candidate_mutation_arrayref (@candidate_mutations) {
	my $candidate_mutation_ref_len = length($candidate_mutation_arrayref->[4]);
	my $candidate_mutation_var_len = length($candidate_mutation_arrayref->[5]);
	if ($candidate_mutation_ref_len >= 1 && $candidate_mutation_var_len == 0) {
	    # simple deletion
	} elsif ($candidate_mutation_ref_len == 0 && $candidate_mutation_var_len >= 1) {
	    # simple insertion
	} else {
	    # complex indel
	    print join("\t", $fields[0], @$candidate_mutation_arrayref, $direction)."\n"  if $rawoutput;
	    my $key = join("|",$fields[2], @{$candidate_mutation_arrayref}[2..5]);
	    if (!exists $detected_mutations{$key}) {
		# forward depth, reverse depth, forward mean qual score, reverse mean qual score
		$detected_mutations{$key} = [0, 0, [], []];
	    }
	    $detected_mutations{$key}->[ $direction eq "-" ? 1 : 0 ]++;
	    push @{$detected_mutations{$key}->[ $direction eq "-" ? 3 : 2 ]}, $candidate_mutation_arrayref->[6];
	}
    }
}
print Dumper(\%detected_mutations) if $debug;
# VCF output by desecnding order of combined depth +&-
if (!$rawoutput) {
    print "##fileformat=VCFv4.1\n";
    print "##source=INDELseek\n";
    print "##reference=file://$refseq\n";
    print "##INFO=<ID=DP2,Number=2,Type=Integer,Description=\"# alt-foward and alt-reverse reads\">\n";
    print "##INFO=<ID=QS2,Number=2,Type=Float,Description=\"Mean quality scores of alt-foward and alt-reverse bases\">\n";
    if (defined $depth_bam && $depth_bam ne "") {
	print "##INFO=<ID=AF,Number=1,Type=Float,Description=\"Allele Frequency\">\n";
	print "##INFO=<ID=RDP,Number=1,Type=Float,Description=\"Mean read depth of REF positions\">\n";
        print "##FILTER=<ID=LowAF,Description=\"AF below $min_af\">\n";
    }
    print "##FILTER=<ID=LowDepth,Description=\"ALT depth below $min_depth\">\n";
    print "##FILTER=<ID=LowQual,Description=\"Mean quality scores below $quality_threshold or ALT contains N\">\n";
    print join("\t", "#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO")."\n";
    MUTATION: foreach my $key (sort {$detected_mutations{$b}->[0] + $detected_mutations{$b}->[1] <=> $detected_mutations{$a}->[0] + $detected_mutations{$a}->[1]} keys %detected_mutations) {
	print "key: $key\n" if $debug;
	my ($chrom, $start, $end, $ref, $alt) = split(/\|/,$key);
	my $qualscore_forward = mean(@{$detected_mutations{$key}->[2]});
	my $qualscore_reverse = mean(@{$detected_mutations{$key}->[3]});
	my $qualscore_combined = mean(@{$detected_mutations{$key}->[2]},@{$detected_mutations{$key}->[3]});
	my $combined_depth = $detected_mutations{$key}->[0] + $detected_mutations{$key}->[1];
	my $filter_score = 0;
	my @filter_tags;
	my $info = "";
	if ($qualscore_combined < $quality_threshold || $alt =~ /N/) {
	    next MUTATION if $skip_lowqual;
	    push @filter_tags, "LowQual";
	}
	if ($combined_depth < $min_depth) {
	    next MUTATION if $skip_lowdepth;
	    push @filter_tags, "LowDepth";
	}
	# prepare info
	$info = sprintf("DP2=%d,%d;QS2=%.2f,%.2f",
		   $detected_mutations{$key}->[0],
		   $detected_mutations{$key}->[1],
		   $qualscore_forward,
		   $qualscore_reverse,
		   );
	if (defined $depth_bam && $depth_bam ne "") {
	    my $region = sprintf("%s:%d-%d",$chrom, $start, $end);
	    my $pos_depth = depth($region);
	    my $af = $combined_depth/$pos_depth;
	    die "ERROR: unexpected pos_depth $pos_depth for $region" if $pos_depth == 0;
	    if ($af < $min_af) {
		next MUTATION if $skip_lowaf;
		push @filter_tags, "LowAF";
	    }
	    $info .= sprintf(";AF=%.3f;RDP=%.1f",$af,$pos_depth);
	    
	}


	print join("\t",
		   $chrom,
		   $start,
		   ".",
		   $ref,
		   $alt,
		   sprintf("%.2f", $qualscore_combined),
		   scalar @filter_tags == 0 ? "PASS" : join(";", sort @filter_tags),
		   $info)
	."\n";
    }
}

sub reconstruct_alignment {
    my ($cigar_total_length, $cigar_arrayref, $readseq, $readqual, $refseq, $refpos_start) = @_;
    my @output;
    my @readseq = split ("", $readseq); # 1-based, first element is not used
    unshift @readseq, undef;
    my @readqual = split ("", $readqual); # 1-based, first element is not used
    unshift @readqual, undef;
    my @refseq = split ("", $refseq); # 1-based, first element is not used
    unshift @refseq, undef;
    
    my @cigar_op; # 1-based, first element is not used
    my $cigar_pos_offset = 1;
    
    my @cigar_readseq;
    my @cigar_readqual;
    my @cigar_refpos;
    my @cigar_refseq;
    
    my $readbase_pos_offset = 0;
    my $refpos_pos_offset = $refpos_start;
    foreach my $cigarop (@$cigar_arrayref) {
	my $len = $cigarop->[0];
	my $op = $cigarop->[1];
	map {$cigar_op[$_] = $op} ($cigar_pos_offset .. $cigar_pos_offset + $len - 1);
	map {$cigar_readseq[$_] = $op eq "D" ? "*" : $readseq[$_ + $readbase_pos_offset]} ($cigar_pos_offset .. $cigar_pos_offset + $len - 1);
	map {$cigar_readqual[$_] = $op eq "D" ? " " : $readqual[$_ + $readbase_pos_offset]} ($cigar_pos_offset .. $cigar_pos_offset + $len - 1);
	map {$cigar_refpos[$_] = ($op eq "I" || $op eq "S") ? "*" : $refpos_pos_offset + $_ - $cigar_pos_offset} ($cigar_pos_offset .. $cigar_pos_offset + $len - 1);
	map {$cigar_refseq[$_] = ($op eq "I" || $op eq "S") ? "*" : $refseq[$refpos_pos_offset + $_ - $cigar_pos_offset - $refpos_start + 1]} ($cigar_pos_offset .. $cigar_pos_offset + $len - 1);
	$readbase_pos_offset -= $len if $op eq "D";
	$refpos_pos_offset += $len if ($op ne "I" && $op ne "S");
	$cigar_pos_offset += $len;
    }
    if ($debug) {
	print join("\t", map {defined $_ ? $_ : "_"} @cigar_refpos[1..$cigar_total_length])."\n";
	print join("\t", map {defined $_ ? $_ : "_"} @cigar_refseq[1..$cigar_total_length])."\n";
	print join("\t", map {defined $_ ? $_ : "_"} @cigar_op[1..$cigar_total_length])."\n";
	print join("\t", map {$cigar_op[$_] eq "M" ? ($cigar_readseq[$_] eq $cigar_refseq[$_] ? "=" : "X") : $cigar_op[$_]} (1..$cigar_total_length))."\n";
	print join("\t", map {defined $_ ? $_ : "_"} @cigar_readseq[1..$cigar_total_length])."\n";
	print join("\t", map {defined $_ ? $_ : "_"} @cigar_readqual[1..$cigar_total_length])."\n";
	print join("\t", map {defined $_ ? ord($_) - $phredoffset : -1} @cigar_readqual[1..$cigar_total_length])."\n";
    }
    my @cigar_cluster;
    # scan if there is cluster of mismatch between cigar_readseq and cigar_refseq
    my $max_distance = 5;
    my $start = -1;
    my $end = -1 * $max_distance;
    
    my $min_window_size = 2;
    my ($window_start, $window_stop, $current_window_score);
    CIGARPOS:for (my $i = 1; $i<=$cigar_total_length; $i++) {
	next CIGARPOS unless (($cigar_op[$i] eq "M" && $cigar_readseq[$i] ne $cigar_refseq[$i]) || $cigar_op[$i] eq "I" || $cigar_op[$i] eq "D");
	# X I D
	if (($i - $end) > $max_distance) {
	    if ($start >= 0) {
		if (($end - $start + 1) >= $min_window_size ) {
		    push @output, [$start, $end, min(grep {$_ ne "*"} @cigar_refpos[$start..$end]), max(grep {$_ ne "*"} @cigar_refpos[$start..$end]), join("", grep {$_ ne "*"} @cigar_refseq[$start..$end]), join("", grep {$_ ne "*"} @cigar_readseq[$start..$end]), mean(map {ord($_) - $phredoffset} grep {$_ ne " "} @cigar_readqual[$start..$end])];
		}
	    }
	    $start = $i;
	    $end = $i;
	} else {
	    if ($i > $end) {
		$end = $i;
	    }
	}
    }
    if ($start >= 0) {
	    if (($end - $start + 1) >= $min_window_size ) {
		push @output, [$start, $end, min(grep {$_ ne "*"} @cigar_refpos[$start..$end]), max(grep {$_ ne "*"} @cigar_refpos[$start..$end]), join("", grep {$_ ne "*"} @cigar_refseq[$start..$end]), join("", grep {$_ ne "*"} @cigar_readseq[$start..$end]), mean(map {ord($_) - $phredoffset} grep {$_ ne " "} @cigar_readqual[$start..$end])];
	    }
    }
    return @output;
}

sub min {
    my @array = sort {$a <=> $b} @_;
    return shift @array;
}

sub max {
    my @array = sort {$b <=> $a} @_;
    return shift @array;
}

sub mean {
    my $i = 0;
    my $sum = 0;
    map {$i++; $sum+=$_} @_;
    return $i == 0 ? -1 : $sum/$i
}

sub parse_cigar {
    my ($cigar_string) = @_;
    my @cigar;
    my $i = 0;
    while ($cigar_string =~ m/([0-9]+)([MIDNSHPX=])/g) {
	# skip hard-clipping
	die "ERROR: Unexpected CIGAR operation $2" if $2 eq "N" || $2 eq "P";
	if ($2 ne "H") {
	    push @cigar, [$1, $2];
	    $i += $1;
	}
    }
    return ($i, @cigar);
}

my %faidxcache;
my $faidxcache_count = 0;
sub faidx {
    my ($region) = @_;
    return $faidxcache{$region} if exists $faidxcache{$region};
    my $output = `$samtools faidx $refseq $region 2>&1`;
    my @outputlines = split("\n",$output);
    die "Error: unexpected output from samtools faidx" unless $outputlines[0] eq ">$region";
    my $seq = join("",@outputlines[1..$#outputlines]);
    $faidxcache{$region} = $seq;
    $faidxcache_count++;
    return $seq;
}

my %depthcache;
my $depthcache_count = 0;
sub depth {
    my ($region) = @_;
    return $depthcache{$region} if exists $depthcache{$region};
    my $output = `$samtools depth -d $max_samtools_depth -r $region -q $quality_threshold $depth_bam`;
    my @outputlines = split("\n",$output);
    my @depth;
    foreach my $line (@outputlines) {
	my @fields = split (/\t/, $line);
	push @depth, $fields[2];
    }
    $depthcache{$region} = mean(@depth);
    $depthcache_count++;
    return $depthcache{$region};
}