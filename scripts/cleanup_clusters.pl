#!/usr/bin/env perl
use strict;
use warnings;
use File::Basename;
use String::Approx;
my $BNAME = basename($ARGV[0], ".contigs.fasta");
my $DNAME = dirname($ARGV[0]);

my $KEEP_COV_RATIO = 0.1;

exit 0 if -s $ARGV[0] == 0;

print "Processing $ARGV[0]\n";

# identify which contigs to keep on the basis of coverage
open(MPILE, "samtools mpileup $DNAME/$BNAME.raw1.pe.sort.bam 2> /dev/null |");
my %contigs;
while(my $line = <MPILE>){
	my @d = split(/\t/, $line);
	$contigs{$d[0]}{cov} = 0 unless defined $contigs{$d[0]}{cov};
	$contigs{$d[0]}{cov} += $d[3];
	$contigs{$d[0]}{sites} = 0 unless defined $contigs{$d[0]}{sites};
	$contigs{$d[0]}{sites} ++;
}
my $max_cov = 0;
foreach my $contig(keys(%contigs)){
	$contigs{$contig}{cov} /= $contigs{$contig}{sites};
	$max_cov = $max_cov > $contigs{$contig}{cov} ? $max_cov : $contigs{$contig}{cov};
}
# print STDERR "Max cov is $max_cov\n";

foreach my $contig(keys(%contigs)){
	 $contigs{$contig}{keep} = ($contigs{$contig}{cov} / $max_cov > $KEEP_COV_RATIO);
#	print STDERR "keep is ".($contigs{$contig}{keep}?"true":"false")."\n";
}

# identify which contigs need to be reversed based on which direction the fill-in
# reads align to the contig
open(VIEW, "samtools view $DNAME/$BNAME.raw1.pe.sort.bam |");
while(my $line = <VIEW>){
	my @d = split(/\t/, $line);
	$contigs{$d[2]}{r1_fwd} = 0 unless defined $contigs{$d[2]}{r1_fwd};
	$contigs{$d[2]}{r1_rev} = 0 unless defined $contigs{$d[2]}{r1_rev};
	$contigs{$d[2]}{r2_fwd} = 0 unless defined $contigs{$d[2]}{r2_fwd};
	$contigs{$d[2]}{r2_rev} = 0 unless defined $contigs{$d[2]}{r2_rev};
	next if $d[2] eq "*";
	if($d[0] =~ /:r1pair:/ && $d[3] < 100){
		$contigs{$d[2]}{r1_fwd} ++ if $d[1] & 0x10;
		$contigs{$d[2]}{r1_rev} ++ unless $d[1] & 0x10;
	}
	next unless defined($contigs{$d[2]}) && defined $contigs{$d[2]}{sites};	# if the coverage was so low it didnt get into the pileup
						# then it will get removed later and this step isnt needed
	if($d[0] =~ /:r2pair:/ && $d[3] > $contigs{$d[2]}{sites} - 100 - length($d[9])){
		$contigs{$d[2]}{r2_fwd} ++ if $d[1] & 0x10;
		$contigs{$d[2]}{r2_rev} ++ unless $d[1] & 0x10;
	}
}

# if more End1 reads are reversed, then need to reverse complement the whole contig
foreach my $contig(keys(%contigs)){
	if($contigs{$contig}{r1_fwd} <= $contigs{$contig}{r1_rev} && $contigs{$contig}{r2_fwd} >= $contigs{$contig}{r2_rev}){
		$contigs{$contig}{revcomp} = 1;
	}
#	print STDERR "r1_fwd is $contigs{$contig}{r1_fwd} r1_rev is $contigs{$contig}{r1_rev} r2_fwd is $contigs{$contig}{r2_fwd} r2_rev is $contigs{$contig}{r2_rev}\n";
}


# now filter contigs
open(FA, "$DNAME/$BNAME.contigs.fasta");
my $name;
my %faseqs;
while(my $line = <FA>){
	chomp $line;
	if($line =~ /^>(.+)/){
		$name = $1;
	}else{
		$faseqs{$name} .= $line;
	}
}

open(FQ, "$DNAME/$BNAME.contigs.fastq");
my %fqseqs;
my $section = 0;
while(my $line = <FQ>){
	chomp $line;
	if($section == 0 && $line =~ /^@(.+)/ || ($section == 2 && defined $fqseqs{$name}{quals} && length($fqseqs{$name}{seq}) == length($fqseqs{$name}{quals}))){
		$line =~ /^@(.+)/;
		$name = $1;
		$section=1;
		$fqseqs{$name}{seq} = "";
	}elsif($section == 1 && $line =~ /^\+/){
		$section = 2;
		$fqseqs{$name}{quals} = "";
	}elsif($section == 1){
		$fqseqs{$name}{seq} .= $line;
	}elsif($section == 2){
		$fqseqs{$name}{quals} .= $line;
	}
}

# reverse complement any sequences that are in the wrong orientation
# then trim off anything that assembled beyond into the adapter region
my %keep_fa;
foreach my $name(keys(%faseqs)){	
	next unless defined($contigs{$name}) && $contigs{$name}{keep};
#	next unless $contigs{$name}{sites} > 1250;
	$faseqs{$name} = reverse($faseqs{$name}) if $contigs{$name}{revcomp};
	$faseqs{$name} =~ tr/ABCDGHMNRSTUVWXYabcdghmnrstuvwxy/TVGHCDKNYSAABWXRtvghcdknysaabwxr/ if $contigs{$name}{revcomp};
	$keep_fa{$name} = $faseqs{$name};
	$fqseqs{$name}{seq} = reverse($fqseqs{$name}{seq}) if $contigs{$name}{revcomp};
	$fqseqs{$name}{seq} =~ tr/ABCDGHMNRSTUVWXYabcdghmnrstuvwxy/TVGHCDKNYSAABWXRtvghcdknysaabwxr/ if $contigs{$name}{revcomp};
	$fqseqs{$name}{quals} = reverse($fqseqs{$name}{quals}) if $contigs{$name}{revcomp};
	my $match_index = String::Approx::aindex("AGAGTTTGATCMTGGCTCAG", [ "I2","D2","S25%" ], $faseqs{$name});
	if($match_index >= 0){
		$faseqs{$name} = substr($faseqs{$name}, $match_index);
		$fqseqs{$name}{seq} = substr($fqseqs{$name}{seq}, $match_index);
		$fqseqs{$name}{quals} = substr($fqseqs{$name}{quals}, $match_index);
#		print STDERR "trimming start at $match_index\n";
	}
	my $rev_annealing = "TGYACACACCGCCCGTC";
	$match_index = String::Approx::aindex($rev_annealing, [ "I2","D2","S25%" ], $faseqs{$name});
	if($match_index >= 0){
		$faseqs{$name} = substr($faseqs{$name}, 0, $match_index + length($rev_annealing));
		$fqseqs{$name}{seq} = substr($fqseqs{$name}{seq}, 0, $match_index + length($rev_annealing));
		$fqseqs{$name}{quals} = substr($fqseqs{$name}{quals}, 0, $match_index + length($rev_annealing));
#		print STDERR "trimming end at $match_index\n";
	}
}

open(FACLEAN, ">$DNAME/$BNAME.contigs.cleaned.fasta");
open(FQCLEAN, ">$DNAME/$BNAME.contigs.cleaned.fastq");
foreach my $name(keys(%keep_fa)){
	print FACLEAN ">$BNAME:$name\n$keep_fa{$name}\n";
	print FQCLEAN "@"."$BNAME:$name\n$fqseqs{$name}{seq}\n+\n$fqseqs{$name}{quals}\n";
}
close FACLEAN;
close FQCLEAN;

