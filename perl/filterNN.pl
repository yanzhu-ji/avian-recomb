#!/usr/bin/perl

# Description: filter sequences with number of consecutive Ns larger than input in fasta sequence; better use countN.pl to check firstly.
#
# USAGE: $0 -s|sequence	[fasta] 
# 	    -n|number	[integer]
# YJ.3.25.2015, modified from filterN.pl

use strict;
use Getopt::Long;
use Bio::SeqIO;
use Bio::PrimarySeq;

my ( $fasta, $n_max, $outfile );

GetOptions (
	's|sequence=s' => \$fasta,
	'n|number=i' => \$n_max,
);

if ( $fasta eq "" || $n_max eq "" ){
	die "Please check the usage!\n";
}

$fasta =~ /(.*)\.fa.*/;

$outfile = $1."_filtered-NN.fasta"; # "NN" to distinguish from the preceding version (filterN.pl)

my @length_info;
my $inseq = Bio::SeqIO->new(-file => $fasta, 
			    -format => 'Fasta',
			   );
my $outseq = Bio::SeqIO->new(-file => ">$outfile",
			     -format => 'Fasta',
			   );
my $match = "";

for my $i (1..$n_max){
	$match .= "N";
}
#t print "$match\n";

while (my $seq = $inseq->next_seq() ){
	my $seq_string = $seq->seq();

	if ( $seq_string !~ $match ){
		$outseq->write_seq($seq);
	}

}

