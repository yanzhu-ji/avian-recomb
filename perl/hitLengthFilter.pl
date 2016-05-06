#!/usr/bin/perl

# Description: Designed as the first step of RepeatModeler library refinement.  
# 	       After RepeatModeler library is blasted against genome sequence (where it was derived from), it reads a fmt6 blast output, and remove the queries with less than 3 hits that covered 80% of the queries.
# 
#d Usage: $0 <consensi.fasta> <genome_seq_blast_output.fmt6> <min_hits(integer)> <min_coverage(floating number)>
#
# Usage: $0 	-s|sequence	<consensi.fasta>
# 		-f|fmt6 	<genome_seq_blast_output.fmt6>
# 		-h|hits 	<min_hits(integer)>
# 		-c|coverage 	<min_coverage(floating number)>
# 		-p|print 	<hits_satiesfied.fmt6>	print hits meet the requirement to a separate file
#		>		<summary of how many hits meet requirement each>
#
# Yanzhu Ji 7.16.2014
# ji20@purdue.edu

use strict;
use Getopt::Long;
use Bio::SeqIO;


my ( $fastaFile, $fmt6, $min_hits, $min_coverage, $out_fmt6 );

GetOptions(
	"s|sequence=s" => \$fastaFile,
	"f|fmt6=s" => \$fmt6,
	"h|hits=i" => \$min_hits,
	"c|coverage=f" => \$min_coverage,
	"p|print=s" => \$out_fmt6,
);

my %info;

my $inseq = Bio::SeqIO->new(-file => $fastaFile,
			-format => 'Fasta',
		);
#my $db = Bio::DB::Fasta->new( $fastaFile );
#my @ids = $db->ids;
# my $n;

while (my $seq = $inseq->next_seq() ){
	my $id = $seq->id();
	my $length	= $seq->length();
	$info{$id}{length} = $length;
}

## old script...
#foreach my $id ( @ids ){
#	my $length = $db->length( $id );
#	if ( !defined $length ){
#		print STDERR "Warning: $id not found! Weird!\n";
#		next;
#	}
#	$info{$id}{length} = $length;
#	print STDOUT "length input: $id\t$info{$id}->{length}\n";
#	$n++;
#}

open (FMT6, "<", $fmt6) or die ("Cannot open file $fmt6: $!\n");

if ( $out_fmt6 ne "" ) {
	open (OUT, ">", $out_fmt6 ) or die ("Cannot create outfile: $!\n");
}


while (my $line = <FMT6> ){
	chomp $line;
	my @fields = split " ", $line;
	my $id = $fields[0];
	my $hit_length = $fields[3];
#	print STDOUT "hits input: $id\t$hit_length\n";
	push  @{$info{$id}{hit_length}}, $hit_length;
	if ( $out_fmt6 ne ""  && $hit_length >= $min_coverage * $info{$id}{length} ){
		print OUT  $line, "\n";
	}
}

for my $id (keys %info){
#	my $temp = join ":", @{$info{$id}{hit_length}};
#	print "$id\t$temp\n";

	my $n = 0;
	my $pass = $min_coverage * $info{$id}{length};
	
	foreach my $hit_length ( @{ $info{$id}{hit_length} } ) {
#		print STDOUT "this hit: $hit_length; this pass: $pass\n";
		if ($hit_length >= $pass ){
#			print "compared: $hit_length\t$pass\n";
			$n ++;
		}
	}
	
#	if ( $n >= $min_hits ){
		print STDOUT "$id\t$info{$id}{length}\t$pass\t$n\n";
#	}

}

# modification log
# 10/29/14, yj
# 	disabled "Bio::DB" since it requires something hard to figure out about the fasta sequence...so switched to Bio::SeqIO instead.
