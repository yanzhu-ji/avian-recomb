#!/usr/bin/perl

# fetches sequences from fasta file.  This is especially useful for protein sequences, since fnafile works for DNA sequences only---plus not every time fnafile is available. -2/14/14
# 6/2/15 updated:
# 	print the full header line instead of just the seq ID.

# copied online: 
# USAGE: $0 [fasta_file] [name_list_with_or_without_">"]

use strict;
use Bio::DB::Fasta;

my $fastaFile = shift;
my $queryFile = shift;

my $db = Bio::DB::Fasta->new( $fastaFile );
#my $test_seq      = $db->seq('E51K_GIEG3IT02HXLZ3');
#print "Testing: \n$test_seq\n";

open (IN, $queryFile);
while (<IN>){
    chomp; 
    my $seqID = $_;
    if ( $seqID =~ /^>/ ){
	$seqID =~ s/^>//;
    }
    my $sequence = $db->seq($seqID);
    my $header = $db->header($seqID);
    if  (!defined( $sequence )) {
            print STDERR "Warning: Sequence $seqID not found. \n";
	    next; 
    }   
    print ">$header\n", "$sequence\n";
}

# Modification Log
# 2.14.14
# accept list file with or without ">".
