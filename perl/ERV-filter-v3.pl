#!/usr/bin/perl
#
# Description: to sort domtblout of hmmscan output (--domtblout) quickly for potential ERVs (or any other suites of pfam domains with one pfam domain as "bait"), just print out matching queries with the pfam domain of interest (as input) and any other co-linear pfam domains...see if this works and if a systematic script (ERV-filter.pl) is needed.
#
# USAGE:  $0	-i/info	<file>	information about pfam domain organization in tab-delimited file with two columns
# 		-e/exclude[optional]	if mentioned, then queries NOT satiesfying required conditions (>=3 domains) will be prnted instead.
# 		-o/outfile <string>	output file names; if "-e" is specified, then two files with suffix ".out1" and ".out2" will be produced; otherwise, only ".out1" will be produced.
# 		<file>		    .hmmscan.domtblout 
#
# Author: YJ
# Date: 10.9.2014 - ?
# 4.30.2015: changed " " delimited output to "," delimited output (within each column).

use strict;
use Getopt::Long;

### Read in arguments
my ( $struc_info, $print_all, $out, $num_dom, $e_value );


GetOptions (
	'i|info=s' => \$struc_info,
	'o|outfile=s' => \$out,
	'p|print-all' => \$print_all,
	'n|number=i' => \$num_dom,
	'e|e-value=f' => \$e_value,
);

$num_dom ||= 0;
$e_value ||= 0.00001;

if ( $struc_info eq "" || $out eq "" || @ARGV != 1 ) {
        die "Error: wrong input!\nPlease check the usage!\n";
}


my $domtblout = $ARGV[0];

### Parse structural information in info.txt
my %structure;

open (STRUC, "<", $struc_info ) or die ( "Cannot open txt file with structual organization: $!\n" );

while ( my $line = <STRUC> ){
	if ( $line !~ /^#/ ){
		chomp $line;
		my @line = split "\t", $line;
		$structure{$line[0]} = $line[1];
	}
}

close STRUC;

### Parse domtblout, select the queries with pfam of interest and print
my ( $queryID, $readingFrame, %hmmsc );

open ( DOMTBLOUT, "<", $domtblout ) or die ( "Cannot open domtblout file $domtblout: $!\n" );

while ( my $line = <DOMTBLOUT> ){
	if ( $line !~ /^#/ ){
		chomp $line;
		my @line = split " ", $line;

		# filter E-value (c-Evalue)
		if ( $line[11] > $e_value ) {	 next; }

		$line[3] =~ /\b(.+)_([1-6])\b/;
		$queryID = $1;
		$readingFrame = $2;
		push @{ $hmmsc{$queryID}{pfam} }, $line[0];

		if ( $structure{$line[0]} eq "add_on" ){
			push @{ $hmmsc{$queryID}{add_on} }, $structure{$line[0]};
		}elsif ( $structure{$line[0]} ){ # non add_on domains, i.e. really informative domains
			push @{ $hmmsc{$queryID}{domain} }, $structure{$line[0]};
		}else {
			push @{ $hmmsc{$queryID}{undefined} }, $line[0];
		}
	}
}

my $total_id;

my $out1 = $out.".out1";
my $out2 = $out.".out2";
my $out3 = $out.".out3";

open (OUT1, ">", $out1) or die ("Cannot create output file $out1: $!\n");
if ( $print_all eq  "1" ) {
	open (OUT2, ">", $out2 ) or die ( "Cannot create output file $out2: $!\n");
	open (OUT3, ">", $out3 ) or die ( "Cannot create output file $out3: $!\n");
}


LOOP:
foreach my $id ( keys %hmmsc   ){
	$total_id ++;
	 @{$hmmsc{$id}{domain}} = get_unique (  @{$hmmsc{$id}{domain}} );
	 @{$hmmsc{$id}{undefined}} = get_unique (  @{$hmmsc{$id}{undefined}} );
	 @{$hmmsc{$id}{pfam}} = get_unique ( @{$hmmsc{$id}{pfam}} );

	 if ( @{$hmmsc{$id}{undefined}} == 0) {
		if ( @{$hmmsc{$id}{domain}} >= $num_dom ) {
			if ( @{$hmmsc{$id}{domain}} == 1 && ${$hmmsc{$id}{pfam}}[0] eq "RVT_1"){			 	
				print OUT2 "RT only\t$id\t", join (",", @{ $hmmsc{$id}{domain} } ), "\t", join (",",  @{ $hmmsc{$id}{pfam} } ), "\n";
			} else { # ( @{$hmmsc{$id}{domain}} == 1 && ${$hmmsc{$id}{pfam}}[0] eq "RVT_1"){
					print OUT1 "prelim-ERV\t$id\t", join (",", @{ $hmmsc{$id}{domain} } ), "\t", join (",",  @{ $hmmsc{$id}{pfam} } ), "\n";
			}
#			delete $hmmsc{$id}{pfam};  # this should also fix the redundant line issue.
#			next LOOP;
		} elsif (  $print_all eq "1" ){
			print OUT2 "unsure\t$id\t", join (",", @{ $hmmsc{$id}{domain} } ), "\t", join (",",  @{ $hmmsc{$id}{pfam} } ), "\n";
#			next LOOP;
		}
	}elsif (  $print_all eq "1" ){
		print OUT3 "others\t$id\tir-pfam\t", join (",", @{ $hmmsc{$id}{pfam} } ), "\n";  # modified 5.13, to make sure info get copied later using column-transfer.pl
	}
}


print "Total ids: $total_id\n";
# print "Number of domains parameter: $num_dom\n";

######## sub get_unique
## returns unique items in the list @array: 
## @array = get_unique ( @array );
sub get_unique{
	my %hash = map { $_, 1 } @_;
		my @unique = keys %hash;
			return @unique;
}


###################### La FIN ##############################


#if ( $print_all ){
#	print "sequences without interested domains:\n";
#	foreach my $id ( keys %hmmsc ){
#		foreach my $pfam_hit ( @{$hmmsc{$id}} ){
#			print "$id:", join ("  ", @{$hmmsc{$id}} ), "\n";
#		}
#	}
#}


# foreach my $id ( keys %hmmsc ){
#	foreach my $pfam_hit ( @{$hmmsc{$id}} ){
#		if ( $pfam_hit eq $pfam ) {
#			print "$id:", join ("  ", @{$hmmsc{$id}} ), "\n";
#		}
#	}
#}

