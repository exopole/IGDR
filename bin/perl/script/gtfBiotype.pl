#!/usr/bin/perl -w


##########################################################################################
# february 2014
# alexandre.nicaise28@gmail.com
# 
# Aims :
#	- Parse a gtf files
#	- Get different biotype
#
##########################################################################################

##################################### Library ############################################

use strict;
use warnings;
use Parser;
use Pod::Usage;
use Getopt::Long;
use Data::Dumper;
use gtf_parse;


#################################### Variables ###########################################

my %h = ();
my %h_count=();
my $count = 0;
my $verbose = 0;
my $infile;
my $man;
my $help;

#################################### Get Options ############################################

# Parsing parameters
GetOptions(
	"infile|i=s"		=> \$infile,
	"verbose|v=i"	=> \$verbose,	
	"man|m"			=> \$man,	
	"help|h"		=> \$help);	


# Print help if needed
pod2usage(-verbose =>2) if $help;
pod2usage(-verbose => 2) if $man;
pod2usage ("-i option: I need a non_empty gff file ... \n") unless (-s $infile);

################################### Main #################################################

%h = gtf_parse::parseLevelGtfhashKey($infile,'exon',"", "KNOWN", $verbose);

for my $k (keys(%h)){
	if exists($h_count{${$h{$k}->{"feature"}}[0]->{"biotype"}}){
		$h_count{${$h{$k}->{"feature"}}[0]->{"biotype"}} = 1;
	}
	else{
		$h_count{${$h{$k}->{"feature"}}[0]->{"biotype"}} += 1;		
	}
	$count++;
}

print "there is $count \n";
for my $b (keys(%h_count)){
	print "$b:  $h_count{$b}\n";
}