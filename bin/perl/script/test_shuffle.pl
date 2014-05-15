#!/usr/bin/perl -w

#
# 
#
#
#
#
#
#
###################################################################

use strict;
use warnings;
use Data::Dumper;
use Getopt::Long;
use Pod::Usage;


	# my $temp_header = `awk 'NR==1' temp.$outfile.V2.dat > temp_header`;
# 	my $temp_corps = `awk 'NR>1' temp.$outfile.V2.dat > temp_corps`;
# 	my $temp_shuffle_corps = `shuffle.pl -s 10 temp_corps > temp_shuffle_corps`;
# 	my $shuffle = `cat  temp_header temp_shuffle_corps > $outfile.V2.dat`;
# 	
	
	
my $in = shift;
my  $c = 0;
my $Header = "";
open IN, $in or die "Cannot open $in\n";
open OUT, ">", "shuffleV2";
open TEMP, "temp_shuffle_corps";
open CORP, ">", "temp_corps";

while (<IN>){
	$c++;
	if ($c==1){
		$Header = $_;
	}
	else{
		print CORP "$_";
	}
}

my $temp_shuffle_corps = `shuffle.pl -s 10 temp_corps > temp_shuffle_corps`;

print OUT "$Header";
while (<TEMP>){
	print OUT "$_";
}