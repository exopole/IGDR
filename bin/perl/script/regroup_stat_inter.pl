#!/usr/bin/perl -w


##########################################################################################
# 
# PROGRAM 	regroup_stat_inter
#
# VERSION 	0.0.1
#
# DATE 		30 april 2014
# 
# AUTHOR 	alexandre.nicaise28@gmail.com
# 
# GOAL
# 	regroup all result do in intergenic type analyse in several folder in one file
#
# INPUT
#	- a temp file 
#	- the final file
#	
# OUTPUT
# 	- a table with result of the statistic test
#
######################################################################################### 			

use strict;
use warnings;


my $file = shift;
my $out = shift;
print `cat intergenic.*/*stat > $file`;
open FILE, $file;
open OUT, ">", $out;
my $cpt =0;
while(<FILE>){
	$cpt++;
	if ($cpt == 1 || $_ =~ /intergenic/){
		print OUT $_;
	}
}

close(FILE);
close(OUT);