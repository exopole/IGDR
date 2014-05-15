#!/usr/bin/perl

use strict; use warnings;


# qw() is a function to make an array.
my @gene_name = qw(unc-10 cyc-1 act-1 let-7 dyf-2);
print "array : @gene_name\n";

# join() is a function uses to make a string from an array.
my $joined_string =join(", ", @gene_name);
print "string : $joined_string\n";

# spit() is a function to make an array from a string
my @gene_name2 = split(", ",$joined_string);
print "array2 : @gene_name2[1]\n";

# Example of split uses 
my $dna = "aaaaaaGAATTCttttttttttGAATTCggggggggg";
my $EcoRI = "GAATTC";
my @digest = split($EcoRI, $dna);
print "@digest\n";


# convert a string into an array and split he string at every possible position.
my @dna = split("",$dna);
print "@dna\n";
