#!/usr/bin/perl -w


use Utils;
use strict; 
use warnings;

my $a = shift;
my $b = shift;
my $c = shift;

my $test = Utils::max3($a,$b,$c);
print ("$test\n");