#!/usr/bin/perl

use warnings;
use strict;
 
my $pi = 3.14;

print "pi = $pi\n";

my $pi2 = 3.141593;
print "pi = $pi2\n";

my $x = 3;
my $y = 2;
my $z = ($x+$y)/2;

print "$z\n";

print "$x plus $y is ", $x + $y, "\n";
print "$x minus $y is ", $x - $y, "\n";
print "$x times $y is ", $x * $y, "\n";
print "$x divided by $y is ", $x / $y, "\n";
print "$x modulo $y is ",$x % $y, "\n";
print "$x to the power of $y is ",$x ** $y, "\n";
print "the abolu value of -$x is ", abs(-$x), "\n";
print "the natural log of $x is ", log($x), "\n";
print "the square root of $x is ", sqrt($x), "\n";
print "the sin of $x is ", sin($x), "\n";
print "a random number up to $y is ", rand($y), "\n";
print "a random integer up to $x x $y is ", int(rand($x * $y)), "\n"; 
