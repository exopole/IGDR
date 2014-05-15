#!/usr/bin/perl

#list.pl

use strict; use warnings;

my ($x,$y,$z) = (1,2,3);

print "x = $x, y=$y, z=$z\n";

($x, $y) = ($y,$x);
print "x = $x, y=$y, z=$z\n";