#!/usr/bin/perl


use strict; use warnings;


my ($n) = (@ARGV);
my $ln_factorial = (0.5*log(2*3.14159265358979))+($n + 0.5) * log($n) - $n + 1 / (12 * $n) - 1 / (360 * ($n ** 3));

print 2.71828 ** $ln_factorial, "\n";
