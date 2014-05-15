#!/usr/bin/perl

#dnastats.pl by alex

use strict; use warnings;

die "usage: dnastats.pl <dna sequence>\n" unless @ARGV == 1;
my ($seq) = @ARGV;

