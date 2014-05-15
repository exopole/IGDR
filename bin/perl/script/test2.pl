#!/usr/bin/perl -w


use strict;
use warnings;
use test;
use shuffle;
# use gtf_parse;
use fasta;
use Data::Dumper;
use List::Util qw(shuffle);
# use List::Tuples qw(:all);
use Utils;
use gtf_parse;

my $gtf = shift;
my $verbose = shift;

my $href = gtf_parse::parseGTFSplitchr75_V2($gtf,'exon', 1,0,"protein_coding",$verbose);

print Dumper $href;