#!/usr/bin/perl -w


use strict;
use warnings;
use test;
use shuffle;
use gtf_parse;
use fasta;
use Data::Dumper;
use Utils;
use List::Util qw(shuffle);
use Bio::DB::Fasta;
use Scalar::Util qw(looks_like_number);
# use Bio::Index::Fasta;

my $gtf = shift;
my $verbose = shift;

gtf_parse::parseGTFSplitchr75($gtf,'exon', $verbose);



