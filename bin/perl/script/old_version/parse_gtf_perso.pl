#!/usr/bin/perl -w

use strict;
use warnings;
use Data::Dumper;
use File::Basename;
# For parsing options
use Getopt::Long;
use Pod::Usage;

# to use liste
use List::Util qw(shuffle);

# to use Bio::perl
use Bio::SeqIO;

# to perso lib
use Utils; 
use gtf_parse;
use parse_cpat;


my $gtf = shift;
my $list = shift;
my $verbose = shift;
$verbose||= 0;
$list||="";

my %h_coord;
%h_coord = gtf_parse::parseLevelGtfhashKey_chr_genecoord($gtf,'exon', $list, $verbose);

print Dumper \%h_coord;