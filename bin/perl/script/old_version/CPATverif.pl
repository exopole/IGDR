#!/usr/bin/perl -w

#
# program : CPATverif.pl
# 
# DATE : 27/01/14
# 
# Goal : do a file with the result of cpat that he presente :
#	- 1 : coding gene
#	- 0 : noncoding gene
#
#######################################
use strict;
use warnings;
use Data::Dumper;

# table contain coding and non coding id
my %hash_coding =();
my %hash_noncoding = ();

# arguments : file of coding gene, non coding gene, result of cpat and the output file
my $infilec_fa = shift;
my $infilenc_fa=shift;
my $infile_cpat = shift;
my $outfile = shift;

# opening of files
open CODING, $infilec_fa or die "Cannot open $infilec_fa\n";
open NONCODING, $infilenc_fa or die "Cannot open $infilenc_fa\n";
open CPAT, $infile_cpat or die "Cannot open $infile_cpat\n";
open OUT,">", $outfile or die "Cannot open $outfile\n";

# reading of the coding gene file 
while(<CODING>){
	chomp; # eliminate the last charactere, here \n
	if ($_=~/^>(\S+)/){ # select line begining by > (^>) and have a sequence before (S+)
		$hash_coding{lc $1}++; # add the id in the hash and the occurence
	}
}

#print Dumper \%hash_coding;
#die;
# reading of the noncoding gene file
while(<NONCODING>){
        chomp;
        if ($_=~/^>(\S+)/){
			$hash_noncoding{lc $1}++;
        }
}

#print Dumper \%hash_noncoding;
#die;

# reading of the cpat file
while(<CPAT>){
	chomp;
	my @fields = split (/\t/); # stock
	my $id =$fields[0];
	$id =~ tr/A-Z/a-z/;
	print "id de cpat:$id\n";
	#print $hash_coding{$id};
	if (exists($hash_coding{lc $id})){
		print OUT "$_\t1\n";
	}
	elsif (exists($hash_noncoding{lc $id})){
		print OUT "$_\t0\n";
	}
	else {
		print OUT "$_\tLabel\n";
	}
}


