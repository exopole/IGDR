#!/usr/bin/perl

use strict;
use warnings;
use Data::Dumper;

my %hash_coding=();
my $infile_fa=shift;
my $infile_cpat=shift;

open FASTA, $infile_fa or die "Cannnot open $infile_fa\n";
open CPAT, $infile_cpat or die "Cannnot open $infile_cpat\n";


while(<FASTA>){
	#next if (/^#/);
	chomp; # retire le dernier caractère qui est un retour chariot ici 
	if ($_ =~ /^>(\S+)/) { # si la ligne commence par un > si on veut recup tout les lignes avec un > on ne metterais pas de ^	
		$hash_coding{$1}++;
	}
}

while (<CPAT>){

	chomp;
	my @fields= split(/\t/);
	my $id = $fields[0];
	print $id,"\n";
	if (exists($hash_coding{$id})){
		print "$_\t1\n";
	}else {
		print ""	
	}
}



