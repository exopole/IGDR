#!/usr/bin/perl -w

use strict; use warnings;
use Text::CSV;
use Data::Dumper;


# stocke le path du fichier dans une variable
my $file = shift;
my $out = shift;
# stocke les options de csv, separateur :  \t
my $csv = Text::CSV->new({eol =>$/, sep_char => "\t"});

# ouvre le fichier csv
open my $fh, "<:encoding(utf8)", $file or die "test.csv: $!";

open my $fh2, ">:encoding(utf8)", $out or die "$out : $!";

#print Dumper \$csv; 

#parcours le fichier csv
while (my $row = $csv->getline($fh)){
	my @field = @$row;
	print @$row,"\n";
	print $field[1], "\n";
	print my $length = @field,"\n";
	$csv -> print ($fh,$_) for 
}


