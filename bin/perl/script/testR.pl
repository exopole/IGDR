#!/usr/bin/perl -w


use Data::Dumper;


my $CommandeR =  `PredictSize.r 1000 800 8`;

$CommandeR	 =~ s/\[\d\]//g;
$CommandeR =~ s/\[\d\d\]//g;
my @list = split (' ', $CommandeR);
print "$CommandeR\n";
print "\nlis: $list[5]\n";
$len = @list;
print "len: $len\n"