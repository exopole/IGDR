#!/usr/bin/perl


use strict; use warnings;

my $sequence = "AACTAGCGGAATTCCGACCGT"; # stocke une séquence dans une variable
my $protein = "MVGGKKKTKICDKVSHEEDRISQLPEPLISEILFHLSTKDLWQSVPGLD";
my $input = "ACNGTARGCCTCACACQ";

# permet de savoir si la séquence contient des occurences de GAATTC
if ($sequence =~ m/GAATTC/) {print "EcoRI site found\n"}
else {print "no EcoRI site found\n"}


print "la séquence = ", $sequence,"\n";

# substitue la première occurence de GAATTC en gaattc
$sequence =~s/GAATTC/gaattc/;
print "la séquence modifié : $sequence\n";

# substitue la première occurence de A en adenine 
$sequence =~s/A/adenine/;
print "$sequence\n";

# substitue toutes les occurences de C en "" grâce à l'ajout de g
$sequence =~s/C//g;
print "$sequence\n";

# print la séquence si elle présente au moins une proline
print "Protein contains proline\n" if ($protein =~m/p/i);

# die va bloquer le code et empéchera la suite du code même si le if passe
die "non-DNA character in input\n" if ($input =~m/[efijlopqxz]/i);
print "we never get here\n";




