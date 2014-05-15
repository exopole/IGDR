#!/usr/bin/perl



use strict; use warnings;

my $text = "these are letters: abcdef, and these are numbers, 123456";
print "text = $text\n";


my $sequence = "AACTAGCGGAATTCCGACCGT";


$text =~tr/a/b/; # changes any occurrence of 'a' to 'b'
print "text1 = $text\n";

$text =~tr/bs/at/; # the letter 'b' become 'a', and 's' become 't'
print "text2 = $text\n";

$text =~tr/123/321/; # 1 become 3, 2 stays as 2, 3 becomes 1
print "text3 = $text\n";

$text =~tr/abc/ABC/; # capitalize the letters a, b and c
print "text4 = $text\n";

$text =~tr/ABC/X/; # any 'A','B', or 'C' will become an X
print "text5 = $text\n";

$text =~tr/d/DE/; # incorrect use, only 'd' will be changed to 'D'
print "text6 = $text\n";

# translate all letter in the first [] in letter in the same place in the second [] 
# first version
$text =~tr [abcdefgh]
	[hgfedcba];
print "text7 = $text\n";

# second version
$text =~tr [abcdefgh][hgfedcba];
print "text8 = $text\n";

# count the number of G in the sequence
my $g_count = ($sequence =~ tr/G/G/);
print "the letter G occurs $g_count times in $sequence\n";

# transformation of uppercase letter in lowercase
$sequence =~ tr/A-Z/a-z/;
print "the sequence $sequence\n"
