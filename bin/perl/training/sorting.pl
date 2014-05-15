#!/usr/bin/perl

use strict; use warnings;

# an unsorted list
my @list = qw( c b a C B A a b c 3 2 1); 
print "liste unsorted: @list\n";
my @list1 = qw( c b a C B A a b c 3 2 1);


# create of a sort list from te first list
my @sorted_list = sort @list;
print "default: @sorted_list\n";

# sort the same name of list of the first list
@list = sort @list;
print "list sort: @list\n";

# sort by numeric but we have some no numerics elements so we have errors product
@sorted_list = sort {$a <=> $b} @list;
print "numeric: @sorted_list\n";


# sort by numeric without errors
my @list2 = qw (2 4 -10000 1.6 8 121 73.2 0);
print "list2 non sorted: @list\n";
@sorted_list = sort {$b <=> $a} @list2;
print "reversed numeric: @sorted_list\n";


# 
@sorted_list = sort{$a <=> $b or uc($a) cmp uc($b)} @list1;
print "combined: @sorted_list\n";