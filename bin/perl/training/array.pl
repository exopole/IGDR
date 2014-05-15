#!/usr/bin/perl

use strict; use warnings;use List::Util qw(shuffle);

my @animals = ('cat','dog','pig');
print "1st animal in array is: $animals[0]\n";
print "2nd animal in array is: $animals[1]\n";
print "Entire animals array contains: @animals\n";

# normally when we use @animals[0], it's an error.
print "$animals[0]\n";

push @animals, "fox"; # the array is now longer.

# when we assign a list to a scalar variable, then the scalar variable becomes the length of the list.
my $length = @animals; 
print "the array now contains $length elements and the last animal is a $animals[-1]\n";


# we create a liste with the first element in animals.
my ($list) = @animals;
print "we have a $list\n";

# take the first and second element on the array.
my ($first, $second) = @animals;
print "First two animals: $first $second\n";

# make a copy of @animals
my @animals2 = @animals;

# assign @animals an empty list -> destroys contents
@animals = (); 

print "Animals array now contains: @animals\n";
print "Animals2 array sill contains: @animals2\n";

# More about Array Indexes
print "Animal at array position 1-2 is @animals2[1.2]\n";
print "Animal at array position 1.8 is $animals2[1.8]\n";
print "Animal at array position -1 is $animals2[-1]\n";
print "array length = ", scalar(@animals2), "\n";
print "Anima ar array position 'foobar' is ", $animals2["foobar"],"\n";

# Take several element in a list
my @numner =  qw(xin fiddle ziggs taric leona vi undertower feeder mid bot top);
print ("@numner\n");

print ("@numner[0..3]\n");
