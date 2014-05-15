package common::sense;

our $VERSION = 3.72;

# overload should be included

sub import {
   local $^W; # work around perl 5.16 spewing out warnings for next statement
   # use warnings
   ${^WARNING_BITS} ^= ${^WARNING_BITS} ^ "\x3c\x3f\x33\x00\x0f\xf0\x0f\xc0\xf0\xfc\x33\x00";
   # use strict, use utf8; use feature;
   $^H |= 0x800700;
}

1
