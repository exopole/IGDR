package R_parse;

$VERSION = v0.0.1;


##########################################################################################
# parse_cpat.pm
#  
# Goal : Parsing of a R result
# 
# ! : some function need a modification of the file with CPATverif.pl
#
#
#
#
#########################################################################################

################################# Library ###############################################
use strict;
use warnings;
use Data::Dumper;


$|=1;




############################# Function ###################################################


sub list{
	my $R_result =  shift;
# 	print "R_result: $R_result";
	$R_result	 =~ s/\[\d*\]//g;	
	my @list = split (' ', $R_result);
	return @list;
}

#########################################################################################
1;
