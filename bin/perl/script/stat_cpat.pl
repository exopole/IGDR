#!/usr/bin/perl -w


##########################################################################################
#
#
# PROGRAM 	PipCPAT
#
# VERSION 	0.0.2
#
# DATE 		20 march 2014
# 
# AUTHOR 	alexandre.nicaise@etu.univ-rennes1.fr 
# 
# GOAL		do some statistic test 
#
##########################################################################################

use strict;
use warnings;
use Data::Dumper;
use test;
use Getopt::Long;
use parse_cpat;

my $infile;
my $intervalle;
my $out="result.V2.stat";
my $c;
my $n =0;
my $end =1;
my $help =0;
my $man =0;
GetOptions(
		'i=s'	=>	\$infile,
		'n=f'	=> \$intervalle,
		'b=f'	=> \$n,
		'e=f'	=> \$end,
		'c=f'	=> \$c;
		'o=s'	=> \$out,
)
	or pod2usage(2);
pod2usage(1) if $help;
pod2usage(-verbose =>2) if $man;
test::testReadableNoEmptyFileOpt(-file => $infile, -option => "i");

if ($c){
	
}


open RES,">", $out;
print RES "file\tTP\tTN\tFP\tFN\tSensibility\tSpecificity\tPrecision\tAccuracy\tcuttof\n";

while ($n<=$end+$intervalle){
	my ($TP,$TN,$FP,$FN,$TPR,$SPC,$PPV,$ACC) = parse_cpat::test_stat($infile, $n);
	print RES "$infile\t$TP\t$TN\t$FP\t$FN\t$TPR\t$SPC\t$PPV\t$ACC\t$n\n";
	$n = $intervalle+$n;

}

close(RES);
# my $command_r = `CPATAn.r $out`;

__END__

=head1 NAME

PipCPAT.pl - run all options for cpat

=head1 AUTHOR

Nicaise Alexandre : alexandre.nicaise28@gmail.com

=head1 SYNOPSIS

perl PipCPATV2.pl -i <file of CPAT result> -n <intervalle of cutoff> [Options...]
	
	    -o                      : output file prefix 
		
=head1 DESCRIPTION

=over 8

=item * use CPAT to know the coding potential of genes

=item * get a table with coding potential

Options:

	-i 	 	:  input file of result test
	
	-o 			: name of output file
	
	-n  		: intervalle of cutoff for statistics test
	
	-help|h		: Help message [default : 0]

	-man|m		: man help [default : 0]

	
=head1 Details

=over 8

=item B<-i>
: input file of CPAT result

=item B<-n>
: intervalle

=item B<-o>
: output file prefix

=item B<-man>
: Print man help and exit

=item B<-help>
: Print help message and exit

=back

=cut

