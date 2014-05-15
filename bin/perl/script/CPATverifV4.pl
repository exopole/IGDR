#!/usr/bin/perl -w

#
# program: CPATverifV3.pl
# 
# DATE: 27/01/14
# 
# Goal: do a file with the result of cpat that he presente :
#	- 1: coding gene
#	- 0: noncoding gene
#
# upgrade:
#	- support by the script R for analysis
#	- put before the prob columns : 
#		+ 1 (contain by a file of coding gene)
#		+ 0 (contain by a file of non coding gene)
#########################################################################################
use strict;
use warnings;
use Data::Dumper;
use Getopt::Long;
use Pod::Usage;
use parse_cpat;
use shuffle;

# table contain coding and non coding id
my %hash_coding =();
my %hash_noncoding = ();

# arguments : file of coding gene, non coding gene, result of cpat and the output file
my $infilec_fa= "";
my $infilenc_fa= "";
my $infile_cpat= "";
my $outfile = "";
my $help =0;
my $man = 0;
my $count =0;

GetOptions(
    'infilec|c=s'            => \$infilec_fa,
	'infilenc|nc=s'		 => \$infilenc_fa,
	'cpat|t=s'		 => \$infile_cpat,
	'output|o=s'		 => \$outfile,
	'help|?|h'		 => \$help,
	'man' 			 => \$man,
	)	
or pod2usage("Try '$0 --help' for more informations" );

pod2usage(-verbose => 1) if $help;
pod2usage(-verbose => 2) if $man;


pod2usage ("-c option: I need a non-empty input file...\n") unless (-s $infilec_fa);
pod2usage ("-nc option: I need a non-empty input file...\n") unless (-s $infilenc_fa);
pod2usage ("-t option: I need a non-empty input file...\n") unless (-s $infile_cpat);


parse_cpat::CPATverif($infilec_fa, $infilenc_fa,$infile_cpat, "$outfile.dat");
shuffle::shuffle_result("$outfile.dat", "$outfile.V2.dat");

__END__

=head1 NAME

CPATverif.pl - make an output file with 0 for noncoding gene and 1 for coding

=head1 AUTHOR

Alexandre Nicaise : alexandre.nicaise28@gmail.com

=head1 SYNOPSIS

CPATverif.pl -c <fasta file of coding gene> -nc <fasta file of non coding gene> -t <cpat result> -o output file [Options ....]

=cut 
