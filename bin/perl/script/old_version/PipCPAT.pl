#!/usr/bin/perl -w

##############################################################################################################
# 
# PROGRAM 	PipCPAT
#
# DATE 		23 january 2014
# 
# AUTHOR 	alexandre.nicaise@etu.univ-rennes1.fr 
# 
# GOAL
# 	run CPAT, make_hexamer, and make_logitModel, 
#
# INPUT
#	- a file of gene coding, 
#		protein coding transcripts either in BED format 
#		or mRNA sequence in FASTA format
#	- a file of gene non coding (optional)
#
# OUTPUT
# 	- a table
#	- a R analysis (optional?)
# 			
##############################################################################################################


# Libraries
use strict;
use warnings;

# For parsing options
use Getopt::Long;
use Pod::Usage;

# file of coding gene
my $Coding;

# file of non coding gene of training
my $NonCoding;

# file of test
my $ftest;

# ouput file name
my $outfile;

my $help =0;
my $man = 0;
my $n = 0.36;
my $stat= "";

GetOptions(
	'c=s'		=> \$Coding,
	'nc=s'		=>  \$NonCoding,
	'test|t=s'		=> \$ftest,
	'outfile|o=s'		=> \$outfile,
	's=s'			=> \$stat,
	'n=f'		=> \$n,
	'help|?|h'	=> \$help,
	'man'		=> \$man,
	)
	or pod2usage(2);
print @ARGV,"\n";
pod2usage(1) if $help;
pod2usage(-verbose =>2) if $man;

print "$Coding $NonCoding $ftest\n";
print "$stat\n";

print "build the hexamere table\n";
# build the hexamere table
#my $hexamere = `make_hexamer_tab.py -c $Coding -n $NonCoding  >"$outfile"_hexamer.table`;

print "build the logit\n";
# build the logit
#my $logit = `make_logitModel.py -x "$outfile"_hexamer.table -c $Coding -n $NonCoding -o $outfile`;

print "use the CPAT\n";
# use CPAT
#my $cpat = `cpat.py -d $outfile.logit.RData -x "$outfile"_hexamer.table -g $ftest  -o $outfile.out`;

print "use R\n";
# R: build a plot and summary
# my $command = `Rcpat.r $outfile.out`;

if ($stat ne ""){ 
	# Verfication 
	print "entrer dans les statistiques\n";
	my $verif = `CPATverif.pl $Coding $stat $outfile.out $outfile.dat`;
	my $verif2 = `CPATverifV2.pl -c $Coding -nc $stat -t $outfile.out -o temp.$outfile.V2.dat`;
	print "debut du shuffle\n";
	# shuffle
	my $temp_header = `awk 'NR==1' temp.$outfile.V2.dat > temp_header`;
	my $temp_corps = `awk 'NR>1' temp.$outfile.V2.dat > temp_corps`;
	my $temp_shuffle_corps = `shuffle.pl -s 10 temp_corps > temp_shuffle_corps`;
	my $shuffle = `cat  temp_header temp_shuffle_corps > $outfile.V2.dat`;
	print "fin du shuffle\n";
	my $rm = `rm temp*`;
	# VP 
	print "calcule des VP\n";
	my $VPCodant = `awk 'NR>1 && \$NF ==1{i++; if (\$(NF-1) > $n){cpt++}}END{print cpt/i}' $outfile.dat`;
	
	my $VPNonCodant = `awk 'NR>1 && \$NF ==0{i++; if (\$(NF-1) < $n){cpt++}}END{print cpt/i}' $outfile.dat`;
	

	# R
	print "R\n";
	my $commandeR = `Human_10fold_crossValidation_allFeatureV2.r  $outfile.V2.dat`;
	
	
	# results
	print "\nRESULTS\n\n";
	print "VP: \n";
	print "VP Codant: $VPCodant";
	print "VP Non codant: $VPNonCodant\n";
	print "cuttofs:\n $commandeR\n";
}
__END__

=head1 NAME

PipCPAT.pl - run all options for cpat

=head1 AUTHOR

Nicaise Alexandre : alexandre.nicaise28@gmail.com

=head1 SYNOPSIS

perl PipCPAT.pl -c <coding fasta file> [Options...]


=head1 DESCRIPTION

=over 8

=item * Parse a gtf files

=item * Get Transcripts Fasta sequences of concatened exons

Options:

	-c 	:  Coding genes fasta file

	-outfile|o	: output file prefix [default : STDOUT]
	
	-t fasta file of test
	
	-nc file of non coding for training
	
	-n cutoff
	
	-s		: file of non coding gene where we can find gene of interest [default: 0]
	
	-help|h		: Help message [default : 0]

	-man|m		: man help [default : 0]

	
=head1 Details

=over 8

=item B<-c>
: file of coding gene 

=item B<-nc>
: file of non coding gene if exist

=item B<-t>
: fasta file of genes test

=item B<-n>
: cutofff

=item B<-s>
: file of non coding gene where we can find gene of interest

=item B<-o>
: output file prefix

=item B<-man>
: Print man help and exit

=item B<-help>
: Print help message and exit

=back




=cut