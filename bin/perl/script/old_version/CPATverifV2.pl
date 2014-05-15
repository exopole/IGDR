#!/usr/bin/perl -w

#
# program: CPATverifV2.pl
# 
# DATE: 27/01/14
# 
# Goal: do a file with the result of cpat that he presente :
#	- 1: coding gene
#	- 0: noncoding gene
#
# upgrade:
#	- support by the script R for analysis
#	- eliminate the prob columns to put 
#		+ 1 (contain by a file of coding gene)
#		+ 0 (contain by a file of non coding gene)
#########################################################################################
use strict;
use warnings;
use Data::Dumper;
use Getopt::Long;
use Pod::Usage;

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


# opening of files
open CODING, $infilec_fa or die "Cannot open $infilec_fa\n";
open NONCODING, $infilenc_fa or die "Cannot open $infilenc_fa\n";
open CPAT, $infile_cpat or die "Cannot open $infile_cpat\n";
open OUT,">", $outfile;

# reading of the coding gene file 
while(<CODING>){
	chomp; # eliminate the last charactere, here \n
	if ($_=~/^>(\S+)/){ # select line begining by > (^>) and have a sequence before (S+)
		$hash_coding{lc $1}++; # add the id in the hash and the occurence
	}
}

#print Dumper \%hash_coding;
#die;
# reading of the noncoding gene file
while(<NONCODING>){
        chomp;
        if ($_=~/^>(\S+)/){
                $hash_noncoding{lc $1}++;
        }
}

#print Dumper \%hash_noncoding;
#die;

# reading of the cpat file
while(<CPAT>){
	chomp;
	my @fields = split (/\t/); # stock
	my $id =$fields[0];
	#print $hash_coding{$id};
	if ($count==0){
		print OUT "ID\tmRNA\tORF\tFickett\tHexamer\tLabel\n";
		$count++;
	}
	elsif (exists($hash_coding{lc $id})){
		print OUT "$fields[0]\t$fields[1]\t$fields[2]\t$fields[3]\t$fields[4]\t1\n";
	}
	elsif (exists($hash_noncoding{lc $id})){
		print OUT "$fields[0]\t$fields[1]\t$fields[2]\t$fields[3]\t$fields[4]\t0\n";
	}
	else {
		print OUT "$fields[0]\t$fields[1]\t$fields[2]\t$fields[3]\t$fields[4]\trien\n";
	}
}


__END__

=head1 NAME

CPATverif.pl - make an output file with 0 for noncoding gene and 1 for coding

=head1 AUTHOR

Alexandre Nicaise : alexandre.nicaise28@gmail.com

=head1 SYNOPSIS

CPATverif.pl -c <fasta file of coding gene> -nc <fasta file of non coding gene> -t <cpat result> -o output file [Options ....]

=cut 
