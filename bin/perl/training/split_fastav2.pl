#!/usr/bin/perl -w


# Feb 2013
# alexandre.nicaise28@gmail.com
# Split a multi fasta file in X fasta 
# run 'split_fasta.pl -h' for help
#####################################

use warnings;
use strict;
use Bio::SeqIO;
use Getopt::Long;
use Pod::Usage;


# Parameters
############

# mfasta file
my $infile	= "";
# output prefix file
my $prefix 	= 	"split";
# number of sequences per output file(s)
my $number	= 	1 ;
# help and man
my $help	=	0;
my $man		=	0;


# GetOptions
GetOptions (
	"infile|i=s" 	=>	\$infile,
	"prefix|p=s"	=>	\$prefix,
	"number|n=i"	=>	\$number,
	"help|h"	=>	\$help,
	"man|m"		=>	\$man
)
	or pod2usage ("Try '$0 --help' for more informations");

# Help and man
pod2usage ( -verbose => 1) if $help;
pod2usage ( -verbose => 2) if $man;

pod2usage ("-i option: I need a non-empty input file...\n") unless (-s $infile);


# Parse input file with bioperl
my $in  = new Bio::SeqIO(-file  => $infile);

my $count = 1;
my $fcount = 1;


my $nb_line = `more $infile | grep ">"| wc -l `;
print "number of file : $number\nnumber of sequence of the file: $nb_line";
my $Div = $nb_line/$number;

my $out = new Bio::SeqIO(-file => ">$prefix.$fcount.fa", -format=>'fasta');
#my $out;
while (my $seq = $in->next_seq) {
        if ($count % $Div == 0 and $fcount != $number) {
				print "file $prefix.$fcount.fa done\n";
                $fcount++;
                $out = new Bio::SeqIO(-file => ">$prefix.$fcount.fa", -format=>'fasta');
        }
        $out->write_seq($seq);
        $count++;
}
print "file $prefix.$fcount.fa done\n";


__END__

=head1 NAME

split_fasta.pl - Split a multi-fasta file in N files

=head1 SYNOPSIS

perl split_fasta.pl -i <MultiFasta_File> [Options]

	Required parameters
		-i <MultiFasta_File>
	
	Optional Parameters
		-p prefix of the output files			[default 'split']
		-n number of soutput files 				[default '1']
		-h help message							[default '0']
		-m man message							[default '0']

=head1 Optional parameters

=item B<--prefix (or -p)>

The string prefix of the output files (default 'split')

=item B<--number (or -n)>

The desired number of output file(s) (default '1')

=item B<--help (or -h)>

This Help message

=item B<--man (or -m)>

Th manual page

=cut