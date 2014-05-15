#!/usr/bin/perl -w


##########################################################################################
# february 2014
# alexandre.nicaise28@gmail.com
# 
# Aims :
#	- Parse a gtf files
#	- Get Fasta sequences of concatened exons
# Note
#	- Updated version of -sav with personal perl packages in : /home/genouest/umr6061/recomgen/tderrien/bin/ThomasPerl/lib/
#	- Option so that $slop is different in 5' and 3'
#	- August test if -f is a dir or file with assumption that if later, it is multifastafile
##########################################################################################


################################### library ##############################################
# Uses
use strict;
use Pod::Usage;
use Getopt::Long;
use Data::Dumper;
# Own package
use StringUtils;
use Utils;
use Parser;
use gtf_parse;

################################## Variables #############################################
# Global Variables
my $gtffile;
my $fastaDirFile    ="/omaha-beach/tderrien/DATA/canFam3/sequence/seq_by_chr/"; # default is canFam3
# my $fastaDirFile    ="/omaha-beach/tderrien/DATA/hg19/"; # for human
# nb of extranucleotide around transcript
my $slop=0;
my $outfile;
my $Coding;

# list all feature of RNA
# antisense, lincRNA
my $list = "protein_coding";

# iniatialisation of predict size list
my @list_size;

my $help = 0;
my $verbose=0;
my $man;
################################# Get Options ############################################
# Parsing parameters
my $result = GetOptions(
	"gtffile|i=s"		=> \$gtffile,
	"outfile|o=s"		=> \$outfile,
	"fastaDirFile|f=s"	=> \$fastaDirFile,
	"coding|c=s"		=> \Codingfile;
	"verbose|v=i"	=> \$verbose,	
	"man|m"			=> \$man,	
	"help|h"		=> \$help);	


# Print help if needed
pod2usage(-verbose =>2) if $help;
pod2usage(-verbose => 2) if $man;
pod2usage("I need a gtf file") unless ($gtffile && -r $gtffile);


pod2usage("I need a valid directory containing fasta or multifasta file ") unless ($fastaDirFile && (-d $fastaDirFile || -r $fastaDirFile));
#################################### Function ############################################




############################## Main ######################################################
print "create the h_transcript\n";
%h_coord = gtf_parse::parseLevelGtfhashKey_chr_genecoord($gtffile,'exon', $list, $verbose);

print "find de la creation du h\n";


#%h_transcript = parse_h(\$h_temp, $list,$st, $verbose);
my %h_fasta= parse_fasta::fasta2hash($fastaDirFile, $verbose);
my ($cin, %fasta) = parse_fasta::seq_ID($Coding);
@list_size = Utils::predict_size(\%fasta);
parse_fasta::random_seq(\%fasta,\@list_size,\%h_coord);








__END__

=head1 NAME

gtf2intergenicv2.pl - Extract fasta sequence from a gtf file

=head1 AUTHOR

Alexandre Nicaise : alexandre.nicaise28@gmail.com

=head1 SYNOPSIS

perl gtf2fasta.pl -i <gtf file> [Options...]
	
	-fastaDirFile|f 
	
	-outfile|o
	
	-slop|s
	
	-help|h
	
	-man|m
	
	-verbose|v

=head1 DESCRIPTION

=over 8

=item * Parse a gtf files

=item * Get Transcripts Fasta sequences of concatened exons

Options:

	-fastaDirFile|f 	:  full path to a multi-fasta file with the genomic sequences for all input mappings, OR a directory with single-fasta files
      (one per genomic sequence, with file names matching sequence_id/chr names)  [default: /omaha-beach/tderrien/DATA/canFam3/sequence/seq_by_chr/]

	-outfile|o	: output file name [default : STDOUT]
	
	-slop|s		: number of nucleotide to be added in 5' and 3' of the transcript sequence [default: 0]
	
	-help|h		: Help message [default : 0]

	-man|m		: man help [default : 0]

	-verbose|v	: level of verbosity [default : 0]

	
=head1 Details

=over 8

=item B<-fastaDirFile>
: Path to the directory that contains fasta sequences of species chromosomes

=item B<-outfile>
: output file name

=item B<-slop>
: number of nucleotide to be added in 5' and 3' of the transcript sequence

=item B<-verbose>
: Level of verbosity to follow process

=item B<-man>
: Print man help and exit

=item B<-help>
: Print help message and exit

=back




=cut
