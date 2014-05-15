#!/usr/bin/perl -w

#
# Separate_fasta.pl
# 
# Date : 24/01/14
# 
#
# Goal :
#	separate 2 file and regroup them in n files
#
# v2 : eliminate list to use hash 
#####################################################

######################################## Libraries ####################################### 

# Libraries
use strict; 
use warnings;
use Bio::SeqIO;
use Getopt::Long;
use Pod::Usage;
use Data::Dumper;


######################################## Variables #######################################

# infile
my $infile = "";

# output prefix file
my $outfile 	= 	"size_id";

# help and man
my $help	=	0;
my $man		=	0;



######################################### Get Options ####################################

# GetOptions
GetOptions (
	"i=s"				=>  \$infile,
	"o=s"				=>  \$outfile,
	"help|h"			=>	\$help,
	"man|m"				=>	\$man
)
	or pod2usage ("Try '$0 --help' for more informations");

# Help and man
pod2usage ( -verbose => 1) if $help;
pod2usage ( -verbose => 2) if $man;

pod2usage ("-i option: I need a non-empty input file...\n") unless (-s $infile);



######################################## Function ########################################

sub parse_fasta{
	my ($fafile) =@_;
	my $in = Bio::SeqIO ->new(-file => $fafile);
	my %fasta=();
	my $cin=0;
	while (my $seq=$in -> next_seq){
		$cin++;
		$fasta{$seq->id} = length($seq->seq);
				
	}
	return ($cin,%fasta);
}



sub table2c{
	my ($href,$out,$col1n,$col2) = @_;
	my %h	= %{$href};
	open(OUT, ">$out");
	print OUT "$col1n\t$col2\n";
	foreach my $k (keys(%h)){
		print OUT "$k\t$h{$k}\n";
	}
}

#################################### Main ################################################

my ($cin, %fasta) = parse_fasta($infile);


table2c(\%fasta,$outfile,"ID","Size");


__END__

=head1 NAME

Table_seqsize.pl - give a table with the id and the size of the sequence

=head1 SYNOPSIS

Table_seqsize.pl -i <fasta_file1> [Options...]

	Required parameters
		-i <fasta_file>
	
	Optional Parameters
		-o prefix of the output files			[default 'group']
		-h help message							[default '0']
		-m man message							[default '0']

=head1 Optional parameters

=item B<--prefix (or -o)>

The string prefix of the output files (default 'size')

=item B<--help (or -h)>

This Help message

=item B<--man (or -m)>

Th manual page

=cut



