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
my $infile;
my $fastaDirFile    ="/omaha-beach/tderrien/DATA/canFam3/sequence/seq_by_chr/"; # default is canFam3
# my $fastaDirFile    ="/omaha-beach/tderrien/DATA/hg19/"; # for human
# nb of extranucleotide around transcript
my $slop=0;
my $outfile;
my $biotype;
my $status;

# list all feature of RNA
# antisense, lincRNA
my $list = "antisense";
my $st = "KNOWN";

my $help = 0;
my $verbose=0;
my $man;
my %h_temp=();
my %h_transcript = ();
################################# Get Options ############################################
# Parsing parameters
my $result = GetOptions(
	"infile|i=s"		=> \$infile,
	"outfile|o=s"		=> \$outfile,
	"slop|s=i"		=> \$slop,
	"fastaDirFile|f=s"	=> \$fastaDirFile,
	"verbose|v=i"	=> \$verbose,	
	"man|m"			=> \$man,	
	"help|h"		=> \$help);	


# Print help if needed
pod2usage(-verbose =>2) if $help;
pod2usage(-verbose => 2) if $man;
pod2usage("I need a gtf file") unless ($infile && -r $infile);
pod2usage("The slop value should be a positive number and lower than 1Mb") unless ($slop >=0 && $slop <1000000);

pod2usage("I need a valid directory containing fasta or multifasta file ") unless ($fastaDirFile && (-d $fastaDirFile || -r $fastaDirFile));
#################################### Function ############################################


# Parse file at the exon level
sub parse_h{
	my ($href, $list,$st, $verbose) = @_;
	$verbose ||= 0;
	%h_temp = %{$href};
	my $j=0;
	foreach my $tx(keys(%h_temp)){
		$biotype = ${$h_temp{$tx}->{"feature"}}[0]->{"transcript_biotype"};
		$status = ${$h_temp{$tx}->{"feature"}}[0]->{"transcript_status"};

		if ($list =~ /$biotype/ and $status eq $st){
			$h_transcript{$tx}=$h_temp{$tx};
		
		}
		if ($verbose> 0) {
			Utils::showProgress(scalar(keys(%h_temp)), $j++, "Selection of lncRNA: ");
		}
	}
	return %h_transcript;
}

						##############################
# Get Sequences

sub gtf2fasta{
	my ($h_transcript, $fastaDirFile, $outfile, $verbose) = @_;
	my $i=0;
	my $h_transcript_size	= keys(%h_transcript); 
	my $seqdata;

	if ($verbose >0){ print STDERR "Get Sequences\n";}


	if (defined $outfile){
		open(OUTPUTFASTA, "> $outfile") || die ("# Error in writing fasta file") ;
	}

	# If Directory
	if (-d $fastaDirFile){

		for my $tr (keys %h_transcript){
	
			#Initalize sequence:
			my $seqstring	=	"";
			my $id_sequence	=	$h_transcript{$tr}->{'chr'};
			my $cpt=0;

	
			foreach my $exon (@{$h_transcript{$tr}->{"feature"}}) {
		
				$cpt++;

				# if first, we add $slop bp to START
				if ( $cpt == 1){
					$seqstring .=	StringUtils::getSubSequence($fastaDirFile."/".$h_transcript{$tr}->{"chr"}.".fa", $id_sequence, $exon->{"start"} - $slop, $exon->{"end"});
		
				# if last exon, we add $slop bp to END
				} elsif ($cpt == scalar(@{$h_transcript{$tr}->{"feature"}})) {
					$seqstring .=   StringUtils::getSubSequence($fastaDirFile."/".$h_transcript{$tr}->{"chr"}.".fa", $id_sequence, $exon->{"start"} , $exon->{"end"} + $slop);
				
				} else{
				$seqstring .=   StringUtils::getSubSequence($fastaDirFile."/".$h_transcript{$tr}->{"chr"}.".fa", $id_sequence, $exon->{"start"} , $exon->{"end"});
				}
			}
	
			#RevComp if strand -
			if ( $h_transcript{$tr}->{"strand"} eq '-'|| $h_transcript{$tr}->{"strand"} eq '-1') {
				$seqstring 	= StringUtils::getRevComp($seqstring);
			}


			# Summarize data e.g >TCONS_00005869 XLOC_001028_-_1:2753268-2784339_Cufflinks
			# and fasta sequence
			# header
			$seqdata 	 =  ">$tr ";
			$seqdata 	.= $h_transcript{$tr}->{"gene_id"}."_".$h_transcript{$tr}->{"strand"}."_".$h_transcript{$tr}->{"chr"}.":".$h_transcript{$tr}->{"startt"}."-".$h_transcript{$tr}->{"endt"}."_".$h_transcript{$tr}->{"biotype"}."\n";
			# sequence
			$seqdata	.= "$seqstring\n";
	
			# Print in a file
			if (defined $outfile){
				print OUTPUTFASTA $seqdata;
			# Or in STDOUT
			}else{
				print $seqdata;	
			}

			if ($verbose > 0){
				Utils::showProgress($h_transcript_size, $i++, "Print ".$tr.": ");
			}

		}

	} elsif (-r $fastaDirFile){
		## If multifasta file
	
		for my $tr (keys %h_transcript){
	
			#Initalize sequence:
			my $seqstring	=	"";
			my $id_sequence	=	$h_transcript{$tr}->{'chr'};
			my $cpt=0;

	
			foreach my $exon (@{$h_transcript{$tr}->{"feature"}}) {
		
				$cpt++;

				# if first, we add $slop bp to START
				if ( $cpt == 1){
					$seqstring .=	StringUtils::getSubSequenceSamtools($fastaDirFile, $id_sequence, $exon->{"start"} - $slop, $exon->{"end"});
		
				# if last exon, we add $slop bp to END
				} elsif ($cpt == scalar(@{$h_transcript{$tr}->{"feature"}})) {
					$seqstring .=   StringUtils::getSubSequenceSamtools($fastaDirFile, $id_sequence, $exon->{"start"} , $exon->{"end"} + $slop);
				
				} else{
					$seqstring .=   StringUtils::getSubSequenceSamtools($fastaDirFile, $id_sequence, $exon->{"start"} , $exon->{"end"});
				}
			}
	
			#RevComp if strand -
			if ( $h_transcript{$tr}->{"strand"} eq '-'|| $h_transcript{$tr}->{"strand"} eq '-1') {
				$seqstring 	= StringUtils::getRevComp($seqstring);
			}


			# Summarize data e.g >TCONS_00005869 XLOC_001028_-_1:2753268-2784339_Cufflinks
			# and fasta sequence
			# header
			$seqdata 	 =  ">$tr ";
			$seqdata 	.= $h_transcript{$tr}->{"gene_id"}."_".$h_transcript{$tr}->{"strand"}."_".$h_transcript{$tr}->{"chr"}.":".$h_transcript{$tr}->{"startt"}."-".$h_transcript{$tr}->{"endt"}."_".$h_transcript{$tr}->{"biotype"}."\n";
			# sequence
			$seqdata	.= "$seqstring\n";
	
			# Print in a file
			if (defined $outfile){
				print OUTPUTFASTA $seqdata;
			# Or in STDOUT
			}else{
				print $seqdata;	
			}

			if ($verbose > 0){
				Utils::showProgress($h_transcript_size, $i++, "Print ".$tr.": ");
			}

		}


	} else {
		pod2usage("Cannot understand defined option -f \"$fastaDirFile\"... ");
	}

	# Close file if defined
	close OUTPUTFASTA unless (!defined $outfile);
}    

############################## Main ######################################################
print "create the h_transcript\n";
%h_transcript = gtf_parse::parseLevelGtfhashKey_transcript($infile,'exon', $list, $st, $verbose);
# print Dumper \%h_transcript;
# die;
print "find de la creation du h\n";
#%h_transcript = parse_h(\$h_temp, $list,$st, $verbose);
gtf2fasta(\%h_transcript, $fastaDirFile, $outfile, $verbose)



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
