#!/usr/bin/perl -w

##############################################################################################################
# 
# PROGRAM 	PipCPAT
#
# VERSION 	0.0.2
#
# DATE 		20 march 2014
# 
# AUTHOR 	alexandre.nicaise@etu.univ-rennes1.fr 
# 
# GOAL
# 	run CPAT, make_hexamer, make_logitModel, when we haven't a file of non coding genes
#	analysis of result
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
# UPGRADE
#	- add intergenic options
#	- transfert of all sub in perl modules
#	- add some test
# 			
##########################################################################################

###################################### Library #######################################

use strict;
use warnings;
use Data::Dumper;
use File::Basename;
# For parsing options
use Getopt::Long;
use Pod::Usage;

# to test is a variable is a number
use Scalar::Util qw(looks_like_number);

# to use liste
use List::Util qw(shuffle);

# to use Bio::perl
use Bio::SeqIO;
use Bio::DB::Fasta;
# to perso lib
use Utils; 
use gtf_parse;
use parse_cpat;
use fasta;
use R_parse;
use shuffle;
use test;

##################################### Variables ###################################

# file of coding gene of ORF
my $CodingORF;

# file of coding gene of mRNA
my $CodingmRNA;

# file of non coding gene to training
my $NonCoding;

# file of test
my $ftest;

# ouput file name
my $outfile;

my $help =0;
my $man = 0;
my $verbose = 0;

# cutoff use to do statistic
my $n = 0;

# number of shuffle file 
my $nfile =1;

# Minimal size to be used when prediced/inferred lncRNAs size is lower than this value
# BNy definition, we assume that no lncRNA are lower in size than $cut 
my $pmin = "100";

# Minimal size to be used when prediced/inferred lncRNAs size is lower than this value
# BNy definition, we assume that no lncRNA are lower in size than $cut 
my $pmax =0; 


# let to do statistic test
my $stat;

# file of coding gene for statistic test
my $statc;

# file of non coding gene for statistic test
my $statnc;

# file of coding gene to predict size
my $PredictC;

# file of non coding gene to predict size
my $PredictNC;

# choice of the analyse
my $choice = 1;

# gtf file
my $gtf;

# list of biotype 
my $list = "";

# fasta file or directory that contain fasta file
my $fastaDirFile;

# use the predict size fonction or the size of ORF
my $predict;

############################# Get Options ######################################

GetOptions(
	'co=s'					=> \$CodingORF,
	'cm=s'					=> \$CodingmRNA,
	'nc=s'					=> \$NonCoding,
	'test|t=s'				=> \$ftest,
	'stat|s=s'				=> \$stat,
	'statc|so=s'			=> \$statc,
	'statnc|sn=s'			=> \$statnc,
	'predict|p=s'			=> \$predict,
	'predictc|pc=s'			=> \$PredictC,
	'predictnc|pnc=s'		=> \$PredictNC,	
	'choice|ch=s'			=> \$choice,
	'pmin=i'				=> \$pmin,
	'pmax=i'				=> \$pmax,
	'n=f'					=> \$n,
	'ns=s'					=> \$nfile,
	'gtf|g=s'				=> \$gtf,
	'fastaDirFile|f=s'		=> \$fastaDirFile,
	'outfile|o=s'			=> \$outfile,
	'verbose|v=i'			=> \$verbose,
	'help|?|h'				=> \$help,
	'man'					=> \$man,
	)
	or pod2usage(2);
pod2usage(1) if $help;
pod2usage(-verbose =>2) if $man;


################################### Test ################################################
# test if the coding file for the training and the file test can be read and not empty 
test::testReadableNoEmptyFileOpt(-file => $CodingORF, -option => "co");
test::testReadableNoEmptyFileOpt(-file => $ftest, -option => "t");
test::testReadableNoEmptyFileOpt(-file => $CodingmRNA, -option => "cm");

# test if $stat if defined.
if ($stat){	
	# test if un file coding is put for statistic analyses
	unless ($statc){
		# if there is not file of coding we take the file of coding for the training
		$statc = $CodingORF;
	}
	else {
		# if there is a file we test if that can be read and not empty
		test::testReadableNoEmptyFileOpt(-file => $statc, -option => "so");
	}
	# same with the non coding gene
	unless ($statnc){
		$statnc = $NonCoding;
		test::testReadableNoEmptyFileOpt(-file => $statnc, -option => "sn");
	}
	else{
		test::testReadableNoEmptyFileOpt(-file => $statnc, -option => "sn");
	}
}

# test if we have not a prefixe for the output file name
if (!$outfile){
	# test if the file has the extension of fasta file
	if ($ftest =~ m/.fa/){
		# extract of the basename without the extension
		$outfile = basename(substr($ftest,0,-3));
	} 
	else {
		$outfile = basename($ftest);
	}	
}

# it's for all choice of analyses type without the normal
if ($choice != 1){
	# test if have the file of Coding gene for the prediction of a distribution of size
	unless ($PredictC){
		# if there is not this file. We take this file and we check if it is readable and no empty
		
		$PredictC = $statc;
		test::testReadableNoEmptyFileOpt(-file => $CodingORF, -option => "pc");
	}
	else{
		# same thing if there is a file
		test::testReadableNoEmptyFileOpt(-file => $PredictC, -option => "pc");
	}
	# we do the same for the file of non coding gene
	unless ($PredictNC){
		$PredictNC = $statnc;
		test::testReadableNoEmptyFileOpt(-file => $CodingORF, -option => "sn");
	}
	else{
		test::testReadableNoEmptyFileOpt(-file => $PredictNC, -option => "pnc");
	}
} 

# if the choice is a normal analyses we check if the file of non coding for the training is readable and no empty
elsif($choice == 1){
	test::testReadableNoEmptyFileOpt(-file => $NonCoding, -option => "nc");
}

# if the choice is an intergenic analyses 
elsif($choice == 3){
	# check if the gtf file is non empty and readable
	test::testReadableNoEmptyFileOpt(-file => $gtf, -option => "g");
	# test if the fastaDirfile (option -f) is a directory or a file
	if (-d $fastaDirFile){
		# if it is a directory we check if we can open and if there is fasta file in
		test::testFastaInDirectory(-dirname => $fastaDirFile, -option => "f", -format => "fasta");
	}
	elsif (-r $fastaDirFile){
		# if it is a file we check if it is readable and non empty
		test::testReadableNoEmptyFileOpt(-file => $fastaDirFile, -option => "f");
		
	}
}

# check if module use can be executable and in the Path 
test::testExecPathModule("cpat.py");
test::testExecPathModule("make_hexamer_tab.py");
test::testExecPathModule("make_logitModel.py");
test::testExecPathModule("PredictSize.r");
test::testExecPathModule("Human_crossValidation_allFeatureV2.r");


################################ Function ######################################

sub extract_size{
	my %hash = @_;
	my @list;
	foreach my $str (keys %hash){
		push(@list, length($hash{$str}->{"sequence"}))
	}
	return @list;
}


################################ Main ##################################################

# When we want run normaly cpat
if($choice==1){
 	parse_cpat::cpat(-ORF => $CodingORF,-mRNA => $CodingmRNA, -noncoding => $NonCoding, -test => $ftest,-out => $outfile);
 	# do some statistic test (optional)
	if ($stat){
		open RES,">", "$outfile.stat2";
		print RES "file\tTP\tTN\tFP\tFN\tSensibility\tSpecificity\tPrecision\tAccuracy\tcuttof_accuracy\tcuttof_performance\n";
		# statistic test to calculate: true positive, true negative, false positive, false negative, precision, accuracy, specificity, sensibility, cutoff of accuracy and performance
		my ($TP,$TN,$FP,$FN,$TPR,$SPC,$PPV,$ACC, $cutf1,$cutf2) = parse_cpat::Stat($statc, $statnc, $outfile, $n); 
		print RES "$outfile\t$TP\t$TN\t$FP\t$FN\t$TPR\t$SPC\t$PPV\t$ACC\t$cutf1\t$cutf2\n";
	}
}

# when we haven't the file of noncoding we can use the shuffle option or intergenic option
# shuffle option [choice 2]
elsif ($choice==2){
	# create of n fasta file of shuffle sequences
   	parse_fasta::fasta_shuffle(-seqA => $PredictC, -seqB => $PredictNC, -seqAprime => $CodingmRNA, -number => $nfile, -pmin => $pmin, pmac => $pmax,-verbose => $verbose);
	if ($stat ne ""){
		open RES,">", $outfile.".stat";
		print RES "file\tTP\tTN\tFP\tFN\tSensibility\tSpecificity\tPrecision\tAccuracy\tcuttof_accuracy\tcuttof_performance\n";
 	}
	for (my $i=1; $i<= $nfile; $i++){
			# do CPAT on 1 shuffle fasta
    		parse_cpat::cpat(-ORF => $CodingORF, -mRNA => $CodingmRNA, -noncoding => "shuffle$i.fa", -test => $ftest, -out => "$outfile.$i");
		if ( $stat){
			print ("shuffle$i analyse\n");
			# statistic test to calculate: true positive, true negative, false positive, false negative, precision, accuracy, specificity, sensibility, cutoff of accuracy and performance
			my ($TP,$TN,$FP,$FN,$TPR,$SPC,$PPV,$ACC, $cutf1,$cutf2) = parse_cpat::Stat($statc, $statnc, "$outfile.$i", $n); 
			print RES "$outfile\t$TP\t$TN\t$FP\t$FN\t$TPR\t$SPC\t$PPV\t$ACC\t$cutf1\t$cutf2";
		}
	}

}	

# intergenic option [choice 3]
elsif ($choice==3){
	print "step1: create the h_transcript\n";
	# create of variables
	my %h_coord;
	my @list_size;
	my %h_fasta;
	my %h_fasta2;
	# store all coord of exon with the biotype in the in a hash{chr}{ID_transcript}{...}
# 	%h_coord = gtf_parse::parseGTFSplitchr75($gtf,'exon', $verbose);
	my $href = gtf_parse::parseGTFSplitchr75_V2($gtf,'exon', 1,0,"protein_coding",$verbose);

	print "step2: parsing of a fasta file\n";
	# store all sequence in a hash{sequenceID} = sequence
 	my $db  = Bio::DB::Fasta->new($fastaDirFile); # create an index if does not exist
	# test if the directory is present and if not make him if not here 
	if(! -d $outfile){print `mkdir $outfile`;}
	print "step3: predict the size\n";
	# choose if we use the predict size fonction or the same size of ORF to predict the size of pseudo lncRNA
	if (!$predict){
		@list_size = Utils::predict_size_DB(-seqfileA => $PredictC, -seqfileB => $PredictNC, -seqfileAprime => $CodingORF, -min => $pmin, -max => $pmax);
		print "predict\n";
	}
	elsif($predict eq "ORF"){
		%h_fasta2= parse_fasta::fasta2hash(-fastaDirFile => $CodingORF, -verbosity => $verbose);
		@list_size = extract_size(%h_fasta2);
		%h_fasta2 = ();
		print "ORF\n";
	}
	print "step4: create a fasta file of random sequence of intergenic\n";
	# Make a fasta file that contain intergenic sequences
	parse_fasta::random_seq_DB(-Index => $db, -list_size => \@list_size, -h_coord => \%h_coord, -outfile => "$outfile/intergenic.fa", -verbosity => $verbose);
	# Do CPAT on this intergenic fasta file
	parse_cpat::cpat(-ORF => $CodingORF, -mRNA => $CodingmRNA, -noncoding => "$outfile/intergenic.fa", -test => $ftest, -out => "$outfile/$outfile");
	if ($stat){
		open RES,">", "$outfile/$outfile.stat";
		print RES "file\tTP\tTN\tFP\tFN\tSensibility\tSpecificity\tPrecision\tAccuracy\tcuttof_accuracy\tcuttof_performance\n";
		# statistic test to calculate: true positive, true negative, false positive, false negative, precision, accuracy, specificity, sensibility, cutoff of accuracy and performance
		my ($TP,$TN,$FP,$FN,$TPR,$SPC,$PPV,$ACC, $cutf1,$cutf2) = parse_cpat::Stat($statc, $statnc,"$outfile/$outfile", $n); 
		print RES "$outfile\t$TP\t$TN\t$FP\t$FN\t$TPR\t$SPC\t$PPV\t$ACC\t$cutf1\t$cutf2\n";
	}
}



else{
	die "You need to choice between\n\t-0(or nothing): normal analyse\n\t-1: analyse by shuffle\n\t-2: analyse by intergenic\n";
}








__END__

=head1 NAME

PipCPAT.pl - run all options for cpat

=head1 AUTHOR

Nicaise Alexandre : alexandre.nicaise28@gmail.com

=head1 SYNOPSIS

perl PipCPATV2.pl -co <coding fasta file of ORF (without)> -co <coding fasta file of mRNA (with UTR)> -t <test fasta file> [Options...]
	
	    -o                      : output file prefix 
		
	        -nc 		    	: file of non coding for training
	
		-ns 			: number of shuffle file 
	
		-n 	 		: cutoff
	
		-s			: put if you want do some test of statictic [default 0]
		
		-so			: file of coding gene for statistic analyses
				
		-sn			: file of non coding gene for statistic analyses
		
		-pc		 	: file of coding gene for predict the size distribution
		
		-pnc			: file of non coding gene for predict the size distribution
		
		-ch			: choice of the analyses type [1: normal, 2: shuffle, 3: intergenic]
		
		-ns 			: number of shuffle fasta create
		
		-g			: gtf file
		
		-f			: fastaDirFile

=head1 DESCRIPTION

=over 8

=item * use CPAT to know the coding potential of genes

=item * get a table with coding potential

Options:

	-co 	 	:  Coding genes fasta file of ORF (without UTR)
	
	-cm 	 	:  Coding genes fasta file of mRNA (with UTR)

	-outfile|o	: output file prefix [default : CPAT]
	
	-t 			: fasta file of test
	
	-nc			: file of non coding for training
	
	-n  		: cutoff for statistics test
	
	-ns 		: number of shuffle file
	
	-s			: put if you want do some test of statictic [default: 0]
	
	-so			: file of coding gene for statistic analyses
			
	-sn			: file of non coding gene for statistic analyses
	
	-pc 		: file of coding gene for predict the size distribution
	
	-pnc		: file of non coding gene for predict the size distribution
	
	-ch			: choice of the analyses type [1: normal, 2: shuffle, 3: intergenic]
		
	-ns 		: number of shuffle fasta create
		
	-g			: gtf file
		
	-f			: fastaDirFile
	
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

=item B<-pc>
: file of coding gene for predict the size distribution

=item B<-pnc>
: file of non coding gene for predict the size distribution

=item B<-ch>
: choice of the analyses type [1: normal, 2: shuffle, 3: intergenic]

=item B<-ns>		
: number of shuffle fasta create
		
=item B<-g>
: gtf file

=item B<-f>		
: fastaDirFile

=item B<-o>
: output file prefix

=item B<-man>
: Print man help and exit

=item B<-help>
: Print help message and exit

=back

=cut

