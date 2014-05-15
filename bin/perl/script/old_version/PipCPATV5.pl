#!/usr/bin/perl -w

##############################################################################################################
# 
# PROGRAM 	PipCPAT
#
# VERSION 	0.0.2
#
# DATE 		23 january 2014
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
#	- Take of the value of performance find with R
#	- Put all sequence under 50 nu at 50
#	- Change of the sub fasta shuffle 
# 			
##############################################################################################################

###################################### Library #######################################

use strict;
use warnings;
use Data::Dumper;
use File::Basename;
# For parsing options
use Getopt::Long;
use Pod::Usage;

# to use liste
use List::Util qw(shuffle);

# to use Bio::perl
use Bio::SeqIO;

# to perso lib
use Utils; 
use parse_cpat;
use fasta;
use R_parse;
use shuffle;

##################################### Variables ###################################

# file of coding gene
my $Coding;

# file of non coding gene to training
my $NonCoding;

# file of test
my $ftest;

# ouput file name
my $outfile = "";

my $help =0;
my $man = 0;
my $verbose = 0;

# cutoff use to do statistic
my $n = 0;

# number of shuffle file 
my $nshuffle =1;

# Minimal size to be used when prediced/inferred lncRNAs size is lower than this value
# BNy definition, we assume that no lncRNA are lower in size than $cut 
my $cut = "100"; 


# let to do statistic test
my $stat;

# file of coding gene for statistic test
my $statc;

# file of non coding gene for statistic test
my $statnc;

############################# Get Options ######################################

GetOptions(
	'c=s'					=> \$Coding,
	'nc=s'					=> \$NonCoding,
	'test|t=s'				=> \$ftest,
	'stat|s=s'				=> \$stat,
	'statc|so=s'			=> \$statc,
	'statnc|sn=s'			=> \$statnc,
	'cs=s'					=> \$cut,
	'n=f'					=> \$n,
	'ns=s'					=> \$nshuffle,
	'outfile|o=s'			=> \$outfile,
	'verbose|v=i'			=> \$verbose,
	'help|?|h'				=> \$help,
	'man'					=> \$man,
	)
	or pod2usage(2);
pod2usage(1) if $help;
pod2usage(-verbose =>2) if $man;

pod2usage ("-c option: I can't read the input file...\n") unless (-r $Coding);
pod2usage ("-c option: I need a non-empty input file...\n") unless (-s $Coding);
pod2usage ("-t option: I can't read the input file...\n") unless (-r $ftest);
pod2usage ("-t option: I need a non-empty input file...\n") unless (-s $ftest);

if ($stat ne ""){	
	unless ($statc){
		$statc = $Coding;
	}
	else {
		pod2usage ("-so option: I can't read the input file...\n") unless (-r $statc);
		pod2usage ("-so option: I need a non-empty input file...\n") unless (-s $statc);
	}
	unless ($statnc){
		$statnc = $NonCoding;
		pod2usage ("-nc option: I can't read the input file...\n") unless (-r $NonCoding);
		pod2usage ("-nc option: I need a non-empty input file...\n") unless (-s$NonCoding);
	}
	else{
		pod2usage ("-sn option: I can't read the input file...\n") unless (-r $statnc);
		pod2usage ("-sn option: I need a non-empty input file...\n") unless (-s $statnc);
	}
}

if ($outfile eq ""){
	if ($ftest =~ m/.fa/){
		$outfile = basename(substr($ftest,0,-3));
	} 
	else {
		$outfile = basename($ftest);
	}	
}

pod2usage ("Can't exec cpat.py...\n") if (-x "cpat.py");
pod2usage ("Can't exec make_hexamer_tab.py...\n") if (-x "make_hexamer_tab.py");
pod2usage ("Can't exec make_logitModel.py...\n") if (-x "make_logitModel.py");
pod2usage ("Can't exec PredictSize.r...\n") if (-x "PredictSize.r");
pod2usage ("Can't exec Human_8fold_crossValidation_allFeatureV2.r...\n") if (-x "Human_8fold_crossValidation_allFeatureV2.r");

################################ Function ######################################

# shuffle sequence of a fasta file
# Input:
#	- $seqfile: 	file of sequence
#	- $number: 		number of shuffle file
# Output:
# 	- a file with the same number of sequence than $seqfile but with a shuffle sequence
sub fasta_shuffle{
	
	my ($seqfile, $number, $cut) = @_;

	print "create of $number file of shuffle sequence fasta \n";
	print "Parsing of $seqfile\n";
	# use of a library to parse $seqfile. 
	# input : $seqfile
	# output : 
	#	- the number of sequence
	# 	- A has that contain identifiant of sequence in keys and sequence in values. 
	my ($cin, $fasta_ref) = parse_fasta::seq_ID($seqfile);
	my %fasta = %{$fasta_ref};
	my $sequence = group_seq_h(\%fasta);
	my $mean = mean_h(\%fasta);
	my $sd = sd_h(\%fasta, $mean);
	for (my $i =1; $i<=$number;$i++){
		print "write the file $i\n";
		# use R to extract some value who follow a normal law
		# input: 
		#	- a mean ($mean)
		#	- a standard deviation ($sd)
		#	- the number of sequence contain by $seqfile
		# output:
		#	- a string 
 		print "R\n";
 		print "$cut\n";
		my $commandR = `PredictSize.r $mean $sd $cin $cut`;
 		# parse he result of R in a list of size
		my @list_size = R_parse::list($commandR);
		print "R done\n";
		# shuffle the list that contain the sequence
		print "shuffle\n"; 
		write_seqfa_rand($sequence, $cin, \@list_size, $i);
		print "file $i done\n";
	}
 }



sub mean_h {
	my $href = shift;
	my %h = %{$href};
	my $cpt = 0;
	my $len = 0;
	foreach my $seq (keys(%h)){
		$cpt++;
		$len += $h{$seq}->{"size"}; 
	}
	return $len/$cpt;
}


sub sd_h {
	my ($href, $mean)= @_;
	my %h = %{$href};
	my $cpt = 0;
	my $sdb = 0;
	foreach my $key (keys(%h)){
		$cpt++;
		$sdb += ($h{$key}->{"size"}-$mean) * ($h{$key}->{"size"}-$mean); 
	}
	return sqrt($sdb/$cpt);
}

sub group_seq_h {
	my $href =shift;
	# initilazation of $sequence like a void string.	
	my %fasta = %{$href};
	my $sequence = "";
	my $cpt = 0;
	# all sequence contain in value of %fasta is put in $sequence 
	foreach my $key (keys(%fasta)){
		$sequence =$sequence.$fasta{$key}->{"sequence"};
		#print "encore ", $cin-$cpt,"sequences to group\r"; 
		$cpt++;
	}
	print "all sequence group\n";

	return $sequence;
}

sub write_seqfa_rand {
	my ($seq,$cin, $ref_list, $i) =@_;
	my @list_size = @{$ref_list};
	my $newid;
	# initialization of a file in fasta format
	my $seqOUT =Bio::SeqIO ->new(-format => 'fasta', -file => ">shuffle".$i.".fa");
	for (my $x = 1; $x <=$cin; $x++){
		# create the identifiant of the sequence
		$newid	= $x."_shuffle_".$i;
		# extract a random number between 0 and the length of the sequence - the random size
		my $rand = int(rand($cin-int($list_size[$x-1])));
		my $len = length($seq);
		# extract the list between the coord give by $rand and the size random
		my @chars = substr($seq, $rand, $rand+$list_size[$x-1]);
		@chars = shuffle @chars;
		# transform the list in string
		my $seq_shuffle = join '',@chars;
# 		print "sequence: ",$seq_shuffle,"\n";
		# create a new sequence with $newid in identifiant and $seq_shuffle in sequence
		my $new_seq = Bio::Seq->new(-id => $newid, -seq => $seq_shuffle);
		# write sequence in the fasta file
		$seqOUT->write_seq($new_seq);
	}
}


# shuffle file of result
# input: file of result
# ouput: shuffle file of result but not the header  
sub shuffle_result{
	my ($in,$out) = @_;
	# a count to know if its a header
	my  $c = 0;
	my @corps;
	# Open the file of result
	open IN, "$in" or die "Cannot open $in\n";
# 	open OUT, ">", "$out" or die "Cannot write in $out\n";
	# Browse the input file
	while (<IN>){
		$c++;
		# the line 1 is write in the ouput file
		chomp;
		push (@corps, $_);
	}
	# shuffle line in the temp file and store in another file
	shuffle::new_shuffle(-seed => 10, -table => \@corps, -output => $out);

	close (IN);
} 


# Uses all commande of cpat and a small script r
# Input: 
#	- file of coding gene($Coding)
#	- file of non coding gene ($NonCoding)
#	- file that contain sequence to test ($ftest)
#	- name of the ouputfile ($outfile) 
# Ouput:
#	- result of cpat
sub cpat{
	my ($Coding, $NonCoding, $ftest, $outfile) = @_;

	print "build the hexamere table\n";
	# build the hexamere table
	my $hexamere = `make_hexamer_tab.py -c $Coding -n $NonCoding  >"$outfile"_hexamer.table`;

	print "build the logit\n";
	# build the logit
	my $logit = `make_logitModel.py -x "$outfile"_hexamer.table -c $Coding -n $NonCoding -o $outfile`;

	print "use the CPAT\n";
	# use CPAT
	my $cpat = `cpat.py -d $outfile.logit.RData -x "$outfile"_hexamer.table -g $ftest  -o $outfile.cpat`;
	
	print "use R\n";
	# R: build a plot and summary
	my $command = `Rcpat.r $outfile.cpat`;
}

# This sub do some statistic test and use the library parse_cpat
# Input: 
#	- a file result of cpat
# 	- a cutoff to delimitate coding rna and non coding rna
# Output:
#	- TP (True Positive): number of real coding RNA detect 
#	- TN (True Negative): number of real non-coding RNA detect
#	- FP (False Positive): number of non-coding RNA detect like coding
#	- FN (False Negative): number of coding RNA detect like non-coding
#	- TPR (Sensitivity): TP/(TP+FN)
#	- SPC (Specificity): TN/(TN+FP)
#	- PPV (Precision): TP/(TP+FP)
#	- ACC (Accuracy): (TP+TN)/(Total)
sub  test_stat{
	my ($infile, $cutoff) =@_;
	open IN, $infile or die "Cannot open $infile\n";
	my %h_parse = parse_cpat::parse_cpat($infile);
	my  ($TP,$TN,$FP,$FN)= parse_cpat::Cond(\%h_parse, $cutoff);
	my $TPR = parse_cpat::Sensitivity($TP,$FN);
	my $SPC = parse_cpat::Specificity($TN,$FP);
	my $PPV = parse_cpat::Precision($TP,$FP);
	my $ACC = parse_cpat::Accuracy($TP,$TN,$FP,$FN);
 	return ($TP,$TN,$FP,$FN,$TPR,$SPC,$PPV,$ACC);
}


# Make some analysis of the result give by CPAT
# Input:
#	- $Coding: file of coding gene in fasta
# 	- $stat: file of non coding gene in fasta
#	- $outfile: name of the output file
#	- $n: cutoff use to know TP, FP, TN and FN
# outut:
#	- ouput of test_stat
#	- 2 cutoff : accuracy and performance 
sub Stat{
	my ($Coding, $stat, $outfile, $n)=@_;

	# Verification 
	print "entrer dans les statistiques\n";
	# Change the file of result to do some analysis
	# Add the label of coding (1) or non-coding (0)
	CPATverif($Coding, $stat,"$outfile.cpat", "$outfile.dat");

	print "debut du shuffle\n";
	# shuffle line of the file result without the header to have a file ready for the R analysis
	shuffle_result("$outfile.dat", "$outfile.V2.dat");
	print "fin du shuffle\n";
	# erase of all file temp
	my $rm = `rm temp*`; 
	
	# R
	print "R\n";
	# A script that going to make some graph and calcul cutoff of accuracy and performance
	my @commandeR = `Human_8fold_crossValidation_allFeatureV2.r  $outfile.V2.dat`;
	my $cutf1 = $commandeR[2];
	my $cutf2 = $commandeR[3];
	chomp ($cutf1);
	chomp ($cutf2);
	 $cutf1 =~ s/\[1\]//;
	 $cutf2 =~ s/\[1\]//;
	 print "$cutf1\n";
	 print "$cutf2\n";
	
	# test stat
	print "Test stat\n";
	if ($n ==0){
		$n = $cutf2;
	}
	my ($TP,$TN,$FP,$FN,$TPR,$SPC,$PPV,$ACC) = test_stat("$outfile.dat", $n);
	return ($TP,$TN,$FP,$FN,$TPR,$SPC,$PPV,$ACC, $cutf1,$cutf2)

}


# Add a column of label in a file of cpat result 
# Input: 
#	- Fasta file of coding gene use to do the test file ($infilec_fa) 
#	- Fasta file of non-coding gene use to do the test file ($infilenc_fa)
# 	- File of cpat result (infile_cpat)
#	- Prefix of the output file ($outfile) 
sub CPATverif {
	
	my ($infilec_fa, $infilenc_fa, $infile_cpat, $outfile) = @_;
	# table contain coding and non coding id
	# create of a hash to id of coding gene and another for non coding gene
	my %hash_coding =();
	my %hash_noncoding = ();
	my $count = 0;
	# opening of files
	open CODING, $infilec_fa or die "Cannot open $infilec_fa\n";
	open NONCODING, $infilenc_fa or die "Cannot open $infilenc_fa\n";
	open CPAT, $infile_cpat or die "Cannot open $infile_cpat\n";
	open OUT,">", $outfile or die "Cannot write in $outfile ";
	
	# reading of the coding gene file 
	while(<CODING>){
		chomp; # eliminate the last charactere, here \n
		if ($_=~/^>(\S+)/){ # select line begining by > (^>) and have a sequence before (S+)
			$hash_coding{lc $1}++; # add the id in the hash and the occurence
		}
	}

	# reading of the noncoding gene file
	while(<NONCODING>){
			chomp;
			if ($_=~/^>(\S+)/){
					$hash_noncoding{lc $1}++;
			}
	}

	# reading of the cpat file
	while(<CPAT>){
		chomp;
		my @fields = split (/\t/); # stock the line in the list, separate is \t
		my $id =$fields[0];
		my $join;
		if ($count==0){
			print OUT "ID\tmRNA\tORF\tFickett\tHexamer\tLabel\tCoding_prob\n";
			$count++;
		}
		elsif (exists($hash_coding{lc $id})){
			$join = join "\t", @fields[0...4],"1",$fields[5]."\n"; 
			print OUT "$join";
		}
		elsif (exists($hash_noncoding{lc $id})){
			$join = join "\t", @fields[0...4],"0",$fields[5]."\n"; 
			print OUT "$join";
		}
		else {
			$join = join '\t', @fields[0...4],"rien",$fields[5]."\n"; 
			print OUT "$join";
		}
	}
	close (CODING);
	close(NONCODING);
	close(CPAT);
	close(OUT);
 }


################################ Main ##################################################
# when we haven't the file of noncoding we can use the shuffle option or intergenic option
if (!$NonCoding){
#     	fasta_shuffle($Coding, $nshuffle, $cut);die;
	if ($stat ne ""){
		open RES,">", "$outfile.stat";
		print RES "file\tTP\tTN\tFP\tFN\tSensibility\tSpecificity\tPrecision\tAccuracy\tcuttof_accuracy\tcuttof_performance\n";
	}
	for (my $i=1; $i<= $nshuffle; $i++){
		#print ("shuffle$i\n");
#    		cpat($Coding, "shuffle$i.fa",$ftest,"$outfile.$i");
		if ($stat ne ""){
			print ("shuffle$i analyse\n");
			#print ("avant entrer dans stat: $outfile$i\n");
			my ($TP,$TN,$FP,$FN,$TPR,$SPC,$PPV,$ACC, $cutf1,$cutf2) = Stat($statc, $statnc, "$outfile.$i", $n); 
			print RES "$outfile\t$TP\t$TN\t$FP\t$FN\t$TPR\t$SPC\t$PPV\t$ACC\t$cutf1\t$cutf2\n";
		}
	}

}

# When we want run normaly cpat
else{
	cpat($Coding, $NonCoding,$ftest,$outfile);
	if ($stat ne ""){
		open RES,">", "$outfile.stat";
		print RES "file\tTP\tTN\tFP\tFN\tSensibility\tSpecificity\tPrecision\tAccuracy\tcuttof_accuracy\tcuttof_performance\n";
		my ($TP,$TN,$FP,$FN,$TPR,$SPC,$PPV,$ACC, $cutf1,$cutf2) = Stat($statc, $statnc, $outfile, $n); 
		print RES "$outfile\t$TP\t$TN\t$FP\t$FN\t$TPR\t$SPC\t$PPV\t$ACC\t$cutf1\t$cutf2\n";
	}
}








__END__

=head1 NAME

PipCPAT.pl - run all options for cpat

=head1 AUTHOR

Nicaise Alexandre : alexandre.nicaise28@gmail.com

=head1 SYNOPSIS

perl PipCPATV2.pl -c <coding fasta file> -t <test fasta file> [Options...]
	
	-o		: output file prefix 
		
	-nc 			: file of non coding for training
	
	-ns 			: number of shuffle file 
	
	-n 	 		: cutoff
	
	-s			: put if you want do some test of statictic
	
	-so			: file of coding gene for statistic analyses
			
	-sn			: file of non coding gene for statistic analyses

=head1 DESCRIPTION

=over 8

=item * use CPAT to know the coding potential of genes

=item * get a table with coding potential

Options:

	-c 			:  Coding genes fasta file

	-outfile|o	: output file prefix [default : CPAT]
	
	-t 			: fasta file of test
	
	-nc			: file of non coding for training
	
	-n  		: cutoff for statistics test
	
	-ns 		: number of shuffle file
	
	-s			: put if you want do some test of statictic
	
	-so			: file of coding gene for statistic analyses
			
	-sn			: file of non coding gene for statistic analyses
	
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

