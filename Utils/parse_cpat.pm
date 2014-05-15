package parse_cpat;

$VERSION = v0.0.1;


##########################################################################################
# parse_cpat.pm
#  
# Goal : Parsing of a CPAT file
# 
# ! : some function need a modification of the file with CPATverif.pl
#
#
#
#
#########################################################################################

################################# Library ###############################################
use strict;
use warnings;
use Data::Dumper;
use Utils;
use Scalar::Util qw(looks_like_number);

$|=1;




############################# Function ###################################################


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
			print OUT "mRNA\tORF\tFickett\tHexamer\tLabel\tCoding_prob\n";
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
 
sub parse_cpat{

	my $cpat_file = shift;
	open(IN, $cpat_file) || die "cannot open zobi!\n";
	my %h =();
	my $cpt =0;
	while (<IN>){
	
		$cpt++;
		next if($cpt == 1); #ignore header
		chomp;
	
		my %attribs = ();

		# parse line
	# 	my @row = split("\t");
	# 	print "Taille tab = ", scalar(@row),"\n";
		my ($id, $mRNA_size,$ORF_size, $Fickett_score,$Hexamer_score,$Label,$coding_prob) = split ("\t");
		next if !$mRNA_size;

		# store SNP data
		$h{$id}->{"mRNA_size"}			=   $mRNA_size;
		$h{$id}->{"ORF_size"}	    	=   $ORF_size;
		$h{$id}->{"Fickett_score"}	    	=   $Fickett_score;
		$h{$id}->{"Hexamer_score"}         =   $Hexamer_score;
		$h{$id}->{"coding_prob"}         =   $coding_prob;
		$h{$id}->{"Label"}         =   $Label;
		
	}
	return %h;

}

sub Cond{

	my ($href, $cut, $verbose) = @_;
	my %h	= %{$href};
	$verbose||=0;
	#$cut 	||= 	0.5;
	#$label 	||= 	0;	
#   	print Dumper \%h;
#   	die;
	my $TP=0;
	my $TN=0;
	my $FN=0;
	my $FP=0;
	my $i =1;
	my $h_size = keys(%h);
	foreach my $tx (keys %h){
		if ($h{$tx}->{'coding_prob'} >= $cut){
			if ($h{$tx}->{'Label'} eq 1){
				$TP++;
			}
			if ($h{$tx}->{'Label'} eq 0){
				$FP++;
			}
		} 
		else{
			if ($h{$tx}->{'Label'} eq 1){
				$FN++;
			}
			if ($h{$tx}->{'Label'} eq 0){
				$TN++;
			}
		}
		if ($verbose > 9){
			Utils::showProgress2($h_size, $i++, "Print transcript: ");
		}
	}
	return ($TP,$TN,$FP,$FN);
	#return ($TP/($TP+$FN));
}

sub Sensitivity {
	my ($TP,$FN) = @_;
	return ($TP/($TP+$FN));
}

sub Specificity {
	my ($TN,$FP) = @_;
	return ($TN/($TN+$FP));
}

sub Precision {
	my ($TP,$FP)=@_;
	return ($TP/($TP+$FP));
}

sub Accuracy {
	my ($TP,$TN,$FP,$FN)=@_;
	return (($TP+$TN)/($TP+$TN+$FP+$FN));
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
	my %h_parse = parse_cpat($infile);
	my  ($TP,$TN,$FP,$FN)= Cond(\%h_parse, $cutoff);
	my $TPR = Sensitivity($TP,$FN);
	my $SPC = Specificity($TN,$FP);
	my $PPV = Precision($TP,$FP);
	my $ACC = Accuracy($TP,$TN,$FP,$FN);
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
	CPATverif($Coding, $stat,"$outfile.cpat", "$outfile.V2.cpat");
	print "debut du shuffle\n";
	# shuffle line of the file result without the header to have a file ready for the R analysis
	shuffle::shuffle_result("$outfile.V2.cpat", "$outfile.V3.cpat");
	print "fin du shuffle\n";	
	# R
	print "R\n";
	# A script that going to make some graph and calcul cutoff of accuracy and performance
	my @commandeR = `crossValidation_allFeatureV2.r  $outfile.V3.cpat $outfile.png`;
	print "crossValidation_allFeatureV2.r  $outfile.V3.cpat $outfile.png";
#  	my @commandeR = `Human_crossValidation_allFeatureV2.r  $outfile.V2.dat $outfile.png`;

	print "commandeR: ",@commandeR,"\n";
	
# 	my $cutf1 = $commandeR[2];
# 	my $cutf2 = $commandeR[3];

 	my $cutf1 = $commandeR[0];
 	my $cutf2 = $commandeR[1];
 	print "$cutf1\n;";
 	print "$cutf2\n;";

	chomp ($cutf1);
	chomp ($cutf2);
	 $cutf1 =~ s/\[1\]//;
	 $cutf2 =~ s/\[1\]//;
	# test stat
	print "Test stat\n";
	if ($n ==0){
		$n = $cutf2;
	}
	die "the cuttof is not a number" if (!looks_like_number($n));
	my ($TP,$TN,$FP,$FN,$TPR,$SPC,$PPV,$ACC) = parse_cpat::test_stat("$outfile.V2.cpat", $n);
	return ($TP,$TN,$FP,$FN,$TPR,$SPC,$PPV,$ACC, $cutf1,$cutf2)

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
	my @args = @_;
	my %opts = @args;
	my $CodingORF = $opts{"-ORF"};
	my $CodingmRNA = $opts{"-mRNA"};
	my $NonCoding = $opts{"-noncoding"};
	my $ftest = $opts{"-test"};
	my $outfile = $opts{"-out"};


	print "build the hexamere table\n";
	
	# build the hexamere table
	my $hexamere = `make_hexamer_tab.py -c $CodingORF -n $NonCoding  >"$outfile"_hexamer.table`;
	
	print "build the logit\n";
	# build the logit
	my $logit = `make_logitModel.py -x "$outfile"_hexamer.table -c $CodingmRNA -n $NonCoding -o $outfile`;
	
	print "use the CPAT\n";
	# use CPAT
	my $cpat = `cpat.py -d $outfile.logit.RData -x "$outfile"_hexamer.table -g $ftest  -o $outfile.cpat`;
	print "use R\n";
	# R: build a plot and summary
	my $command = `Rcpat.r $outfile.cpat`;
}

1;
