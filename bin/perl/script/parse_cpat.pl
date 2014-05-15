#!/usr/bin/perl -w

# parse_cpat.pl
#  
# Goal : calcul of the sensibility and specificiy for a cpatverifV2
# 
#
#
#
#
#
#########################################################################################

################################# Library ###############################################
use strict;
use warnings;
use Data::Dumper;
use Getopt::Long;
use Pod::Usage;


################################# Scalar ################################################


# Input file
my $infile = "";

# incrementation of the cuttof
my $n = 1;

# start
my $begin =0;

# end 
my $end = 1;

# outputfile
my $output = "parsing_cpat";

# man, help
my $help = 0;
my $man = 0;



############################## Get Options ###############################################

GetOptions(
    'infilec|i=s'            => \$infile,
	'n=f'		 			 => \$n,
	'b=f'					 => \$begin,
	'e=f'					 => \$end,
	'output|o=s'			 => \$output,
	'help|?|h'			     => \$help,
	'man' 			 		 => \$man,
	)	
or pod2usage("Try '$0 --help' for more informations" );

pod2usage(-verbose => 1) if $help;
pod2usage(-verbose => 2) if $man;


pod2usage ("-i option: I need a non-empty input file...\n") unless (-s $infile);





############################# Function ###################################################



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
		my ($id, $mRNA_size,$ORF_size, $Fickett_score,$Hexamer_score,$coding_prob,$Label) = split ("\t");



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

	my ($href, $cut) = @_;
	my %h	= %{$href};
	#$cut 	||= 	0.5;
	#$label 	||= 	0;	
# 	print Dumper \%h;
# 	die;
	my $TP=0;
	my $TN=0;
	my $FN=0;
	my $FP=0;
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

####################################### Main ############################################

my %h_parse = parse_cpat($infile);
open OUT, ">", $output;
print OUT "Cutoff\tTP\tFP\tTN\tFN\tSensibility\tSpecificity\tPrecision\tAccuracy\n";
for ($begin; $begin<=$end; $begin=$begin+$n  ){
# 	print Dumper \%h_parse;
# 	die;
	my  ($TP,$TN,$FP,$FN)= Cond(\%h_parse, $begin);
	my $TPR = Sensitivity($TP,$FN);
	my $SPC = Specificity($TN,$FP);
	my $PPV = Precision($TP,$FP);
	my $ACC = Accuracy($TP,$TN,$FP,$FN);
	print OUT "$begin\t$TP\t$FP\t$TN\t$FN\t$TPR\t$SPC\t$PPV\t$ACC\n"
}
close(OUT);
# print Dumper \%h_parse;