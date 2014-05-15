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

###################################### Library #######################################

# Libraries
use strict;
use warnings;

# For parsing options
use Getopt::Long;
use Pod::Usage;

# to use liste
use List::Util qw(shuffle);

# to use Bio::perl
use Bio::SeqIO;

use Utils;
##################################### Variables ###################################

# file of coding gene
my $Coding;

# file of non coding gene to training
my $NonCoding = "";

# file of test
my $ftest;

# ouput file name
my $outfile;

my $help =0;
my $man = 0;
my $n = 0.36;

# number of shuffle file 
my $nshuffle =1;
my $shufflesize = 0;
my $shufflediv = 1;

# number of sequence shuffle put in the file
my $Number_seq = "INF";

# file of non coding gene to stat
my $stat= "";



############################# Get Options ######################################

GetOptions(
	'c=s'				=> \$Coding,
	'nc=s'				=>  \$NonCoding,
	'test|t=s'			=> \$ftest,
	'outfile|o=s'		=> \$outfile,
	's=s'				=> \$stat,
	'n=f'				=> \$n,
	'ns=s'				=> \$nshuffle,
	'size|z=i'			=> \$shufflesize,
	'divide|d=i'		=> \$shufflediv,
	'nseq=i'			=> \$Number_seq,
	'help|?|h'			=> \$help,
	'man'				=> \$man,
	)
	or pod2usage(2);
pod2usage(1) if $help;
pod2usage(-verbose =>2) if $man;


################################ Function ######################################

# shuffle sequence of a fasta file
sub fasta_shuffle{
	print ("file of non coding gene undefined so: shuffle of the fasta file of coding gene\n");
	my ($seqfile, $number, $num,$size, $div) = @_;
	my $in = Bio::SeqIO->new(-file => $seqfile);
	my @OUT = [];
	my $c = 0;
	my @chars = ();
	my @chars2 = ();
	my $size2 = $size;
	while (my $seq = $in->next_seq){
		@chars = split '', $seq->seq;
		if ($size == 0){
			$size2 = int(@chars/$div) ;
# 			print "size: $size\n";
# 			print "length:", length(@chars),"\n";
# 			my $len = @chars;
# 			print "len: $len\n";
# 			print "seq: @chars[1..10]\n";
		}
# 		if ($size < 200){
# 			$size =200 ;
# 		}

		$c++;
		for my $i (1 .. $number) {
			@chars		= shuffle @chars;
# 			print "size: $size2\n";
# 			my $len = @chars;
# 			print "len: $len\n";
# 			print "div: $div\n" ;
# 			print @chars[0..$size2],"\n";
			
			my $joint = join '', @chars[0..$size2];
# 		print "$size\n";
# 		print "joint: $joint\n";
 		
			my $newid	= $seq->primary_id."_shuffle".$i;
			my $new_seq = Bio::Seq->new(-id => $newid, -seq => $joint);
			#$hash_seq{$newid} = join '', @chars;
			#print ("coucou$i\n")
			if ($c == 1){
				my $seqOUT =new Bio::SeqIO ->new(-format => 'fasta', -file => ">shuffle$i.fa");
				push (@OUT, $seqOUT);
			}
			
			$OUT[$i]->write_seq($new_seq);
		}
		if ($c eq $num){last;}
	}
	
 }


# shuffle file of result
sub shuffle_result{
		
	my $in = shift;
	my  $c = 0;
	my $Header = "";
	open IN, "temp.$in.V2.dat" or die "Cannot open $in\n";
	
	open CORP, ">", "temp_corps";
	open OUT, ">", "$in.V2.dat";

	while (<IN>){
		$c++;
		if ($c==1){
			$Header = $_;
		}
		else{
			print CORP $_;
		}
	}

	my $temp_shuffle_corps = `shuffle.pl -s 10 temp_corps > temp_shuffle_corps`;
	open TEMP, "temp_shuffle_corps" or die;
	print OUT "$Header";
	while (<TEMP>){
		#print "$_";
		print OUT $_;
	}
	close (OUT);
	close (TEMP);
	close (CORP);
	close (IN);
} 


# uses all commande of cpat and a small script r
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
	my $cpat = `cpat.py -d $outfile.logit.RData -x "$outfile"_hexamer.table -g $ftest  -o $outfile.out`;
	
	print "use R\n";
	# R: build a plot and summary
	my $command = `Rcpat.r $outfile.out`;
}

sub VP {
	my ($in, $cutoff) =@_;
	my $cpt1 =0;
	my $cpt0 =0;
	my $cpt_codant = 0;
	my $cpt_NonCodant = 0;
	my @fields;
	my $Cod_prob;
	my $label;
	my $count =0;
	open IN, $in or die "Cannot open $in\n";
	while (<IN>){
		#print ("coucou: $count\n");
		if ($count ne 0){ 
			@fields = split (/\t/);
			#print ("coucou\n");
			$Cod_prob = $fields[-2];
			$label = $fields[-1];
			#print ("$label\n");
			if ($label == 0){ 
				$cpt0++;
				if ($Cod_prob<$cutoff){
					$cpt_NonCodant++;
				}
			}
			if ($label == 1){
				#print ("$cpt1\n");
				$cpt1++;
				if ($Cod_prob>$cutoff){
					$cpt_codant++;
				}
			}
		}
		$count++;
	}
	my $VPC = $cpt_codant/$cpt1;
	my $VPNC = $cpt_NonCodant/$cpt0;
	return ($VPC, $VPNC);	
		
}


# Make some analysis of the result if we have the file where we know if coding gene and noncoding 
sub Stat{
	my ($Coding, $stat, $outfile, $n)=@_;
	#print ("apres entrer dans stat: $outfile\n"); 
	# Verification 
	print "entrer dans les statistiques\n";
	my $verif = `CPATverif.pl $Coding $stat $outfile.out $outfile.dat`;
	my $verif2 = `CPATverifV2.pl -c $Coding -nc $stat -t $outfile.out -o temp.$outfile.V2.dat`;
	print "debut du shuffle\n";
	# shuffle
	shuffle_result($outfile);
	print "fin du shuffle\n";
	my $rm = `rm temp*`;
	# VP 
	print "calcule des VP\n";
	#print "awk 'NR>1 && \$NF ==1{i++; if (\$(NF-1) > $n){cpt++}}END{print cpt/i}' $outfile.dat\n";
# 	my $VPCodant = `awk 'NR>1 && \$NF ==1{i++; if (\$(NF-1) > $n){cpt++}}END{print cpt/i}' $outfile.dat`;
# 	my $VPNonCodant = `awk 'NR>1 && \$NF ==0{i++; if (\$(NF-1) < $n){cpt++}}END{print cpt/i}' $outfile.dat`;
# 	chomp ($VPCodant);
# 	chomp ($VPNonCodant);
	my @VP= VP("$outfile.dat", $n);
	my $VPCo = $VP[0];
	my $VPNCo = $VP[1];
	# R
	print "R\n";
	my @commandeR = `Human_8fold_crossValidation_allFeatureV2.r  $outfile.V2.dat`;
	my $cutf1 = $commandeR[1];
	my $cutf2 = $commandeR[2];
	chomp ($cutf1);
	chomp ($cutf2);
	 $cutf1 =~ s/\[1\]//;
	 $cutf2 =~ s/\[1\]//;
	return ($VPCo, $VPNCo, $cutf1,$cutf2)

}





################################ Main ##################################################

#cpat: my ($Coding, $NonCoding, $ftest, $outfile)
#fasta_shuffle: my ($start_file, $Number_outfile, $Number_seq)
#  
if ($NonCoding eq "" ){
	fasta_shuffle($Coding, $nshuffle, $Number_seq, $shufflesize, $shufflediv);
	open RES,">", "$outfile.stat";
	print RES "file\tVP_coding\tVP_NonCoding\tcuttof_accuracy\tcuttof_performance\n";
 	for (my $i=1; $i<= $nshuffle; $i++){
 		#print ("shuffle$i\n");
 		#cpat($Coding, "shuffle$i.fa",$ftest,"$outfile$i");
 		if ($stat ne ""){
 			print ("shuffle$i analyse\n");
 			#print ("avant entrer dans stat: $outfile$i\n");
 			my ($VPCodant, $VPNonCodant, $cut1,$cut2) = Stat($Coding, $stat, "$outfile$i", $n); 
 			print RES "$outfile$i\t$VPCodant\t$VPNonCodant\t$cut1\t$cut2\n";
 		}
 	}

}
else{
	cpat($Coding, $NonCoding,$ftest,$outfile);
	if ($stat ne ""){
		open RES,">", "$outfile.stat";
		print RES "file\tVP_coding\tVP_NonCoding\tcuttof_accuracy\tcuttof_performance\n";
		my ($VPCodant, $VPNonCodant, $cut1,$cut2) = Stat($Coding, $stat, $outfile, $n); 
		print RES "$outfile\t$VPCodant\t$VPNonCodant\t$cut1\t$cut2\n";
	}
}








__END__

=head1 NAME

PipCPAT.pl - run all options for cpat

=head1 AUTHOR

Nicaise Alexandre : alexandre.nicaise28@gmail.com

=head1 SYNOPSIS

perl PipCPATV2.pl -c <coding fasta file> [Options...]


=head1 DESCRIPTION

=over 8

=item * use CPAT to know the coding potential of genes

=item * get a table with coding potential

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

