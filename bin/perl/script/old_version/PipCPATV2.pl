#!/usr/bin/perl -w

##############################################################################################################
# 
# PROGRAM 	PipCPAT
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
my $Number_seq = "INF";

# file of non coding gene to stat
my $stat= "";

# table hash
my %hash_stat = ();

############################# Get Options ######################################

GetOptions(
	'c=s'				=> \$Coding,
	'nc=s'				=>  \$NonCoding,
	'test|t=s'			=> \$ftest,
	'outfile|o=s'		=> \$outfile,
	's=s'				=> \$stat,
	'n=f'				=> \$n,
	'ns=s'				=> \$nshuffle,
	'nseq=i'			=> \$Number_seq,
	'help|?|h'			=> \$help,
	'man'				=> \$man,
	)
	or pod2usage(2);
print @ARGV,"\n";
pod2usage(1) if $help;
pod2usage(-verbose =>2) if $man;


################################ Function ######################################

# shuffle sequence of a fasta file
sub fasta_shuffle{
	my ($seqfile, $number, $num) = @_;
	my $in = Bio::SeqIO->new(-file => $seqfile);
	my %hash_seq;
	#print Dumper \$in;
	#die;
	my $c = 0;

	while (my $seq = $in->next_seq){
		my @chars = split '', $seq->seq;
		$c++;
		#print "$c\n";
		#print Dumper @chars;
		for my $i (1 .. $number) {
			@chars 		= shuffle @chars;
			my $newid	= $seq->primary_id."_shuffle".$i;
			#my $new_seq = Bio::Seq->new(-id => $newid, -seq => join '', @chars);
			$hash_seq{$newid} = join '', @chars;
			#print ("coucou$i\n")
		}
		if ($c eq $num){last;}
	}
	
	for (my $i=1; $i <= $number;$i++){
		print ("file $i\n");
		my $seqOUT =new Bio::SeqIO ->new(-format => 'fasta', -file => ">shuffle$i.fa");# create a file in fasta format
		foreach my $ind (keys(%hash_seq)){
			my $new_seq = Bio::Seq -> new(-id => $ind, -seq => $hash_seq{$ind}); # create a new sequence
			$seqOUT->write_seq($new_seq);# write in the file
		}
	}
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


# Make some analysis of the result if we have the file where we know if coding gene and noncoding 
sub Stat{
	my ($Coding, $stat, $outfile)=@_; 
	# Verification 
	print "entrer dans les statistiques\n";
	my $verif = `CPATverif.pl $Coding $stat $outfile.out $outfile.dat`;
	my $verif2 = `CPATverifV2.pl -c $Coding -nc $stat -t $outfile.out -o temp.$outfile.V2.dat`;
	print "debut du shuffle\n";
	# shuffle
	my $temp_header = `awk 'NR==1' temp.$outfile.V2.dat > temp_header`;
	my $temp_corps = `awk 'NR>1' temp.$outfile.V2.dat > temp_corps`;
	my $temp_shuffle_corps = `shuffle.pl -s 10 temp_corps > temp_shuffle_corps`;
	my $shuffle = `cat  temp_header temp_shuffle_corps > $outfile.V2.dat`;
	print "fin du shuffle\n";
	my $rm = `rm temp*`;
	# VP 
	print "calcule des VP\n";
	my $VPCodant = `awk 'NR>1 && \$NF ==1{i++; if (\$(NF-1) > $n){cpt++}}END{print cpt/i}' $outfile.dat`;
	my $VPNonCodant = `awk 'NR>1 && \$NF ==0{i++; if (\$(NF-1) < $n){cpt++}}END{print cpt/i}' $outfile.dat`;

	# R
	print "R\n";
	my $commandeR = `Human_10fold_crossValidation_allFeatureV2.r  $outfile.V2.dat`;
	return ($VPCodant, $VPNonCodant, $commandeR)
	# results
# 	print "\nRESULTS\n\n";
# 	print "VP: \n";
# 	print "VP Codant: $VPCodant";
# 	print "VP Non codant: $VPNonCodant\n";
# 	print "cuttofs:\n $commandeR\n";
}





################################ Main ##################################################

#cpat: my ($Coding, $NonCoding, $ftest, $outfile)
# fasta_shuffle: my ($start_file, $Number_outfile, $Number_seq)
#  
if ($NonCoding eq "" ){
	fasta_shuffle($Coding, $nshuffle, $Number_seq);
	for (my $i=0; $i< $nshuffle; $i++){
		cpat($Coding, "shuffle$i",$ftest,"$outfile$i");
		if (Stat ne ""){
			@hash_stat{"Stat$i"} = Stat($Coding, $stat, "$outfile$i")
		}
	}
}
else{
	cpat($Coding, $NonCoding,$ftest,$outfile$i)
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

