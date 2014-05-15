#!/usr/bin/perl -w

#
# Separate_fasta.pl
# 
# Date : 24/01/14
# 
#
# Goal :
#	group 2 file in 1 
#
# update : 
#	- eliminate list to use hash 
#	- use sub
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

# mfasta file
my $fafile1	= "";
my $fafile2 = "";

# output prefix file
my $outfile 	= 	"group";
# number of sequences per output file(s)
my $number	= 	1 ;

# number of seq
my $seq_number = 0;

# help and man
my $help	=	0;
my $man		=	0;

# number of sequence by files
my $nd1;
my $nd2;
my %fasta2;


######################################### Get Options ####################################

# GetOptions
GetOptions (
	"coding|c=s" 		=>	\$fafile1,
	"noncoding|nc=s"	=>  \$fafile2,
	"prefix|p=s"		=>	\$outfile,
	"number|n=i"		=>	\$number,
	"s=i"	   			=>  \$seq_number, 
	"help|h"			=>	\$help,
	"man|m"				=>	\$man
)
	or pod2usage ("Try '$0 --help' for more informations");

# Help and man
pod2usage ( -verbose => 1) if $help;
pod2usage ( -verbose => 2) if $man;

pod2usage ("-c option: I need a non-empty input file...\n") unless (-s $fafile1);
pod2usage ("-nc option: I need a non-empty input file...\n") unless (-s $fafile2);


######################################## Function ########################################

sub parse_fasta{
	my ($fafile) =@_;
	my $in = Bio::SeqIO ->new(-file => $fafile);
	my %h=();
	my $cin=0;
	my $id;
	while (my $seq=$in -> next_seq){
		$cin++;
		$h{$seq->id} = $seq->seq;
	}
	return ($cin,\%h);
}


sub separate{
	my @args = @_;
	my %opts = @args;
	my $outfile;
	my %h1;
	my %h2;
	my $nd1;
	my $nd2;
	if ($opts{"-fasta1"}){
		%h1 = %{$opts{"-fasta1"}};
	}
	if ($opts{"-fasta2"}){
		%h2 = %{$opts{"-fasta2"}};
	}
	if (!$opts{"-fasta1"} and !$opts{"-fasta2"}){
		die "You don't give reference of hash";
	}
	my $number = $opts{"-number"};
	if ($outfile){	
		$outfile = $opts{"-outfile"}; 
	}
	else {
		$outfile = "group";
	}
	$nd1 = $opts{"-number_of_seq_f1"};
	$nd2 = $opts{"-number_of_seq_f2"};
	# Count for navigate in @seq and @id
	my $c1 = 0;
	my $c2 = 0;
	my @file=();
	my $i1 = 0;
	my $i2 = 0;

	
	# Creation of the output file with coding and non coding gene
	for (my $i=1; $i <= $number;$i++){
		open(FILE, ">$outfile.$i.fa");
		close(FILE);
			my $seqOUT =new Bio::SeqIO ->newFh(-file => ">$outfile.$i.fa", -format => "fasta");# create a file in fasta format
			push (@file, $seqOUT);
	}    	

	foreach my $k (keys(%h1)){
		my $OUT;
		if ($c1%$nd1==0 and $c1!=0){
			$i1++;
			print "coucou\n";
		}
		if (exists($file[$i1])){
			$OUT= $file[$i1];
		}
		else{
			print "last sequence of the first file\n";
			last;
		}
		print "c1: $c1 and nd1: $nd1\n";
		print "difference between c1 et nd1",$c1%$nd1, "\n";
		my $file = @file;
		print "nombre de file ",$file, "\n";
		my $new_seq = Bio::Seq -> new(-id => $k, -seq => $h1{$k},); # create a new sequence
		$OUT->write_seq($new_seq); 
		$c1++;
	}
	if (%h2){
		foreach my $k (keys(%h2)){
			my $OUT;
			if ($c2%$nd1==0 and $c2!=0){
				$i2++;
			}
			if(exists($file[$i2])){		
				$OUT= $file[$i2];
			}
			else{
				print "last sequence\n";
				last;
			}
			my $new_seq = Bio::Seq -> new(-id => $k, -seq => $h2{$k},); # create a new sequence
			$OUT->write_seq($new_seq);
			$c2++; 
		}
	}
}


#################################### Main ################################################
my $cin2;
my $ref;

my ($cin1, $ref2) = parse_fasta($fafile1);
my %fasta1 = %{$ref2};
if ($fafile2 ne ""){
	($cin2, $ref) = parse_fasta($fafile2);
	%fasta2 = %{$ref};
}
if ($seq_number == 0){
	# divide the file in $number
	$nd1 = int($cin1/$number);
	print ("We have $nd1 sequence of the second file in the output fasta files\n");
	if ($fafile2 ne ""){
		$nd2 = int($cin2/$number);
		print ("We have $nd2 sequence of the second file in the output fasta files\n");
	}
}
else {
	$nd1 = $seq_number;
	print ("We have $nd1 sequence of the second file in the output fasta files\n");
	if ($fafile2 ne ""){
		$nd2 = $seq_number;
		print ("We have $nd2 sequence of the second file in the output fasta files\n");
	}
}

separate(-fatsa1 => \%fasta1,-fasta2 => \%fasta2,-number => $number,-output =>  $outfile, -number_of_seq_f1 => $nd1,-number_of_seq_f2 => $nd2);



__END__

=head1 NAME

separate_fastav2.pl - regroup 2 file in n files

=head1 SYNOPSIS

separate_fastav2.pl -c <fasta_file1> -nc <fasta_file2>  [Options]

	Required parameters
		-c <fasta_file1>
		-nc <fasta_file2>
	
	Optional Parameters
		-p prefix of the output files			[default 'group']
		-n number of output files 				[default '1']
		-s number of sequence
		-h help message							[default '0']
		-m man message							[default '0']

=head1 Optional parameters

=item B<--prefix (or -p)>

The string prefix of the output files (default 'group')

=item B<--number (or -n)>

The desired number of output file(s) (default '1')

=item B<--help (or -h)>

This Help message

=item B<--man (or -m)>

Th manual page

=cut



