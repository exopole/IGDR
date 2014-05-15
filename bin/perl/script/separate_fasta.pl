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


# Libraries
use strict; 
use warnings;
use Bio::SeqIO; 


# file fasta, number
my ($fafile1,$fafile2, $outfile,$number, $ns)= @ARGV;
my $nd1;
my $nd2;

# open fasta file with SeqIO
my $in1 = Bio::SeqIO ->new(-file => $fafile1);
my $in2 = Bio::SeqIO ->new(-file => $fafile2);

$ns |= 0;
$number |= 1;	
# number of sequence in the file 1
my $cin1 = 0;
my @seq1= ();
my @id1 = ();
while (my $seq= $in1 -> next_seq){
	$cin1++;
	push (@seq1, $seq->seq);
	push (@id1, $seq ->id);
}
print ("number of sequence in the file 1: $cin1\n");
#print  $id1[0], "\n" , $seq1[0],"\n";

# number of sequences in the file 2
my $cin2;
my @seq2 = ();
my @id2 = ();
my %id_seq=();
while (my $seq=$in2 -> next_seq){
	$cin2++;
#	$id_seq{$seq->id}=$seq->seq;
	push (@id2,$seq->id); 
	push (@seq2,$seq->seq); 
}
print ("number of sequence in the file 2: $cin2\n");


# divide the file in $number
if ($ns == 0){
	$nd1 = int($cin1/$number);
	print ("We have $nd1 sequence of of the first file in the output fasta files\n");
	$nd2 = int($cin2/$number);
	print ("We have $nd2 sequence of the second file in the output fasta files\n");
}
else {
	$nd1 = $ns;
	$nd2 = $ns;
}
# Count for navigate in @seq and @id
my $c1 = 0;
my $c2 = 0;

# Creation of the output file with coding and non coding gene
for (my $i=1; $i <= $number;$i++){
	open(FILE,">$outfile$i.fa");
	close(FILE);
    my $seqOUT =new Bio::SeqIO ->new(-format => 'fasta', -file => ">>$outfile$i.fa");# create a file in fasta format
	for (my $x=0; $x<$nd1; $x++){
		my $new_seq = Bio::Seq -> new(-id => $id1[$c1], -seq => $seq1[$c1]); # create a new sequence
		$seqOUT->write_seq($new_seq);# write in the file
		$c1++;
	}
	for (my $x=0; $x<$nd2; $x++){
                my $new_seq = Bio::Seq -> new(-id => $id2[$c2], -seq => $seq2[$c2]); # create a new sequence
                $seqOUT->write_seq($new_seq); # write in the file
                $c2++;
    	
	
	}

print "file $i done/$number\r";
close(FILE);
}
print "file $number done/$number\n";