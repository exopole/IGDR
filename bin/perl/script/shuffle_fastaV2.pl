#!/usr/bin/perl -w

# $ shuffle.pl input_sequence.fasta 20 > random_sequences.fasta
# Produce $number (20) shuffled sequences from a (multi)fasta file
# ie the resulting sequences will have
#	- same length as input_sequence.fasta
#	- same composition
###################################################################

use List::Util qw(shuffle);
use Bio::SeqIO;
use Data::Dumper;


my ($seqfile, $number, $num) = @ARGV;
 
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
	foreach $ind (keys(%hash_seq)){
		# print ("$ind\n");
# 		print ("$hash_seq{$ind}\n");
		my $new_seq = Bio::Seq -> new(-id => $ind, -seq => $hash_seq{$ind}); # create a new sequence
		$seqOUT->write_seq($new_seq);# write in the file
	}
}