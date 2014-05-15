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


my ($seqfile, $number) = @ARGV;
 
my $in = Bio::SeqIO->new(-file => $seqfile);
my $fh = Bio::SeqIO->newFh(-format => 'fasta');

# print Dumper \$in;

# die;

while (my $seq = $in->next_seq ){

	
	my @chars = split '', $seq->seq;
 
 	print Dumper @chars;
 	
	for my $i (1 .. $number) {
    	@chars 		= shuffle @chars;
	    my $newid	= $seq->primary_id."_shuffle".$i;
    	my $new_seq = Bio::Seq->new(-id => $newid, -seq => join '', @chars);
	    print $fh $new_seq;
	}
}
