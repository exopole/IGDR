#!/usr/bin/perl

## nombre de fois ou on pipcpat (nombre de fichier de sequence et de shuffle 
my $number = shift;


for (my $i = 1; $i<=$number; $i++){
	for (my $x =1; $x<=$number;$x++){
		my $command = `time PipCPAT.pl -c Human_ORF.fa -nc shuffle$i.fa -t sequence$x.fa -o test$x_sh$i`;
	}
}
