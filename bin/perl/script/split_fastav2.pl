#!/usr/bin/perl -w

#
# February 2014
#
# alexandre.nicaise28@gmail.com
#
# Split a multi fasta file in X fasta 
#
# run 'split_fasta.pl -h' for help
#
##########################################################################################


############################ library #####################################################
use warnings;
use strict;
use Bio::SeqIO;
use Getopt::Long;
use Pod::Usage;
use Data::Dumper;

############################# Parameters #################################################


# mfasta file
my $infile	= "";
# output prefix file
my $prefix 	= 	"split";
# number of sequences per output file(s)
my $number	= 	1 ;
# number of line
my $nb_line = 0;
# Identifiant
my $id = "";
# help and man
my $help	=	0;
my $man		=	0;


################################ GetOptions###############################################
GetOptions (
	"infile|i=s" 	=>	\$infile,
	"prefix|p=s"	=>	\$prefix,
	"number|n=i"	=>	\$number,
	"l=s"			=> \$nb_line,
	"id=s"				=> \$id,
	"help|h"	=>	\$help,
	"man|m"		=>	\$man
)
	or pod2usage ("Try '$0 --help' for more informations");

# Help and man
pod2usage ( -verbose => 1) if $help;
pod2usage ( -verbose => 2) if $man;

pod2usage ("-i option: I need a non-empty input file...\n") unless (-s $infile);

################################## function ##############################################


sub write_fasta{
	my ($in, $number, $Div, $nb_line) =@_;
	my $count = 0;
	my $fcount = 0;
	my $out;
	
	while ($count <=$nb_line and my $seq = $in->next_seq) {

        if ($count % $Div == 0 and $fcount != $number) {
        		if ($fcount != 0){
        			$out->write_seq($seq);
        			$count++;
				}
				print "file $prefix.$fcount.$count.fa done\n";
                $fcount++;
                $out = new Bio::SeqIO(-file => ">$prefix.$fcount.$Div.fa", -format=>'fasta');
                print "$prefix.$fcount.$Div.fa\n";
                $count = 1;
        }
        $out->write_seq($seq);
        $count++;
        if ($count == $Div+1 and $fcount == $number){
        	last;
        }
	}

}

sub write_seq{
	my ($in, $id)=@_;
	my $out = new Bio::SeqIO(-file => ">$id.fa", -format=>'fasta');
	while (my $seq = $in->next_seq){
		if ($seq->id =~ $id){
			$out->write_seq($seq);
			print "coucou\n";
			last;
		}
	}
}

################################# Main ###################################################
# Parse input file with bioperl
my $in  = new Bio::SeqIO(-file  => $infile);


if ($id ne ""){
	write_seq($in, $id);
}
elsif ($nb_line == 0){
	my $nb_line = `more $infile | grep ">"| wc -l `;
	print "number of file : $number\nnumber of sequence of the file: $nb_line";

	my $Div = $nb_line/$number;
	write_fasta($in,$number,$Div, $nb_line);
}


else {
	write_fasta($in, $number, $nb_line);
}




__END__

=head1 NAME

split_fasta.pl - Split a multi-fasta file in N files

=head1 SYNOPSIS

perl split_fasta.pl -i <MultiFasta_File> [Options]

	Required parameters
		-i <MultiFasta_File>
	
	Optional Parameters
		-p prefix of the output files			[default 'split']
		-n number of soutput files 				[default '1']
		-h help message							[default '0']
		-m man message							[default '0']

=head1 Optional parameters

=item B<--prefix (or -p)>

The string prefix of the output files (default 'split')

=item B<--number (or -n)>

The desired number of output file(s) (default '1')

=item B<--help (or -h)>

This Help message

=item B<--man (or -m)>

Th manual page

=cut
