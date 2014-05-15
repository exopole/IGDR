#!/usr/bin/perl -w


#
# VP.pl
#
# Goal : calcul of True positiv of 
#
#
#
#
##########################################################################################

############################### Library ##################################################
use strict;
use warnings;
use Data::Dumper;
use Getopt::Long;
use Pod::Usage;

############################### Variables ################################################

# Input file
my $infile = "";

# incrementation of the cuttof
my $n = 1;

# start
my $begin =0;

# end 
my $end = 1;

# outputfile
my $output = "VP";

# man, help
my $help = 0;
my $man = 0;

############################### Get options ##############################################

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



############################### Function #################################################

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
	print ("$VPC\n");
	return ($VPC, $VPNC);	
		
}

################################ Main #################################################

open OUT, ">", $output;
print OUT "Cutoff\tVPCodant\tVPNonCodant\n";
for ($begin; $begin<=$end	; $begin=$begin+$n  ){
	my ($VPC,$VPNC) = VP($infile,$begin);
	#print scalar(@VP),"\n";
# 	print "$VPC\n";
# 	print "$VPNC\n";
	print OUT "$begin\t", $VPC,"\t",$VPNC,"\n";
}


########################################################################################



__END__

=head1 NAME

VP.pl -i file, give true positive

=head1 AUTHOR

Nicaise Alexandre : alexandre.nicaise28@gmail.com

=head1 SYNOPSIS

perl VP -i <Result CPATverif.pl> [Options...]


=head1 DESCRIPTION

=over 8

=item * get a table with True postive of non coding and coding

Options:

	-i 			:  Input file of result of CPATverif

	-outfile|o	: output file prefix [default : VP]
	
	-n			: step of cutoff [default : 1]
	
	-help|h		: Help message [default : 0]

	-man|m		: man help [default : 0]

	
=head1 Details

=over 8

=item B<-i>
: file of result of CPATverif.pl 

=item B<-n>
: step of the cutoff

=item B<-o>
: output file prefix

=item B<-man>
: Print man help and exit

=item B<-help>
: Print help message and exit

=back

=cut
