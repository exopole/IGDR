#!/usr/bin/perl -w

# AIM
# Exemple de programme 

# Libraries
use strict;
use warnings;

# For parsing options
use Getopt::Long;
use Pod::Usage;

# Very useful to print a data structure
use Data::Dumper;

# input  file 
my $vcffile;
my $gtffile;

# verbosity, help and man
my $verbosity	=	0;
my $help		=	0;
my $man			=	0;

# Parsing parameters
GetOptions (
    "vcf|v=s"           =>	\$vcffile,
	"gtffile|g=s"		=>	\$gtffile,
	"verbosity" 		=> \$verbosity,
	"help|h"			=> \$help,
	"man|m"				=>	\$man,
)
    or pod2usage ("Try '$0 --help' for more information");


# Parsing parameters 
if ($verbosity>0){ print STDERR "Parsing parameters...\n";}

pod2usage ( -verbose => 1) if $help;
pod2usage ( -verbose => 2) if $man;
pod2usage ("I need an input .vcf file (option -v)...\n") unless (defined $vcffile);
pod2usage ("I need an input .gtf file (option -g)...\n") unless (defined $gtffile);

print STDERR "Checking options end ...\n";

# Intersect command from BEDTOOLS
my $command="intersectBed -a $vcffile -b $gtffile -wa -wb -loj |";
print STDERR "$command\n";
open(IN, $command) || die "can't run $command";

my %h =();
while (<IN>){
    
    next if(/^##/); #ignore header
    chomp;
    
    my %attribs = ();

    # parse line
    my ($chr, $pos, $snpid, $ref, $alt, $qual, $filter, $info, 
    $chrgtf, $source, $type, $start, $end, $score, $strand, $phase, $attributes) = split("\t");

    # store SNP data
    $h{$snpid}->{"chr"}			=   $chr;
	$h{$snpid}->{"pos"}	    	=   $pos;
	$h{$snpid}->{"REF"}	    	=   $ref;
	$h{$snpid}->{"ALT"}         =   $alt;
			
    # parse gtf attributes
    my @add_attributes = split(";", $attributes);
    # store ids and additional information in second hash
    foreach my $attr ( @add_attributes ) {
        next unless $attr =~ /^\s*(.+)\s(.+)$/;
        my $c_type  = $1;
        my $c_value = $2;
        $c_value=~ s/\"//g;
       if(!exists($attribs{$c_type})){
         warn "$attribs{$c_type} already exists...\n";
       }
               if($c_type  && $c_value){
            $attribs{$c_type} = $c_value;
        }
    }
    # if (%attribs) tests whether %attrib is not empty
    push (@{$h{$snpid}->{"OVERLAP"}}, \%attribs) if (%attribs);

}
# To print data structure after parsing
# print Dumper %h;

# Extract unique gene data overlap in hash '%genes'
#   - gene_id
#   - gene_biotype
#   - gene_name
for my $snp (keys %h){
#     print $snp, scalar @{$h{$snp}->{"OVERLAP"}}, "\n";
    
    my %genes;
    
    foreach my $feat (@{$h{$snp}->{"OVERLAP"}}) {
        if (exists ($feat->{'gene_id'})){
            my $name;
            my $biotype;

            if (!exists($feat->{'gene_name'})) {$name = "NA"} else {$name = $feat->{'gene_name'}}
            $genes{$feat->{'gene_id'}}->{"gene_name"}    = $name;
            
            if (!exists($feat->{'gene_biotype'})) {$biotype = "NA"} else {$biotype = $feat->{'gene_biotype'}}
            $genes{$feat->{'gene_id'}}->{"gene_biotype"} = $biotype;        
        }
    }
    push (@{$h{$snp}->{"GENES"}}, \%genes) if (%genes);

}

# To print data structure after data extraction 
# print Dumper %h;
