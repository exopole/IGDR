#!/usr/bin/perl -w


# Aim : adapted from intersedBed from the Bedtools
###############################################################

=head1 NAME

intersectRetouch.pl - (inspired by intersectBed from Bedtools) : Intersect 2 .gtf files at the transcript level.

=head1 SYNOPSIS

intersectRetouch.pl  -a infileA.gtf -b infileB.gtf [options]

=head1 OPTIONS

=over 8

=item B<-a|infileA>

Input .gtf fileA. [mandatory]

=item B<-b|infileB>

Input .gtf fileB. [mandatory]
	
=item B<-s|strand>

Require same strandedness for overlap. [default no strand required]

=item B<-f|fraction>

Minimum overlap required as a fraction of A [default : overlap > 0 ].

=item B<-r|reciprocal>

Require that the fraction overlap be reciprocal for A and B.

=item B<-p|proc>

Number of threads  [default 4].

=item B<-v|verbosity>

Level of verbosity [default 0].

=item B<-help>

Print a brief help message and exits.

=item B<-man>

Prints the manual page and exits.

=back

=head1 DESCRIPTION

B<This program> will parse 2 .gtf files and return the fraction at the transcript level

=cut

# Perl libs
use warnings;
use strict;
use Getopt::Long;
use Pod::Usage;
use File::Basename;
use Data::Dumper;
use Cwd;
use File::Temp qw/ tempfile tempdir /;
use Parallel::ForkManager; # in /local/perl/5.12.2/lib/site_perl/5.12.2/Parallel/ForkManager.pm


# Own lib : /home/genouest/umr6061/recomgen/tderrien/bin/ThomasPerl/lib/
use Parser;
use ExtractFromHash;
use Intersect;
use Utils;

my $progname=basename($0);

# Variables
my $infileA='';
my $infileB='';
my $strand=0;
my $fraction=0;
my $reciprocal=0;
my $man = 0;
my $help = 0;
my $verbosity=0;
my $processus = 4;

## Parse options and print usage if there is a syntax error,
## or if usage was explicitly requested.
GetOptions('a|infilea=s' => \$infileA,
   'b|infileb=s' => \$infileB,
   's|strand' => \$strand,
   'f|fraction=f' => \$fraction,
   'r|reciprocal' => \$reciprocal,
   'p|proc=i' => \$processus,   
   'v|verbosity=i' => \$verbosity,
   'help|?' => \$help,
   'man' => \$man
) or pod2usage(2);

# Help and man
pod2usage(1) if $help;
pod2usage(-verbose => 2) if $man;

# Test parameters
pod2usage("Error: Cannot read your input .gtf file A '$infileA'...\nFor help, see:\n$progname --help\n") unless( -r $infileA);
pod2usage("Error: Cannot read your input .gtf file B '$infileB'...\nFor help, see:\n$progname --help\n") unless( -r $infileB);

pod2usage ("- Error: \$fraction option '$fraction' should be a float between 0 and 1 [0-1] (e.g 0.5 if 50% overlap is required)\n") unless ($fraction >= 0 and $fraction <= 1);

# Create a temp directory to store all chr files for both infiles
my $templateA = "A_XXXX";
my $dirnameA = tempdir($templateA, CLEANUP => 1);
if ($verbosity > 5){
	print STDERR "Splitting '$infileA'...\n";
}
splitFilebyChr($infileA,$dirnameA);

my $templateB = "B_XXXX";
my $dirnameB = tempdir($templateB, CLEANUP => 1);
if ($verbosity > 5){
	print STDERR "Splitting '$infileB'...\n";
}
splitFilebyChr($infileB,$dirnameB);


######## get all chr files and sort by descending size
my @filesA = glob "$dirnameA/*.out";
@filesA = sort { -s $b <=> -s $a } @filesA;

my @filesB = glob "$dirnameB/*.out";
@filesB = sort { -s $b <=> -s $a } @filesB;

# Temp dir to write all output
my $templateOUT = "OUT_XXXX";
my $dirnameOUT = tempdir($templateOUT, CLEANUP => 1);

# Launch forks
my $pm = Parallel::ForkManager->new($processus);

# Process 2 directories
foreach my $fileA (@filesA ) {

	
	# start fork
    my $pid = $pm->start and next;

    foreach my $fileB (@filesB ) {
    
	    my $fileA_pref = basename($fileA);
        my $fileB_pref = basename($fileB);

        next unless ($fileB_pref eq $fileA_pref);
        
        my %h1  = Parser::parseLevelGtfhashKey($fileA, 'exon', $verbosity);
#         if ($verbosity>0){print STDERR "\n"}
        my %h2  = Parser::parseLevelGtfhashKey($fileB, 'exon', $verbosity);
#         if ($verbosity>0){print STDERR "\n"}

		my $outputfile = $fileA_pref."_".$fileB_pref;
        # Launch Intersect
        Intersect::Intersect2HspliteqCHR(\%h1, \%h2, $strand, $fraction, $reciprocal, $verbosity, $dirnameOUT, $outputfile);
    }
    # end fork
    $pm->finish;
}
# wait all sub process
$pm->wait_all_children;

###################
# Print Results
my @filesOUT = glob "$dirnameOUT/*.out";

foreach my $fileOUT (@filesOUT ) {
	open  my $outFile, '<', $fileOUT or die "Unable to open $fileOUT for reading : $!";
	while (<$outFile>){print;}
	close $outFile;
}



#-------------------
my %file_handles = ();
sub get_file_handle {
   my $chr = shift;
   my $dir = shift;
   
   my $key = $dir."_".$chr;
   
   # if the file handle already exists, then we don't need to create a new one
   unless (exists ($file_handles{$key})) {
       # open a file handle (hostname.out) to write to and store it in the %file_handles hash
       open $file_handles{$key}, '>', $dir.'/'.$chr.'.out' or die "ERROR: Unable to open $dir/$chr.out: $! \n";
   }

   # return the file handle
   return $file_handles{$key};
}

sub splitFilebyChr {

	my $infile	=	shift; # infile o be splited
	my $dir 	=	shift;	# directory to put output splited file # default cwd
	
	die unless (-s $infile);
	
	open (IN, $infile) or die "Cannot open '$infile'\n";
	while (my $line = <IN>) {
	   # remove new line from end of $line
	   chomp $line;

	   # split the line on pipes but only return the first field
	   my $chr = ( split(/\t/, $line) )[0];

	   # retrieve the file handle to write to from the hash based on the chr name
	   my $fh = get_file_handle($chr, $dir);

	   # write the line to the file handle (files)
	   printf $fh "%s\n", $line;
	}
}