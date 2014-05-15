package shuffle;
$| = 1;
 

# Shuffle performs a random permutation on the lines from the standard input,
# with possible constraints on successive lines.

# Version 1.0 (perl)
# Author: Christophe Pallier (pallier@lscp.ehess.fr)
# Date: 13 Mar. 1999 (awk) 20 July 2000 (perl)
#
# cf. www.pallier.org for a paper describing "shuffle" in details.
#
# This program is copyrighted under the terms of the GNU 
# license (see http://www.gnu.org).
#
# Options 
#   -constraint 'n1 n2 ... '  constraints
#   -seed xxx           seed
#   -number_line xxx           number of output lines 
#   -iteration xxx           maximum number of iterations
#   -e               turns on the algorithm "equiprobable permutations" 

use Getopt::Std;
use Data::Dumper;
use strict;
use warnings;




sub new_shuffle {
	$[ = 1;			# set array base to 1
	$, = ' ';		# set output field separator
	$\ = "\n";		# set output record separator

	my $ExitValue = 0;
    my (@args) = @_;
    my %opts = @args;
    @opts{ map { lc $_ } keys %opts } = values %opts;

	# set randomizer's seed 
	if ($opts{'-seed'}) {
		srand($opts{'-seed'});
	} 

	##############################################
	# reads the input lines and store them in 'table'
	# eliminating empty lines

	my @table=();
	my @table2;
	my $header;
	
	if ($opts{'-file'}){
		open FILE, $opts{'-file'}, or die "cannot open ".$opts{'-file'}."\n";
		my $cpt = 0;
		while (<FILE>) {
		  chomp;
		  if ($opts{'-header'} && $cpt == 0){
		  	$header = $_
		  }
		  else{
		  		push(@table, $_) unless ($_ eq ""); 
			}	
		}
		close(FILE);
	}
	
	elsif ($opts{'-table'}){
		@table = @{$opts{'-table'}};
		if (defined $opts{'-header'}){
			$header = $table[0];
			my $len = @table;
# 			print "longueur de la liste: $len\n";

			splice @table, 1,1;
# 			print "header: $header\n";
#  			print Dumper \@table;
		}
# 		print Dumper \%opts;
# 	print Dumper \@table;

	 		
	}

	if ($opts{'-output'}){
		open OUT,">" ,$opts{'-output'} or die "Cannot write in $opts{'-output'} file\n";
	}

	my $nlines=$#table;
	my $n=$nlines;

	if ($opts{'-number_line'}) { $n=$opts{'-number_line'}; }

	my $iter=$nlines;
	if ($opts{'-iteration'}) { $iter=$opts{'-iteration'}; }
	if (!$opts{'-constraint'}) { 
# 		print Dumper \@table; 
		permute(\@table,$n); # generate an unconstrained, random, permutation 
# 		print Dumper \@table;

		if ($opts{'-header'}){
			splice @table, 1, 0, $header; 
		}
		
# 		print Dumper \@table;
		if ($opts{'-output'}){
			for (my $i=$nlines-$n+1;$i<=$nlines;$i++){
				if ($table[$i]){ 
					print OUT $table[$i];
				}
			}
			close(OUT);
		}
		return @table;
	
	}
	else
	{
		my $constr=$opts{'-constraint'};

		if ($opts{'-equiprobable'})
		{ # generate random permutations until one fullfills the constraints
		  # or the maximum number of iterations is reached.
			my $loop=1; 
			do { permute(\@table,$#table); } 
			until (filter(\@table,$constr,$n) or ($loop++>$iter));
			if ($loop>$iter) { $ExitValue=1; }
		}
		else 
		{ # build a permutation that fullfills the constraints
			shuffle(\@table,$constr,$n,$iter, $nlines) || ($ExitValue=1);
		}


		if ($ExitValue==0) {
		   for (my $i = 1; $i <= $n; $i++) 
		   {
				if (!$table[$i]){
					splice @table, $i ,1;
				}
				elsif ($opts{'-output'} && $table[$i]){
					print OUT $table[$i];
				}
					print $table[$i];
		   	}
		   return @table;
		   
		if ($opts{'-output'}){
			close(OUT);
			print "fermeture\n";
		}
		}
		else { print STDERR "Shuffle could not find a solution. Try again"; }
	}
}


########################## subroutine shuffle ###############
sub shuffle {
    my @table = shift;
    my $constr= shift;
    my $n = shift;
    my $iter = shift;
    my $nlines = shift;
# 	print Dumper \@table;
	# return 0 if no solution, 1 else. 
	#  print @table;

	  if ($constr eq "") {
		  permute(\@table,$n);
		  return 1;
	  }

	my @constraint = split(' ', $constr);

	my $loop = 0;
	my $MaxLoops = $iter ? $iter : $nlines;
	my $output_nlines = $n ? $n : $nlines;

	my @rep;
	my @prev;
	my @current;
	my @tmp;
	my ($nf,$pass_line,$fail,$bad_permut);
	my ($i,$j,$k,$back,$ll);

	do {
		permute(\@table, $#table);
		$bad_permut = 0;    # will be set to 1 if the permutation is bad...
		$loop++;

		# scan the table line by line
		# swapping lines to try and respect the constraints

		$nf = (@prev = split(' ', $table[1]));

		for ($k = 1; $k <= $nf; $k++) {
		$rep[$k] = 1;	# rep[k]=number of repetitions in column [k]
		}
		for ($i = 2; (!$bad_permut) and ($i <= $output_nlines); $i++) {	
		$pass_line = 0;
		for ($j = $i; (!$pass_line) and ($j <= $nlines); $j++) {	
			$nf = (@current = split(' ', $table[$j]));
			$fail = 0;
			for ($k = 1; (!$fail) and ($k <= $nf); $k++) {	

			if ($constraint[$k]>0) { 
				if (($prev[$k] eq $current[$k]) and 
				($rep[$k] + 1 > $constraint[$k])) {
				$fail = 1;
				}
			}

			if ($constraint[$k]<0) { 
				$back=$i-(-$constraint[$k]);
				if ($back<1) { $back=1; }; 
				for ($ll = $back; ($ll < $i) and (!$fail);$ll++) {
				@tmp=split(' ',$table[$ll]);
				if ($tmp[$k] eq $current[$k]) {
					$fail=1;
				}
				}
			}

			} # next k
			if ($fail == 0) {
			$pass_line = 1;
			}
		} # next j
		# now, if pass_line==1, line j fulfills the constraints
		if ($pass_line) {
			$j--;
			if ($j > $i) { 
			swap(\@table, $i, $j);
			}
			for ($k = 1; $k <= $nf; $k++) {
			if ($prev[$k] eq $current[$k]) {       
				$rep[$k]++;
			}
			else {
				$rep[$k] = 1;
			}
			$prev[$k] = $current[$k];
			}
		}
		else {
			# it is not possible to rearrange this permutation;
			# let's try another one.
			$bad_permut = 1;
		}
		}
		# next i
	} while ( $bad_permut & ($loop < $MaxLoops));  

	if ($bad_permut) { 0 } 
	else { 1 } 

}

sub swap {
    my $array = shift;
    my $i = shift;
    my $j = shift;
    @$array[$i,$j]=@$array[$j,$i];
}


sub permute {
# random permutation of array's elements
# only the "size" last elements are randomized

    my $array = shift;
    my $size = shift;
    my $i;
    my $j;
    for ($i = $size-1; $i > 1; $i--) {
		$j = 1 + int(rand(1) * $i);
		@$array[$i,$j] = @$array[$j,$i];
    }
#     print Dumper \$array;
}


sub filter  {
# this routine just checks if a permutation respects the constraints
# returning 1 if pass or 0 if failure
  my @table = shift;
  my $constr = shift;
  my $nlines = shift;


  my @constraint = split(' ', $constr);

  my @prev;
  my  $nf = (@prev = split(' ', $table[1]));

  my @current;
  my @tmp;
  my @rep;
  my $back;

    for (my $k = 1; $k <= $nf; $k++) {
	$rep[$k] = 1;	# rep[k]=number of repetitions in column [k]
    }

    for (my $i = 2; $i <= $nlines; $i++) {	
	    $nf = (@current = split(' ', $table[$i]));
	    my $fail = 0;
	    for (my $k = 1; (!$fail) and ($k <= $nf); $k++) {	

		if ($constraint[$k]>0) { 
		    if (($prev[$k] eq $current[$k]) and 
			($rep[$k] + 1 > $constraint[$k])) {
			$fail = 1;
			return 0;
		    }
		}

		if ($constraint[$k]<0) { 
		    $back=$i-(-$constraint[$k]);
		    if ($back<1) { $back=1; }; 
		    for (my $ll = $back; ($ll < $i) and (!$fail);$ll++) {
			@tmp=split(' ',$table[$ll]);
			if ($tmp[$k] eq $current[$k]) {
			    $fail=1; 
			    return 0;
			}
		    }
		}

	    } # next k

	    for (my $k = 1; $k <= $nf; $k++) {
		if ($prev[$k] eq $current[$k]) {       
		    $rep[$k]++;
		}
		else {
		    $rep[$k] = 1;
		}
		$prev[$k] = $current[$k];
	    }

    }
    # next i
return 1; 

}
# shuffle file of result
# input: file of result
# ouput: shuffle file of result but not the header 
sub shuffle_result{
	my ($in,$out) = @_;
	# a count to know if its a header
	my  $c = 0;
	my @corps;
	# Open the file of result
	open IN, "$in" or die "Cannot open $in\n";
# 	open OUT, ">", "$out" or die "Cannot write in $out\n";
	# Browse the input file
	while (<IN>){
		$c++;
		# the line 1 is write in the ouput file
		chomp;
		push (@corps, $_);
	}
	# shuffle line in the temp file and store in another file
	shuffle::new_shuffle(-seed => 10, -table => \@corps, -output => $out, -header => "true");

	close (IN);

} 
1;

