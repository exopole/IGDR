package Utils;

$VERSION = v0.0.1;

use strict;
use warnings;
use Data::Dumper;
use Bio::DB::Fasta;

$| = 1;

# get the minimum of two elements
sub min2 {
  my ($a, $b) = @_;

  if (! defined $a) {
    return $b;
  }
  if (! defined $b) {
    return $a;
  }
  return (($a < $b)? $a: $b);
}

# get the minimum of three elements
sub min3 {
  my ($a, $b, $c) = @_;

  return min2(min2($a, $b), $c);
}

# get the maximum of two elements
sub max2 {
  my ($a, $b) = @_;

  if (! defined $a) {
    return $b;
  }
  if (! defined $b) {
    return $a;
  }
  return (($a > $b)? $a: $b);
}

# get the maximum of three elements
sub max3 {
  my ($a, $b, $c) = @_;

  return max2(max2($a, $b), $c);
}

# get the average of elements
sub average {
  my (@elements) = @_;
  my $sum = 0;
  my $nb = 0;
  foreach my $element (@elements) {
    if (defined $element) {
      $sum += $element;
      $nb++;
    }
  }
  if ($nb == 0) {
    return 0;
  }

  return ($sum / $nb);
}

# write to a file
sub writeNewFile {
  my ($fileName, $content) = @_;

  if (-e $fileName) {
    unlink $fileName;
  }
  open (TMPFILE, "+> $fileName") or die $!;
  print TMPFILE $content;
  close TMPFILE;
}

# show a progress bar
sub showProgress {
  my ($aim, $current, $message) = @_;

  if (! defined $message) {
    $message = "";
  }
  my $percent = int ((($current+1) / ($aim+1)) * 100);
  print STDERR $message . " " x (30 - (length $message)) . "[" . "-" x $percent . " " x (100 - $percent) . "]\r";
}

sub showProgress2 {
  my ($aim, $current, $message) = @_;

  if (! defined $message) {
    $message = "";
  }
  my $percent = int ((($current+1) / ($aim+1)) * 100);
  print STDERR $message . " ".$percent." %\r";
}

# compute (G+C)%
sub getGC {
  my ($sequence) = @_;
  my $size = 0;
  my $gc = 0;

  for (my $i = 0; $i < length $sequence; $i++) {
    if ((substr $sequence, $i, 1) =~ /[GCgc]/) {
      $gc++;
    }
    $size++;
  }
  return (($gc / $size) * 100);
}

# remove leading or trailing whitespace characters
sub trim {
  my ($string) = @_;

  $string =~ s/^\s+//;
  $string =~ s/\s+$//;

  return $string;
}

# Sub routines for computing overlap between 2 ranges
# return 1 if overlap, 0 otherwise
sub foverlap
{
    my ($beg1,$end1,$beg2,$end2) = @_;
    die "beg1 :$beg1\n" if (!$beg1);
    die "beg2 :$beg2\n" if (!$beg2);
	die "end1 :$end1\n" if (!$end1);    
	die "end2 :$end2\n" if (!$end2);
    return (($end1>=$beg2)&&($beg1<=$end2)) ? 1 : 0;
}


# Sub routines for computing overlap between 2 ranges
# return size of the overlap if overlap, 0 otherwise
# Note : to be tested extensively
sub foverlapsizerange
{
    my ($b1,$e1,$b2,$e2) = @_;
	
	if ( foverlap ($b1,$e1,$b2,$e2) ){
	
		if ($b1<=$b2 && $e1<=$e2) {return ($e1-$b2+1)}
		if ($b1>=$b2 && $e1>=$e2) {return ($e2-$b1+1)}
		if ($b2>=$b1 && $e1>=$e2) {return ($e2-$b2+1)}
		if ($b2<=$b1 && $e1<=$e2) {return ($e1-$b1+1)}						
	} else { return 0}
}

# sub routine for computing overlap between a list of coord in a hash and a range
# return 0 if do not overlap and 1 if overlap
sub overlapHash_chr{

	#my @args = @_;
	my %opts= @_;
	my %hash_coord = %{$opts{"-hash_coord"}};
	my $coord_begin = $opts{"-coord_begin"};
	my $coord_end = $opts{"-coord_end"};
	my $overlap =0;
	my $i=0;
	foreach my $tr (
		sort {$hash_coord{$a}->{'startt'} <=> $hash_coord{$b}->{'startt'} }
		keys %hash_coord
        ){
        last if ($coord_begin > $hash_coord{$tr} -> {'endt'});
        if ($hash_coord{$tr}->{'biotype'} eq "protein_coding"){
        	next if ($coord_begin >  $hash_coord{$tr} -> {'endt'});
        	$overlap = Utils::foverlap($hash_coord{$tr} -> {'startt'},$hash_coord{$tr} -> {'endt'},$coord_begin,$coord_end);
			last if ($overlap);	
		}
		
	}
	return $overlap;

}

# build a file that contain a table
# input : a hash, the name of the ouptut file, and names of columns
# output: a table with 2 columns
sub table2c{
	my ($href,$out,$col1n,$col2) = @_;
	my %h	= %{$href};
	open(OUT, ">$out");
	print OUT "$col1n\t$col2\n";
	foreach my $k (keys(%h)){
		print OUT "$k\t$h{$k}\n";
	}
}

# to have the standard deviation
sub stdev{
        my($data) = @_;
        if(@$data == 1){
                return 0;
        }
        my $average = &average($data);
        my $sqtotal = 0;
        foreach(@$data) {
                $sqtotal += ($average-$_) ** 2;
        }
        my $std = ($sqtotal / (@$data-1)) ** 0.5;
        return $std;
}
# calculate the mean from a hash{identifier} -> {"size"}-> {size}
sub mean_h {
	my $href = shift;
	my %h = %{$href};
	my $cpt = 0;
	my $len = 0;
	foreach my $seq (keys(%h)){
		$cpt++;
		$len += $h{$seq}->{"size"}; 
	}
	return $len/$cpt;
}

sub mean_DB{
	my $db = shift;
	my @ids = $db->ids;
	my $count =@ids;
	my $total =0;
	foreach my $id (@ids){
		$total += $db->length($id);
	}
	return $total/$count;
}

sub mean_max_DB{
	my $db = shift;
	my @ids = $db->ids;
	my $count =@ids;
	my $total =0;
	my $max =0;
	foreach my $id (@ids){
		$total += $db->length($id);
		$max = length($id) if (length($id)>$max);
	}
	my $mean = $total/$count;
	return ($mean, $max);
}

sub sd_h {
	my ($href, $mean)= @_;
	my %h = %{$href};
	my $cpt = 0;
	my $sdb = 0;
	foreach my $key (keys(%h)){
		$cpt++;
		$sdb += ($h{$key}->{"size"}-$mean) * ($h{$key}->{"size"}-$mean); 
	}
	return sqrt($sdb/$cpt);
}

sub sd_DB {
	my ($db, $mean)= @_;
	my @ids = $db->ids;
	my $count =@ids;
	my $sdb =0;
	foreach my $id (@ids){
		$sdb += ($db->length($id) - $mean) * ($db->length($id) - $mean);
	}
	return sqrt($sdb/$count);
}

sub division2xd1n{
	my ($den1, $den2, $num) =@_;
	return ($den1*$den2)/$num;
}

# take random size of sequence that follow a normal law
sub predict_size{
	my @args = @_;
	my %opt = @args;
	my $seqfileA;
	my $seqfileB;
	my $seqfileAprime;
	my $cut; 
	if (! $opt{"-seqfileA"}){
		die "you forgot to put a name file in -seqfileA\n";
	}
	elsif (! $opt{"-seqfileB"}){
		die "you forgot to put a name file in -seqfileB\n";
	}
	elsif (!$opt{"-seqfileAprime"}){
		die "you forgot to put a name file in -seqfileAprime\n";
	}
	elsif (!$opt{"-cuttof"}){
		die "you forgot to put an integer in the cuttof option\n";
	}
	else{
		$seqfileA = $opt{"-seqfileA"};
		$seqfileB = $opt{"-seqfileB"};
		$seqfileAprime = $opt{"-seqfileAprime"};
		$cut = $opt{"-cuttof"}; 
	}
	my @mean;
	my @sd;
	my $cin;
	my $sequence;
	# stock all name of file in a list
	my @seqfile = ($seqfileA, $seqfileB, $seqfileAprime);
	# use of a library to parse $seqfile. 
	# input : $seqfile
	# output : 
	#	- the number of sequence
	# 	- A has that contain identifiant of sequence in keys and sequence in values. 
	for (my $i =0;$i<$#seqfile+1; $i++){
		print "Parsing of ".$seqfile[$i]."\n";
		my ($cin_temp, $fasta_ref) = parse_fasta::seq_ID($seqfile[$i]);
		my %fasta = %{$fasta_ref};
		if ($i ==2){
			$sequence = parse_fasta::group_seq_h(\%fasta);
			$cin = $cin_temp;
		}
		push (@mean, Utils::mean_h(\%fasta));
		push (@sd, Utils::sd_h(\%fasta, $mean[$i]));
	}
	my $mean = Utils::division2xd1n($mean[2], $mean[1], $mean[0]);
	my $sd = Utils::division2xd1n($sd[2], $sd[1], $sd[0]);
	my $commandR = `PredictSize.r $mean $sd $cin $cut`;
	# parse he result of R in a list of size
	my @list_size = R_parse::list($commandR);
	return @list_size;

}

# take random size of sequence that follow a normal law
sub predict_size_DB{
	my @args = @_;
	my %opt = @args;
	my $seqfileA;
	my $seqfileB;
	my $seqfileAprime;
	my $pmin;
	my $pmax; 
	if (! $opt{"-seqfileA"}){
		die "you forgot to put a name file in -seqfileA\n";
	}
	elsif (! $opt{"-seqfileB"}){
		die "you forgot to put a name file in -seqfileB\n";
	}
	elsif (!$opt{"-seqfileAprime"}){
		die "you forgot to put a name file in -seqfileAprime\n";
	}
	elsif (!$opt{"-min"}){
		die "you forgot to put an integer in the option pmin?\n";
	}
	else{
		$seqfileA = $opt{"-seqfileA"};
		$seqfileB = $opt{"-seqfileB"};
		$seqfileAprime = $opt{"-seqfileAprime"};
		$pmin = $opt{"-min"};
		$pmax = $opt{"-max"}; 
	}
	my $cin;
	my $sequence;
	my %mean;
	my %sd;
	my @seqfile = ($seqfileA, $seqfileB, $seqfileAprime);
	my $max;
	my @result;
	# use of a library to parse $seqfile. 
	# input : $seqfile
	# output : 
	#	- the number of sequence
	# 	- A has that contain identifiant of sequence in keys and sequence in values. 

	foreach my $file (@seqfile){
		print "Parsing of ".$file."\n";
	    my $db      = Bio::DB::Fasta->new($file);
		if ($file eq $seqfileB){
			$cin = $db->ids;
			@result = mean_max_DB($db);
			$mean{$file} =$result[0];
			$pmax = $result[1] if ($pmax ==0);
			$sd{$file} = sd_DB($db, $mean{$file});
		}
		else{
			$mean{$file} = mean_DB($db);
			$sd{$file} = sd_DB($db, $mean{$file});
		}
	}
	my $mean = Utils::division2xd1n($mean{$seqfileA}, $mean{$seqfileB}, $mean{$seqfileAprime});
	my $sd = Utils::division2xd1n($sd{$seqfileA}, $sd{$seqfileB}, $sd{$seqfileAprime});
	my $commandR = `PredictSize.r $mean $sd $cin $pmin $pmax`;
	# parse he result of R in a list of size
	my @list_size = R_parse::list($commandR);
	return @list_size;

}

# check if the string sequence contain a N
sub CheckNinSeq {
	my @args = @_;
	my %opts= @args;
	my $subseq;
	if ($opts{"-begin"} && $opts{"-end"}){
		$subseq = substr $opts{"-seq"}, int($opts{"-begin"}), int($opts{"-end"});
	}
	else{
		$subseq = $opts{"-seq"};
	}
	if ($subseq =~ /N/){
		return $subseq, "N";
	}	
	else{
		return $subseq, "";
	}

}


1;
