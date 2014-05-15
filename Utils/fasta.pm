package parse_fasta;

$VERSION = v0.0.1;


##########################################################################################
# parse_fasta.pm
#  
# Goal : Parsing of a fasta file
# 
#
#
#
#########################################################################################

################################# Library ###############################################

use Bio::SeqIO;
use strict;
use warnings;
use Data::Dumper;
# to use liste
use List::Util qw(shuffle);
use Pod::Usage;
use Bio::DB::Fasta;



$|=1;

################################# function ##############################################


# let to know the sizes of sequences fasta
# input : name of file
# ouput : number of sequences and hash  with ID in keys and length in values
sub lengthseq{
	my ($fafile) =@_;
	# allow to read a fasta file 
	my $in = Bio::SeqIO ->new(-file => $fafile);
	my %fasta=();
	my $cin=0;
	# read the file sequence by sequence 
	while (my $seq=$in -> next_seq){
		$cin++;
		# stocke sequence ID and the length of the sequence in a hash: has{ID} -> sequence
		$fasta{$seq->id} = length($seq->seq);
		print "$cin\n";
				
	}
	return ($cin,%fasta);
}



# let to know the sizes of sequences fasta
# input : name of fasta file [optional: a hash of reference and a verbosity]
# hash: if present =>hash{ID} = ""
# ouput : number of sequences and hash  with ID in keys and sequences in values
sub seq_ID{
	my ($fafile, $hash_ref, $verbose) =@_;
	$verbose ||= 0;
	# use of the library SeqIo from Bioperl to stocke a file in a variable
	my $in = Bio::SeqIO ->new(-file => $fafile);
	my %fasta=();
	my %genID; # this is to know if there is a duplicate of ID
	if ($hash_ref) {
		%genID = %{$hash_ref};
	}
	my $cin =0;
	# while there is a sequence stock in the variable
	while (my $seq=$in -> next_seq){
		$fasta{$seq->id}->{"sequence"} = $seq->seq;
		$fasta{$seq->id}-> {"size"} = length($seq->seq);
		$cin++;
		if (!exists $genID{$seq->id}){
			$genID{$seq->id}="";
		}
		else{
			die "The programme found 2 similar gene ID:", $seq->id,"\n"; 
		}
		if ($verbose >9){
			  print "sequence: ",$cin,"/",length $in ," with a size of ", length($seq->seq), " done\r";
		}
			
	}
	print "gene $cin\n";
	return ($cin,\%fasta,\%genID);
}

# put all sequence and identifiant of the sequence in a hash
# input: 
#	- fasta file or directory that contain genome
#	- a hash of reference where we take identifiant of sequence
# output:
#	- a hash that contain sequence and id sequence
sub fasta2hash{
	# variable
	my @args = @_;
	my %opts = @args;
	my $fastaDirFile;
	my $i=0;
	my $list_id;
	my %hash_coord;
	my %hash_chr;
	# Let to stop if we have two time the same identifiant
	my %hash_genID;
	my $cin;
	# test on option
	if (!$opts{"-fastaDirFile"}){
		die "you have not specify a name of directory or file in the optsion -fastaDirFile in the sub fasta2hash. \ncheck if you have well write -fastaDirFile\n";
	}
	elsif ($opts{"-fastaDirFile"}){
		$fastaDirFile =$opts{"-fastaDirFile"};
	}
	my $verbose = $opts{"-verbosity"};
	$verbose ||=0;
	if ($opts{"-hash_coord"}){
		my $href = $opts{"-hash_coord"};
		%hash_coord = %{$href};
	}
	if ($verbose >0){ print STDERR "Get Sequences\n";}

	# If fastaDirFile is a Directory
	if (-d $fastaDirFile){
		my %hash_chr_temp;
		my $cpt =0;
		foreach my $chr (keys %hash_coord){
			#Initalize sequence:
			my $file= $fastaDirFile."/".$chr.".fa";
			my ($cin, $hash_chr_ref,$hash_genID_ref) = parse_fasta::seq_ID($file,\%hash_genID, $verbose);		
			%hash_chr_temp = %{$hash_chr_ref};
			foreach my $key (keys %hash_chr_temp){
				$hash_chr{$key}= $hash_chr_temp{$key};
			}
			%hash_genID = %{$hash_genID_ref};
			if ($verbose >9){
			   Utils::showProgress2(scalar (keys %hash_coord), $i++, "Parse input file: ");
			}
			
		}
		return %hash_chr;
	}
	
	# if fastaDirfFile is a repertory
	 elsif (-r $fastaDirFile){
		## If multifasta file
	
		#Initalize sequence:
			my ($cin, $hash_chr_ref,$hash_genID_ref) = parse_fasta::seq_ID($fastaDirFile,\%hash_genID, $verbose);		
			%hash_chr = %{$hash_chr_ref};
			%hash_genID = %{$hash_genID_ref};
		return %hash_chr;
		}


	 else {
		pod2usage("Cannot understand defined optsion -f \"$fastaDirFile\"... ");
	}

}  


# get a random sequence
# input: 
#	- list_size: list of size
#	- h_fasta: a hash that contain sequence of the genome 
#	- h_coord: a hash result of the parsing of a gtf
# output:
#	- a fasta file that contain random sequence 
sub random_seq{
	# variable
	my (@args) = @_;
    my %opts = @args;
	my @list_size;
	my %hash_fasta;
	my %hash_coord;
	my $outfile;
	my $rand_coord;
	my $new_seq;
	my $newid;
	my $verbose;
	
	# tests
	if (! $opts{"-list_size"}){
		die "you forgot to put the reference to the list of size in the optsion -list_size\n";
	}
	elsif (!$opts{"-h_fasta"}){
		die "you forgot to put the reference to the hash of fasta in the optsion -h_fasta\n";
	}
	elsif (!$opts{"-h_coord"}){
		die "you forgot to put the reference to the hash of coord in the optsion h_coord\n";
	}
	elsif (!$opts{"-outfile"}){
		warn "You have not given a name of the output file then this name will be: random.fa\n";
		$outfile = "random".($#list_size+1).".fa";	
	}
	else{
		@list_size = @{$opts{'-list_size'}};
		%hash_fasta = %{$opts{'-h_fasta'}};
		if (!keys %hash_fasta){
			die "the hash from the reference put in the option -h_fasta is empty\n";
		}
		%hash_coord = %{$opts{'-h_coord'}};
		if (!keys %hash_coord){
			die "the hash from the reference put in the option -h_coord is empty\n";
		}
		$outfile = $opts{'-outfile'};
	}
	if ($opts{'-verbosity'}){
		$verbose = $opts{'-verbosity'};
	}
	else{$verbose =0}
	print $verbose;
	# main
	my $seqOUT =Bio::SeqIO ->new(-format => 'fasta', -file => ">".$outfile);
	my @list_chr =grep /^[\d+XY]/, (keys %hash_fasta);
# 	foreach my $chr (keys %hash_fasta){
# 		if ($chr =~ /^[\d+XY]/){		
# 			push (@list_chr, $chr);
# 		}
# 	}
	my $i=0;
	foreach my $size (@list_size){
		my $chromosome = @list_chr[rand @list_chr];
		my $seq = "";
		my $N = "N";
		my $overlap=1;
		my $href;

		while ($overlap ==1 or $N eq "N"){
			$rand_coord = int(rand($hash_fasta{$chromosome}->{"size"} - int($size)));

			$overlap = Utils::overlapHash_chr(-hash_coord => $hash_coord{$chromosome} , -coord_begin => $rand_coord, -coord_end => int($rand_coord+$size));
			if ($overlap == 0){
				($seq, $N) = Utils::CheckNinSeq(-seq => $hash_fasta{$chromosome}->{"sequence"}, -begin => $rand_coord, -end => ($size-1));
			}

		}		
		
		$newid = "intergenic_".$chromosome."_".$rand_coord."-".($rand_coord+$size);
		$new_seq = Bio::Seq->new(-id => $newid, -seq => $seq);
		$seqOUT->write_seq($new_seq);
 		if ($verbose>9){
			Utils::showProgress2(scalar @list_size, $i++, "Sequence done: ");
		}
	}
}



# same of random_seq but try to use Bio::SB::Fasta
sub random_seq_DB{
	# variable
	my (@args) = @_;
    my %opts = @args;
	my @list_size;
	my $db;
	my %hash_coord;
	my $outfile;
	my $rand_coord;
	my $new_seq;
	my $newid;
	my $verbose;
	
	# tests
	if (! $opts{"-list_size"}){
		die "you forgot to put the reference to the list of size in the optsion -list_size\n";
	}
	elsif (!$opts{"-Index"}){
		die "you forgot to put the reference to the hash of fasta in the optsion -h_fasta\n";
	}
	elsif (!$opts{"-h_coord"}){
		die "you forgot to put the reference to the hash of coord in the optsion h_coord\n";
	}
	elsif (!$opts{"-outfile"}){
		warn "You have not given a name of the output file then this name will be: random.fa\n";
		$outfile = "random".($#list_size+1).".fa";	
	}
	else{
		@list_size = @{$opts{'-list_size'}};
		$db = $opts{'-Index'};
		if (!$db){
			die "The index is not present\n";
		}
		%hash_coord = %{$opts{'-h_coord'}};
		if (!keys %hash_coord){
			die "the hash from the reference put in the option -h_coord is empty\n";
		}
		$outfile = $opts{'-outfile'};
	}
	if ($opts{'-verbosity'}){
		$verbose = $opts{'-verbosity'};
	}
	else{$verbose =0}
	print $verbose;
	# main
	my $seqOUT =Bio::SeqIO ->new(-format => 'fasta', -file => ">".$outfile);
	my @list_chr =grep /^[\d+XY]/, $db->ids;

	my $i=0;
	foreach my $size (@list_size){
		my $chromosome = @list_chr[rand @list_chr];
		my $seq = "";
		my $N = "N";
		my $overlap=1;
		my $href;

		while ($overlap ==1 or $N eq "N"){
			$rand_coord = int(rand($db->length($chromosome) - int($size)));

			$overlap = Utils::overlapHash_chr(-hash_coord => $hash_coord{$chromosome} , -coord_begin => $rand_coord, -coord_end => int($rand_coord+$size));
			if ($overlap == 0){
				($seq, $N) = Utils::CheckNinSeq(-seq => $db->seq($chromosome, $rand_coord => int($rand_coord+$size)));
			}

		}		
		
		$newid = "intergenic_".$chromosome."_".$rand_coord."-".($rand_coord+$size);
		$new_seq = Bio::Seq->new(-id => $newid, -seq => $seq);
		$seqOUT->write_seq($new_seq);
 		if ($verbose>9){
			Utils::showProgress2(scalar @list_size, $i++, "Sequence done: ");
		}
	}
}


# group all sequence of a hash in one sequence
# input: 
#	- a hash : keys => ID, value => hash2
#		+ hash2: keys => "sequence" , value: string of the sequence
#		+ hash2: keys => "size", 	value: size of the sequence 
# output:
#	a sting that regroup all sequence 
sub group_seq_h {
	my $href =shift;
	# initilazation of $sequence like a void string.	
	my %fasta = %{$href};
	my $sequence = "";
	# all sequence contain in value of %fasta is put in $sequence 
	foreach my $key (keys(%fasta)){
		$sequence =$sequence.$fasta{$key}->{"sequence"};
	}
	print "all sequence group\n";

	return $sequence;
}

# write random sequence in a fasta file
# input: 
#	- $seq: string of a sequence
#	- $cin: number of sequence
#	- $ref_list: reference of a list that contain a list of size
#	- $i number of the file
# output:
#	- a file of random sequence
sub write_seqfa_rand {
	my ($seq,$cin, $ref_list, $i) =@_;
	my @list_size = @{$ref_list};
	my $newid;
	# initialization of a file in fasta format
	my $seqOUT =Bio::SeqIO ->new(-format => 'fasta', -file => ">shuffle".$i.".fa");
	for (my $x = 1; $x <=$cin; $x++){
		# create the identifiant of the sequence
		$newid	= $x."_shuffle_".$i;
		# extract a random number between 0 and the length of the sequence - the random size
		my $rand = int(rand($cin-int($list_size[$x-1])));
		my $len = length($seq);
		# extract the list between the coord give by $rand and the size random
		my @chars =split '', substr($seq, $rand, $rand+$list_size[$x-1]);
		@chars = shuffle @chars;
		# transform the list in string
		my $seq_shuffle = join '',@chars;
# 		print "sequence: ",$seq_shuffle,"\n";
		# create a new sequence with $newid in identifiant and $seq_shuffle in sequence
		my $new_seq = Bio::Seq->new(-id => $newid, -seq => $seq_shuffle);
		# write sequence in the fasta file
		$seqOUT->write_seq($new_seq);
   	    Utils::showProgress2($cin, $x++, "write: ");
	}

}

# shuffle sequence of a fasta file
# Input:
#	- $seqfileA: 		file of sequence coding 
#	- $seqfileB: 		file of sequence non coding
#	- $seqfileAprime: 	file of sequence coding of the interest
#	- $number: 		number of shuffle file
# Output:
# 	- a file with the same number of sequence than $seqfile but with a shuffle sequence
sub fasta_shuffle{
	my @args = @_;
	my %opts = @args;
	my $seqfileA = $opts{"-seqA"};
	my $seqfileB = $opts{"-seqB"};
	my $seqfileAprime = $opts{"-seqAprime"};
	my $number = $opts{"-number"};
	my $cut = $opts{"-cut"};
	my $verbosity = $opts{"-verbosity"};
	my @mean;
	my @sd;
	my $cin;
	my $sequence;
	print "create of $number file of shuffle sequence fasta \n";
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
	for (my $i =1; $i<=$number;$i++){
		print "write the file $i\n";
		# use R to extract some value who follow a normal law
		# input: 
		#	- a mean ($mean)
		#	- a standard deviation ($sd)
		#	- the number of sequence contain by $seqfile
		# output:
		#	- a string 
 		print "R\n";
 		print "$cut\n";
		my $commandR = `PredictSize.r $mean $sd $cin $cut`;
 		# parse he result of R in a list of size
		my @list_size = R_parse::list($commandR);
		print "R done\n";
		# shuffle the list that contain the sequence
		print "shuffle\n"; 
		parse_fasta::write_seqfa_rand($sequence, $cin, \@list_size, $i);
		print "file $i done\n";
	}
 }



1;
