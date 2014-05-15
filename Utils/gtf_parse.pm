package gtf_parse;

$VERSION = v0.0.1;


##########################################################################################
# gtf_parse.pm
#  
# Goal : Parsing of a GTF file
# 
# 
#
#########################################################################################

################################# Library ###############################################
use strict;
use warnings;
use Data::Dumper;
use File::Basename;


$|=1;




############################# Function ###################################################
sub parseLevelGtfhashKey_transcript{

	my ($infile, $levels, $list, $st, $verbosity) = @_;
	$infile ||="";
	$list ||="";	
 	my $baseName = basename($infile);
	my @all_level = split (",", $levels);
	
  	$verbosity = 0 unless (defined $verbosity);

	open GTFFILE, "$infile" or die "Error! Cannot open GTF File ". $infile . ": ".$!;
	my @lines	=	<GTFFILE>;
	
	# Store data
	my %h_transcript;
	
	# Store data temp
	my %h_transcript_temp;
	
	# counter line
	my $i = 0;
	
	# hash that contain gene id in key and transcript id in values
	my %geneid;
	my $biotype;
	my $status;
	foreach my $line  (@lines){
		chomp $line;
		if ($line =~ /^(\S+)\s+(\S+)\s+(\S+)\s+(\d+)\s+(\d+)\s+(\S+)\s+(\S+)\s+(\S+)\s+gene_id "(\S+)"; transcript_id "(\S+)"; (\w.*)$/ || 	
		$line =~ /^(\S+)\s+(\S+)\s+(\S+)\s+(\d+)\s+(\d+)\s+(\S+)\s+(\S+)\s+(\S+)\s+gene_id "(\S+)"; transcript_id "(\S+)";$/){
			my $chr     	=   $1;
			my $biotype		=   $2; 
			my $level		=   $3; 
			my $beg_feat    =   $4; 
			my $end_feat    =   $5;  
			my $score	    =   $6; 
			my $strand      =   $7; 
			my $frame		=   $8;
			my $gene_id		=   $9;
			my $transcript	=   $10;
			my $extrafield	=	$11;
# 			print "level: $level\n";
			for my $lvl (@all_level) {
# 				print "lvl: $lvl\n";
				if ($lvl eq $level){
# 					print "entrer dans la boucle \n";
			# 		Transcript levels
					$chr =~ s/chr//;
					$h_transcript_temp{$transcript}->{"chr"}			=   $chr;
					$h_transcript_temp{$transcript}->{"biotype"}		=   $biotype;
					$h_transcript_temp{$transcript}->{"startt"}		=   Utils::min2($h_transcript_temp{$transcript}->{"startt"}, $beg_feat);
					$h_transcript_temp{$transcript}->{"endt"}		=   Utils::max2($h_transcript_temp{$transcript}->{"endt"}, $end_feat);
					$h_transcript_temp{$transcript}->{"score"}		=   $score;
					$h_transcript_temp{$transcript}->{"strand"}      =   $strand;
					$h_transcript_temp{$transcript}->{"frame"}		=   $frame;
					$h_transcript_temp{$transcript}->{"gene_id"}		=   $gene_id;

					# parse gtf attributes after transcript_id
					my %attribs;
					if ($extrafield){
						my @extrafields = split(";", $extrafield);
						# store ids and additional information in second hash
						foreach my $attr ( @extrafields ) {
							next unless $attr =~ /^\s*(.+)\s(.+)$/;
							my $c_type  = $1;
							my $c_value = $2;
							$c_value=~ s/\"//g;
							if(exists($attribs{$c_type}) and $c_type ne 'tag' and $c_type ne 'ont'){
								warn "$c_type key already exists with",$attribs{$c_type}," value...\n";	
								print "$line\n";					
							
							}
							
							# NOTE : CORRECT nomenclature is:
							# - transcript_biotype (not transcript_type as in Ensembl)
                            # - gene_biotype (not gene_type as in Ensembl)
                            $c_type =~ s/transcript_type/transcript_biotype/g;
                            $c_type =~ s/gene_type/gene_biotype/g;                            
							$attribs{$c_type} = $c_value if($c_type  && $c_value);
						}
					}
					# Feature infos are stored in N hashes from a array reference
					# $h{$tr}->{"feature"} is the reference on a array
					# @{$h{$tr}->{"feature"}} to dereference it
					# To obtain a particular element in the array : ${ $h{$tr}->{"feature"} }[0]
					my %feature = ( "start" => $beg_feat, "end" => $end_feat, "feat_level" => $level);
					
					my %merge = (%attribs, %feature);
					push (@{$h_transcript_temp{$transcript}->{"feature"}}, \%merge);
					$biotype = 
			${$h_transcript_temp{$transcript}->{"feature"}}[0]->{"transcript_biotype"};
					$status = 
				${$h_transcript_temp{$transcript}->{"feature"}}[0]->{"transcript_status"};
					if ($list eq "" or $list =~ /$biotype/ and $status eq $st){

						if (exists($geneid{$gene_id})){
							my $size	=
							cumulSize($h_transcript_temp{$transcript}->{"feature"});
							my $size2 	=  
							cumulSize($h_transcript{$geneid{$gene_id}}->{"feature"});
						# 	print "size: $size\n";
	# 						print "size2: $size2\n";
							if ($size > $size2){
								delete $h_transcript{$geneid{$gene_id}};
								delete $geneid{$gene_id};
								# stocke of gene_id
								$geneid{$gene_id}= $transcript;
								$h_transcript{$transcript} = 
									$h_transcript_temp{$transcript} ;
								@{$h_transcript{$transcript}->{"feature"}} = sort 
									{ $a->{"start"} <=> $b->{"start"}}
									@{$h_transcript{$transcript}->{"feature"}};
							}
						}
						else{
							$h_transcript{$transcript} = $h_transcript_temp{$transcript} ;
							$geneid{$gene_id}= $transcript;
							@{$h_transcript{$transcript}->{"feature"}} = sort 
								{ $a->{"start"} <=> $b->{"start"}} 
								@{$h_transcript{$transcript}->{"feature"}};
						}				
					}
				}
			}

		
		}	 
		elsif (($line =~ /^#/) || ($line =~ /^\s*$/) || $line =~ /^track/)  {
     	 # skip
   		} else {
			die "\nError! Cannot parse this line in '$infile':\n", $line,"\n! Aborting.\nCheck the format of your input file\n";
		}

		
# 		verbose
		if ($verbosity > 0){
    	    Utils::showProgress(scalar @lines, $i++, "Parse input file: ");
    	}
    	

	}
	
	# close handle
	close GTFFILE;
	
	# Test parsing i.e empty hash
	die " Error! Invalid parsing step (empty hash):\n- Check that your file is in .GTF format.\n- Or check that level '$levels' are present in field 3 of '$infile'... " unless (keys(%h_transcript)>0);



	# return hash
	return %h_transcript;
}

# parse a gtf file to have coord of gene in the genome
# input: 
#	- $infile: gtf file
#	- $levels: expression level (exon intron) 
#	- $list: list of biotype don't accept in the parsing
# output:
#	- a hash:
#		+ keys: chromosome identifiant
#		+ value: list of tuple that contain gene ID, begin and end
sub parseLevelGtfhashKey_chr_genecoord{

	my ($infile, $levels,$list , $verbosity) = @_;
	$infile ||="";
	$list ||="";
 	my $baseName = basename($infile);
	my @all_level = split (",", $levels);
	
  	$verbosity = 0 unless (defined $verbosity);

	open GTFFILE, "$infile" or die "Error! Cannot open GTF File ". $infile . ": ".$!;
	my @lines	=	<GTFFILE>;
	
	# Store data
	my %h_coord;
	
	# initialization of coord list
	my @list_coord;
	
	# counter line
	my $i = 0;
	
	# hash that contain gene id in key and transcript id in values
	my %geneid;
	
	my $biotype;
	my $status;
	foreach my $line  (@lines){
		chomp $line;
		if ($line =~ /^(\S+)\s+(\S+)\s+(\S+)\s+(\d+)\s+(\d+)\s+(\S+)\s+(\S+)\s+(\S+)\s+gene_id "(\S+)"; transcript_id "(\S+)"; (\w.*)$/ || 	
		$line =~ /^(\S+)\s+(\S+)\s+(\S+)\s+(\d+)\s+(\d+)\s+(\S+)\s+(\S+)\s+(\S+)\s+gene_id "(\S+)"; transcript_id "(\S+)";$/ || 
	$line =~ /^(\S+)\s+(\S+)\s+(\S+)\s+(\d+)\s+(\d+)\s+(\S+)\s+(\S+)\s+(\S+)\s+gene_id "(\S+)";(\w.*)$/){
			my $chr     	=   $1;
			my $level		=   $3; 
			my $beg_feat    =   $4; 
			my $end_feat    =   $5;  
			my $strand      =   $7;
			my $gene_id		=   $9; 
			for my $lvl (@all_level) {
				if ($list eq "" or $list =~ /$biotype/ and $lvl eq $level){
			# 		Transcript levels
					$chr =~ s/chr//;
					# parse gtf attributes after transcript_id
					@list_coord = [$gene_id, Utils::min2($end_feat, $beg_feat), Utils::max2($beg_feat, $end_feat) ];
					push (@{$h_coord{$chr}}, @list_coord);
				}
			}
			foreach my $chromosome (keys %h_coord){
				@{$h_coord{$chromosome}} = sort { $a->[1] <=> $b->[1]} @{$h_coord{$chromosome}};
			}
		}	 
		elsif (($line =~ /^#/) || ($line =~ /^\s*$/) || $line =~ /^track/)  {
     	 # skip
   		} 
   		
   		else {
			die "\nError! Cannot parse this line in '$infile':\n", $line,"\n! Aborting.\nCheck the format of your input file\n";
		}

		
# 		verbose
		if ($verbosity > 9){
    	    Utils::showProgress2(scalar @lines, $i++, "Parse input file: ");
    	}
    	

	}
	
	# close handle
	close GTFFILE;
	
	# Test parsing i.e empty hash
	die " Error! Invalid parsing step (empty hash):\n- Check that your file is in .GTF format.\n- Or check that level '$levels' are present in field 3 of '$infile'... " unless (keys(%h_coord)>0);



	# return hash
	return %h_coord;
}

sub getCumulSizeFromGtfHash{

	# parsing a hash in sub need dereference in shift
	my %h 			=	%{shift()};
	my $longest 	=   shift;
	my $verbosity 	=   shift;
	$longest		= 	0 unless defined ($longest);
	$verbosity		= 	0 unless defined ($verbosity);
	
# 	print Dumper %h;
	# print if verbosity	
	print STDERR "- Extract size from GtfHash ...\n" if ($verbosity > 0);	
	
	# Mode : only compute size of feature 
	if ($longest == 0) {
	
		# Get size of the array ref and assign to transcript
		for my $tr (keys %h){
		
			#compute size
			my $size	=	cumulSize($h{$tr}->{"feature"});
			$h{$tr}->{"size"}	= $size;
		}
		return %h;


	} else { # MODE: get longest transcript per locus

		# print if verbosity	
		print STDERR "- Extract Longest feature per gene_id from GtfHash ...\n" if ($verbosity > 0);	
		
		# cp the hash into one that will only contain longest tx per locus 
		my %h_filtered = %h;
		
		# Hash for storing $tx/isoforms of equal sizes that should not be seen more than once 
		# if not used, 2 txs of equal size will be removed themselves (line
		my %notpassagain;
		
		# counter for progression
		my $i=0;	
		
		# 1ST: Get size of the array ref and assign to transcript
		#-----
		for my $tr (keys %h){
		
			if ($verbosity> 10) {
				Utils::showProgress(scalar(keys(%h)), $i++, "Get longest Tx per gene: ");
			}
			
			# For case of n isoforms of a same gene of equal size
			next if ($notpassagain{$tr});
			
			#compute size
			my $size	=	cumulSize($h{$tr}->{"feature"});
			$h{$tr}->{"size"}	= $size;

			# counter for same size isoform
			my $cpt =0;	

			# 2ND : Loop on the cpied hash
			#-----
			for my $tr2 (keys %h_filtered){				
								
				# next if not same locus since we want biggest Tx *per* locus
				next if ($h{$tr}->{"gene_id"} ne $h_filtered{$tr2}->{"gene_id"});
				
				# next if we are on the same tx!
				next if ($tr eq $tr2);				
				

				#compute size
				my $size2	=	cumulSize($h_filtered{$tr2}->{"feature"});
				$h_filtered{$tr2}->{"size"}	= $size2;
				
				# if the 2 transcripts $tr and $tr2 have the same gene_id (isoforms of the same locus)
				# if size of $tr > size $tr2
				if ($h{$tr}->{"size"} >= $h_filtered{$tr2}->{"size"} ){
 	 				# print "$tr $size >= $tr2 $size2 -- delete $tr2\n";
					delete $h_filtered{$tr2};
					
					# Store $tr2 in hash in order that $tr when encounter in $h loop will not remove current $tr
					$notpassagain{$tr2}++;
				}
			}
		} 
		return %h_filtered;		
	}
}


# parse a GTF in a hash
# example of a transcript: 
#				{ '32' => {
#                     'ENSCAFT00000015937' => {
#                                               'gene_id' => 'ENSCAFG00000010016',
#                                               'chr' => '32',
#                                               'frame' => '.',
#                                               'score' => '.',
#                                               'endt' => '16808406',
#                                               'biotype' => 'protein_coding',
#                                               'feature' => [
#                                                              {
#                                                                'transcript_type' => 'protein_coding',
#                                                                'feat_level' => 'exon',
#                                                                'gene_type' => 'protein_coding',
#                                                                'start' => '16807336',
#                                                                'end' => '16808406'
#                                                              }
#                                                            ],
#                                               'startt' => '16807336',
#                                               'strand' => '+'
#                                             },
#					...}

sub parseGTFSplitchr73{

	my ($infile, $levelfilter, $verbosity) = @_;
	$infile ||=	"";	
	$levelfilter ||=	"exon";	
	my $baseName = basename($infile);
	  $verbosity ||= 0;

	open GTFFILE, "$infile" or die "Error! Cannot open GTF File ". $infile . ": ".$!;
	my @lines	=	<GTFFILE>;
	# Store data
	my %h_transcript;
	# Store data by chr
	my %h_chr = ();
	# counter line
	my $i = 0;

	foreach my $line  (@lines){
		chomp $line;
		if ($line =~ /^(\S+)\s+(\S+)\s+(\S+)\s+(\d+)\s+(\d+)\s+(\S+)\s+(\S+)\s+(\S+)\s+gene_id "(\S+)"; transcript_id "(\S+)"; (\w.*)$/ ||
			$line =~ /^(\S+)\s+(\S+)\s+(\S+)\s+(\d+)\s+(\d+)\s+(\S+)\s+(\S+)\s+(\S+)\s+gene_id "(\S+)"; transcript_id "(\S+)";$/){
			my $chr     =   $1; 
			my $biotype	=   $2; 
			my $level	=   $3; 
			my $beg_feat    =   $4; 
			my $end_feat    =   $5;  
			my $score	   =   $6; 
			my $strand      =   $7; 
			my $frame	=   $8;
			my $gene_id	=   $9;
			my $transcript	=   $10;
			my $extrafield	=	$11;

			if ($level eq $levelfilter) {
				# Transcript levels
				$h_transcript{$transcript}->{"chr"}	=   $chr;
				$h_transcript{$transcript}->{"biotype"}	=   $biotype;
				$h_transcript{$transcript}->{"startt"}	=   Utils::min2($h_transcript{$transcript}->{"startt"}, $beg_feat);
				$h_transcript{$transcript}->{"endt"}	=   Utils::max2($h_transcript{$transcript}->{"endt"}, $end_feat);
				$h_transcript{$transcript}->{"score"}	=   $score;
				$h_transcript{$transcript}->{"strand"}      =   $strand;
				$h_transcript{$transcript}->{"frame"}	=   $frame;
				$h_transcript{$transcript}->{"gene_id"}	=   $gene_id;

				# parse gtf attributes after transcript_id
				my $ref_attrib = parseExtrafield($extrafield, $verbosity,'transcript_type,gene_type');
				# Feature infos are stored in N hashes from a array reference
				# $h{$tr}->{"feature"} is the reference on a array
				# @{$h{$tr}->{"feature"}} to dereference it
				# To obtain a particular element in the array : ${ $h{$tr}->{"feature"} }[0]
				my %feature = ( "start" => $beg_feat, "end" => $end_feat, "feat_level" => $level);
				my %merge = (%{$ref_attrib}, %feature);
				push (@{$h_transcript{$transcript}->{"feature"}}, \%merge);	
			}	
		} elsif (($line =~ /^#/) || ($line =~ /^\s*$/) || $line =~ /^track/)  {
				# skip
		} else {
			die "\nError! Cannot parse this line in '$infile':\n", $line,"\n! Aborting.\nCheck the format of your input file\n";
		}
		# verbose
		if ($verbosity > 0){
			   Utils::showProgress2(scalar @lines, $i++, "Parse input file: ");
			}


		}
	# close handle
	close GTFFILE;
	# Test parsing i.e empty hash
	die " Error! Invalid parsing step (empty hash):\n- Check that your file is in .GTF format.\n" unless (keys(%h_transcript)>0);

	################################################################################
	print STDERR "Sort coordinates....\n" if ($verbosity > 5);
	# Sort features coordinates wrt to tr
	for my $t (values %h_transcript) {

		# Sort Exons from the transcripts
	   @{$t->{"feature"}} = sort { $a->{"start"} <=> $b->{"start"}} @{$t->{"feature"}};
	}

	## Add a upper level (hash key) with chromosome information for parralellization by chr
	print STDERR "Populate h_chr....\n" if ($verbosity > 5);
	foreach my $tr (keys %h_transcript){

		my $chr = $h_transcript{$tr}->{"chr"};
		$h_chr{$chr}->{$tr} = $h_transcript{$tr};
	}
	undef %h_transcript;
	# return hash
	return %h_chr;
}

sub parseGTFSplitchr75{

	my ($infile, $levelfilter, $verbosity) = @_;
	$infile ||=	"";	
	$levelfilter ||=	"exon";	
	my $baseName = basename($infile);
	  $verbosity ||= 0;

	open GTFFILE, "$infile" or die "Error! Cannot open GTF File ". $infile . ": ".$!;
	my @lines	=	<GTFFILE>;
	# Store data
	my %h_transcript;
	# Store data by chr
	my %h_chr = ();
	# counter line
	my $i = 0;

	foreach my $line  (@lines){
		chomp $line;
		if ($line =~ /^(\S+)\s+(\S+)\s+(\S+)\s+(\d+)\s+(\d+)\s+(\S+)\s+(\S+)\s+(\S+)\s+gene_id "(\S+)"; transcript_id "(\S+)"; (\w.*)$/ 
	||
			$line =~ /^(\S+)\s+(\S+)\s+(\S+)\s+(\d+)\s+(\d+)\s+(\S+)\s+(\S+)\s+(\S+)\s+gene_id "(\S+)"; transcript_id "(\S+)";$/
	|| 
	$line =~ /^(\S+)\s+(\S+)\s+(\S+)\s+(\d+)\s+(\d+)\s+(\S+)\s+(\S+)\s+(\S+)\s+gene_id "(\S+)"; (\w.*)$/ 
	){
			my $chr     =   $1; 
			my $biotype	=   $2; 
			my $level	=   $3; 
			my $beg_feat    =   $4; 
			my $end_feat    =   $5;  
			my $score	   =   $6; 
			my $strand      =   $7; 
			my $frame	=   $8;
			my $gene_id	=   $9;
			my $transcript	=   $10;
			my $extrafield	=	$11;

			if ($level eq $levelfilter) {
				# Transcript levels
				$h_transcript{$transcript}->{"chr"}	=   $chr;
				$h_transcript{$transcript}->{"biotype"}	=   $biotype;
				$h_transcript{$transcript}->{"startt"}	=   Utils::min2($h_transcript{$transcript}->{"startt"}, $beg_feat);
				$h_transcript{$transcript}->{"endt"}	=   Utils::max2($h_transcript{$transcript}->{"endt"}, $end_feat);
				$h_transcript{$transcript}->{"score"}	=   $score;
				$h_transcript{$transcript}->{"strand"}      =   $strand;
				$h_transcript{$transcript}->{"frame"}	=   $frame;
				$h_transcript{$transcript}->{"gene_id"}	=   $gene_id;

				# parse gtf attributes after transcript_id
				my $ref_attrib = parseExtrafield($extrafield, $verbosity,'transcript_type,gene_type');
				# Feature infos are stored in N hashes from a array reference
				# $h{$tr}->{"feature"} is the reference on a array
				# @{$h{$tr}->{"feature"}} to dereference it
				# To obtain a particular element in the array : ${ $h{$tr}->{"feature"} }[0]
				my %feature = ( "start" => $beg_feat, "end" => $end_feat, "feat_level" => $level);
				my %merge = (%{$ref_attrib}, %feature);
				push (@{$h_transcript{$transcript}->{"feature"}}, \%merge);	
			}	
		} elsif (($line =~ /^#/) || ($line =~ /^\s*$/) || $line =~ /^track/)  {
				# skip
		} else {
			die "\nError! Cannot parse this line in '$infile':\n", $line,"\n! Aborting.\nCheck the format of your input file\n";
		}
		# verbose
		if ($verbosity > 0){
			   Utils::showProgress2(scalar @lines, $i++, "Parse input file: ");
			}


		}
	# close handle
	close GTFFILE;
	# Test parsing i.e empty hash
	die " Error! Invalid parsing step (empty hash):\n- Check that your file is in .GTF format.\n" unless (keys(%h_transcript)>0);

	################################################################################
	print STDERR "Sort coordinates....\n" if ($verbosity > 5);
	# Sort features coordinates wrt to tr
	for my $t (values %h_transcript) {

		# Sort Exons from the transcripts
	   @{$t->{"feature"}} = sort { $a->{"start"} <=> $b->{"start"}} @{$t->{"feature"}};
	}

	## Add a upper level (hash key) with chromosome information for parralellization by chr
	print STDERR "Populate h_chr....\n" if ($verbosity > 5);
	foreach my $tr (keys %h_transcript){

		my $chr = $h_transcript{$tr}->{"chr"};
		$h_chr{$chr}->{$tr} = $h_transcript{$tr};
	}
	undef %h_transcript;
	# return hash
	return %h_chr;
}


sub parseGTFSplitchr75_V2{

	 my ($infile, $levels, $split, $filtertag,$biotype_filter, $verbosity,) = @_;
     $infile ||='';
     $levels ||='exon';
     $biotype_filter||='';
     $split ||= 0; # split data structure by chr
     $filtertag ||= '';  #  filter some interesting (beyond gene_id and transcript_id) for instance 'transcript_biotype' : use transcript_id,gene_id for just 12th fields
     $verbosity ||= 10;     
     my $baseName = basename($infile);     
     my @all_level = split (",", $levels);
     # Open file
     open GTFFILE, "$infile" or die "Error! Cannot open GTF File ". $infile . ": ".$!;
     my @lines = <GTFFILE>;
     #     Store data
     my %h_transcript;
     # counter line     
     my $i = 0;
     foreach my $line  (@lines){
          chomp $line;
          next if (($line =~ /^#/) || ($line =~ /^\s*$/) || $line =~ /^track/); # remove weird line
          # split by tab
          my ($chr, $biotype, $level, $beg_feat, $end_feat, $score, $strand, $frame,  $attributes) = split(/\t+/, $line);
        # verbose
		Utils::showProgress(scalar @lines, $i++, "Parse input file: ") if ($verbosity > 0);

		  next if ($biotype ne "" and $biotype ne  $biotype_filter);
		  
          # check number of columns
          if ( !defined $attributes ) {
               die "[$line] does not look like GTF... not 9 columns\n";
          }
          # extract interesting information
          my $refattrib = parseAttributes($attributes, $filtertag, $verbosity);
          # check for mandatory tags : gene_id and transcript_id

          die "[$line] does not contain 'gene_id' attribute...\n" unless (defined $refattrib->{'gene_id'});
          my $transcript = $refattrib->{'transcript_id'};
          my $gene_id = $refattrib->{'gene_id'};
          if (!defined $refattrib->{'transcript_biotype'}){
               ### check for protein_Coding
          } else{ 
               warn "Overlap will be done on all annotations\n";
          }

          # delete from ref since we already check they exist and we dont need them anymore
          delete $refattrib->{'transcript_id'};
          delete $refattrib->{'gene_id'};
          for my $lvl (@all_level) {
               if ($lvl eq $level){
                    $h_transcript{$gene_id}->{"chr"} =   $chr;
                    $h_transcript{$gene_id}->{"biotype"} =   $biotype;
                    $h_transcript{$gene_id}->{"startt"} =   Utils::min2($h_transcript{$gene_id}->{"startt"}, $beg_feat);
                    $h_transcript{$gene_id}->{"endt"} =   Utils::max2($h_transcript{$gene_id}->{"endt"}, $end_feat);
#                     $h_transcript{$gene_id}->{"transcript id"} =   $transcript;
                    # Feature infos are stored in N hashes from a array reference
                    # $h{$tr}->{"feature"} is the reference on a array
                    # @{$h{$tr}->{"feature"}} to dereference it
                    # To obtain a particular element in the array : ${ $h{$tr}->{"feature"} }[0]
#                     my %feature = ( "start" => $beg_feat, "end" => $end_feat, "feat_level" => $level, 'strand' => $strand, 'frame' => $frame);
                    # add refeature
#                     my %merge = (%{$refattrib}, %feature);

#                     push (@{$h_transcript{$gene_id}->{"feature"}}, \%merge);
               }
            }
            # verbose
     }
     # close handle
     close GTFFILE;
     # Test parsing i.e empty hash
     die " Error! Invalid parsing step (empty hash):\n- Check that your file is in .GTF format.\n- Or check that level '$levels' are present in field 3 of '$infile'... " unless (keys(%h_transcript)>0);

     # Sort trx wrt to chr and start coordinates (and exons start)
  #    print STDERR "\nSorting transcript data..." if ($verbosity>0);
#      my $refhordered = ExtractFromHash::sortTranscript(\%h_transcript);
# 
#      # clear unordered hash
#      undef %h_transcript;
		my $refh = \%h_transcript;
        # Do we split by chr?
        if ($split){
                my $refhchr = splitbyChr($refh, $verbosity);
                return $refhchr;

        } else {
                return $refh;
        }
}


# split hash by chr
sub splitbyChr{

        my ($hreftx, $verbosity) = @_;
        $verbosity ||= 0 ;

        die "splitbyChr: hash not defined " unless (ref($hreftx));
        my %h_chr;
        foreach my $tr (keys %{$hreftx}){

                my $chr = $hreftx->{$tr}->{"chr"};
                $h_chr{$chr}->{$tr} = $hreftx->{$tr};

        }
        # return ref hash split
        return \%h_chr;
}


sub parseAttributes{

     my ($string,  $filtertag, $verbosity) = @_;
     $verbosity ||= 10 ;

     my %attribs;
     die "parseAttributes: string '$string' not defined\n" unless defined $string;     
     my @extrafields = split(";", $string);;
     # store ids and additional information in second hash
     foreach my $attr ( @extrafields ) {
          next unless $attr =~ /^\s*(.+)\s(.+)$/;
          my $c_type  = $1;
          my $c_value = $2;
          $c_value=~ s/\"//g;
          if (exists($attribs{$c_type})){
               warn "Warning: key '$c_type' already exists with ",$attribs{$c_type}," value...\n" if ($verbosity > 15);
          }
          # nomenclature is to use : 'biotype' instead of 'type'
          $c_type =~ s/transcript_type/transcript_biotype/g;
          $c_type =~ s/gene_type/gene_biotype/g;
          $attribs{$c_type} = $c_value;
     }

     # filtering tag if defined
     if ($filtertag ne ""){
          # tag to be parse
          my @tag = split (",", $filtertag);
          # by default_add transcript and gene_id
          push @tag, "transcript_id", "gene_id";
          # perlish way of extracting a sub hash from an array
          my %filterattribs = map { $_ => $attribs{$_}} @tag;

          return \%filterattribs;
     } else {
          return \%attribs;
     }
}

# Sort Tx according to chr and start coordinates
# and sort exons forach transcript wrt to start coordinates
sub sortTranscript{

     my ($hreftx, $verbosity) = @_;
     $verbosity ||= 0 ;
     tie my %hash_sorted , 'Tie::IxHash';
     die "sortTranscript: hash not defined " unless (ref($hreftx));
     # sort exons
     $hreftx = sortExons($hreftx);
     # Sort features coordinates wrt to tr
     print STDERR "\nSorting Transcript data structure...\n" if ($verbosity > 0);
     # return an array of sorted transcript id
     my @keys_sorted =  sort {
                $hreftx->{$a}->{'chr'} cmp $hreftx->{$b}->{'chr'} or # sort by chr
                $hreftx->{$a}->{'startt'} <=> $hreftx->{$b}->{'startt'}  # and start coordinates
     } keys %{$hreftx} ;
     # map keys sorted to hash
     %hash_sorted = map { $_ => $hreftx->{$_} } @keys_sorted;
     return \%hash_sorted;

} 

sub cumulSize{
	my $refonarray 	=	shift;	
	my $verbosity 	= 	shift;
	$verbosity 		= 0 unless defined ($verbosity);
	
	# print if verbosity	
	print STDERR "Getting cumulative size...\n" if ($verbosity > 0);
	
	# cumulsize
	my $cumulsize	=	0;
	
	# Parse gtfHash to be printed
	foreach my $feat1 (@{$refonarray}) {

		$cumulsize += $feat1->{"end"} - $feat1->{"start"} +1;

	}
	return $cumulsize;
}

# Parse file at the exon level
sub parse_h{
	my ($href, $list,$st, $verbose) = @_;
	$verbose ||= 0;
	my %h_temp = %{$href};
	my $j=0;
	my $biotype;
	my $status;
	my %h_transcript;
	foreach my $tx(keys(%h_temp)){
		$biotype = ${$h_temp{$tx}->{"feature"}}[0]->{"transcript_biotype"};
		$status = ${$h_temp{$tx}->{"feature"}}[0]->{"transcript_status"};

		if ($list =~ /$biotype/ and $status eq $st){
			$h_transcript{$tx}=$h_temp{$tx};
		
		}
		if ($verbose> 0) {
			Utils::showProgress(scalar(keys(%h_temp)), $j++, "Selection of lncRNA: ");
		}
	}
	return %h_transcript;
}


# write seq in $outfile
sub gtf2fasta{
	my ($href, $fastaDirFile, $outfile, $verbose, $slop) = @_;
	my $i=0;
	$slop ||= 0;
	my %h_transcript = %{$href};
	my $h_transcript_size	= keys(%h_transcript); 
	my $seqdata;

	if ($verbose >0){ print STDERR "Get Sequences\n";}


	if (defined $outfile){
		open(OUTPUTFASTA, "> $outfile") || die ("# Error in writing fasta file") ;
	}

	# If Directory
	if (-d $fastaDirFile){

		for my $tr (keys %h_transcript){
	
			#Initalize sequence:
			my $seqstring	=	"";
			my $id_sequence	=	$h_transcript{$tr}->{'chr'};
			my $cpt=0;

	
			foreach my $exon (@{$h_transcript{$tr}->{"feature"}}) {
		
				$cpt++;

				# if first, we add $slop bp to START
				if ( $cpt == 1){
					$seqstring .=	StringUtils::getSubSequence($fastaDirFile."/".$h_transcript{$tr}->{"chr"}.".fa", $id_sequence, $exon->{"start"} - $slop, $exon->{"end"});
		
				# if last exon, we add $slop bp to END
				} elsif ($cpt == scalar(@{$h_transcript{$tr}->{"feature"}})) {
					$seqstring .=   StringUtils::getSubSequence($fastaDirFile."/".$h_transcript{$tr}->{"chr"}.".fa", $id_sequence, $exon->{"start"} , $exon->{"end"} + $slop);
				
				} else{
				$seqstring .=   StringUtils::getSubSequence($fastaDirFile."/".$h_transcript{$tr}->{"chr"}.".fa", $id_sequence, $exon->{"start"} , $exon->{"end"});
				}
			}
	
			#RevComp if strand -
			if ( $h_transcript{$tr}->{"strand"} eq '-'|| $h_transcript{$tr}->{"strand"} eq '-1') {
				$seqstring 	= StringUtils::getRevComp($seqstring);
			}


			# Summarize data e.g >TCONS_00005869 XLOC_001028_-_1:2753268-2784339_Cufflinks
			# and fasta sequence
			# header
			$seqdata 	 =  ">$tr ";
			$seqdata 	.= $h_transcript{$tr}->{"gene_id"}."_".$h_transcript{$tr}->{"strand"}."_".$h_transcript{$tr}->{"chr"}.":".$h_transcript{$tr}->{"startt"}."-".$h_transcript{$tr}->{"endt"}."_".$h_transcript{$tr}->{"biotype"}."\n";
			# sequence
			$seqdata	.= "$seqstring\n";
	
			# Print in a file
			if (defined $outfile){
				print OUTPUTFASTA $seqdata;
			# Or in STDOUT
			}else{
				print $seqdata;	
			}

			if ($verbose > 0){
				Utils::showProgress($h_transcript_size, $i++, "Print ".$tr.": ");
			}

		}

	} elsif (-r $fastaDirFile){
		## If multifasta file
	
		for my $tr (keys %h_transcript){
	
			#Initalize sequence:
			my $seqstring	=	"";
			my $id_sequence	=	$h_transcript{$tr}->{'chr'};
			my $cpt=0;

	
			foreach my $exon (@{$h_transcript{$tr}->{"feature"}}) {
		
				$cpt++;

				# if first, we add $slop bp to START
				if ( $cpt == 1){
					$seqstring .=	StringUtils::getSubSequenceSamtools($fastaDirFile, $id_sequence, $exon->{"start"} - $slop, $exon->{"end"});
		
				# if last exon, we add $slop bp to END
				} elsif ($cpt == scalar(@{$h_transcript{$tr}->{"feature"}})) {
					$seqstring .=   StringUtils::getSubSequenceSamtools($fastaDirFile, $id_sequence, $exon->{"start"} , $exon->{"end"} + $slop);
				
				} else{
					$seqstring .=   StringUtils::getSubSequenceSamtools($fastaDirFile, $id_sequence, $exon->{"start"} , $exon->{"end"});
				}
			}
	
			#RevComp if strand -
			if ( $h_transcript{$tr}->{"strand"} eq '-'|| $h_transcript{$tr}->{"strand"} eq '-1') {
				$seqstring 	= StringUtils::getRevComp($seqstring);
			}


			# Summarize data e.g >TCONS_00005869 XLOC_001028_-_1:2753268-2784339_Cufflinks
			# and fasta sequence
			# header
			$seqdata 	 =  ">$tr ";
			$seqdata 	.= $h_transcript{$tr}->{"gene_id"}."_".$h_transcript{$tr}->{"strand"}."_".$h_transcript{$tr}->{"chr"}.":".$h_transcript{$tr}->{"startt"}."-".$h_transcript{$tr}->{"endt"}."_".$h_transcript{$tr}->{"biotype"}."\n";
			# sequence
			$seqdata	.= "$seqstring\n";
	
			# Print in a file
			if (defined $outfile){
				print OUTPUTFASTA $seqdata;
			# Or in STDOUT
			}else{
				print $seqdata;	
			}

			if ($verbose > 0){
				Utils::showProgress($h_transcript_size, $i++, "Print ".$tr.": ");
			}

		}


	} else {
		pod2usage("Cannot understand defined option -f \"$fastaDirFile\"... ");
	}

	# Close file if defined
	close OUTPUTFASTA unless (!defined $outfile);
}    

sub parseExtrafield{

	my ($string, $verbosity, $filtertag) = @_;
	$filtertag ||= undef;
	$verbosity ||= 0 ;

	my %attribs;
	if ($string){
		my @extrafields = split(";", $string);

		# store ids and additional information in second hash
		foreach my $attr ( @extrafields ) {
			next unless $attr =~ /^\s*(.+)\s(.+)$/;
			my $c_type  = $1;
			my $c_value = $2;
			$c_value=~ s/\"//g;
			if (exists($attribs{$c_type}) && $verbosity > 10){
				warn "Warning: key '$c_type' already exists with ",$attribs{$c_type}," value...\n";
			}
			# filtering tag if defined
			if (defined $filtertag){
				# tag to be parse
				my @tag = split (",", $filtertag);
				if (grep /^$c_type$/, @tag) { $attribs{$c_type} = $c_value if($c_type  && $c_value)};
			} else {
				$attribs{$c_type} = $c_value if($c_type  && $c_value);
			}
		}
	}
	return \%attribs;
}	