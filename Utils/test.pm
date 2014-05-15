package test;

$VERSION = v0.0.1;

use strict;
use warnings;
use Data::Dumper;
use Pod::Usage;



$| = 1;

# test if a Directory can be open and have one or more fasta file
# input:
#	- dirname: directory name
#	- option: name of the option
#	- format: name of the format
sub  testFastaInDirectory {
	my (@args)= @_;
	my %opt = @args;
	my $count =0;
	my $dirname = $opt{"-dirname"};
	if ($opt{"-option"} && $opt{"-dirname"} ){
		opendir(DIR, $dirname) or pod2usage("-".$opt{"-option"}." option: I can't open the directory".$dirname."\n");
	}
	if ($opt{"-option"} && !$opt{"-dirname"} ){
		pod2usage ("-".$opt{"-option"}." option: I have not directory\n");
	}
	else{
		opendir(DIR, $dirname) or die "I can't open the directory".$dirname."\n";
	}
	while (my $file = readdir(DIR)){
		if ($opt{"-format"}){
			if ($opt{"-format"} == "fasta" && $opt{"-file"} =~ /.fa/){
				$count++;
			}
		}
		else {
			$count++;
		}
	}
	closedir(DIR);
	if ($opt{"-format"}&& $opt{"-option"}){
		pod2usage( "-".$opt{"-option"}." option: There is not ".$opt{"-format"}." file in your directory $dirname") if $count == 0;
	}
	elsif ($opt{"-format"}){
		pod2usage("There is not ".$opt{"-format"}." file in your directory $dirname") if $count == 0;				
	}
	else{
		die "There is not file in your directory $dirname" if $count == 0;
	}
}

sub testReadableNoEmptyFileOpt {
	my (@args) = @_;
	my %opt = @args;
	if (!$opt{"-file"} && $opt{"-option"}){
		pod2usage ("-".$opt{"-option"}." option: you have not given file ...\n");
	}
	
	elsif (!$opt{"-file"}){
		pod2usage ("You have not given file ...\n");
	}
	
	elsif ($opt{"-option"} && $opt{"-file"}){
		pod2usage ("-".$opt{"-option"}." option: I can't read the input file: ".$opt{"-file"}."...\n") unless (-r $opt{"-file"});
		pod2usage ("-".$opt{"-option"}." option: I need a non-empty input file:".$opt{"-file"}."...\n") unless (-s $opt{"-file"});
	}
	

	else{
		pod2usage ("I can't read the input file: ".$opt{"-file"}."...\n") unless (-r $opt{"-file"});
		pod2usage ("I need a non-empty input file:".$opt{"-file"}."...\n") unless (-s $opt{"-file"});
	}
}


sub testExecPathModule {
	my $module =shift;
	if ($module){
		pod2usage ("Can't exec ".$module."...\n") if (-x $module);
		pod2usage ($module." not in your PATH") unless (grep { -x "$_/".$module}split /:/,$ENV{PATH});
	}
	elsif(!$module){
		die ("You don't give module name...\n") ;
	}
}

1;