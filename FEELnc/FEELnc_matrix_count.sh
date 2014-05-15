#!/bin/bash

# Aim : 
#	- Read a matrix file (0 - 1)
#	- count the sum of each fields
######################################

if [ ! -r "$1" ];then
	echo "ERROR: cannot read your input file">&2
	echo "USAGE:
	`basename $0` <matrix_file>">&2
	exit ;
else
	infile=$1
fi

awk 'NR>1{sum=0;for (i=2; i<=NF;i++){sum+=$(i)}; print $1, sum}' $infile
