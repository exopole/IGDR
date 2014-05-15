#!/bin/bash


##############################################################################################################
# 
# PROGRAM 	FEELnc_Classifier2.sh
#
# DATE 		June 2013
# 
# AUTHOR 	tderrien@univ-rennes1.fr 
# 
# GOAL
#			Classify a set of input transcript in file1 wrt to a set of second transcript file2
#			Categories are based on file1 txs
#			- Intergenic 	: Divergent - Convergent - SameOrient
#			- Genic			: Exonic	- Intronic	- Encompassing (with sense and antisense)
#
# USAGE
#			./FEELnc.sh -i test.gtf -g /omaha-beach/tderrien/DATA/canFam3/annotation/Ensembl70_RefSeq/canFam3_v70.RefSeq.3lvl.format.protein_coding.gtf 
##############################################################################################################

############
# PARAMETERS
############

# Prog name
PROGRAM=`basename $0`

# Bin directory
# SCRIPTPATH
SCRIPTPATH=$( cd $(dirname $0) ; pwd -P )

# Bin directory
BIN_DIR="$SCRIPTPATH/bin/"

# Out directory name default
CLASSDIR="Classes"


# Verbosity
VERBOSITY=10


##############################################################################################################
# FONCTIONS
##############################################################################################################

########
# USAGE
########
function usage() {
    echo " " >&2
    echo "# USAGE	: ./$PROGRAM  -i <file1.gtf> -g <file2.gtf>">&2
    echo "# HELP	: ./$PROGRAM -h ">&2
    echo "# OPTIONS :">&2          
        echo "#    -o  <CLASSDIRectory>		: Output directory for classes [default 'Classes']">&2                  
        echo "#    -v <verbosity>			: level of verbosity [default '10']">&2
        echo "# ">&2
        echo "# EXAMPLES :">&2
        echo "# ./$PROGRAM -i transcripts.gtf -g Canis_familiaris.gtf">&2
        echo "#		Classify transcripts from file transcripts.gtf wrt to input reference annotation (-g option)">&2
	    echo " " >&2
}


##############################################################################################################
# Programs path
##############################################################################################################
# OVERLAP_BEDTOOLS
if [ -s "$BIN_DIR/overlap" ];then
	OVERLAP_BEDTOOLS="$BIN_DIR/overlap"
elif [ -s "$BIN_DIR/getOverlap" ];then
	OVERLAP_BEDTOOLS="$BIN_DIR/getOverlap"
elif type getOverlap >/dev/null 2>&1 ;then 
	OVERLAP_BEDTOOLS="getOverlap"
elif type  overlap >/dev/null 2>&1 ;then 
	OVERLAP_BEDTOOLS="overlap"
else
	echo >&2 "I require \"overlap\" (or getOverlap) from BedTools but it's not in your PATH nore in curent dir '`pwd`'.  Aborting."; exit 1; 
fi

# MAKEINTRONS
if [ -s "$BIN_DIR/make_introns.awk" ];then
	MAKEINTRONS="$BIN_DIR/make_introns.awk"
else
	type make_introns.awk >/dev/null 2>&1 || { echo >&2 "I require \"make_introns.awk\" but it's not in your PATH nore in curent dir '`pwd`'.  Aborting."; exit 1; }
	MAKEINTRONS=`type make_introns.awk | awk '{print $3}'`
fi


# INTERSECT
if [ -s "$BIN_DIR/intersectBed" ];then
	INTERSECT="$BIN_DIR/intersectBed"
else
	type intersectBed >/dev/null 2>&1 || { echo >&2 "I require \"$INTERSECT\" from BedTools but it's not in your PATH nore in curent dir '`pwd`'.  Aborting."; exit 1; }
	INTERSECT="intersectBed"
fi

# WINDOWBED
if [ -s "$BIN_DIR/windowBed" ];then
	WINDOWBED="$BIN_DIR/windowBed"
else
	type windowBed >/dev/null 2>&1 || { echo >&2 "I require \"$INTERSECT\" from BedTools but it's not in your PATH nore in curent dir '`pwd`'.  Aborting."; exit 1; }
	WINDOWBED="windowBed"
fi

# closestBed
if [ -s "$BIN_DIR/closestBed" ];then
	CLOSEST="$BIN_DIR/closestBed"
else
	type closestBed >/dev/null 2>&1 || { echo >&2 "I require \"closestBed\" from BedTools but it's not in your PATH nore in curent dir '`pwd`'.  Aborting."; exit 1; }
	CLOSEST="closestBed"
fi

# exonGtf2TxGene.pl
if [ -s "$BIN_DIR/exonGtf2TxGene.pl" ];then
	EXONGTF2TXGENE="$BIN_DIR/exonGtf2TxGene.pl"
else
	type exonGtf2TxGene.pl >/dev/null 2>&1 || { echo >&2 "I require \"$EXONGTF2TXGENE\" but it's not in your PATH nore in $BIN_DIR dir...  Aborting."; exit 1; }
	EXONGTF2TXGENE="exonGtf2TxGene.pl"
fi



##############################################################################################################
# Parse command line
##############################################################################################################
if [ $# -lt 4 ];then
	echo "ERROR: Wrong number of arguments...">&2
	usage
	exit;
elif [ "$1" == "-h" ] || [ "$1" == "-help" ] || [ "$1" == "--help" ] || [ "$1" == "" ];then
	usage
	exit;
fi
# Command line
# -------------
while [ "$1" ];do 
	if [ "$1" == "-i" ];then
		INFILE=$2
	elif [ "$1" == "-o" ];then
		CLASSDIR=$2
	elif [ "$1" == "-g" ];then
		ANNOTATIONFILE=$2
	 elif [ "$1" == "-v" ];then
		VERBOSITY=$2 
	else
		echo " ERROR: in Command_line ">&2
		echo "	Parameter option \""$1"\" with \""$2"\" value not known! ">&2
		echo " ">&2
		usage
		exit;
	fi
	shift 2
done


##############################################################################################################
# Test input
##############################################################################################################
### INFILE test
if [ "$VERBOSITY" -gt 0 ];then
	echo "=> Parsing input file 1: '$INFILE'...">&2 
fi
#  infile not GTF
test_gtf=`awk '{cpt++; if ($1 !~/^#/ && $1 !~/^track/){ if (($0 !~/gene_id/ || $0 !~/transcript_id/) || ($7 !="+" && $7 !="-" && $7 !=".") ){ print cpt; exit}}}' $INFILE`
if [ -n "$test_gtf" ] ;then
		wrongline=`awk -v test_gtf=$test_gtf 'NR==test_gtf' $INFILE`
		echo " ERROR:  Your file `basename $INFILE` does not seem to be .gtf file !">&2
		echo "Line $test_gtf :
		$wrongline ">&2
		echo " ">&2
		usage;
		exit 1;
fi


#  ANNOTATIONFILE not GTF
if [ "$VERBOSITY" -gt 0 ];then
	echo "=> Parsing input file 2: '$ANNOTATIONFILE'...">&2 
fi
test_gtf=`awk '{cpt++; if ($1 !~/^#/ && $1 !~/^track/){ if (($0 !~/gene_id/ || $0 !~/transcript_id/) || ($7 !="+" && $7 !="-" && $7 !=".") ){ print cpt; exit}}}' $ANNOTATIONFILE`
if [ -n "$test_gtf" ] ;then
		wrongline=`awk -v test_gtf=$test_gtf 'NR==test_gtf' $ANNOTATIONFILE`
		echo " ERROR:  Your file `basename $ANNOTATIONFILE` does not seem to be .gtf file !">&2
		echo "Line $test_gtf :
		$wrongline ">&2
		echo " ">&2
		usage;
		exit 1;
fi

## TEST IF EXON LEVEL
text_exon=`awk '$3=="exon"' $INFILE | wc -l`
if [ "$text_exon" -eq 0 ];then
		echo " ERROR:  Your file `basename $INFILE` does not have exon level in field 3... mandatory !">&2
		echo " ">&2
		usage;
		exit 1;    
fi

text_exon=`awk '$3=="exon"' $ANNOTATIONFILE | wc -l`
if [ "$text_exon" -eq 0 ];then
		echo " ERROR:  Your file `basename $ANNOTATIONFILE` does not have exon level in field 3... mandatory !">&2
		echo " ">&2
		usage;
		exit 1;    
fi

## END TEST ##################################################################################################

INFILEPREF=`basename $INFILE .gtf`
ANNOTATIONFILEPREF=`basename $ANNOTATIONFILE .gtf`

# Input Files
##############################################################################################################
if [ "$VERBOSITY" -gt 0 ];then
	echo "=> Basic Statistcis FILE 1: '$INFILE'">&2 
# 	gtfstat.sh 	$INFILE
	echo "=> Basic Statistcis FILE 2: '$ANNOTATIONFILE'">&2 
	# gtfstat.sh 	$ANNOTATIONFILE
fi




###############################
# Create Class dir
###############################
if [ -d "$CLASSDIR" ];then
	echo " WARNING: Output directory for FEELnc_filter: \""$CLASSDIR"\" already exists... ">&2
	echo " ">&2
else
	mkdir $CLASSDIR
fi

###############################
# Create exon and transcript levels 
###############################
if [ "$VERBOSITY" -ge 10 ];then
	echo "=> Create exon and transcripts levels for file : '$INFILE' ">&2
	$EXONGTF2TXGENE -i $INFILE  | gawk -v CLASSDIR=$CLASSDIR -v INFILEPREF=$INFILEPREF '{if ($3=="exon") {print $0 > CLASSDIR"/"INFILEPREF".exon.gtf"}else if ($3=="transcript") {print $0 > CLASSDIR"/"INFILEPREF".transcript.gtf"}}'

	echo "=> Create exon and transcripts levels for file '$ANNOTATIONFILE' ">&2
	$EXONGTF2TXGENE -i $ANNOTATIONFILE  | gawk -v CLASSDIR=$CLASSDIR -v ANNOTATIONFILEPREF=$ANNOTATIONFILEPREF '{if ($3=="exon") {print $0 > CLASSDIR"/"ANNOTATIONFILEPREF".exon.gtf"}else if ($3=="transcript") {print $0 > CLASSDIR"/"ANNOTATIONFILEPREF".transcript.gtf"}}'
fi

echo "";

# Important Files:
INFILEPREF_TX="$CLASSDIR/${INFILEPREF}.transcript.gtf"
INFILEPREF_EX="$CLASSDIR/${INFILEPREF}.exon.gtf"

echo "Stats"
gtfstat.sh 	$CLASSDIR/${INFILEPREF}.exon.gtf
echo ""

# Coding files
ToBeCompared_file_t="$CLASSDIR/${ANNOTATIONFILEPREF}.transcript.gtf"
ToBeCompared_file_e="$CLASSDIR/${ANNOTATIONFILEPREF}.exon.gtf"

# Merge both file to retrieve exon coordinates in genic class
cat $INFILEPREF_TX $ToBeCompared_file_t > file1file2

##############################################################################################################

if [ "$VERBOSITY" -ge 5 ];then
	echo "# Create Overlapping non-coding categories (in \"$CLASSDIR\")" >&2
fi
##############################################################################################################
# 1- INTERGENIC
####################

# Create empty class files
touch ${CLASSDIR}/${INFILEPREF}.transcript.noover.gtf
touch ${CLASSDIR}/${INFILEPREF}.transcript.noover.closestSS.gtf2
touch ${CLASSDIR}/${INFILEPREF}.transcript.noover.closestDiv.gtf2
touch ${CLASSDIR}/${INFILEPREF}.transcript.noover.closestConv.gtf2
# If a lncRNA is on chr without any annotation (often unknown chromosome)
touch ${CLASSDIR}/${INFILEPREF}.transcript.noover.closestNone.gtf2

$INTERSECT -a ${INFILEPREF_TX} -b ${ToBeCompared_file_t} -wa  -v  >  $CLASSDIR/${INFILEPREF}.transcript.noover.gtf
if [ "$VERBOSITY" -ge 5 ];then
	echo "- Non overlapping (Intergenic): `wc -l $CLASSDIR/${INFILEPREF}.transcript.noover.gtf`" >&2
fi

## Unknow chrosomosome:
# ---------------------
$CLOSEST -a $CLASSDIR/${INFILEPREF}.transcript.noover.gtf -b ${ToBeCompared_file_t} -t first |\
awk '{delete (arr); cpt=0;for (i=1;i<=NF;i++){ if ( $(i) == "gene_id"){cpt++; arr[cpt]=i}}
if (length(arr)==1){print $0}
}' > ${CLASSDIR}/${INFILEPREF}.transcript.noover.closestNone.gtf2
## Same strand      : 
# -------------------
# closestBed with -t option => report only the first one
$CLOSEST -a $CLASSDIR/${INFILEPREF}.transcript.noover.gtf -b ${ToBeCompared_file_t} -t first |\
awk '{delete (arr); cpt=0;for (i=1;i<=NF;i++){ if ( $(i) == "gene_id"){cpt++; arr[cpt]=i}}; 
if (length (arr)>1){
strd1=$(arr[1]-2)
strd2=$(arr[2]-2);
if (strd1 == strd2){print $0}}}' > ${CLASSDIR}/${INFILEPREF}.transcript.noover.closestSS.gtf2

## Divergent        : 
#--------------------
$CLOSEST -a $CLASSDIR/${INFILEPREF}.transcript.noover.gtf -b ${ToBeCompared_file_t} -t first| \
awk  '{delete (arr); cpt=0;for (i=1;i<=NF;i++){ if ( $(i) == "gene_id"){cpt++; arr[cpt]=i}}; 
if (length (arr)>1){
strd1=$(arr[1]-2);
strd2=$(arr[2]-2);
start1=$(arr[1]-5);
start2=$(arr[2]-5);
if (strd1 != strd2 && ( (start1<start2 && strd1 =="-") || (start1>start2 && strd1=="+")) ){print $0}}}' > ${CLASSDIR}/${INFILEPREF}.transcript.noover.closestDiv.gtf2

# Convergent        : 
#--------------------
$CLOSEST -a $CLASSDIR/${INFILEPREF}.transcript.noover.gtf -b ${ToBeCompared_file_t} -t first | \
awk  '{delete (arr); cpt=0;for (i=1;i<=NF;i++){ if ( $(i) == "gene_id"){cpt++; arr[cpt]=i}}; 
if (length (arr)>1){
strd1=$(arr[1]-2);
strd2=$(arr[2]-2);
start1=$(arr[1]-5);
start2=$(arr[2]-5);
if (strd1 != strd2 && ( (start1<start2 && strd1 =="+") || (start1>start2 && strd1=="-")) ){print $0}}}' > ${CLASSDIR}/${INFILEPREF}.transcript.noover.closestConv.gtf2

if [ "$VERBOSITY" -ge 5 ];then
	echo "-		Non-Overlapping SS   	:		`wc -l ${CLASSDIR}/${INFILEPREF}.transcript.noover.closestSS.gtf2`" >&2
	echo "-		Non-Overlapping Div 	:		`wc -l ${CLASSDIR}/${INFILEPREF}.transcript.noover.closestDiv.gtf2`" >&2
	echo "-		Non-Overlapping Conv	:		`wc -l ${CLASSDIR}/${INFILEPREF}.transcript.noover.closestConv.gtf2`" >&2
	echo "-		Non-Overlapping None	:		`wc -l ${CLASSDIR}/${INFILEPREF}.transcript.noover.closestNone.gtf2`" >&2
fi



##############################################################################################################
# 2- GENIC
##############################################################################################################

# Tx
$INTERSECT -a ${INFILEPREF_TX} -b ${ToBeCompared_file_t} -wa  -u > ${CLASSDIR}/${INFILEPREF}.transcript.over.gtf

# Get the Exons coordinates of the overlapping tx in FILE1
 awk -v fileRef=${CLASSDIR}/${INFILEPREF}.transcript.over.gtf 'BEGIN{while (getline<fileRef>0){tr[$12]++}} {if ($12 in tr){print $0}}' ${INFILEPREF_EX}  > ${CLASSDIR}/${INFILEPREF}.exon.over.gtf

# Redefinition of variables : FILE1 levels overlapping FILE2
INFILEPREF_TX_OVER="${CLASSDIR}/${INFILEPREF}.transcript.over.gtf"
INFILEPREF_EX_OVER="${CLASSDIR}/${INFILEPREF}.exon.over.gtf"

if [ "$VERBOSITY" -ge 5 ];then
	echo "-  Overlapping (Genic) '$ANNOTATIONFILEPREF'  : `wc -l ${CLASSDIR}/${INFILEPREF}.transcript.over.gtf`" >&2
fi

    ############
    # EXON 
    ############
    
    # -1 Get exon coordinates
    # 2 -- We intersect FILE1 exons with FILE2 exons 
    # (3 -- Then format in gtf2 )
    # 4 -- Then compute the nb bp of Tx FILE1 overlapped all exons of FILE2 exons : 
    # (overlap stdin => therefore sort | uniq (if UTR and exon... twice count) ) => remove because exon filtered before
    # 5 the line  sort -k4,4rn -k1,1 ==> sort by overlap size and tx1tx2 IDs to get the best tx2 partner
    # awk '{split($1,pair,"_"); seen[pair[1]]++; if (seen[pair[1]]==1){print $0}}' ===> display only the best (first) one
 	
	$INTERSECT -a  ${INFILEPREF_EX_OVER} -b ${ToBeCompared_file_e} -wa -wb | \
    # format just to get ~bed positions as bedtools overlap make a seg fault 11 if native gtf is given in stdin
       awk '{delete (arr); cpt=0;for (i=1;i<=NF;i++){ if ( $(i) == "gene_id"){cpt++; arr[cpt]=i}}; 
chr1=$(arr[1]-8);
chr2=$(arr[2]-8);
start1=$(arr[1]-5);
start2=$(arr[2]-5);
end1=$(arr[1]-4);
end2=$(arr[2]-4);
txid1=$(arr[1]+3);
txid2=$(arr[2]+3); OFS="\t"; print chr1,start1,end1,txid1,chr2,start2,end2,txid2}' | $OVERLAP_BEDTOOLS -i stdin -cols 2,3,6,7 | awk '{pair_sizeOver[$4"#"$8]+=$NF+1; 
	
	    # 	then compute the proportion of mRNAs overlaped i.e (sum of overlapped exons) / (sum of all mRNA exons)
    		txoverlaped_size[$8]+=$7-$6+1
		}END{ for (i in pair_sizeOver){
				split (i, ar, "#");
				
				printf ("%s\t%.2f\n", i, pair_sizeOver[i]*100/txoverlaped_size[ar[2]])
			}
			}' | sort -k2,2rn -k1,1 | awk '{split($1,pair,"#"); seen[pair[1]]++; if (seen[pair[1]]==1){print $0}}' > ${CLASSDIR}/${INFILEPREF}.bestpartnerOverlap.exon.lst 
    
    
# 	Extract trancript level information from FILE2 for the 2 partners and write in Sense or AntiSense files
	awk -v INFILE=$INFILE -v CLASSDIR=$CLASSDIR -v INFILEPREF=$INFILEPREF 'BEGIN{
		while (getline<"file1file2">0){ 
			if ($3 =="transcript"){
				file1file2[$12]=$0
			}
		}}{
		split ($1,pairid,"#");
		if (pairid[1] in file1file2){tx1=file1file2[pairid[1]]; split(tx1,f,"\t"); strand_lnc=f[7]}
		if (pairid[2] in file1file2){tx2=file1file2[pairid[2]];split(tx2,f,"\t"); strand_partner=f[7]}			

		# Write in Sense and AntiSense files		
		if (strand_lnc == strand_partner){
			filenameS=CLASSDIR"/"INFILEPREF".transcript.over.BEST_exon.S.gtf2";
			print tx1, tx2, $2 > filenameS
		} else {
			filenameAS=CLASSDIR"/"INFILEPREF".transcript.over.BEST_exon.AS.gtf2"
			print tx1, tx2, $2 > filenameAS		
		}
		
	}' ${CLASSDIR}/${INFILEPREF}.bestpartnerOverlap.exon.lst 

	# Test if exists TXS in files 1 and create an empty one if no TX1 in this category
	if [ ! -s "${CLASSDIR}/${INFILEPREF}.transcript.over.BEST_exon.S.gtf2" ];then
		touch ${CLASSDIR}/${INFILEPREF}.transcript.over.BEST_exon.S.gtf2
	fi
	if [ ! -s "${CLASSDIR}/${INFILEPREF}.transcript.over.BEST_exon.AS.gtf2" ];then
		touch ${CLASSDIR}/${INFILEPREF}.transcript.over.BEST_exon.AS.gtf2
	fi

	# Combine S et AS
	cat ${CLASSDIR}/${INFILEPREF}.transcript.over.BEST_exon.S.gtf2 ${CLASSDIR}/${INFILEPREF}.transcript.over.BEST_exon.AS.gtf2 > ${CLASSDIR}/${INFILEPREF}.transcript.over.BEST_exon.gtf2
	#Redefinition of 
	INFILEPREF_EX_OVER_TX="${CLASSDIR}/${INFILEPREF}.transcript.over.BEST_exon.gtf2"

if [ "$VERBOSITY" -ge 5 ];then

	echo "-		Overlapping Exon   	:		`wc -l ${CLASSDIR}/${INFILEPREF}.transcript.over.BEST_exon.gtf2`" >&2
	echo "-		Overlapping Exon S 	:		`wc -l ${CLASSDIR}/${INFILEPREF}.transcript.over.BEST_exon.S.gtf2`" >&2
	echo "-		Overlapping Exon AS	:		`wc -l ${CLASSDIR}/${INFILEPREF}.transcript.over.BEST_exon.AS.gtf2`" >&2
fi	

# Cleaning
# if [ -s "${CLASSDIR}/${INFILEPREF}.bestpartnerOverlap.exon.lst" ];then
# 	rm ${CLASSDIR}/${INFILEPREF}.bestpartnerOverlap.exon.lst
# fi


	############
    # INTRONS : 
    ############
    # 1 --  COnstruct Introns from  FILE exons : TODO Check why l(introns) <0 
    # Note: Maybe to do with only coding transcipt
    sort -k12,12 -k4,4n -k5,5n $ToBeCompared_file_e | \
    awk -v fldgn=10 -v fldtr=12 -f $MAKEINTRONS | \
    awk '$5>$4'   > ${CLASSDIR}/${ANNOTATIONFILEPREF}.intron.gtf


	# Overlapping LnRNAs left: : Get exons
	awk -v fileRef=$INFILEPREF_EX_OVER_TX 'BEGIN{while (getline<fileRef>0){tr[$12]=$0}} {if ($12 in tr) {} else {print $0}}' $INFILEPREF_EX_OVER 	> .temp
	$INTERSECT -a  .temp -b ${CLASSDIR}/${ANNOTATIONFILEPREF}.intron.gtf -wa -wb | \
    # format just to get ~bed positions as bedtools overlap make a seg fault 11 if native gtf is given in stdin
    awk '{delete (arr); cpt=0;for (i=1;i<=NF;i++){ if ( $(i) == "gene_id"){cpt++; arr[cpt]=i}}; 
chr1=$(arr[1]-8);
chr2=$(arr[2]-8);
start1=$(arr[1]-5);
start2=$(arr[2]-5);
end1=$(arr[1]-4);
end2=$(arr[2]-4);
txid1=$(arr[1]+3);
txid2=$(arr[2]+3); OFS="\t"; print chr1,start1,end1,txid1,chr2,start2,end2,txid2}'  | $OVERLAP_BEDTOOLS -i stdin -cols 2,3,6,7 | awk '{pair_sizeOver[$4"#"$8]+=$NF+1; 
	
	    # 	then compute the proportion of mRNAs overlaped i.e (sum of overlapped exons) / (sum of all mRNA exons)
    		txoverlaped_size[$8]+=$7-$6+1
		}END{ for (i in pair_sizeOver){
				split (i, ar, "#"); 
				printf ("%s\t%.2f\n", i, pair_sizeOver[i]*100/txoverlaped_size[ar[2]])
			}
			}' | sort -k2,2rn -k1,1 | awk '{split($1,pair,"#"); seen[pair[1]]++; if (seen[pair[1]]==1){print $0}}' > ${CLASSDIR}/${INFILEPREF}.bestpartnerOverlap.intron.lst
		
			
	# 	Extract trancript level information from the gencode file for the 2 partners and write in Sense or AntiSense files
	awk -v INFILE=$INFILE -v CLASSDIR=$CLASSDIR -v INFILEPREF=$INFILEPREF 'BEGIN{
		while (getline<"file1file2">0){ 
			if ($3 =="transcript"){
				file1file2[$12]=$0
			}
		}}{
		split ($1,pairid,"#");
		if (pairid[1] in file1file2){tx1=file1file2[pairid[1]]; split(tx1,f,"\t"); strand_lnc=f[7]}
		if (pairid[2] in file1file2){tx2=file1file2[pairid[2]];split(tx2,f,"\t"); strand_partner=f[7]}			

		# Write in Sense and AntiSense files from AWK		
		if (strand_lnc == strand_partner){
			filenameS=CLASSDIR"/"INFILEPREF".transcript.over.BEST_intron.S.gtf2";
			print tx1, tx2, $2 > filenameS
		} else {
			filenameAS=CLASSDIR"/"INFILEPREF".transcript.over.BEST_intron.AS.gtf2"
			print tx1, tx2, $2 > filenameAS		
		}
		
	}' ${CLASSDIR}/${INFILEPREF}.bestpartnerOverlap.intron.lst 

	# Test if exists Tx1 in these files and create an empty one if none in this category
	if [ ! -s "${CLASSDIR}/${INFILEPREF}.transcript.over.BEST_intron.S.gtf2" ];then
		 touch ${CLASSDIR}/${INFILEPREF}.transcript.over.BEST_intron.S.gtf2
	fi
	if [ ! -s "${CLASSDIR}/${INFILEPREF}.transcript.over.BEST_intron.AS.gtf2" ];then
		 touch ${CLASSDIR}/${INFILEPREF}.transcript.over.BEST_intron.AS.gtf2
	fi

	# Combine S et AS
	cat ${CLASSDIR}/${INFILEPREF}.transcript.over.BEST_intron.S.gtf2 ${CLASSDIR}/${INFILEPREF}.transcript.over.BEST_intron.AS.gtf2 > ${CLASSDIR}/${INFILEPREF}.transcript.over.BEST_intron.gtf2
	#Redefinition of 
	INFILEPREF_INT_OVER_TX="${CLASSDIR}/${INFILEPREF}.transcript.over.BEST_intron.gtf2"


if [ "$VERBOSITY" -ge 5 ];then
	echo "-		Overlapping Intron   	:		`wc -l ${CLASSDIR}/${INFILEPREF}.transcript.over.BEST_intron.gtf2`" >&2
	echo "-		Overlapping Intron S 	:		`wc -l ${CLASSDIR}/${INFILEPREF}.transcript.over.BEST_intron.S.gtf2`" >&2
	echo "-		Overlapping Intron AS	:		`wc -l ${CLASSDIR}/${INFILEPREF}.transcript.over.BEST_intron.AS.gtf2`" >&2
fi

# Cleaning
rm ${CLASSDIR}/${ANNOTATIONFILEPREF}.intron.gtf


    ###############
    # ENCOMPASSING
    ###############
    
    # 1 : Gather Overlapping Tx1 from previous EXON + INTRON 
	cat ${INFILEPREF_INT_OVER_TX} ${INFILEPREF_EX_OVER_TX} > .temp
	# Get the TRANSCRIPT level .gtf of the left Tx1s (why Tx level? because at this stage it should not have overlap ar exon level with coding)
	awk 'BEGIN{while (getline<".temp">0){seen_over[$12]++}}{if ($12 in seen_over){}else{print $0}}' $INFILEPREF_TX_OVER > ${CLASSDIR}/${INFILEPREF}.transcript.over.encompass.gtf


    #  Overlap at the transcript level and get the transcript with maximum cov (as previously)
    $INTERSECT -a  ${CLASSDIR}/${INFILEPREF}.transcript.over.encompass.gtf -b ${ToBeCompared_file_t} -wa -wb  | \
    # format just to get ~bed positions as bedtools overlap make a seg fault 11 if native gtf is given in stdin
    awk '{delete (arr); cpt=0;for (i=1;i<=NF;i++){ if ( $(i) == "gene_id"){cpt++; arr[cpt]=i}}; 
chr1=$(arr[1]-8);
chr2=$(arr[2]-8);
start1=$(arr[1]-5);
start2=$(arr[2]-5);
end1=$(arr[1]-4);
end2=$(arr[2]-4);
txid1=$(arr[1]+3);
txid2=$(arr[2]+3); OFS="\t"; print chr1,start1,end1,txid1,chr2,start2,end2,txid2}' | $OVERLAP_BEDTOOLS -i stdin -cols 2,3,6,7 | awk '{pair_sizeOver[$4"#"$8]+=$NF+1; 
	
	    # 	then compute the proportion of mRNAs overlaped i.e (sum of overlapped exons) / (sum of all mRNA exons)
    		txoverlaped_size[$8]+=$7-$6+1
		}END{ for (i in pair_sizeOver){
				split (i, ar, "#"); 
				printf ("%s\t%.2f\n", i, pair_sizeOver[i]*100/txoverlaped_size[ar[2]])
			}
			}' | sort -k2,2rn -k1,1 | awk '{split($1,pair,"#"); seen[pair[1]]++; if (seen[pair[1]]==1){print $0}}' > ${CLASSDIR}/${INFILEPREF}.bestpartnerOverlap.encompass.lst
    
    
# 	Extract trancript level information from the gencode file for the 2 partners and write in Sense or AntiSense files
	awk -v INFILE=$INFILE -v CLASSDIR=$CLASSDIR -v INFILEPREF=$INFILEPREF 'BEGIN{
		while (getline<"file1file2">0){ 
			if ($3 =="transcript"){
				file1file2[$12]=$0
			}
		}}{
		split ($1,pairid,"#");
		if (pairid[1] in file1file2){tx1=file1file2[pairid[1]]; split(tx1,f,"\t"); strand_lnc=f[7]}
		if (pairid[2] in file1file2){tx2=file1file2[pairid[2]];split(tx2,f,"\t"); strand_partner=f[7]}			

		# Write in Sense and AntiSense files		
		if (strand_lnc == strand_partner){
			filenameS=CLASSDIR"/"INFILEPREF".transcript.over.BEST_encompass.S.gtf2";
			print tx1, tx2, $2 > filenameS
		} else {
			filenameAS=CLASSDIR"/"INFILEPREF".transcript.over.BEST_encompass.AS.gtf2"
			print tx1, tx2, $2 > filenameAS		
		}
		
	}' ${CLASSDIR}/${INFILEPREF}.bestpartnerOverlap.encompass.lst 

	# Test if exists Tx1s in these files and create an empty one if none in this category
	if [ ! -s "${CLASSDIR}/${INFILEPREF}.transcript.over.BEST_encompass.S.gtf2" ];then
		touch ${CLASSDIR}/${INFILEPREF}.transcript.over.BEST_encompass.S.gtf2
	fi
	if [ ! -s "${CLASSDIR}/${INFILEPREF}.transcript.over.BEST_encompass.AS.gtf2" ];then
		touch ${CLASSDIR}/${INFILEPREF}.transcript.over.BEST_encompass.AS.gtf2
	fi

	# Combine S et AS
	cat ${CLASSDIR}/${INFILEPREF}.transcript.over.BEST_encompass.S.gtf2 ${CLASSDIR}/${INFILEPREF}.transcript.over.BEST_encompass.AS.gtf2 > ${CLASSDIR}/${INFILEPREF}.transcript.over.BEST_encompass.gtf2

if [ "$VERBOSITY" -ge 5 ];then

	echo "-		Overlapping Encompass   :		`wc -l ${CLASSDIR}/${INFILEPREF}.transcript.over.BEST_encompass.gtf2`" >&2
	echo "-		Overlapping Encompass S :		`wc -l ${CLASSDIR}/${INFILEPREF}.transcript.over.BEST_encompass.S.gtf2`" >&2
	echo "-		Overlapping Encompass AS:		`wc -l ${CLASSDIR}/${INFILEPREF}.transcript.over.BEST_encompass.AS.gtf2`" >&2
fi	

# cleaning
rm ${CLASSDIR}/${INFILEPREF}.transcript.over.encompass.gtf

##############################################################################################################
# cleaning
rm file1file2
rm ${INFILEPREF_TX}
rm ${INFILEPREF_EX}
rm ${ToBeCompared_file_t}
rm ${ToBeCompared_file_e}
rm ${INFILEPREF_TX_OVER}
rm ${INFILEPREF_EX_OVER}
rm  ${INFILEPREF_INT_OVER_TX} ${INFILEPREF_EX_OVER_TX}

#############	END ################

