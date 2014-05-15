#!/bin/bash


##############################################################################################################
# 
# PROGRAM 	FEELnc.sh
#
# DATE 		22 May 2013
# 
# AUTHOR 	tderrien@univ-rennes1.fr 
# 
# GOAL
# 	From an input .gtf file (1) and a reference annotation file (2), the program encapsulates 3 sub programs
# 			- FEELnc_Filter.sh 					: Filtering steps from a .gtf file
#			- FEELnc_CodPot.sh (independant)	: Compute coding potential of .gtf file
# 			- FEELnc_Classifier.sh 				: Classify lncRNAs according to genomic position relative to reference annotation
##############################################################################################################

############
# PARAMETERS
############

# Prog name
PROGRAM=`basename $0`

# SCRIPTPATH
SCRIPTPATH=$( cd $(dirname $0) ; pwd -P )

# Bin directory
BIN_DIR="$SCRIPTPATH/bin/"

# Out directory name default
CLASSDIR="allClasses"

# min lncRNAs size
FILTER_SIZELNC=200

# Filter out Biexonic lncRNAs with min size exon
FILTER_BIEXONIC_SIZE=25

# Filter monoexonic putative lncRNAs
FILTER_MONOEXONIC="Y"

# Filter sense mRNAs overlaping lncRNAs
FILTER_SENSEOVERLAP="Y"

# 5' filtering
FPRIMEmRNA_FILTER=0

# 3' filtering
TPRIMEmRNA_FILTER=0

# min lncRNAs size
VERBOSITY=10


##############################################################################################################
# FONCTIONS
##############################################################################################################

########
# USAGE
########
function usage() {
    echo " " >&2
    echo "# USAGE	: ./$PROGRAM  -i <infile.gtf> -g <Known.annotation.gtf>">&2
    echo "# HELP	: ./$PROGRAM -h ">&2
    echo "# OPTIONS :">&2
		echo "#    -e <exonicCodingSense_FILTER>: filter out lncRNAs that overlap any exon feature in the annotation.gtf in Sense [default 'Y'] - 'N' to activate filtering ">&2
		echo "#    -m <monoexonic_FILTER>		: filter out mono-exonic transcript [default 'Y'] - 'N' to desactivate filtering - 'NEAS' to remove all monoexonic but NotExonicAntiSense ">&2
		echo "#    -b <biexonic_sizeFILTER>		: filter out bi-exonic transcript(s) having one exon < X nt size [default X='25']  ">&2
        echo "#    -s <size_lncrna_FILTER>		: filter out transcripts with size < X [default X=200]">&2
        echo "#    -5p <5'mRNA_size_FILTER>		: filter out transcripts overlapping Known.annotation with 5' extra region of X nt [default X=0] (strand specific)">&2        
        echo "#    -3p <3'mRNA_size_FILTER>		: filter out transcripts overlapping Known.annotation with 3' extra region of X nt [default X=0] (strand specific)">&2      
        echo "#    -o  <classdirectory>		: Output directory for temporary class [default 'allClasses']">&2                  
        echo "#    -v <verbosity>			: level of verbosity [default '10']">&2
        echo "# ">&2
        echo "# EXAMPLES :">&2
        echo "# ./$PROGRAM -i transcripts.gtf -g Canis_familiaris.gtf -e 1  ">&2
        echo "#		Filter lncRNAs (>=200 nt) from file 'transcripts.gtf' in directory 'lncRNAsClasses' and keep onl those not overlapping <annotation.gtf>">&2
	    echo " " >&2
}

# Programs path
# -------------
function test_bin(){

	# BEDTOOLS INTERSECT
	if [ -s "$BIN_DIR/intersectBed" ];then
		INTERSECT="$BIN_DIR/intersectBed"
	else
		type intersectBed >/dev/null 2>&1 || { echo >&2 "I require \"$INTERSECT\" from BedTools but it's not in your PATH nore in $BIN_DIR dir...  Aborting."; exit 1; }
		INTERSECT="intersectBed"
	fi

	# BEDTOOLS WINDOWBED
	if [ -s "./windowBed" ];then
		WINDOWBED="./windowBed"
	else
		type windowBed >/dev/null 2>&1 || { echo >&2 "I require \"$WINDOWBED\" from BedTools but it's not in your PATH nore in $BIN_DIR dir...  Aborting."; exit 1; }
		WINDOWBED="windowBed"
	fi

	# closestBed
	if [ -s "$BIN_DIR/closestBed" ];then
		CLOSEST="$BIN_DIR/closestBed"
	else
		type closestBed >/dev/null 2>&1 || { echo >&2 "I require \"closestBed\" from BedTools but it's not in your PATH nore in curent dir '`pwd`'.  Aborting."; exit 1; }
		CLOSEST="$BIN_DIR/closestBed"
	fi

	# exonGtf2TxGene.pl
	if [ -s "$BIN_DIR/exonGtf2TxGene.pl" ];then
		EXONGTF2TXGENE="$BIN_DIR/exonGtf2TxGene.pl"
	else
		type exonGtf2TxGene.pl >/dev/null 2>&1 || { echo >&2 "I require \"$EXONGTF2TXGENE\" but it's not in your PATH nore in $BIN_DIR dir...  Aborting."; exit 1; }
		EXONGTF2TXGENE="exonGtf2TxGene.pl"
	fi

	# gtfstat.sh
	if [ -s "$BIN_DIR/gtfstat.sh" ];then
		GTFSTAT="$BIN_DIR/gtfstat.sh"
	else
		type gtfstat.sh >/dev/null 2>&1 || { echo >&2 "I require \"$GTFSTAT\" but it's not in your PATH nore in $BIN_DIR dir...  Aborting."; exit 1; }
		GTFSTAT="gtfstat.sh"
	fi
	
	#make_introns
	if [ -s "$BIN_DIR/make_introns.awk" ];then
		MAKEINTRON="$BIN_DIR/make_introns.awk"
	else
		type make_introns.awk >/dev/null 2>&1 || { echo >&2 "I require \"$GTFSTAT\" but it's not in your PATH nore in $BIN_DIR dir...  Aborting."; exit 1; }
		MAKEINTRON="make_introns.awk"
	fi
	
	# Venn?
	if [ -s "$BIN_DIR/CreateVenn.sh" ];then
		CREATEVENN="$BIN_DIR/CreateVenn.sh"
	else
		type CreateVenn.sh >/dev/null 2>&1 || { echo >&2 "I require \"$CREATEVENN\" but it's not in your PATH nore in $BIN_DIR dir...  Aborting."; exit 1; }
		CREATEVENN="CreateVenn.sh"
	fi
}

# 1- FILTER
##############################################################################################################
function FEELnc_filter(){

	if [ -d "${CLASSDIR}" ];then
		echo " WARNING: Output directory for FEELnc_filter: \""${CLASSDIR}"\" already exists... ">&2
		echo " ">&2
	else
		mkdir ${CLASSDIR}
	fi

	###############################
	# Create exon and transcript levels 
	###############################
	if [ "$VERBOSITY" -ge 10 ];then
		echo "=> Create exon and transcripts levels for file : '$INFILE' ">&2
		$EXONGTF2TXGENE -i $INFILE  | awk -v CLASSDIR=${CLASSDIR} -v INFILEPREF=$INFILEPREF '{if ($3=="exon") {print $0 > CLASSDIR"/"INFILEPREF".filter.exon.gtf"}else if ($3=="transcript") {print $0 > CLASSDIR"/"INFILEPREF".filter.transcript.gtf"}}'

		echo "=> Create exon and transcripts levels for file '$ANNOTATIONFILE' ">&2
		$EXONGTF2TXGENE -i $ANNOTATIONFILE  | awk -v CLASSDIR=${CLASSDIR} -v ANNOTATIONFILEPREF=$ANNOTATIONFILEPREF '{if ($3=="exon") {print $0 > CLASSDIR"/"ANNOTATIONFILEPREF".coding.exon.gtf"}else if ($3=="transcript") {print $0 > CLASSDIR"/"ANNOTATIONFILEPREF".coding.transcript.gtf"}}'
	fi

	echo " "

	###############################
	# EXON Sense mRNA over filtering 
	###############################
	if [ "$FILTER_SENSEOVERLAP" == "Y" ] || [ "$FILTER_SENSEOVERLAP" == "y" ];then
	
		if [ "$VERBOSITY" -ge 10 ];then
			echo "=> Filtering step:	- Removing potential lncRNAs Overlaping exonic mRNA in sense ...">&2
		fi
	
		$INTERSECT -a ${CLASSDIR}/${INFILEPREF}.filter.exon.gtf -b  ${CLASSDIR}/${ANNOTATIONFILEPREF}.coding.exon.gtf -s -wa  | awk '{print $12}' | sort | uniq > ${CLASSDIR}//lncRNAs_FILTER_SENSEOVERLAP.txt
	
		# Get back fitlered transcript
		awk 'BEGIN{while (getline<"'${CLASSDIR}/'lncRNAs_FILTER_SENSEOVERLAP.txt">0){id[$1]++}}{if ($12 in id){}else{print $0}}'  ${CLASSDIR}/${INFILEPREF}.filter.exon.gtf > t
		mv t ${CLASSDIR}/${INFILEPREF}.filter.exon.gtf
	
		if [ "$VERBOSITY" -ge 0 ];then
			echo "=> 			:	- `wc -l ${CLASSDIR}/lncRNAs_FILTER_SENSEOVERLAP.txt` ... Check in file: '${CLASSDIR}/lncRNAs_FILTER_SENSEOVERLAP.txt'...">&2
		fi
	fi




	###############################
	# Size: cumulated exon filter 
	# the size is computed on exon levels of transcript
	# but to get correct ids for gene level (which does not hav transcript id)
	# need to print gene_id from correct transcript_id
	###############################
	if [ "$FILTER_SIZELNC" -gt 0 ];then

		if [ "$VERBOSITY" -ge 10 ];then
			echo "=> Filtering step:	- Removing LncRNAs < $FILTER_SIZELNC nt ...">&2
		fi

		# Compute size on exon level
		awk -v FILTER_SIZELNC=$FILTER_SIZELNC '{id[$12]+=$5-($4+1)}END{for (i in id){if (id[i] >= FILTER_SIZELNC){print i}}}' ${CLASSDIR}/${INFILEPREF}.filter.exon.gtf > .correct_id_sizeall.lst
		awk -v FILTER_SIZELNC=$FILTER_SIZELNC '{id[$12]+=$5-($4+1)}END{for (i in id){if (id[i] < FILTER_SIZELNC){print i, id[i]}}}' ${CLASSDIR}/${INFILEPREF}.filter.exon.gtf > ${CLASSDIR}/lncRNAs_FILTER_SIZELNC${FILTER_SIZELNC}.txt
	
		awk  'BEGIN{while (getline<".correct_id_sizeall.lst">0){id[$1]++}}{if ($12 in id){print $0}}' ${CLASSDIR}/${INFILEPREF}.filter.exon.gtf > t
		mv t ${CLASSDIR}/${INFILEPREF}.filter.exon.gtf
		
		# cleaning
		rm .correct_id_sizeall.lst

		if [ "$VERBOSITY" -ge 0 ];then
			echo "=> 			:	- `wc -l ${CLASSDIR}/lncRNAs_FILTER_SIZELNC${FILTER_SIZELNC}.txt` ... Check in file: '${CLASSDIR}/lncRNAs_FILTER_SIZELNC${FILTER_SIZELNC}.txt'...">&2
		fi
	fi


	###############################
	# FILTER_MONOEXONIC  filtering 
	###############################
	if [ "$FILTER_MONOEXONIC" == "Y" ] || [ "$FILTER_MONOEXONIC" == "y" ] || [ "$FILTER_MONOEXONIC" == "NEAS" ] || [ "$FILTER_MONOEXONIC" == "neas" ];then
	
		if [ "$VERBOSITY" -ge 10 ];then
			echo "=> Filtering step:	- Removing monoexonic LncRNAs option '$FILTER_MONOEXONIC'...">&2
		fi	

		awk  -v CLASSDIR=${CLASSDIR} '{id[$12]++}END{for (i in id){if (id[i]==1){print i > CLASSDIR"/lncRNAs_FILTER_MONOEXONIC.txt"  }else{print i}}}' ${CLASSDIR}/${INFILEPREF}.filter.exon.gtf  > .correct_id_mono.lst
		
		# get ALL monoexonic in .gtf 
		awk 'BEGIN{while (getline<"'${CLASSDIR}/lncRNAs_FILTER_MONOEXONIC.txt'">0){id[$1]++}}{if ($12 in id){print $0}}' ${CLASSDIR}/${INFILEPREF}.filter.exon.gtf > t
		mv t ${CLASSDIR}/lncRNAs_FILTER_MONOEXONIC.gtf		
		
		# if we keep mono exonic overlaping mRNA exon antisense , we add to the list ".correct_id_mono.lst"
		if [ "$FILTER_MONOEXONIC" == "NEAS" ] || [ "$FILTER_MONOEXONIC" == "neas" ];then
			$INTERSECT -a ${CLASSDIR}/lncRNAs_FILTER_MONOEXONIC.gtf -b  ${CLASSDIR}/${ANNOTATIONFILEPREF}.coding.exon.gtf -wa -wb | awk '{for (j=12;j<=NF;j++){if ($j == "gene_id"){field_geneid=j}}; strand2=$(field_geneid-2); if ($7!= strand2){print $12}}' | sort | uniq > .correct_mono_exonicAS.lst
		

			# add to other list of multi-exonic
			cat .correct_mono_exonicAS.lst .correct_id_mono.lst > t
			
			mv t .correct_id_mono.lst
			if [ "$VERBOSITY" -ge 0 ];then
				echo "=> 			:	- `wc -l .correct_mono_exonicAS.lst`  mono-exonic AS kept... ls .correct_mono_exonicAS.lst">&2
			fi
			
			# modify monoexonic list
			cat ${CLASSDIR}/lncRNAs_FILTER_MONOEXONIC.gtf .correct_mono_exonicAS.lst | awk 'NF>1{print $12}NF==1{print $1}' | sort | uniq -c | awk '$1==1{print $2}' > .temp
			awk  'BEGIN{while (getline<"'.temp'">0){id[$1]++}}{if ($12 in id){print $0}}' ${CLASSDIR}/${INFILEPREF}.filter.exon.gtf > t
			mv t ${CLASSDIR}/lncRNAs_FILTER_MONOEXONIC.gtf
		fi 
		
		awk  'BEGIN{while (getline<"'.correct_id_mono.lst'">0){id[$1]++}}{if ($12 in id){print $0}}' ${CLASSDIR}/${INFILEPREF}.filter.exon.gtf > t
		mv t ${CLASSDIR}/${INFILEPREF}.filter.exon.gtf
				
		#cleaning
		rm .correct_id_mono.lst

		# Create file if it does not exist i.e 0 monoexonic lncRNAs_FILTER_MONOEXONIC.txt
		# touch will not RE-recreate it if it already exists
		touch ${CLASSDIR}/lncRNAs_FILTER_MONOEXONIC.gtf
		if [ "$VERBOSITY" -ge 0 ];then
			echo "=> 			:	- `wc -l ${CLASSDIR}/lncRNAs_FILTER_MONOEXONIC.gtf` ... Check in file: '${CLASSDIR}/lncRNAs_FILTER_MONOEXONIC.gtf'...">&2
		fi

	fi


	###############################
	# FITLER bi exonic dubious exon
	# 	In order to remove bi-exonic lncRNA that have an exon < XX nucleotide (dubious annotation) [default XX= 10]
	###############################
	if [ "$FILTER_BIEXONIC_SIZE" -gt 0 ];then

		if [ "$VERBOSITY" -ge 10 ];then
			echo "=> Filtering step:	- Removing biexonic LncRNAs with exon < $FILTER_BIEXONIC_SIZE nt ...">&2
		fi

		awk '{id[$12]++}END{for (i in id){if (id[i]==2){print i }}}' ${CLASSDIR}/${INFILEPREF}.filter.exon.gtf > .tmp.biexonic 
	
		# Compute size on exon level
		awk -v FILTER_BIEXONIC_SIZE=$FILTER_BIEXONIC_SIZE 'BEGIN{while (getline<".tmp.biexonic">0){id[$1]++}}{if ($12 in id){if ( ($5-$4+1) < FILTER_BIEXONIC_SIZE){rmbiexonic[$12]++}}}END{for (i in rmbiexonic){print i}}' ${CLASSDIR}/${INFILEPREF}.filter.exon.gtf > ${CLASSDIR}/lncRNAs_FILTER_BIEXONIC_SIZE_$FILTER_BIEXONIC_SIZE.txt
	
		awk  'BEGIN{while (getline<"'${CLASSDIR}/lncRNAs_FILTER_BIEXONIC_SIZE_$FILTER_BIEXONIC_SIZE.txt'">0){id[$1]++}}{if ($12 in id){}else{print $0}}' ${CLASSDIR}/${INFILEPREF}.filter.exon.gtf > t
		mv t ${CLASSDIR}/${INFILEPREF}.filter.exon.gtf
		
		# cleaning
		rm .tmp.biexonic 

		if [ "$VERBOSITY" -ge 0 ];then
			echo "=> 			:	- `wc -l ${CLASSDIR}/lncRNAs_FILTER_BIEXONIC_SIZE_${FILTER_BIEXONIC_SIZE}.txt` ... Check in file: '${CLASSDIR}/lncRNAs_FILTER_BIEXONIC_SIZE_${FILTER_BIEXONIC_SIZE}.txt'...">&2
		fi
	fi



	###############################
	# 5': lncRNA filter
	###############################
	if [ "$FPRIMEmRNA_FILTER" -gt 0 ];then

		if [ "$VERBOSITY" -ge 10 ];then
			echo "=> Filtering step:	- Removing LncRNAs in a window of $FPRIMEmRNA_FILTER nt before mRNAs...">&2
		fi

		# windowBed options
		# -l left size region*
		# -r right size region *
		# -sw define according to strand
		# * both needed at the same time
		# -v	Only report those entries in A that have _no overlaps_ with B. 	- Similar to "grep -v."	
	
		#
		$WINDOWBED -a ${CLASSDIR}/${INFILEPREF}.filter.transcript.gtf -b  ${CLASSDIR}/${INFILEPREF}.coding.transcript.gtf -l $FPRIMEmRNA_FILTER -r $TPRIMEmRNA_FILTER -u -sw  | awk '{print $12}'  > ${CLASSDIR}/lncRNAs_FILTER_FPRIMEmRNA_FILTER${FPRIMEmRNA_FILTER}.txt

		awk 'BEGIN{while (getline<"'${CLASSDIR}/lncRNAs_FILTER_FPRIMEmRNA_FILTER${FPRIMEmRNA_FILTER}.txt'">0){id[$1]++}}{if ($12 in id){}else{print $0}}'  ${CLASSDIR}/${INFILEPREF}.filter.exon.gtf > t
		mv t ${CLASSDIR}/${INFILEPREF}.filter.exon.gtf
	
	
		if [ "$VERBOSITY" -ge 0 ];then
			echo "=> 			:	- `wc -l ${CLASSDIR}/lncRNAs_FILTER_FPRIMEmRNA_FILTER${FPRIMEmRNA_FILTER}.txt` Removed lncRNAs in mRNAs 5' are stored in file: '${CLASSDIR}/lncRNAs_FILTER_FPRIMEmRNA_FILTER${FPRIMEmRNA_FILTER}.txt'...">&2
		fi
	fi

	###############################
	# 3': lncRNA filter
	###############################
	if [ "$TPRIMEmRNA_FILTER" -gt 0 ];then

		if [ "$VERBOSITY" -ge 10 ];then
			echo "=> Filtering step:	- Removing LncRNAs in a window of $TPRIMEmRNA_FILTER nt after mRNAs...">&2
		fi

		$WINDOWBED -a ${CLASSDIR}/${INFILEPREF}.filter.transcript.gtf -b  ${CLASSDIR}/${INFILEPREF}.coding.transcript.gtf -l $FPRIMEmRNA_FILTER -r $TPRIMEmRNA_FILTER -u -sw  | awk '{print $12}'  > ${CLASSDIR}/lncRNAs_FILTER_FPRIMEmRNA_FILTER${FPRIMEmRNA_FILTER}.txt

		awk 'BEGIN{while (getline<"'${CLASSDIR}/lncRNAs_FILTER_TPRIMEmRNA_FILTER${TPRIMEmRNA_FILTER}.txt'">0){id[$1]++}}{if ($12 in id){}else{print $0}}'  ${CLASSDIR}/${INFILEPREF}.filter.exon.gtf > t
		mv t ${CLASSDIR}/${INFILEPREF}.filter.exon.gtf
	
	
		if [ "$VERBOSITY" -ge 0 ];then
			echo "=> 			:	- `wc -l ${CLASSDIR}/lncRNAs_FILTER_TPRIMEmRNA_FILTER${TPRIMEmRNA_FILTER}.txt` Removed lncRNAs in mRNAs 3' are stored in file: '${CLASSDIR}/lncRNAs_FILTER_TPRIMEmRNA_FILTER${TPRIMEmRNA_FILTER}.txt'...">&2
		fi
	fi


	##########################################
	# Check if exist potential lncRNAs
	if [ ! -s "${CLASSDIR}/${INFILEPREF}.filter.exon.gtf" ];then
		echo " Warning: There is no non-coding category in '$INFILE' after filtering! Try modifying your parameters... ">&2
		exit 0;
	
	else
		echo "- End Filters! ">&2
		echo "- Check files: '${CLASSDIR}/${INFILEPREF}.filter.exon.gtf' ">&2
		echo "- Stats on filtered file : ">&2
		gtfstat.sh 	${CLASSDIR}/${INFILEPREF}.filter.exon.gtf
	fi


}
##############################################################################################################
test_bin;
##############################################################################################################



#####################
# Parse command line
#####################
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
		if [ "$2" != "" ];then
			CLASSDIR=$2
		else
			echo " ERROR: Cannot read your output dirrectory FEELnc_filter... ">&2
			echo " ">&2		
			exit
		fi
	elif [ "$1" == "-g" ];then
		ANNOTATIONFILE=$2
	elif [ "$1" == "-s" ];then
		FILTER_SIZELNC=$2
	elif [ "$1" == "-m" ];then
		FILTER_MONOEXONIC=$2
	elif [ "$1" == "-b" ];then
		FILTER_BIEXONIC_SIZE=$2
	 elif [ "$1" == "-e" ];then
		FILTER_SENSEOVERLAP=$2       
	 elif [ "$1" == "-5p" ];then
		FPRIMEmRNA_FILTER=$2  
	 elif [ "$1" == "-3p" ];then
		TPRIMEmRNA_FILTER=$2  
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


# Test Correctness of arguments
#### Test parameters

if [[ "$FILTER_BIEXONIC_SIZE" == *[!0-9]* ]] || [[ "$FILTER_BIEXONIC_SIZE" -lt 0 ]];then
	echo " ERROR: Min Exonic size   for biExonic filtering \""$FILTER_BIEXONIC_SIZE"\" should not be a non-negative integer... ">&2
	echo " ">&2
	usage
	exit;
fi

if [[ "$FILTER_SIZELNC" == *[!0-9]* ]] || [[ "$FILTER_SIZELNC" -lt 0 ]];then
	echo " ERROR: LncRNAs size cutoff \""$FILTER_SIZELNC"\" should not be a non-negative integer... ">&2
	echo " ">&2
	usage
	exit;
fi
if [[ "$FPRIMEmRNA_FILTER" == *[!0-9]* ]] || [[ "$FPRIMEmRNA_FILTER" -lt 0 ]];then
	echo " ERROR: 5' window size  \""$FPRIMEmRNA_FILTER"\" should not be a non-negative integer... ">&2
	echo " ">&2
	usage
	exit;
fi
if [[ "$TPRIMEmRNA_FILTER" == *[!0-9]* ]] || [[ "$TPRIMEmRNA_FILTER" -lt 0 ]];then
	echo " ERROR: 3' window size  \""$FPRIMEmRNA_FILTER"\" should not be a non-negative integer... ">&2
	echo " ">&2
	usage
	exit;
fi


if [ "$FILTER_SENSEOVERLAP" != "Y" ] && [ "$FILTER_SENSEOVERLAP" != "N" ] && [ "$FILTER_SENSEOVERLAP" != "y" ] && [ "$FILTER_SENSEOVERLAP" != "n" ];then
	echo " ERROR: Filtering Sense lncRNAs overlapping exonic mRNAs in Sense = \""$FILTER_SENSEOVERLAP"\" whereas it should be ither \"Y\" (or \"y\") to activate filter or \"N\" (or \"n\") de desactivate... ">&2
	echo " ">&2
	usage
	exit;
fi

if [ "$FILTER_MONOEXONIC" != "Y" ] && [ "$FILTER_MONOEXONIC" != "N" ] && [ "$FILTER_MONOEXONIC" != "y" ] && [ "$FILTER_MONOEXONIC" != "n" ] && [ "$FILTER_MONOEXONIC" != "NEAS" ] && [ "$FILTER_MONOEXONIC" != "neas" ];then
	echo " ERROR: Filtering monoexonic parazmeter = \""$FILTER_MONOEXONIC"\" whereas it should be either Y or N or NEAS (or y|n|neas)... ">&2
	echo " ">&2
	usage
	exit;
fi

if [[ "$VERBOSITY" == *[!0-9]* ]] || [[ "$VERBOSITY" -lt 0 ]];then
	echo " ERROR: Level of verbosity \""$VERBOSITY"\" should not be a non-negative integer... ">&2
	echo " ">&2
	usage
	exit;
fi

### INFILE test
if [ "$VERBOSITY" -gt 0 ];then
	echo "=> Parsing input file : '$INFILE' and $ANNOTATIONFILE...">&2 
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

INFILEPREF=`basename $INFILE .gtf`
ANNOTATIONFILEPREF=`basename $ANNOTATIONFILE .gtf`

# Input File
if [ "$VERBOSITY" -gt 0 ];then
	echo "=> Basic Statistcis: '$INFILE'">&2 
	gtfstat.sh 	$INFILE
fi


echo " 
=> Summarizing Input Parameters:
- INFILE=$INFILE
- ANNOTATIONFILE=$ANNOTATIONFILE
- FILTER_SENSEOVERLAP=$FILTER_SENSEOVERLAP
- FILTER_BIEXONIC_SIZE=$FILTER_BIEXONIC_SIZE
- FILTER_MONOEXONIC=$FILTER_MONOEXONIC
- FILTER_SIZELNC=$FILTER_SIZELNC
- FPRIMEmRNA_FILTER=$FPRIMEmRNA_FILTER
- TPRIMEmRNA_FILTER=$TPRIMEmRNA_FILTER
- CLASSDIR=${CLASSDIR}
- VERBOSITY=$VERBOSITY
">&2

##############################################################################################################
FEELnc_filter;
##############################################################################################################

# Cleaning
rm ${CLASSDIR}/${ANNOTATIONFILEPREF}*
rm ${CLASSDIR}/${INFILEPREF}.filter.transcript.gtf