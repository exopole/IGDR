#!/bin/bash


# May 2013
# tderrien@univ-rennes1.fr
# AIM:
#       - Compute coding potential of transcript.gtf based on a combination of3 methods
#
# INPUT:
#		- a .gtf file with exon coordinates
#		- a parameter file "FEELnc_CodPot.param"
# OUTPUT
#		- a .gtf file with scores for each method and flag lncRNA "1" or "0"
############################################################################################

# Dependecoes
##############
# VennDiagram R: library('VennDiagram') 
# Coding potential tools

############
# PARAMETERS
############

# Prog name
PROGRAM=`basename $0`
# SCRIPTPATH
SCRIPTPATH=$( cd $(dirname $0) ; pwd -P )
# Bin directory
BIN_DIR="$SCRIPTPATH/bin/"

# PARAMETER FILE
PARAMFILE="${BIN_DIR}/FEELnc_CodPot.param"

# Genome Sequence dir
GENOMEDIR="/omaha-beach/tderrien/DATA/canFam3/sequence/seq_by_chr/"
# Methods defaults
GENEID=0;
CPAT=0
TXCDSPREDICT=0;
# Update June 14th 2013
# Add non coding independently found by CPC program
CPC=0;

# Verbose
VERBOSITY=10


##############################################################################################################
########
# USAGE
########
function usage() {
    echo " " >&2
    echo "# USAGE	: ./$PROGRAM  -i <transcript.gtf> -p <parameter_file> [OPTIONS]">&2
    echo "# HELP	: ./$PROGRAM -h ">&2
    echo "# OPTIONS :">&2
        echo "#    -c <cpat>			: compute cpat coding potential -- See 'FEELnc_CodPot.param' for config [default 'no' ]">&2
        echo "#    -g <geneid>		: compute gene_id coding potential -- See 'FEELnc_CodPot.param' for config [default 'no' ]">&2
        echo "#    -t <txcdspredict>		: compute TxCdsPredict coding potential -- See 'FEELnc_CodPot.param' for config [default 'no' ]">&2
		echo "#    -f <cpc_fileNonCoding>: Filename with noncoding txs analyzed by CPC program [default 'no'] ">&2
		echo "#    -s <seqs_fasta_genomic>	: full path to a multi-fasta file with the genomic sequences for all input mappings, OR a directory with single-fasta files (one per genomic sequence, with file names matching sequence names) [default : /omaha-beach/tderrien/DATA/canFam3/sequence/seq_by_chr/]">&2
		echo "#    -o <outputfile>		: name of the output file [default 'infile.feelnc_codpot'] ">&2
		echo "#    -r <outdirectory>		: name of the output directory [default in current directory ] ">&2		
		echo "#    -d <diagram_venn>		: Directory name where venn data will be stored [default 'no'] ">&2
        echo "#    -v <verbosity>		: level of verbosity [default '10']">&2
        echo "# ">&2
        echo "# EXAMPLES :">&2
        echo "# ./$PROGRAM -i transcripts.gtf -p FEELnc_CodPot.param -cgt -d myVennDir # also work with ./$PROGRAM -i transcripts.gtf -p FEELnc_CodPot.param -c -g -t -d myVennDir">&2
	    echo " Get Coding potential for transcript in 'transcripts.gtf' by 3 methods: cpat (-c option), geneid(-g option) and txcdspredict (-t) and compute a Venn diagram of the intersection/union of each method" >&2
	    echo " " >&2
        echo "# ./$PROGRAM -i transcripts.gtf -p FEELnc_CodPot.param -c">&2
	    echo " Get Coding potential for transcript in 'transcripts.gtf' by only cpat method (-c option)" >&2
	    echo " " >&2
}

##############################################################################################################
# HELP
#####################
if [ "$1" == "-h" ] || [ "$1" == "-help" ] || [ "$1" == "--help" ] || [ "$1" == "" ];then
	echo "Needing help...">&2
	usage
	exit;
elif [ $# -lt 4 ];then
	echo "ERROR: Wrong number of arguments...">&2
	usage
	exit;
fi

##############################################################################################################
# PARSING command line
while getopts ":i:p:o:s:v:d:r:cpcf:gct" opt; do
	case $opt in
   i) # GTF File 
    	# is not readable
    	if [ ! -r "$OPTARG" ];then
    		echo "Option -i (input file) = '$OPTARG'... cannot read input file $OPTARG" >&2
    		exit;
    	else 
    		INFILE=$OPTARG;
    	fi
      	;;
    p)	# PARAM FILE
    	# File is not readable
    	if [ ! -r "$OPTARG" ];then
    		echo "Option -p (parameter file) = '$OPTARG'... cannot read parameter_file $OPTARG" >&2
    		exit;
    	else 
    		PARAMFILE=$OPTARG;
    	fi
      	;;

	o)	# OUTFILE 
		# Not writable
    	if [ -r "$OPTARG" ] && [ ! -w "$OPTARG" ];then
    		echo "Option -o (output file) = '$OPTARG'... file is not writable" >&2
    		exit;
    	else 
    		OUTFILE=$OPTARG;
    	fi
      	;;
    s)	# FASTA SEQUENCE 
		# Dir File is not readable
    	if [ ! -d "$OPTARG" ] && [ ! -r "$OPTARG" ];then
    		echo "Option -s cannot read $OPTARG" >&2
    		exit;
    	else 
    		GENOMEDIR=$OPTARG;
    	fi
      	;;      	
    v)	# verbose
    	if [[ "$OPTARG" == *[!0-9]* ]] ;then
    		echo "Option -v (verbosity) = '$OPTARG'... should be an integer" >&2
    		exit;
    	else 
			VERBOSITY=$OPTARG;
		fi
      	;;
    d)	# Venn diagram directory
    		DIAGVENN=$OPTARG;
      	;;       	
    r)	# Venn diagram directory
    		OUTDIR=$OPTARG;
      	;;  
    g)	# Launch GENEID
		GENEID=1;
      	;;
    c)	# Launch CPAT
		CPAT=1;
      	;;
    t)	# Launch TXCDSPREDICT
		TXCDSPREDICT=1;
      	;;
   f) # CPC file with 
    	# is not readable
    	if [ ! -r "$OPTARG" ];then
    		echo "Option -f (cpc input file) = '$OPTARG'... cannot read input file $OPTARG" >&2
    		exit;
    	else 
    		CPCFILE=$OPTARG;
    		# We turn CPC to 1 to be added to the 
    		CPC=1;
    	fi
      	;;      	
    \?)	echo "Invalid option: -$OPTARG" >&2
		usage
		exit 1
		;;
    :) echo "Option -$OPTARG requires an argument." >&2
    	usage;
		exit 1
		;;
	*) echo "Unimplimented option: -$OPTARG" >&2; exit 1;;
	esac
done

##############################################################################################################
# Test parameters
### test Methods
if [ -z "$CPAT" ]&& [ -z "$GENEID" ]&&[ -z "$TXCDSPREDICT" ] &&[ -z "$CPC" ];then

	echo " ERROR:  None of the 3 methods (i.e cpat, geneid or txcdspredict) (option -c, -g or -t) are flagged...">&2
	echo " $PROGRAM requires at least one!">&2
	exit 1;
else
	let "NUMBER_METHOD = $CPAT + $GENEID + $TXCDSPREDICT +$CPC"
fi
if [ "$NUMBER_METHOD" == "0" ] ;then
		echo " ERROR:  You need to precise at least one method i.e -c or -g or -t options...">&2
		echo " ">&2
		usage;
		exit 1;
fi



### test Outdir
if [ -d "${OUTDIR}" ];then
	echo " WARNING: Output directory for FEELnc_CodPot: \""${OUTDIR}"\" already exists... ">&2
	echo " ">&2
elif [ ! -z "$OUTDIR" ];then
	mkdir ${OUTDIR}
else
	OUTDIR=`pwd`
fi
	
	
### INFILE test
if [ "$VERBOSITY" -gt 0 ];then
	echo "- Parsing input file : '$INFILE' ...">&2 
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
## TEST IF EXON LEVEL
text_exon=`awk '$3=="exon"' $INFILE | wc -l`
if [ "$text_exon" -eq 0 ];then
		echo " ERROR:  Your file `basename $INFILE` does not have exon level in field 3... mandatory !">&2
		echo " ">&2
		usage;
		exit 1;    
fi



##############################################################################################################
# Outfile init
if [ -z "$OUTFILE" ];then
	OUTFILE=$(basename $INFILE).feelnc_codpot
fi
cp $INFILE ${OUTDIR}/${OUTFILE}
OUTFILE=${OUTDIR}/${OUTFILE}
OUTFILE_PREF=$(basename $OUTFILE .feelnc_codpot)


##############################################################################################################
# Get Fasta sequence with cufflink gffread (fast)
#   -w  write a fasta file with spliced exons for each GFF transcript
if [ "$VERBOSITY" -gt 0 ];then
	echo "- Get FASTA sequence for ${INFILE}">&2

fi
# I used gffread but I had to add extra nucleotide sequence around tx to improve coding potential calcul (at least for GeneId)
# gffread ${INFILE} -g ${GENOMEDIR} -w ${INFILE_FASTA}
# So now I used a own perl prog : gtf2fasta.pl
INFILE_FASTA=$(basename $INFILE .gtf).fa
SLOP=25 # number of extra nucleotide around transcript
echo "perl ${BIN_DIR}/gtf2fasta.pl -i ${INFILE} -f ${GENOMEDIR} -o ${OUTDIR}/${INFILE_FASTA} -s ${SLOP} -v ${VERBOSITY}"
perl ${BIN_DIR}/gtf2fasta.pl -i ${INFILE} -f ${GENOMEDIR} -o ${OUTDIR}/${INFILE_FASTA} -s ${SLOP} -v ${VERBOSITY}


if [ "$VERBOSITY" -gt 0 ];then
	echo "- FASTA sequence are stored in '${INFILE_FASTA}'">&2

fi

##############################################################################################################
# Run paramfile according to SHELL
shel=`ps -p $$ | awk 'NR==2{print $NF}'`
if [ "$shel" == "bash" ];then 
	source ${PARAMFILE}
else
	. ${PARAMFILE}
fi

##############################################################################################################
# Compute CodPot
# Launch Coding Potential Program
if [[ "$CPAT" == "1" ]];then

	if [ "$VERBOSITY" -gt 0 ];then
		echo "- Launch CPAT program...">&2
	fi
	
	METHOD_DIR="${OUTDIR}/DIR_CPAT";
	INFILE_METHOD=$(basename $INFILE .gtf).CPAT
	
	# create method dir
	if [ -d "$METHOD_DIR" ];then
		echo " WARNING: ${METHOD_DIR} directory  already exists... ">&2
	else
		mkdir $METHOD_DIR
	fi
	
	# Launch Mathod prog
	${CPAT_BIN} -d ${CPAT_DAT}/Human_train.RData -x ${CPAT_DAT}/Human_Hexamer.tab -g ${OUTDIR}/${INFILE_FASTA} -o ${METHOD_DIR}/${INFILE_METHOD}
	
	# Add cpat score and orf size to .gtf infile
	awk 'BEGIN{while(getline<"'${METHOD_DIR}/${INFILE_METHOD}'">0){id["\""$1"\";"]=$NF"#"$3}}{
		if ($12 in id){
			split (id[$12],ar,"#");
			print $0,"cpat_score \""ar[1]"\"; cpat_orf \""ar[2]"\";"
		}
	}' $OUTFILE > .tmp
	mv .tmp $OUTFILE

	if [ "$VERBOSITY" -gt 0 ];then
		echo "- Results are stored in : '${METHOD_DIR}' ...">&2
	fi	 
	
fi

#### GENEID ##################################################################################################
if [[ "$GENEID" == "1" ]];then

	if [ "$VERBOSITY" -gt 0 ];then
		echo "- Launch GeneId program...">&2
	fi
	
	METHOD_DIR="${OUTDIR}/DIR_GENEID";
	INFILE_METHOD=$(basename $INFILE .gtf).GENEID
	
	# create method dir
	if [ -d "$METHOD_DIR" ];then
		echo " Warning: ${METHOD_DIR} directory  already exists... ">&2
	else
		mkdir $METHOD_DIR
	fi
	
	# Launch Method prog
	${GENEID_BIN} -i ${OUTDIR}/${INFILE_FASTA} -P ${GENEID_DAT}  > ${METHOD_DIR}/${INFILE_METHOD}
	
	# Add cpat score to infile
	awk 'BEGIN{while(getline<"'${METHOD_DIR}/${INFILE_METHOD}'">0){id["\""$1"\";"]=$4"#"$3}}{
		
		if ($12 in id){
			split (id[$12],ar,"#");
			print $0,"geneid_score \""ar[1]"\"; geneid_orf \""ar[2]"\";"
		}
	}' $OUTFILE > .tmp
	mv .tmp $OUTFILE

	if [ "$VERBOSITY" -gt 0 ];then
		echo "- Results are stored in : '${METHOD_DIR}' ...">&2
	fi	 
fi

#### TxCdsPred ##################################################################################################
if [[ "$TXCDSPREDICT" == "1" ]];then

	if [ "$VERBOSITY" -gt 0 ];then
		echo "- Launch TxCdsPredict program...">&2
	fi
	
	METHOD_DIR="${OUTDIR}/DIR_TXCDSPREDICT";
	INFILE_METHOD=$(basename $INFILE .gtf).TXCDSPREDICT
	
	# create method dir
	if [ -d "$METHOD_DIR" ];then
		echo " Warning: ${METHOD_DIR} directory  already exists... ">&2
	else
		mkdir $METHOD_DIR
	fi
	
	# Launch Method prog
	${TXCDSPREDICT_BIN} ${OUTDIR}/${INFILE_FASTA} ${METHOD_DIR}/${INFILE_METHOD}
	 
	 # Add TXCDSPREDICT score to infile
	awk 'BEGIN{while(getline<"'${METHOD_DIR}/${INFILE_METHOD}'">0){id["\""$1"\";"]=$6}}{
		if ($12 in id){
			print $0,"txcdspredict_score \""id[$12]"\";";
		}
	}' $OUTFILE > .tmp
	mv .tmp $OUTFILE

	if [ "$VERBOSITY" -gt 0 ];then
		echo "- Results are stored in : '${METHOD_DIR}' ...">&2
	fi	 
fi

#### TxCdsPred ##################################################################################################
if [[ "$CPC" == "1" ]];then

	if [ "$VERBOSITY" -gt 0 ];then
		echo "- Get CPC noncoding...">&2
	fi
	
	METHOD_DIR="${OUTDIR}/DIR_CPC";
	INFILE_METHOD=$(basename $INFILE .gtf).CPC
	
	# create method dir
	if [ -d "$METHOD_DIR" ];then
		echo " Warning: ${METHOD_DIR} directory  already exists... ">&2
	else
		mkdir $METHOD_DIR
	fi
	 
	 # cp file
	 cp $CPCFILE ${METHOD_DIR}/$(basename $INFILE .gtf).CPC
	 
	# Add CPC score to infile
	awk 'BEGIN{while(getline<"'${METHOD_DIR}/${INFILE_METHOD}'">0){id["\""$1"\";"]=$3}}{
    	if ($12 in id){
        	print $0,"cpc_score \""id[$12]"\";";
        }
    }' $OUTFILE > .tmp
	mv .tmp $OUTFILE

	if [ "$VERBOSITY" -gt 0 ];then
		echo "- Results are stored in : '${METHOD_DIR}' ...">&2
	fi	 
fi




#############################################
# Parse output to filter wrt to METHOD CUTOFF
#############################################
if [ "$VERBOSITY" -gt 0 ];then
	echo "- Parse Output Cod.Pot files...">&2
fi
	
awk -v CPAT_CUT=$CPAT_CUT -v CPC_CUT=$CPC_CUT -v GENEID_CUT=$GENEID_CUT -v TXCDSPREDICT_CUT=$TXCDSPREDICT_CUT -v CPAT=$CPAT -v GENEID=$GENEID -v TXCDSPREDICT=$TXCDSPREDICT -v CPC=$CPC '{
	if (CPAT  && GENEID  && TXCDSPREDICT && CPC ){
		for (j=13;j<=NF;j++){
			bool=0;if ($j == "cpat_score"){   		var=$(j+1);gsub(/\"|;/,"",var);   if (var < CPAT_CUT)          {bool=1}; cpat="cpat_NCstatus \""bool"\";"} 
			bool=0;if ($j == "geneid_score"){		var=$(j+1);gsub(/\"|;/,"",var);   if (var < GENEID_CUT)        {bool=1}; geneid="geneid_NCstatus \""bool"\";"} 
			bool=0;if ($j == "txcdspredict_score"){ var=$(j+1);gsub(/\"|;/,"",var);   if (var < TXCDSPREDICT_CUT) {bool=1}; txcdspredict="txcdspredict_NCstatus \""bool"\";"} 		
			bool=0;if ($j == "cpc_score")         { var=$(j+1);gsub(/\"|;/,"",var);   if (var == "noncoding")     {bool=1}; cpc="cpc_NCstatus \""bool"\";"} 		
			
		}
	}else if (CPAT  && GENEID  && CPC && TXCDSPREDICT=="0" ){
		for (j=13;j<=NF;j++){
			bool=0;if ($j == "cpat_score"){   		var=$(j+1);gsub(/\"|;/,"",var);   if (var < CPAT_CUT)          {bool=1}; cpat="cpat_NCstatus \""bool"\";"} 
			bool=0;if ($j == "geneid_score"){		var=$(j+1);gsub(/\"|;/,"",var);   if (var < GENEID_CUT)        {bool=1}; geneid="geneid_NCstatus \""bool"\";"} 
			bool=0;if ($j == "cpc_score")         { var=$(j+1);gsub(/\"|;/,"",var);   if (var == "noncoding")     {bool=1}; cpc="cpc_NCstatus \""bool"\";"} 		
		}
	}else if (CPAT  && TXCDSPREDICT && CPC && GENEID=="0"){
		for (j=13;j<=NF;j++){
			bool=0;if ($j == "cpat_score"){   		var=$(j+1);gsub(/\"|;/,"",var);   if (var < CPAT_CUT)          {bool=1}; cpat="cpat_NCstatus \""bool"\";"} 
			bool=0;if ($j == "txcdspredict_score"){		var=$(j+1);gsub(/\"|;/,"",var);   if (var < TXCDSPREDICT_CUT)        {bool=1}; txcdspredict="txcdspredict_NCstatus \""bool"\";"} 
			bool=0;if ($j == "cpc_score")         { var=$(j+1);gsub(/\"|;/,"",var);   if (var == "noncoding")     {bool=1}; cpc="cpc_NCstatus \""bool"\";"} 		
		}
	}else if (GENEID  && TXCDSPREDICT && CPC && CPAT=="0"){
		for (j=13;j<=NF;j++){
			bool=0;if ($j == "geneid_score"){		var=$(j+1);gsub(/\"|;/,"",var);   if (var < GENEID_CUT)        {bool=1}; geneid="geneid_NCstatus \""bool"\";"} 
			bool=0;if ($j == "txcdspredict_score"){ var=$(j+1);gsub(/\"|;/,"",var);   if (var < TXCDSPREDICT_CUT) {bool=1}; txcdspredict="txcdspredict_NCstatus \""bool"\";"} 		
			bool=0;if ($j == "cpc_score")         { var=$(j+1);gsub(/\"|;/,"",var);   if (var == "noncoding")     {bool=1}; cpc="cpc_NCstatus \""bool"\";"} 		
		}	
	}else if (GENEID  && TXCDSPREDICT && CPAT && CPC=="0"){
		for (j=13;j<=NF;j++){
			bool=0;if ($j == "geneid_score"){		var=$(j+1);gsub(/\"|;/,"",var);   if (var < GENEID_CUT)        {bool=1}; geneid="geneid_NCstatus \""bool"\";"} 
			bool=0;if ($j == "txcdspredict_score"){ var=$(j+1);gsub(/\"|;/,"",var);   if (var < TXCDSPREDICT_CUT) {bool=1}; txcdspredict="txcdspredict_NCstatus \""bool"\";"} 		
			bool=0;if ($j == "cpat_score"){   		var=$(j+1);gsub(/\"|;/,"",var);   if (var < CPAT_CUT)          {bool=1}; cpat="cpat_NCstatus \""bool"\";"} 
		}	
# 2 methods CPC centered		
	}else if (CPAT && CPC && TXCDSPREDICT=="0" && GENEID=="0"){
		for (j=13;j<=NF;j++){
			bool=0;if ($j == "cpat_score"){   		var=$(j+1);gsub(/\"|;/,"",var);   if (var < CPAT_CUT)          {bool=1}; cpat="cpat_NCstatus \""bool"\";"} 
			bool=0;if ($j == "cpc_score")         { var=$(j+1);gsub(/\"|;/,"",var);   if (var == "noncoding")     {bool=1}; cpc="cpc_NCstatus \""bool"\";"} 		
		}
	}else if (CPAT=="0" && CPC && TXCDSPREDICT=="0"  && GENEID  ){
		for (j=13;j<=NF;j++){
			bool=0;if ($j == "geneid_score"){		var=$(j+1);gsub(/\"|;/,"",var);   if (var < GENEID_CUT)        {bool=1}; geneid="geneid_NCstatus \""bool"\";"} 
			bool=0;if ($j == "cpc_score")         { var=$(j+1);gsub(/\"|;/,"",var);   if (var == "noncoding")     {bool=1}; cpc="cpc_NCstatus \""bool"\";"} 		
		}
	}else if (CPAT =="0" && CPC && TXCDSPREDICT && GENEID=="0"){
		for (j=13;j<=NF;j++){
			bool=0;if ($j == "txcdspredict_score"){ var=$(j+1);gsub(/\"|;/,"",var);   if (var < TXCDSPREDICT_CUT) {bool=1}; txcdspredict="txcdspredict_NCstatus \""bool"\";"} 		
			bool=0;if ($j == "cpc_score")         { var=$(j+1);gsub(/\"|;/,"",var);   if (var == "noncoding")     {bool=1}; cpc="cpc_NCstatus \""bool"\";"} 		
		}
# 2 methods CPAT centered		
	} else if (CPAT && CPC=="0" && TXCDSPREDICT=="0"  && GENEID  ){
		for (j=13;j<=NF;j++){
			bool=0;if ($j == "geneid_score"){		var=$(j+1);gsub(/\"|;/,"",var);   if (var < GENEID_CUT)        {bool=1}; geneid="geneid_NCstatus \""bool"\";"} 
			bool=0;if ($j == "cpat_score"){   		var=$(j+1);gsub(/\"|;/,"",var);   if (var < CPAT_CUT)          {bool=1}; cpat="cpat_NCstatus \""bool"\";"} 
		}
	}else if (CPAT && CPC=="0" && TXCDSPREDICT && GENEID=="0"){
		for (j=13;j<=NF;j++){
			bool=0;if ($j == "txcdspredict_score"){ var=$(j+1);gsub(/\"|;/,"",var);   if (var < TXCDSPREDICT_CUT) {bool=1}; txcdspredict="txcdspredict_NCstatus \""bool"\";"} 		
			bool=0;if ($j == "cpat_score"){   		var=$(j+1);gsub(/\"|;/,"",var);   if (var < CPAT_CUT)          {bool=1}; cpat="cpat_NCstatus \""bool"\";"} 
		}
# 2 methods GENEID centered		
	}else if (CPAT=="0" && CPC=="0" && TXCDSPREDICT && GENEID){
		for (j=13;j<=NF;j++){
			bool=0;if ($j == "geneid_score"){		var=$(j+1);gsub(/\"|;/,"",var);   if (var < GENEID_CUT)        {bool=1}; geneid="geneid_NCstatus \""bool"\";"} 
			bool=0;if ($j == "txcdspredict_score"){ var=$(j+1);gsub(/\"|;/,"",var);   if (var < TXCDSPREDICT_CUT) {bool=1}; txcdspredict="txcdspredict_NCstatus \""bool"\";"} 		
		}
	}
	# print
	print $0,cpat,geneid,txcdspredict, cpc;
}' $OUTFILE >  .tmp
mv .tmp $OUTFILE

#### VennDiagram ##################################################################################################
# Only return filtered lncRNAs of for all method selected wrt to cutoff defined .param
###################################################################################################################
awk -v NUMBER_METHOD=$NUMBER_METHOD '{ 	som=0 ;for (j=13;j<=NF;j++){

	if ($j == "cpat_NCstatus"){if ($(j+1) == "\"1\";"){som+=1}} 
	if ($j == "geneid_NCstatus"){if ($(j+1) == "\"1\";"){som+=1}} 
	if ($j == "txcdspredict_NCstatus"){if ($(j+1) == "\"1\";"){som+=1}} 
	if ($j == "cpc_NCstatus"){if ($(j+1) == "\"1\";"){som+=1}} 
	
	} 
	if (som == NUMBER_METHOD) {print $0}
}' $OUTFILE > ${OUTDIR}/$(basename $OUTFILE .gtf).noncoding.gtf

echo "--> Check : ${OUTDIR}/$OUTFILE : for coding potential score of all your input transcripts">&2
echo "--> Check : ${OUTDIR}/$(basename $OUTFILE .gtf).noncoding.gtf : for transcripts identified as lncRNAs in the $NUMBER_METHOD methods...">&2


#### VennDiagram ##################################################################################################
if [ "$DIAGVENN" ];then

	#No Venn for only One method
	if [ "$NUMBER_METHOD" == "1" ];then
		echo " ERROR: Cannot compute a Venn Diagram for only one method...">&2
		exit 1;	
	
	else
		if [ "$VERBOSITY" -gt 0 ];then
			echo "- Create Venn Diagram...">&2
		fi 
		# Results dir
		if [ -d "${OUTDIR}/$DIAGVENN" ];then		
			echo " Warning: ${OUTDIR}/${DIAGVENN} directory  already exists... ">&2	
		else	
			mkdir ${OUTDIR}/$DIAGVENN;
		fi

		# Format OUTFILE file such as
		# #Ids 	CPC_31490	TxCdsPredict_39355	CPAT_35397
		# TCONS_00206799	1	1	1
		# TCONS_00207213	1	1	1
		awk -v CPAT=$CPAT -v GENEID=$GENEID -v TXCDSPREDICT=$TXCDSPREDICT -v CPC=$CPC '{ 
			id[$12]++; #keep all transcript_id
			for (j=13;j<=NF;j++){
				if ($j == "cpat_NCstatus"){		ar_cpat[$12]=$(j+1); if ($(j+1) == "\"1\";"){ar_cpat_pos[$12]++} } 
				else if ($j == "geneid_NCstatus"){		ar_geneid[$12]=$(j+1); if ($(j+1) == "\"1\";"){ar_geneid_pos[$12]++}} 
				else if ($j == "txcdspredict_NCstatus"){ar_txcdspredict[$12]=$(j+1); if ($(j+1) == "\"1\";"){ar_txcdspredict_pos[$12]++}} 
				else if ($j == "cpc_NCstatus"){         ar_cpc[$12]=$(j+1); if ($(j+1) == "\"1\";"){ar_cpc_pos[$12]++}} 
				
			} 	
		}END{
			# print length array
	if (CPAT  && GENEID  && TXCDSPREDICT && CPC){
		print "#Ids	CPAT_"length(ar_cpat_pos)"	GeneId_"length(ar_geneid_pos)"	TxCdsPredict_"length(ar_txcdspredict_pos)"	CPC_"length(ar_cpc_pos);
		for (i in id){
			print i"	"ar_cpat[i]"	"ar_geneid[i]"	"ar_txcdspredict[i]"	"ar_cpc[i];
		}
	}else if (CPAT  && GENEID  && TXCDSPREDICT && CPC=="0"){
		print "#Ids	CPAT_"length(ar_cpat_pos)"	GeneId_"length(ar_geneid_pos)"	TxCdsPredict_"length(ar_txcdspredict_pos);
		for (i in id){
			print i"	"ar_cpat[i]"	"ar_geneid[i]"	"ar_txcdspredict[i];
		}
	} else if (CPAT  && GENEID  && TXCDSPREDICT=="0" && CPC){
		print "#Ids	CPAT_"length(ar_cpat_pos)"	GeneId_"length(ar_geneid_pos)"	CPC_"length(ar_cpc_pos);
		for (i in id){
			print i"	"ar_cpat[i]"	"ar_geneid[i]"	"ar_cpc[i];
		}
	} else if (CPAT  && GENEID=="0"  && TXCDSPREDICT && CPC){
		print "#Ids	CPAT_"length(ar_cpat_pos)"	CPC_"length(ar_cpc_pos)"	TxCdsPredict_"length(ar_txcdspredict_pos);
		for (i in id){
			print i"	"ar_cpat[i]"	"ar_cpc[i]"	"ar_txcdspredict[i];
		}
	} else if (CPAT=="0"  && GENEID  && TXCDSPREDICT && CPC){
		print "#Ids	CPC_"length(ar_cpc_pos)"	GeneId_"length(ar_geneid_pos)"	TxCdsPredict_"length(ar_txcdspredict_pos);
		for (i in id){
			print i"	"ar_cpc[i]"	"ar_geneid[i]"	"ar_txcdspredict[i];
		}

	}else if (CPAT && CPC && TXCDSPREDICT=="0" && GENEID=="0"){
		print "#Ids	CPAT_"length(ar_cpat_pos)"	CPC_"length(ar_cpc_pos);
		for (i in id){
			print i"	"ar_cpat[i]"	"ar_cpc[i];
		}
	}else if (CPAT=="0" && CPC && TXCDSPREDICT=="0"  && GENEID  ){
		print "#Ids	GeneId_"length(ar_geneid_pos)"	CPC_"length(ar_cpc_pos);
		for (i in id){
			print i"	"ar_geneid[i]"	"ar_cpc[i];
		}
	}else if (CPAT =="0" && CPC && TXCDSPREDICT && GENEID=="0"){
		print "#Ids	TxCdsPredict_"length(ar_txcdspredict_pos)"	CPC_"length(ar_cpc_pos);
		for (i in id){
			print i"	"ar_txcdspredict[i]"	"ar_cpc[i];
		}
# 2 methods CPAT centered		
	} else if (CPAT && CPC=="0" && TXCDSPREDICT=="0"  && GENEID  ){
		print "#Ids	CPAT_"length(ar_cpat_pos)"	GeneId_"length(ar_geneid_pos);
		for (i in id){
			print i"	"ar_cpat[i]"	"ar_geneid[i];
		}
	}else if (CPAT && CPC=="0" && TXCDSPREDICT && GENEID=="0"){
		print "#Ids	CPAT_"length(ar_cpat_pos)"	TxCdsPredict_"length(ar_txcdspredict_pos);
		for (i in id){
			print i"	"ar_cpat[i]"	"ar_txcdspredict[i];
		}
# 2 methods GENEID centered		
	}else if (CPAT=="0" && CPC=="0" && TXCDSPREDICT && GENEID){
		print "#Ids	TxCdsPredict_"length(ar_txcdspredict_pos)"	GeneId_"length(ar_geneid_pos);
		for (i in id){
			print i"	"ar_txcdspredict[i]"	"ar_geneid[i];
		}
	}
}' $OUTFILE | tr -d "\";" > ${OUTDIR}/${DIAGVENN}/${NUMBER_METHOD}_Way_VennDiagram_${OUTFILE_PREF}.txt

	# Create Venne png
	${BIN_DIR}/CreateVenn.sh ${OUTDIR}/${DIAGVENN}/${NUMBER_METHOD}_Way_VennDiagram_${OUTFILE_PREF}.txt
	
	fi # Fin du if du nb de method

	echo "--> Check : ${OUTDIR}/${DIAGVENN}/${NUMBER_METHOD}_Way_VennDiagram_${OUTFILE_PREF}.txt : Venn matrix">&2
	echo "--> Check : ${OUTDIR}/${DIAGVENN}/${NUMBER_METHOD}_Way_VennDiagram_${OUTFILE_PREF}.png : Venn image file">&2
fi


