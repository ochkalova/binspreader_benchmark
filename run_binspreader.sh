#!/bin/bash

set -e               # terminate the script if any command exited with a nonzero exit status
set -u               # terminate the script if any command contains a reference to an unset variable name
set -o pipefail      # any program that returns a nonzero exit status in the pipe will cause the entire pipe to return a nonzero status

####################################################################################
############  PLEASE PROVIDE config.sh AS AN ARGUMENT FOR THIS SCRIPT  #############
####################################################################################

# Tools
TOOLS_DIR=/Bmo/sochkalova/helper_scripts
RENAME_BINS_SCRIPT=$TOOLS_DIR/rename_bins.py
TSV_CREATION_SCRIPT=$TOOLS_DIR/convert_fasta_bins_to_biobox_format.py
FASTA_BINS_EXTRACTION_SCRIPT=$TOOLS_DIR/extract_fasta_bins.py
WIDE2LONG_SCRIPT=$TOOLS_DIR/wide2long.py

# Varables will be exported from config.sh
config=$1
assembly_type=scaffold

if [ -e $config ]
then
	source $config
else
	echo 'You should provide config.sh as an argument for this script!'
	exit
fi

#####################################################################################


for binner in metabat2 vamb  metawrap # bin3c 
do
	echo ' '
	echo "Currently working with $binner..."
	
	binner_output_dir=${WORKING_DIR}/$binner
	binspreader_output_dir=${WORKING_DIR}/binspreader_$binner

	if [ ! -e $binspreader_output_dir ]
	then
		echo ' '
		echo 'Step 1. Run binning refining using BinSPreader'
		echo "bin-refine $ASSEMBLY_GRAPH $binner_output_dir/scaffolds2bin.tsv $binspreader_output_dir \
		-t $THREADS -Rcorr -Smle" #-e  3.4e-4  
		bin-refine $ASSEMBLY_GRAPH $binner_output_dir/scaffolds2bin.tsv $binspreader_output_dir \
		-t $THREADS -Rcorr -Smle  #-e  #3.4e-4 
		mv $binspreader_output_dir/binning.tsv $binspreader_output_dir/scaffolds2bin.tsv
	else
		echo ' '
		echo 'Seems like binning is done. Skipping Step 1!'
	fi 

	
	if [ ! -e $binspreader_output_dir/bins ]
		then
			echo ' '
			echo 'Step 2. Create bins/ directory with fasta files'
			echo "$FASTA_BINS_EXTRACTION_SCRIPT -b $binspreader_output_dir/scaffolds2bin.tsv -i $ASSEMBLY -o $binspreader_output_dir/bins"
			mkdir $binspreader_output_dir/bins
			$FASTA_BINS_EXTRACTION_SCRIPT -b $binspreader_output_dir/scaffolds2bin.tsv -i $ASSEMBLY -o $binspreader_output_dir/bins
		
		else
			echo ' '
			echo 'bins/ directory exists. Skipping Step 2!'
	fi

#####################################################################################

done

