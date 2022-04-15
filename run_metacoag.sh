#!/bin/bash

set -e               # terminate the script if any command exited with a nonzero exit status
set -u               # terminate the script if any command contains a reference to an unset variable name
set -o pipefail      # any program that returns a nonzero exit status in the pipe will cause the entire pipe to return a nonzero status

####################################################################################
############  PLEASE PROVIDE config.sh AS AN ARGUMENT FOR THIS SCRIPT  #############
####################################################################################

# average running time - 

# Tools
TOOLS_DIR=/Bmo/sochkalova/helper_scripts
RENAME_BINS_SCRIPT=$TOOLS_DIR/rename_bins.py
TSV_CREATION_SCRIPT=$TOOLS_DIR/convert_fasta_bins_to_biobox_format.py
FASTA_BINS_EXTRACTION_SCRIPT=$TOOLS_DIR/extract_fasta_bins.py
WIDE2LONG_SCRIPT=$TOOLS_DIR/wide2long.py
CONVERT_BINNACLE=$TOOLS_DIR/convert_binnacle_scaffolds.py

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
set -x               # print every command with expanded variables
#####################################################################################

metacoag_output_dir=$WORKING_DIR/metacoag
jgi_abundancies=${WORKING_DIR}/${PREFIX}.${assembly_type}.depth
abundance_file=$WORKING_DIR/abundance.${assembly_type}.tsv

if [ ! -e $abundance_file ] 
then
	echo ''
	echo 'Step 1. Create abundance.tsv for MetaCoAG input'
	cut -f 1,3 $jgi_abundancies | tail -n +2 > $abundance_file
else
	echo ''
	echo 'File abundance.tsv exists. Skipping Step 1!'
fi

if [ ! -e $metacoag_output_dir ]
then
	echo ''
	echo 'Step 2. Bin contigs with MetaCoAG'
	~/tools/MetaCoAG/MetaCoAG --assembler SPAdes --graph $ASSEMBLY_GRAPH --abundance $abundance_file \
	--contigs $ASSEMBLY --output $metacoag_output_dir --nthreads $THREADS --paths $ASSEMBLY_PATHS	
else
	echo ''
	echo 'MetaCoAG output directory exists. Skipping Step 2!'
fi

if [ ! -e $metacoag_output_dir/scaffolds2bin.tsv ]
	then
		echo ' '
		echo 'Step 3. Create one-style formatted binning results'
		echo 'Creating binning file scaffolds2bin.tsv...'
		echo "$TSV_CREATION_SCRIPT $metacoag_output_dir/bins/* -o $metacoag_output_dir/scaffolds2bin.tsv"
		$TSV_CREATION_SCRIPT $metacoag_output_dir/bins/* -o $metacoag_output_dir/scaffolds2bin.tsv 

	else
		echo ' '
		echo 'File scaffolds2bin.tsv exists. Skipping Step 3!'
	fi

