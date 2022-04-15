#!/bin/bash

set -e               # terminate the script if any command exited with a nonzero exit status
set -u               # terminate the script if any command contains a reference to an unset variable name
set -o pipefail      # any program that returns a nonzero exit status in the pipe will cause the entire pipe to return a nonzero status

####################################################################################
############  PLEASE PROVIDE config.sh AS AN ARGUMENT FOR THIS SCRIPT  #############
####################################################################################

# average running time - 5-10 min

# Tools
TOOLS_DIR=/Bmo/sochkalova/helper_scripts
RENAME_BINS_SCRIPT=$TOOLS_DIR/rename_bins.py
TSV_CREATION_SCRIPT=$TOOLS_DIR/convert_fasta_bins_to_biobox_format.py
FASTA_BINS_EXTRACTION_SCRIPT=$TOOLS_DIR/extract_fasta_bins.py
WIDE2LONG_SCRIPT=$TOOLS_DIR/wide2long.py

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

for binner in metabat2  vamb  metawrap # bin3c
do
	echo ' '
	echo "Currently working with $binner..."
	
	binner_output_dir=${WORKING_DIR}/$binner
	das_tool_output_dir=${WORKING_DIR}/das_tool_$binner

	if [ ! -e $das_tool_output_dir ]
	then
		echo ' '
		echo 'Step 1. Run DAS_Tool on generated bins'
		echo ""
		mkdir $das_tool_output_dir
		grep -v '^@' $binner_output_dir/scaffolds2bin.tsv > $das_tool_output_dir/initial_binning.tsv
		~/tools/DAS_Tool-1.1.3/DAS_Tool -i $das_tool_output_dir/initial_binning.tsv -l $binner \
		-o $das_tool_output_dir/$PREFIX --search_engine diamond -t $THREADS -c $ASSEMBLY --write_bins
		mv $das_tool_output_dir/${PREFIX}_DASTool_bins/ $das_tool_output_dir/bins
	else
		echo ' '
		echo 'DAS Tool output dir exists. Skipping Step 1!'
	fi 

	#####################################################################################

	if [ ! -e $das_tool_output_dir/scaffolds2bin.tsv ]
	then
		echo ' '
		echo 'Step 2. Create one-style formatted binning results'
		echo 'Creating binning file scaffolds2bin.tsv...'
		echo "$TSV_CREATION_SCRIPT $das_tool_output_dir/bins/* -o $das_tool_output_dir/scaffolds2bin.tsv"
		$TSV_CREATION_SCRIPT $das_tool_output_dir/bins/* -o $das_tool_output_dir/scaffolds2bin.tsv 

	else
		echo ' '
		echo 'File scaffolds2bin.tsv exists. Skipping Step 2!'
	fi
done

