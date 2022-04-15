#!/bin/bash

set -e               # terminate the script if any command exited with a nonzero exit status
set -u               # terminate the script if any command contains a reference to an unset variable name
set -o pipefail      # any program that returns a nonzero exit status in the pipe will cause the entire pipe to return a nonzero status

####################################################################################
############  PLEASE PROVIDE config.sh AS AN ARGUMENT FOR THIS SCRIPT  #############
####################################################################################

# average running time - 4-5h (without read mapping)

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

# Parameters
minimap2_params='-ax sr'
bam_file=${WORKING_DIR}/${PREFIX}.${assembly_type}.bam

if [ ! -e $bam_file ]
then
	echo ' '
	echo 'Step 1. Align reads to assembly, filter and sort bam file'

	echo "minimap2 -t $THREADS -ax sr $ASSEMBLY $FORWARD_READS $REVERSE_READS | \
	samtools view -F 3584 -b -o - | \
	samtools sort -@ $THREADS -o $bam_file"

	minimap2 -t $THREADS $minimap2_params $ASSEMBLY $FORWARD_READS $REVERSE_READS | \
	samtools view -F 3584 -b -o - | \
	samtools sort -@ $THREADS -o $bam_file
else
	echo ' '
	echo 'Reads are already aligned to assembly. Skipping Step 1!'
fi

#####################################################################################

# Parameters
metabat_params="--seed 42"
metabat_output_dir=${WORKING_DIR}/metabat2
jgi_abundancies=${WORKING_DIR}/${PREFIX}.${assembly_type}.depth

if [ ! -e $metabat_output_dir ]
then
	echo ' '
	echo 'Step 2. Make initial binning with MetaBAT2'
	
	if [ ! -e $jgi_abundancies ]
	then
		echo "jgi_summarize_bam_contig_depths --outputDepth $jgi_abundancies $bam_file "
		jgi_summarize_bam_contig_depths --outputDepth $jgi_abundancies $bam_file
	fi
	echo "metabat -i $ASSEMBLY -o ${metabat_output_dir}/bins/bin -a $jgi_abundancies --numThreads $THREADS $metabat_params"
	mkdir -p ${metabat_output_dir}
	metabat -i $ASSEMBLY -o ${metabat_output_dir}/bins/bin -a $jgi_abundancies --numThreads $THREADS $metabat_params
	
else
	echo ' '
	echo 'Seems like binning is done. Skipping Step 2!'
fi

if [ ! -e $metabat_output_dir/scaffolds2bin.tsv ]
	then
		echo ' '
		echo 'Step 3. Create one-style formatted binning results'
		echo 'Creating binning file scaffolds2bin.tsv...'
		rename -f 's/bin\./bin_/' $metabat_output_dir/bins/bin.*.fa 
		rename -f 's/\.fa/.fasta/' $metabat_output_dir/bins/bin_*.fa
		echo "$TSV_CREATION_SCRIPT $metabat_output_dir/bins/* -o $metabat_output_dir/scaffolds2bin.tsv"
		$TSV_CREATION_SCRIPT $metabat_output_dir/bins/* -o $metabat_output_dir/scaffolds2bin.tsv 

	else
		echo ' '
		echo 'File scaffolds2bin.tsv exists. Skipping Step 3!'
	fi
