#!/bin/bash

set -e               # terminate the script if any command exited with a nonzero exit status
set -u               # terminate the script if any command contains a reference to an unset variable name
set -o pipefail      # any program that returns a nonzero exit status in the pipe will cause the entire pipe to return a nonzero status

####################################################################################
############  PLEASE PROVIDE config.sh AS AN ARGUMENT FOR THIS SCRIPT  #############
####################################################################################

# average running time - 6-9h (metacarvel) + 2-3h (binnacle)

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

# Parameters
minimap2_params='-ax sr'
bam_file=${WORKING_DIR}/${PREFIX}.${assembly_type}.bam

if [ ! -e $bam_file ]
then
	echo ' '
	echo 'Step 1. Align reads to assembly, filter and sort bam file'
	echo "minimap2 -t $THREADS $minimap2_params $ASSEMBLY $FORWARD_READS $REVERSE_READS | \
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
binnacle_output_dir=${WORKING_DIR}/binnacle_metabat2
metacarvel_output_dir=${binnacle_output_dir}/metacarvel

if [ ! -e $metacarvel_output_dir ]
then
	echo ' '
	echo 'Step 2. Make graph from alignment using MetaCarvel'
	echo "python ~/tools/MetaCarvel/run.py -a $ASSEMBLY -m $bam_file -d $metacarvel_output_dir -k true -r true"
	python ~/tools/MetaCarvel/run.py -a $ASSEMBLY -m $bam_file -d $metacarvel_output_dir -k true -r true 
	
else
	echo ' '
	echo 'MetaCarvel output directory exists. Skipping Step 2!'
fi

#####################################################################################

# Parameters


if [ ! -e $binnacle_output_dir/Feature-Matrix-metabat.txt ]
then
	echo ' '
	echo 'Step 3. Estimate scaffold coverage using Binnacle'
	echo "python ~/tools/binnacle/src/Estimate_Abundances.py -g ${metacarvel_output_dir}/oriented.gml \
	-bam $bam_file -c $ASSEMBLY -d $binnacle_output_dir"
	python ~/tools/binnacle/src/Estimate_Abundances.py -g ${metacarvel_output_dir}/oriented.gml \
	-bam $bam_file -c $ASSEMBLY -d $binnacle_output_dir

	echo "python ~/tools/binnacle/src/Collate.py -d $binnacle_output_dir -m metabat"
	python ~/tools/binnacle/src/Collate.py -d $binnacle_output_dir -m metabat
	
else
	echo ' '
	echo 'Binnacle output directory exists. Skipping Step 3!'
fi

#####################################################################################

# Parameters
metabat_params="--seed 42"
metabat_output_dir=${binnacle_output_dir}/bins

if [ ! -e $metabat_output_dir ]
then
	echo ' '
	echo 'Step 4. Use abundances estimated by Binnacle for binning with MetaBAT2'
	echo "metabat2 -i ${binnacle_output_dir}/Scaffolds.fasta -o $metabat_output_dir/bin -a $binnacle_output_dir/Feature-Matrix-metabat.txt \
	--numThreads $THREADS $metabat_params"
	metabat2 -i ${binnacle_output_dir}/Scaffolds.fasta -o $metabat_output_dir/bin -a $binnacle_output_dir/Feature-Matrix-metabat.txt \
	--numThreads $THREADS $metabat_params

else
	echo ' '
	echo 'Seems like binning is done. Skipping Step 4!'
fi

#####################################################################################

if [ ! -e $metabat_output_dir/scaffolds2bin.tsv ]
then
	echo ' '
	echo 'Step 5. Create one-style formatted binning results'
	
	rename -f 's/bin\./bin_/' $metabat_output_dir/bin.*.fa 
	rename -f 's/\.fa/.fasta/' $metabat_output_dir/bin_*.fa

	$TSV_CREATION_SCRIPT $metabat_output_dir/* -o $metabat_output_dir/binnacle_binning.tsv
	$CONVERT_BINNACLE --coordinates ${binnacle_output_dir}/Coords_After_Delinking.txt \
	--binning $metabat_output_dir/binnacle_binning.tsv -o ${binnacle_output_dir}/scaffolds2bin.tsv

else
	echo ' '
	echo 'File scaffolds2bin.tsv exists. Skipping Step 5!'
fi