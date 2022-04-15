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
metabat_params="--seed 42"
metabat_output_dir=${WORKING_DIR}/metabat2
jgi_abundancies=${WORKING_DIR}/${PREFIX}.depth

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

#####################################################################################


# code to run vamb will be here


#####################################################################################

for binner in metabat2 vamb # metawrap bin3c
do
	echo ' '
	echo "Currently working with $binner..."
	
	binner_output_dir=${WORKING_DIR}/$binner/bins
	metamvgl_output_dir=${WORKING_DIR}/metamvgl_$binner

	if [ ! -e $binner/bins ]
	then
		echo ' '
		echo 'Step 3.5. Create bins/ directory with fasta files for initial binning'
		echo "$FASTA_BINS_EXTRACTION_SCRIPT -b $binner/scaffolds2bin.tsv -i $ASSEMBLY -o $binner/bins"
		mkdir $binner/bins
		$FASTA_BINS_EXTRACTION_SCRIPT -b $binner/scaffolds2bin.tsv -i $ASSEMBLY -o $binner/bins
	else
		echo ' '
		echo 'bins/ directory exists. Skipping Step 3.5!'
	fi

	if [ ! -e $metamvgl_output_dir/initial_contig_bins.csv ]
	then
		echo ' '
		echo 'Step 4. Convert initial binning to sutable form for METAMVGL'
		echo "python /Bmo/akorobeynikov/METAMVGL/prepResult.py --binned $binner_output_dir \
		--assembler SPAdes --output $metamvgl_output_dir"
		python /Bmo/akorobeynikov/METAMVGL/prepResult.py --binned $binner_output_dir \
		--assembler SPAdes --output $metamvgl_output_dir
	else
		echo ' '
		echo 'File initial_contig_bins.csv exists. Skipping Step 4!'
	fi 

	#####################################################################################

	pe_graph=$WORKING_DIR/$PREFIX.pe
	ag_graph=$WORKING_DIR/$PREFIX.ag

	if [ ! -e $pe_graph ] && [ ! -e $ag_graph ]
	then 
		echo ' '
		echo 'Step 5. Compute pair-end and assembly graph with METAMVGL'
		echo "/Bmo/akorobeynikov/METAMVGL/prep_graph --assembly-graph $ASSEMBLY_GRAPH_FASTG \
		--bam $bam_file --assembler metaSPAdes --output $WORKING_DIR/$PREFIX --paths $ASSEMBLY_PATHS"
		/Bmo/akorobeynikov/METAMVGL/prep_graph --assembly-graph $ASSEMBLY_GRAPH_FASTG \
		--bam $bam_file --assembler metaSPAdes --output $WORKING_DIR/$PREFIX --paths $ASSEMBLY_PATHS
	else
		echo ' '
		echo 'Seems like graphs are already present. Skipping Step 5!'
	fi

	#####################################################################################

	if [ ! -e $metamvgl_output_dir/bins ]
	then
		echo ' '
		echo 'Step 6. Bin contigs using two types of graphs using METAMVGL'
		echo "python /Bmo/akorobeynikov/METAMVGL/METAMVGL.py --contigs $ASSEMBLY --assembler metaSPAdes \
		--assembly_graph $ag_graph --PE_graph $pe_graph \
		--binned $metamvgl_output_dir/initial_contig_bins.csv --output $metamvgl_output_dir/bins"
		python /Bmo/akorobeynikov/METAMVGL/METAMVGL.py --contigs $ASSEMBLY --assembler metaSPAdes \
		--assembly_graph $ag_graph --PE_graph $pe_graph \
		--binned $metamvgl_output_dir/initial_contig_bins.csv --output $metamvgl_output_dir/bins
	else
		echo ' '
		echo 'METAMVGL already binned the contigs. Skipping Step 6!'
	fi

	#####################################################################################

	if [ ! -e $metamvgl_output_dir/scaffolds2bin.tsv ]
	then
		echo ' '
		echo 'Step 7. Create one-style formatted binning results'

		echo 'Renaming cluster.*.fasta to bin_*.fasta...'
		rename -f 's/cluster\./bin_/' $metamvgl_output_dir/bins/cluster.*.fasta 

		echo 'Creating binning file scaffolds2bin.tsv...'
		echo "$TSV_CREATION_SCRIPT $metamvgl_output_dir/bins/* -o $metamvgl_output_dir/bins/scaffolds2bin.tsv"
		$TSV_CREATION_SCRIPT $metamvgl_output_dir/bins/* -o $metamvgl_output_dir/bins/scaffolds2bin.tsv

		echo 'Retrieving the original contig names...'

		# take original header from assembly fasta file
		for replace_string in $(awk 'sub(/^>/, "")' $ASSEMBLY)
		do 
			# strip original header to METAMVGL short name
		    search_string=$(echo $replace_string | cut -d "_" -f 1,2 --output-delimiter='_')

		    # if this contig was binned
		    if grep -Pq "$search_string\t" $metamvgl_output_dir/bins/scaffolds2bin.tsv
		    then
		    	# find its bin file using scaffolds2bin.tsv
		    	bin_id=$(grep -P "$search_string\t" $metamvgl_output_dir/bins/scaffolds2bin.tsv | cut -f 2)
		    	file=$metamvgl_output_dir/bins/${bin_id}.fasta
		    	# replace shortend contig name to original
		    	sed -i "s/$search_string$/$replace_string/" $file
		    fi
		done

		echo 'Creating new binning file scaffolds2bin.tsv with original contig names...'
		rm $metamvgl_output_dir/bins/scaffolds2bin.tsv
		echo "$TSV_CREATION_SCRIPT $metamvgl_output_dir/bins/* -o $metamvgl_output_dir/scaffolds2bin.tsv"
		$TSV_CREATION_SCRIPT $metamvgl_output_dir/bins/* -o $metamvgl_output_dir/scaffolds2bin.tsv

	else
		echo ' '
		echo 'File scaffolds2bin.tsv exists. Skipping Step 7!'
	fi
done