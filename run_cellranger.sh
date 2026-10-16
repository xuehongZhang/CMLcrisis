#!/usr/bin/env bash

# ============================================================
###### Run Cell Ranger (v7.1.0) to generate the single-cell gene expression count matrix, which was then used for downstream quality control and clustering analysis  #######
# ============================================================

####Processing #######

me=whoami

# ---------- Input files ----------

# Input: directory for this project
ProjectDir="/public/$me/projects/CMLblast";

# Input: directory containing paired-end FASTQ files
rawDataDir=${ProjectDir}'/rawData';


# Input: Cell Ranger transcriptome reference (contains the index)
# This is the reference package, e.g., refdata-gex-GRCh38-2020-A
Index="/public/$me/utilities/reference/singlecell/refdata-gex-GRCh38-2020-A/";

# ---------- Output files ----------

# Output: parent directory for Cell Ranger results
outDir=${ProjectDir}'/results/cellranger';


# if the Cell Ranger Output Dir does not exist, then create it #
if [ ! -d "${outDir}" ]; then mkdir ${outDir}; fi


Allfiles=`ls ${rawDataDir}`
Sub="";
for f in ${Allfiles};do
## Get the sequences that needed to be mapped #
	SUBSTRING=`echo ${f} | cut -d '_' -f 1`;
	result=$(echo $Sub | grep "${SUBSTRING}")
	if [[ "$result" == "" ]]; then Sub=${Sub}" "${SUBSTRING};
	fi
done

num=0;
for f in ${Sub};do
	echo $f;
	let "num=$num +1 "
done

echo "The total numer of seq file = "$num;


	# Using cellranger-7.1.0 to map the sequence #	
for f in ${Sub};do
		#### step 1: STAR mapping ###
	cd ${outDir}
	/public/zxh/utilities/cellranger-7.1.0/bin/cellranger count --id ${f} --transcriptome=${Index} --create-bam true --fastqs=${rawDataDir} --sample=${f} --nosecondary --localcores=20 --localmem=256
		  
done
######################################


