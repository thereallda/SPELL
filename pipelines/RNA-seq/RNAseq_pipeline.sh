#!/usr/bin/bash

help() {
		echo ""
		echo -e "Usage: \n\t bash $0 [options] -d <project_dir> -i <input_csv> --ref <reference_genome> --gtf <GTF_flie>"
		echo ""
		echo -e "  -d \n\t\t Project directory."
		echo -e "  -i \n\t\t Input csv file"
		echo -e "  --ref \n\t\t Reference for STAR mapping."
		echo -e "  --gtf \n\t\t GTF for gene quantification."
		echo ""
		echo "Options:"
		echo -e "  -t, --thread \n\t\t Numbers of threads, default: 1"
		echo -e "  --pair \n\t\t Flag to turn on pair-end mode."
		echo -e "  --syn \n\t\t Using synthetic spike-ins reference 1 or 2."
		echo -e "  -h, --help \n\t\t show help message."
		echo ""
}

if [ $# -eq 0 ]; then
        echo "Run 'bash $0 -h' to see more information."
        exit 0
fi

# NOTE: This requires GNU getopt.  On Mac OS X and FreeBSD, you have to install this
# separately; see below.
TEMP=`getopt -o d:t:i:h --long help,ref:,threads:,pair,gtf:,syn: \
             -n 'RNAseq_pipeline.sh' -- "$@"`

if [ $? != 0 ] ; then echo "Terminating..." >&2 ; exit 1 ; fi

# Note the quotes around `$TEMP': they are essential!
eval set -- "$TEMP"

THREADS=1
PAIR=false
FQ1=
FQ2=
INPUTCSV=
SYNTHETIC=
# Illumina universial adapters
#FORADAPT=AGATCGGAAGAG
#READAPT=AGATCGGAAGAG

while true; do
  case "$1" in
				-d ) PROJECTDIR="$2"; shift 2;;
				-i ) INPUTCSV="$2"; shift 2;;
				-t | --threads ) THREADS="$2"; shift 2;;
				--pair ) PAIR=true; shift ;;
				--ref ) REFERENCE="$2"; shift 2;;
				--gtf ) GTF="$2"; shift 2;;
				--syn ) SYNTHETIC="$2"; shift 2;;
				-h | --help)
                        help
                        exit 0;;
                -- ) shift; break ;;
                * ) echo "Invalid option: ${optionName}" ;
				   echo "Try 'bash `basename $0` -h' for more information" ; 
				   break ;;
        esac
done

if [[ ! -d ${PROJECTDIR}/data ]] || [[ ! -d ${PROJECTDIR}/results ]]; then
	echo " The 'data' and 'results' folders must be in your project directory"
	exit 1
fi

if [[ -z ${INPUTCSV} ]]; then
   echo "Input csv file not specified"
   exit 1
fi

if [[ $REFERENCE = "" ]]; then
	echo " The reference genome is needed "
	exit 1
fi

# Make sure the EOL of input csv is "\n" and the last line have line ending
# Add \n to each line in CSV file
sed 's/\r$//g' $INPUTCSV > tempfile.csv
sed -i -e '$a\' tempfile.csv
mv tempfile.csv $INPUTCSV
rm -f tempfile.csv

## mkdir
# data directory
DATADIR=${PROJECTDIR}/data
if [[ ! -d ${DATADIR}/clean ]]; then mkdir -p ${DATADIR}/clean; fi
FASTQLOC=${DATADIR}/fastq
CLEANLOC=${DATADIR}/clean

# output directory
OUTDIR=${PROJECTDIR}/results
if [[ ! -d ${OUTDIR}/align ]]; then mkdir -p ${OUTDIR}/align ; fi
ALIGNLOC=${OUTDIR}/align
if [[ ! -d ${OUTDIR}/featurecounts ]]; then mkdir -p ${OUTDIR}/featurecounts ; fi
COUNTLOC=${OUTDIR}/featurecounts
if [[ ! -d ${OUTDIR}/QC ]]; then mkdir -p ${OUTDIR}/QC; fi

## Printing 
echo ""
echo " Project directory: ${PROJECTDIR} "
echo " Fastq directory: ${FASTQLOC}"
echo " Reference path: ${REFERENCE}"
echo " Numbers of threads: ${THREADS} "
echo ""

# Read the csv file
while IFS=',' read -r sample fastq1 fastq2; do
    # Create a prefix variable
    PREFIX=${sample}
    if $PAIR; then
	    # Pair-end mode
		FQ1=$(echo "${FASTQLOC}/${fastq1}" | sed -e 's/\r\?$//') # file path of fastq1
		FQ2=$(echo "${FASTQLOC}/${fastq2}" | sed -e 's/\r\?$//') # file path of fastq2
		echo "Processing "${fastq1} "&" ${fastq2}

		## Trimming -- ${DATADIR}/clean/
		echo ""
		echo ">Step01 Trimming & QC"
		echo ""
		trim_galore --nextseq 30 \
		--phred33 \
		--gzip \
		-o ${CLEANLOC}/ \
		--cores ${THREADS} \
		--basename ${PREFIX} \
		--fastqc_args "-o ${OUTDIR}/QC -t ${THREADS}" \
		--paired ${FQ1} ${FQ2}

		## Mapping -- ${OUTDIR}/align/${PREFIX}_align/
		echo ""
		echo ">Step02 Mapping" 
		echo ""
		if [[ ! -d ${ALIGNLOC}/${PREFIX}_align ]]; then mkdir -p ${ALIGNLOC}/${PREFIX}_align; fi
		STAR --genomeDir ${REFERENCE} \
		--readFilesCommand zcat \
		--readFilesIn ${CLEANLOC}/${PREFIX}_val_1.fq.gz ${CLEANLOC}/${PREFIX}_val_2.fq.gz \
		--runThreadN ${THREADS} \
		--outFileNamePrefix ${ALIGNLOC}/${PREFIX}_align/ 

		## samtools -- ${OUTDIR}/align/${PREFIX}_align/
		echo ""
		echo ">Step03 SAMTOOLS" 
		echo ""
		samtools view -@ ${THREADS} -q 30 -f 2 -hSb ${ALIGNLOC}/${PREFIX}_align/Aligned.out.sam |samtools sort - -@ ${THREADS} -o ${ALIGNLOC}/${PREFIX}_align/${PREFIX}.sorted.bam
		samtools index ${ALIGNLOC}/${PREFIX}_align/${PREFIX}.sorted.bam
		rm -f ${ALIGNLOC}/${PREFIX}_align/Aligned.out.sam
		
		## featureCounts -- ${OUTDIR}/featurecounts/
		echo ""
		echo ">Step04 featureCounts" 
		echo ""
		featureCounts \
		-p -B -C \
		-t exon \
		-g gene_id \
		-T ${THREADS} \
		-a $GTF \
		-o ${COUNTLOC}/${PREFIX}_counts.txt \
		${ALIGNLOC}/${PREFIX}_align/${PREFIX}.sorted.bam
	
	else
	    # Single-end mode
		FQ1=$(echo "${FASTQLOC}/${fastq1}" | sed -e 's/\r\?$//') # file path of fastq1
		echo "Processing " ${fastq1}
        
		## Trimming -- ${DATADIR}/clean/
		echo ""
		echo ">Step01 Trimming & QC"
		echo ""
		trim_galore --nextseq 30 \
		--phred33 \
		--gzip \
		-o ${CLEANLOC}/ \
		--cores ${THREADS} \
		--basename ${PREFIX} \
		--fastqc_args "-o ${OUTDIR}/QC -t ${THREADS}" \
		${FQ1} 

		## Mapping -- ${OUTDIR}/align/${PREFIX}_align/
		echo ""
		echo ">Step02 Mapping" 
		echo ""
		if [[ ! -d ${ALIGNLOC}/${PREFIX}_align ]]; then mkdir -p ${ALIGNLOC}/${PREFIX}_align; fi
		STAR --genomeDir ${REFERENCE} \
		--readFilesCommand zcat \
		--readFilesIn ${CLEANLOC}/${PREFIX}_trimmed.fq.gz \
		--runThreadN ${THREADS} \
		--outFileNamePrefix ${ALIGNLOC}/${PREFIX}_align/  

		## samtools -- ${OUTDIR}/align/${PREFIX}_align/
		echo ""
		echo ">Step03 SAMTOOLS" 
		echo ""
		samtools view -@ ${THREADS} -q 30 -F 4 -hSb ${ALIGNLOC}/${PREFIX}_align/Aligned.out.sam |samtools sort - -@ ${THREADS} -o ${ALIGNLOC}/${PREFIX}_align/${PREFIX}.sorted.bam
		samtools index ${ALIGNLOC}/${PREFIX}_align/${PREFIX}.sorted.bam
		rm -f ${ALIGNLOC}/${PREFIX}_align/Aligned.out.sam

		## featureCounts -- ${OUTDIR}/featurecounts/
		echo ""
		echo ">Step04 featureCounts" 
		echo ""
		featureCounts \
		-t exon \
		-g gene_id \
		-T ${THREADS} \
		-a $GTF \
		-o ${COUNTLOC}/${PREFIX}_counts.txt \
		${ALIGNLOC}/${PREFIX}_align/${PREFIX}.sorted.bam

	fi
done < <(tail -n +2 "$INPUTCSV")

# post multiqc
## trimmed reads qc 
multiqc -o ${OUTDIR}/QC/  ${OUTDIR}/QC/*zip
## star alignment reports
multiqc -o ${OUTDIR}/align/  ${OUTDIR}/align/*/Log.final.out

# merge output featureCount files into single csv file
mergecounts=/public/publicUse/script/RNA-seq_downstream/R/mergeCounts.r
Rscript ${mergecounts} --wd ${PROJECTDIR}/

# perform synthetic spike-ins analysis if --syn option provided with 1 or 2
if [[ "$SYNTHETIC" == "1" ]] || [[ "$SYNTHETIC" == "2" ]]; then
    synthetic=/public/publicUse/script/SyntheticAll.sh
	bash ${synthetic} --syn ${SYNTHETIC} -t ${THREADS} 
fi
