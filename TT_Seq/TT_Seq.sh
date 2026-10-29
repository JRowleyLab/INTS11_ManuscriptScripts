#!/bin/bash

# Directory for raw FASTQ files
RAW_DIR=...

# Directory for processed files
PROC_DIR=...

# Reference genome
REFERENCE=...

# Adapter sequences
#ADAPTERS=/path/to/your/adapters.fa

# Make the directory structure
mkdir -p ${PROC_DIR}/fastqc
mkdir -p ${PROC_DIR}/fastp

# Loop over each sample
# for SAMPLE in ${RAW_DIR}/*.fastq.gz
# do
#     # Get the base of the file name
#     BASE=$(basename $SAMPLE .fastq.gz)

#     echo Processing sample pair ${BASE}

#     # Perform FastQC analysis
#     #fastqc ${RAW_DIR}/${BASE}_read1.fastq.gz ${RAW_DIR}/${BASE}_read2.fastq.gz -o ${PROC_DIR}/fastqc
#     fastqc ${RAW_DIR}/${BASE}.fastq.gz -o ${PROC_DIR}/fastqc

#     # Trim adapters (piaired end)
#     # fastp \
# 	#    -i ${RAW_DIR}/${BASE}_read1.fastq.gz \
# 	#    -I ${RAW_DIR}/${BASE}_read2.fastq.gz \
# 	#    -o ${PROC_DIR}/${BASE}_1_trimmed.fq.gz \
# 	#    -O ${PROC_DIR}/${BASE}_2_trimmed.fq.gz \
# 	#    --detect_adapter_for_pe \
# 	#    --thread 5 \
# 	#    --length_required 3 \
# 	#    --html ${PROC_DIR}/fastp/${BASE}_fastp_report.html \
# 	#    --json ${PROC_DIR}/fastp/${BASE}_fastp_report.json
#     # trimmomatic SE -phred33 ${RAW_DIR}/${BASE}.fastq ${PROC_DIR}/${BASE}_trimmed.fastq ILLUMINACLIP:${ADAPTERS}:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

#     # single end
#     fastp \
#             -i ${RAW_DIR}/${BASE}.fastq.gz \
#             -o ${PROC_DIR}/${BASE}_trimmed.fq.gz \
#             --detect_adapter_for_pe \
#             --thread 5 \
#             --length_required 3 \
#             --html ${PROC_DIR}/fastp/${BASE}_fastp_report.html \
#             --json ${PROC_DIR}/fastp/${BASE}_fastp_report.json

#     # Align with Bowtie2
#     # bowtie2 -x ${REFERENCE} -1 ${PROC_DIR}/${BASE}_1_trimmed.fq.gz -2 ${PROC_DIR}/${BASE}_2_trimmed.fq.gz -S ${PROC_DIR}/${BASE}.sam -p 10
#     bowtie2 -x ${REFERENCE} -U ${PROC_DIR}/${BASE}_trimmed.fq.gz -S ${PROC_DIR}/${BASE}.sam -p 10

#     # Convert SAM to BAM, sort and index
#     samtools view -S -b ${PROC_DIR}/${BASE}.sam > ${PROC_DIR}/${BASE}.bam
#     samtools view -b -q 30 ${PROC_DIR}/${BASE}.bam > ${PROC_DIR}/${BASE}_filtered.bam
#     samtools sort ${PROC_DIR}/${BASE}_filtered.bam -o ${PROC_DIR}/${BASE}_sorted.bam
#     samtools index ${PROC_DIR}/${BASE}_sorted.bam

#     # bamCoverage -b ${PROC_DIR}/${BASE}_sorted.bam -o ${PROC_DIR}/${BASE}_ForwardStrand.bigwig --normalizeUsing BPM --binSize 10 --filterRNAstrand forward
#     # bamCoverage -b ${PROC_DIR}/${BASE}_sorted.bam -o ${PROC_DIR}/${BASE}_ReverseStrand.bigwig --normalizeUsing BPM --binSize 10 --filterRNAstrand reverse
#     # bamCoverage -b ${PROC_DIR}/${BASE}_sorted.bam -o ${PROC_DIR}/${BASE}.bigwig --normalizeUsing BPM --binSize 10

#     # Count/quality control of reads with HTSeq
#     #htseq-count -f bam -s no ${PROC_DIR}/${BASE}_sorted.bam /home/genomefiles/human/hg38/hg38_ensemble.sorted.gtf > ${PROC_DIR}/${BASE}_counts.txt
# done

samtools merge ${PROC_DIR}/All_Samples_merged.bam ${PROC_DIR}/*_sorted.bam
samtools sort ${PROC_DIR}/All_Samples_merged.bam -o ${PROC_DIR}/All_Samples_sorted.bam
samtools index ${PROC_DIR}/All_Samples_sorted.bam
bamCoverage -b ${PROC_DIR}/All_Samples_sorted.bam -o ${PROC_DIR}/All_Samples_ForwardStrand.bigwig --normalizeUsing BPM --binSize 10 --filterRNAstrand forward
bamCoverage -b ${PROC_DIR}/All_Samples_sorted.bam -o ${PROC_DIR}/All_Samples_ReverseStrand.bigwig --normalizeUsing BPM --binSize 10 --filterRNAstrand reverse
bamCoverage -b ${PROC_DIR}/All_Samples_sorted.bam -o ${PROC_DIR}/All_Samples.bigwig --normalizeUsing BPM --binSize 10


# At this point, you have to use R/Bioconductor package to perform data normalization and generate a BEDGraph

echo "All Done!"
