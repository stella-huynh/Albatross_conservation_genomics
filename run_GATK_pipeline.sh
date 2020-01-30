#!/usr/bin/bash

module load bwa/0.7.13
module load samtools/1.4

machDIR1="/group/sbs_ssin/anaconda3/share/picard-2.20.6-0"
machDIR2="/group/sbs_ssin/applications/GenomeAnalysisTK-3.8-1-0"
JAVA8="/group/sbs_ssin/applications/jdk1.8.0_221/bin"

list=$1

for SAMPLE in `cat ${list}`
do
  for REP in 1st 2nd 3rd 4th
  do
  
    # Trim PE reads
    NGmerge -1 reseq/${SAMPLE}_${REP}_R1.fastq.gz -2 reseq/${SAMPLE}_${REP}_R2.fastq.gz \
			      -a -v -n ${NTHREADS} -o reseq/${SAMPLE}_${REP}.trimmed

    # Convert FastQ file to unaligned SAM file
    java -Xmx8g -XX:ParallelGCThreads=1 -jar ${machDIR1}/picard.jar FastqToSam \
	        FASTQ=reseq/${SAMPLE}_${REP}.trimmed_1.fastq.gz \
	        FASTQ2=reseq/${SAMPLE}_${REP}.trimmed_2.fastq.gz \
	        OUTPUT=reseq/${SAMPLE}_${REP}.fastqtosam.bam \
	        READ_GROUP_NAME=${FLOWCELL}.${LANE} \
	        SAMPLE_NAME=${SAMPLE} \
	        LIBRARY_NAME=${SAMPLE} \
        	PLATFORM_UNIT=${FLOWCELL}.${LANE}.${SAMPLE} \
        	PLATFORM=ILLUMINA \
        	SEQUENCING_CENTER=BAUER \
        	RUN_DATE=${DATE}

    # Rewrite SAM file with new adapter-trimming tags
    java -Xmx8g -XX:ParallelGCThreads=1 -jar ${machDIR1}/picard.jar MarkIlluminaAdapters \
	        I=reseq/${SAMPLE}_${REP}.fastqtosam.bam \
	        O=reseq/${SAMPLE}_${REP}_markilluminaadapters.bam \
	        M=reseq/${SAMPLE}_${REP}_markilluminaadapters_metrics.txt \
	        TMP_DIR=./tmp

    # Convert SAM file to FastQ.
    java -Xmx8g -XX:ParallelGCThreads=4 -jar ${machDIR1}/picard.jar SamToFastq \
	        I=reseq/${SAMPLE}_${REP}_markilluminaadapters.bam \
	        FASTQ=reseq/${SAMPLE}_${REP}_samtofastq_interleaved.fq \
	        CLIPPING_ATTRIBUTE=XT CLIPPING_ACTION=2 INTERLEAVE=true NON_PF=true \
	        TMP_DIR=./tmp

    # Align raw reads onto REF genome using BWA
    bwa mem -M -t 1 -p WGS/BFAL_genome.fasta \
		                reseq/${SAMPLE}_${REP}_samtofastq_interleaved.fq \
		                > reseq/${SAMPLE}_${REP}_bwa_mem.sam

    # Merge info from aligned SAM file and unmapped BAM file in a new BAM file.
    java -Xmx8g -XX:ParallelGCThreads=4 -jar ${machDIR1}/picard.jar MergeBamAlignment \
          ALIGNED_BAM=reseq/${SAMPLE}_${REP}_bwa_mem.sam \
	        UNMAPPED_BAM=reseq/${SAMPLE}_${REP}.fastqtosam.bam \
          OUTPUT=reseq/${SAMPLE}_${REP}_mergebamalign.bam \
          R=WGS/BFAL_genome.fasta \
	        CREATE_INDEX=true \
          ADD_MATE_CIGAR=true \
          CLIP_ADAPTERS=false \
          CLIP_OVERLAPPING_READS=true \
          INCLUDE_SECONDARY_ALIGNMENTS=true \
          MAX_INSERTIONS_OR_DELETIONS=-1 \
          PRIMARY_ALIGNMENT_STRATEGY=MostDistant \
          ATTRIBUTES_TO_RETAIN=XS \
          TMP_DIR=./tmp
  done

done


for SAMPLE in `cat list_reseq_BFAL_LAAL`
do
    # Mark duplicates
    java -Xmx8g -XX:ParallelGCThreads=1 -jar ${machDIR1}/picard.jar MarkDuplicates
          TMP_DIR=${oDIR}/tmp \
    		  I=reseq/${SAMPLE}_1st_mergebamalign.bam \
    		  I=reseq/${SAMPLE}_2nd_mergebamalign.bam \
    	  	I=reseq/${SAMPLE}_3rd_mergebamalign.bam \
    	  	O=reseq/${SAMPLE}.dedup.bam \
    	  	METRICS_FILE=reseq/${SAMPLE}.dedup.metrics.txt \
      		REMOVE_DUPLICATES=false TAGGING_POLICY=All

    # Index BAM file
    java -Xmx8g -XX:ParallelGCThreads=2 -jar ${machDIR1}/picard.jar BuildBamIndex \
          I=reseq/${SAMPLE}.dedup.sorted.bam

    # Generate summary of alignment metrics from BAM files
    java -Xmx8g -XX:ParallelGCThreads=2 -jar ${machDIR1}/picard.jar CollectAlignmentSummaryMetrics \
          I=reseq/${SAMPLE}.dedup.sorted.bam \
        	R=WGS/BFAL_genome.fasta \
	        METRIC_ACCUMULATION_LEVEL=SAMPLE \
	        METRIC_ACCUMULATION_LEVEL=READ_GROUP \
	        O=reseq/${SAMPLE}.alignment_metrics.txt

    # Report on the validity of BAM files relative to the SMA format specification
    java -Xmx8g -XX:ParallelGCThreads=1 -jar ${machDIR}/picard.jar ValidateSamFile \
          I=reseq/${SAMPLE}.dedup.sorted.bam \
          O=reseq/${SAMPLE}.validate.txt \
          MODE=SUMMARY

done


# Assess indels from all BAM files
${JAVA_8}/java -Xmx16g -XX:ParallelGCThreads=1 -jar ${machDIR2}/GenomeAnalysisTK.jar \
              -T RealignerTargetCreator \
              -nt ${NTHREADS} \
              -R WGS/BFAL_genome.fasta \
              -I reseq/${SAMPLE}_1.dedup.sorted.bam \
              -I reseq/${SAMPLE}_2.dedup.sorted.bam \
              -I ...
              -I reseq/${SAMPLE}_30.dedup.sorted.bam \
              -o reseq/test_BFAL_LAAL_indel.intervals

# Realign BAM files around indels
for SAMPLE in `cat list_reseq_BFAL_LAAL`
do

    ${JAVA_8}/java -Xmx8g -XX:ParallelGCThreads=2 -jar ${machDIR2}/GenomeAnalysisTK.jar \
                          -T IndelRealigner \
                          -R WGS/BFAL_genome.fasta \
                          -I reseq/${SAMPLE}.dedup.sorted.bam \
                          -targetIntervals reseq/test_BFAL_LAAL_indel.intervals \
                          -o reseq/${SAMPLE}.realigned.bam &
done


