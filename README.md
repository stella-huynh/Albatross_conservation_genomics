## **Albatross_conservation_genomics**

Commands and scripts from albatross genomics project

### **Genome assembly**

#### **Satsuma alignment with chicken genome**

Alignment
```
for sp in BFAL LAAL STAL WAAL WAL
do
    SatsumaSynteny2 -q WGS/${sp}_genome.fasta \
                    -t REF/Gallus_gallus.GRCg6a.dna.toplevel_chrom.fa \
                    -o Satsuma/${sp} \
                    -sl_mem 100 -slaves 4 -threads ${NTHREADS} \
                    -km_mem 100 -km_sync 1 -q_chunk 1000 \
                    -min_prob 0.9999
done
```

Pseudochromome
```
for sp in BFAL LAAL STAL WAAL WAL
do
    Chromosemble -q WGS/${sp}_genome.fasta \
                 -t REF/Gallus_gallus.GRCg6a.dna.toplevel_chrom.fa \
                 -o WGS/${sp}_superscaffolds.fasta \
                 -n ${NTHREADS} -thorough 0 -pseudochr 1 -s 0
done
```


### **Read mapping and processing of resequenced data**

Map raw reads 

Build genome index
```
bwa index WGS/BFAL_genome.fasta
samtools faidx WGS/BFAL_genome.fasta
picard CreateSequenceDictionary REFERENCE=REF/BFAL_genome.fasta OUTPUT=REF/BFAL_genome.dict
```

GATK pipeline
```
bash run_GATK_pipeline.sh list_sample_BFAL_LAAL.txt
```


### **Genetic diversity**

#### **Identify sex-linked scaffolds**

1. From read coverage difference between male and females samples
```
samtools idxstats reseq/${SAMPLE}.realigned.bam > reseq/idxstats_${SAMPLE}.txt

Rscript Read_cov_male_females.R 
```
2. From SatsumaSynteny2 mapping results

Combine list of scaffolds mapped onto galGal6 Z/W chromosomes (Satsuma results) and list of scaffolds with 2x lower coverage in females than males ("Read_cov_male_females.R" results).

#### **Infer SFS for BFAL and LAAL populations**

Use short-tailed albatross *Phoebastria albatrus* (STAL) to infer ancestral states of alleles.
```
# Map raw read of STAL onto BFAL following GATK pipeline
bash run_GATK_pipeline.sh list_sample_STAL.txt

# Convert to BAM to FASTA file
angsd -i WGS/STAL_220bp.realigned.bam -doFasta 1 -doCounts 1 -out WGS/STAL_ANGSDgenome.fasta
gunzip WGS/STAL_ANGSDgenome.fasta.gz
```

Estimate allelic site frequencies (SAF) in ANGSD
```
angsd   -ref WGS/BFAL_genome.fasta -anc WGS/STAL_ANGSDgenome.fasta -out reseq/${OUTFILE} \
        -bam reseq/${BAMLIST} -rf reseq/${REGION}.txt \
        -nThreads ${NTHREADS} \
        -remove_bads 1 -uniqueOnly 1 -only_proper_pairs 0 -minMapQ 20 -minQ 20 -trim 0 \
        -doMajorMinor 1 -skipTriallelic 1 -GL 1 -doSaf 1
```
Estimate unfolded 1D-SFS
```
realSFS reseq/${OUTFILE}.saf.idx -P ${NTHREADS} > reseq/${OUTFILE}_1DSFS.sfs #without bootstrap
realSFS reseq/${OUTFILE}.saf.idx -P ${NTHREADS} -bootstrap 20 > reseq/${OUTFILE}_1DSFS.sfs #with 20 bootstraps
```
Estimate unfolded 2D-SFS
```
realSFS reseq/${F1}.saf.idx reseq/${F2}.saf.idx -P ${NTHREADS} > reseq/${OUTFILE}_1DSFS.sfs #without bootstrap
realSFS reseq/${F1}.saf.idx reseq/${F2}.saf.idx -P ${NTHREADS} -bootstrap 20 > reseq/${OUTFILE}_1DSFS.sfs #with 20 bootstraps
```



