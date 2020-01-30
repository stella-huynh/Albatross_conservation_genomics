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
                 -o ${oDIR}/${sp}_scaffoldssyn.fasta \
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
bash run_GATK_pipeline.sh list_reseq_BFAL_LAAL
```


