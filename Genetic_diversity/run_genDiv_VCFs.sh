#!/bin/bash

java8="/group/sbs_ssin/applications/jdk1.8.0_221/bin"

REF="/group/sbs_ssin/stella/albatross/00_genomes/01_WGS/01_dna/BFAL_genome.fasta"
wDIR="/group/sbs_ssin/stella/albatross/05_reseq/initial_pipeline/06f_ANGSD"
oDIR="${wDIR}/01_VCF" ; mkdir -m 775 ${oDIR} ${oDIR}/tmp


############################################################

# A0. GATK HaplotypeCaller #

for b in `cat ${BAMLIST}`
do
       python run_hapCaller.py -b $b \
                               -r $REF \
                               --outDir $oDIR --runName HC \
                               --gatkDir $machDIR --tmpDir $oDIR/tmp \
                               --javaRam 12g --threads ${NCORES} \
                               -ploidy 2 --heterozygosity 0.001 \
                               --outputMode EMIT_ALL_CONFIDENT_SITES
done


#############################################################################

# A1. Combine GVCFs by population (=species) #

for SP in BFAL LAAL
do
        SAMPLES=$(ls ${oDIR}/${SP}_*.g.vcf.gz | sed -e 's/^/--variant /')

        ${java8}/java -Xmx12g -jar \
                ${machDIR}/GenomeAnalysisTK.jar \
                        -T CombineGVCFs \
                        -R $REF \
                        ${SAMPLES} \
                        -o ${oDIR}/${SP}.vcf.gz \
                        -nt ${NCORES} &
done
wait


# A2. Genotype GVCFs by population (=species) #

for SP in BFAL LAAL
do
        ${java8}/java -Xmx12g -jar \
                ${machDIR}/GenomeAnalysisTK.jar \
                        -T GenotypeGVCFs \
                        -R $REF \
                        -V ${oDIR}/${SP}.vcf.gz \
                        -o ${oDIR}/${SP}.HC.vcf.gz \
                        -nt ${NCORES} \
                        --includeNonVariantSites &
done
wait


#############################################################################

# A3. Filter VCFs #

for SP in BFAL LAAL
do
        bcftools view -U ${oDIR}/${SP}.HC.vcf.gz | bgzip > ${oDIR}/${SP}.HC.clean.vcf.gz &
done
wait

#rm ${oDIR}/${SP}.HC.vcf.gz.tbi
mv ${oDIR}/${SP}.HC.clean.vcf.gz ${oDIR}/${SP}.HC.vcf.gz


############################################################

# A4. Make Geno files required as input #


for SP in BFAL LAAL
do
        bcftools filter -e 'FORMAT/DP < 5' --set-GTs "." -O u ${oDIR}/${SP}.HC.vcf.gz | \
        python parseVCF.py --skipIndels | bgzip \
        > ${oDIR}/${SP}.HC.DP5.geno.gz &

        python popgenWindows.py -g ${oDIR}/${SP}.HC.DP5.geno.gz \
                                -f phased --analysis indHet popDist \
                                -o ${oDIR}/${SP}.HC.DP5.het.w100m50s50.csv \
                                -w 100000 -m 50 -s 50000 --roundTo 5 -T 10 &
done







