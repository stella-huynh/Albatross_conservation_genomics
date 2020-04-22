######################################################################################################
### ----   TUTORIAL I : reconstruction of a mitochondrial genome using a two step procedure   ---- ###

# Mapping assembly

echo -e "\n#manifest file for basic mapping assembly with illumina data using MIRA 4\n\nproject = initial-mapping-testpool-to-Salpinus-mt\n\njob=genome,mapping,accurate\n\nparameters = -NW:mrnl=0 -AS:nop=1 SOLEXA_SETTINGS -CO:msr=no\n\nreadgroup\nis_reference\ndata = reference.fa\nstrain = Salpinus-mt-genome\n\nreadgroup = reads\ndata = reads.fastq\ntechnology = solexa\nstrain = testpool\n" > manifest.conf

mira manifest.conf





./MITObim.pl -sample testpool \
			 -ref Salpinus_mt_genome \
			 -readpool reads.fastq \
			 -maf initial-mapping-testpool-to-Salpinus-mt_assembly/initial-mapping-testpool-to-Salpinus-mt_d_results/initial-mapping-testpool-to-Salpinus-mt_out.maf \
			 -start 1 -end 10 \
			 --redirect_tmp /data/huynhs/ \
			 &> log


########################################################################################
### ----   TUTORIAL II : direct reconstruction without prior mapping assembly   ---- ###


./MITObim.pl -sample testpool \
			 -ref Salpinus_mt_genome \
			 -readpool ../testdata1/Tthymallus-150bp-300sd50-interleaved.fastq \
			 --quick ../testdata1/Salpinus-mt-genome-NC_000861.fasta \
			 -start 1 -end 30 \
			 --redirect_tmp /data/huynhs/ \
			 &> log






######################################################################################################
### ----   TUTORIAL III : reconstructing mt genomes from mt barcode seeds   ---- ###

# This tutorial reconstructs the mt genome of T. thymallus solely using a partial mitochondrial COI sequence as starting seed.
# MITObim reconstructs the mitchondrial genome in 82 (or less) iterations.


./MITObim.pl -sample testpool \
			 -ref Tthymallus-COI \
			 -readpool ../testdata1/Tthymallus-150bp-300sd50-interleaved.fastq \
			 --quick ../testdata1/Tthymallus-COI-partial-HQ961018.fasta \
			 -end 100 --clean \
			 --redirect_tmp /data/huynhs/ \
			 &> log



# For "well behaved" datasets the standard mapping assembly can be substituted by a de novo assembly (--denovo flag).
# Utilizing read pair information (--paired flag) can further speed up the reconstruction if run in de novo mode.
# This strategy reconstructs the correct mitochondrial genome in only 31 (or less) iterations.


./MITObim.pl -sample testpool \
			 -ref Tthymallus-COI \
			 -readpool ../testdata1/Tthymallus-150bp-300sd50-interleaved.fastq \
			 --quick ../testdata1/Tthymallus-COI-partial-HQ961018.fasta \
			 -end 50 --clean \
			 --denovo --paired \
			 --redirect_tmp /data/huynhs/ \
			 &> log


