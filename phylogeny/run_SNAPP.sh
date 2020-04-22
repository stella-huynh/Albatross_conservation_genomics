###################################################################
#############                  SNAPP                  #############
###################################################################

plink -tfile BFAL_autosome_1e6snp_05maf_ghwe_LD01 --allow-extra-chr --recode 12 -out BFAL_autosome_1e6snp_05maf_ghwe_LD01_recode

plink --vcf ../SNAPP/ALB5_autosome_05maf_1e6snp_mxDepth10f_baQ_DEPTHpersamp3.fil.SNPs.thin.recode.vcf --make-bed --allow-extra-chr --double-id --bp-space 1000 --out ALB5_autosome_05maf_1e6snp_mxDepth10f_baQ_DEPTHpersamp3.fil.SNPs.thin.recode

plink --bfile ../TreeMIX/ALB5_autosome_05maf_1e6snp_mxDepth10f_baQ_DEPTHpersamp3.fil.SNPs.thin.recode --freq --missing --allow-extra-chr --out ../TreeMIX/ALB5_autosome_05maf_1e6snp_mxDepth10f_baQ_DEPTHpersamp3.fil.SNPs.thin.recode



conda activate /home/huynhs/programs/envs/ruby
ruby ../../scripts/snapp_prep.rb -v ALB5_autosome_05maf_1e6snp_mxDepth10f_baQ_DEPTHpersamp3.fil.SNPs.vcf -t ALB5_samples.txt -c ALB5_constraints.txt -s ALB5_start_tree.newick -l 500000 -q 1000 -x ALB5_autosome_05maf_1e6snp_mxDepth10f_baQ_DEPTHpersamp3.xml -o ALB5_autosome_05maf_1e6snp_mxDepth10f_baQ_DEPTHpersamp3



####################################################################
##  RESULTS  #######################################################
####################################################################

#Found 5586688 SNP records in file '../SNAPP/ALB5_autosome_05maf_1e6snp_mxDepth10f_baQ_DEPTHpersamp3.fil.SNPs.thin.recode.vcf'. (Skipped 0 already filtered-out SNPs and 0 non-SNP records ; more with --verbose.)
#Removed 0 loci that did not pass sample/population constraints from 5586688 loci.
#Kept 5586688 loci, composed of 5586688 sites; 0 of those sites were filtered, 5586688 variant sites remained.
#    5586688 genomic sites, of which 0 were covered by multiple loci (0.0%).
#Mean genotyped sites per locus: 1.00bp (stderr 0.00).
#Population summary statistics (more detail in populations.sumstats_summary.tsv):
#  BFAL: 1 samples per locus; pi: 0.061811; all/variant/polymorphic sites: 5586688/5586688/345321; private alleles: 545461
#  LAAL: 1 samples per locus; pi: 0.064216; all/variant/polymorphic sites: 5586688/5586688/358757; private alleles: 496167
#  STAL: 1 samples per locus; pi: 0.096823; all/variant/polymorphic sites: 5586688/5586688/540919; private alleles: 687756
#  WAAL: 1 samples per locus; pi: 0.060751; all/variant/polymorphic sites: 5586688/5586688/339396; private alleles: 1004144
#  WAL: 1 samples per locus; pi: 0.04738; all/variant/polymorphic sites: 5586688/5586688/264700; private alleles: 2119242
#Populations is done.

#######################

#Found 1754491 SNP records in file '../SNAPP/ALB_autosome_05maf_1e6snp_mxDepth10f_baQ_DEPTHpersamp3.fil.SNPs.thin.recode.vcf'. (Skipped 0 already filtered-out SNPs and 0 non-SNP records ; more with --verbose.)
#Removed 0 loci that did not pass sample/population constraints from 1754491 loci.
#Kept 1754491 loci, composed of 1754491 sites; 0 of those sites were filtered, 1754491 variant sites remained.
#    1754491 genomic sites, of which 0 were covered by multiple loci (0.0%).
#Mean genotyped sites per locus: 1.00bp (stderr 0.00).
#Population summary statistics (more detail in populations.sumstats_summary.tsv):
#  BFAL: 15.531 samples per locus; pi: 0.1174; all/variant/polymorphic sites: 1754491/1754491/697909; private alleles: 488915
#  LAAL: 14.469 samples per locus; pi: 0.15218; all/variant/polymorphic sites: 1754491/1754491/838442; private alleles: 514579
#  STAL: 1 samples per locus; pi: 0.065893; all/variant/polymorphic sites: 1744278/1744278/114935; private alleles: 33
#  WAAL: 1 samples per locus; pi: 0.041245; all/variant/polymorphic sites: 1732127/1732127/71442; private alleles: 142
#  WAL: 1 samples per locus; pi: 0.053991; all/variant/polymorphic sites: 1694072/1694072/91465; private alleles: 208
#Populations is done.
