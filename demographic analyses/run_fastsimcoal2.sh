#!/bin/bash

SPECIES=$1
REGION=$2
FSC_model=$3

for SP in ${SPECIES}
do

    PREFIX="${SP}_${REGION}"
    mkdir -m 775 ${FSC_model}/${PREFIX}
    
    ## Run bootstrapped analysis ##
    
    for i in {1..100}
    do
    
        mkdir -m 775 ${FSC_model}/${PREFIX}/run$i
        cp ${PREFIX}.tpl ${PREFIX}.est ${PREFIX}_DAFpop0.obs ${PREFIX}/run$i"/"
        cd ${PREFIX}/run$i
        
        echo -e "\nRunning fastsimcoal for ${PREFIX}.\n"
        fsc26 -t ${PREFIX}.tpl -e ${PREFIX}.est -d -C 10 -n 500000 -L 40 -s 0 -M -c ${nt} -B ${nt} -b10 -j i-q &
        cd ../..
        
     done
        
        
     ## Summarize likelihoods over bootstraps ##  
     
     cd ${PREFIX}
     fsc-selectbestrun.sh ${PREFIX}
     cd bestrun/
     calculateAIC.sh ${PREFIX}
     
     
     ## Visualize model ##
     
     SFStools.R -t print2D -i ${PREFIX}
     
        
        
        
        
