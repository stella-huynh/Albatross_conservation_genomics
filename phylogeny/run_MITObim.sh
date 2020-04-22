#!/bin/bash

conda install -p /group/sbs_ssin/anaconda3 -c bioconda mira
git clone --recursive git://github.com/chrishah/MITObim.git


#Use the biggest run library for each species
 #-fragment, paired-end library (220bp)
 #-jump, mate-pair library (3kb)

# Copy the raw reads (.FASTQ) in the project subdirectory called 'data'

# Create and modify the 'manifest.conf' for MIRA

# Run sequence mapping with MIRA

mira manifest.conf &> log

