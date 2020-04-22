#!/bin/bash


#==============================================================================================================
### Download UCEs from Avian Phylogenomics Project
#==============================================================================================================
#http://gigadb.org/dataset/101041

tar -xf FASTA_files_of_loci_datasets.tar.gz -C ./
tar -xf uce-filtered-alignments-without-gator.tar.gz -C ./   <--------- this is the 3679 UCEs to be used  

### Reference UCEs species
### 1. Gallus gallus  (chicken)              --> GALGA
### 2. Gavia stellata (red-throated loon)    --> GAVST
### 3. Fulmarus glacialis (northern fulmar)  --> FULGL
### 4. Eurypygia helias (little egret)       --> EURHE



