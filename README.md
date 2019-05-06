# Tools  
## 1. 2pops/fasta2ABC_2pops.py  
uses pypy as python interpreter    
### location on laptop  
/usr/bin/pypy
### location on core.cluster.france-bioinformatique.fr
/shared/mfs/data/software/miniconda/envs/pypy-2.7-5.10.0/bin/pypy  


## 2. 2pops/priorgen_2pop.py  
### imports
import sys  
from numpy.random import uniform  
from numpy.random import binomial  
from numpy.random import beta  
from random import shuffle  


## 3. msnsam
### info  
C code, compiled by executing the command ./clms (calling gcc) in the msnsam/ directory  


## 4. 2pops/mscalc_2pop.py and 2pops/mscalc_2pop_observedDataset.py
### info  
The two versions are required.  
They both need pypy as python interpreter  
 

## 5. RNAseqFGT_src  
### info  
C code compiled by: cc -Wall -o RNAseqFGT RNAseqFGT.c RNAseqFGT_seq_reading.c RNAseqFGT_analysis.c -I RNAseqFGT.h  


## 6. 2pops/model_comp_2pop.R  
### info  
R script (#!/usr/bin/Rscript)  
### dependencies  
only requires the R package 'abcrf' (library('abcrf'))  


# ABConline
ABC.py infile=all_loci.fasta nspecies=2 nameA=txn nameB=ama region=coding Lmin=30 max_N_tolerated=0.2 nMin=15 Nref=100000 mu=0.000000003 rho_over_theta=1 nSimulations=5000

