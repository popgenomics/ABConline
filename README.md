# 1 - snakemake  
**The entire workflow is based on snakemake, which is essential for analysis.**  
https://snakemake.readthedocs.io/en/stable/  

# 2 - python  
2pops/fasta2ABC_2pops.py  
2pops/mscalc_2pop_observedDataset.py  
2pops/mscalc_2pop.py  
2pops/priorgen_2pop.py  
2pops/priorgen_gof_2pop.py  
2pops/submit_simulations_2pop.py  
2pops/submit_simulations_gof_2pop.py  

## dependencies  
**uses pypy as python interpreter**    
import sys  
import os  
from math import ceil  
import random  
from numpy.random import uniform  
from numpy.random import binomial  
from numpy.random import beta  
from random import shuffle  
from numpy.random import randint  
import time  
  
# 3 - R  
**uses Rscript from /usr/bin or elsewhere**  
2pops/collaborative_plot.R  
2pops/estimates_2pop_best.R  
2pops/estimates_2pop.R  
2pops/get_parameters.R  
2pops/gof_2pop.R  
2pops/model_comp_2pop_allModels.R  
2pops/model_comp_2pop.R  
  
## dependencies  
library(plotly)  
library(viridis)  
library(abcrf)  
library('nnet')  
library(ggplot2)  
library(ggpubr)  
  
# 4 - C
## **msnsam**  
### info  
C code, compiled by executing the command ./clms (calling gcc) in the msnsam/ directory  
   
##**RNAseqFGT_src**  
### info  
C code compiled by: cc -Wall -o RNAseqFGT RNAseqFGT.c RNAseqFGT_seq_reading.c RNAseqFGT_analysis.c -I RNAseqFGT.h  
  
# 5 - external  
**pandoc** (https://pandoc.org/index.html)  

