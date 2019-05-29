# Table of contents
1. [Snakemake](#1---snakemake)  
2. [Scripts in python](#2---python)  
	1. [scripts](#scripts)  
	2. [dependencies](#dependencies)  
3. [Scripts in R](#3---r)  
	1. [scripts](#scripts)  
	2. [dependencies](#dependencies)  
4. [Codes in C](#4---c)  
	1. [msnsam (by Jeffrey Ross-Ibarra)](#msnsam)  
	2. [RNAseqFGT (by Laurent Duret)](#RNAseqFGT)  
5. [External codes](#5---external)  

# 1 - snakemake  
**The entire workflow is based on snakemake, which is essential for analysis.**  
https://snakemake.readthedocs.io/en/stable/  

# 2 - python  
**Executables have to be linked to a bin directory**  
## scripts  
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
import random  
import time  
from math import ceil  
from numpy.random import uniform  
from numpy.random import binomial  
from numpy.random import beta  
from numpy.random import randint  
  
# 3 - R  
**Executables have to be linked to a bin directory**  
## scripts  
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
**Executables have to be linked to a bin directory**  
## msnsam  
### info  
C code, compiled by executing the command ./clms (calling gcc) in the msnsam/ directory  
   
## RNAseqFGT  
### info  
C code compiled by: cc -Wall -o RNAseqFGT RNAseqFGT.c RNAseqFGT_seq_reading.c RNAseqFGT_analysis.c -I RNAseqFGT.h  
  
# 5 - external  
**pandoc** (https://pandoc.org/index.html)  

