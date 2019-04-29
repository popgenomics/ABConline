import glob
import re
import sys
from os.path import join 

nmultilocus = 200 # number of multilocus simulations per iteration
nIterations_model_comp = 5 # number of subdirectories for the simulations used in the RF model comparison
nCPU_modelComp = 8 # number of CPUs for the model comp for the model forest R functions
ntree = 1000 # number of tree for the random forest (RF) model comparison

infile = '../../all_loci.fasta'
nspecies = 2
nameA = 'txn'
nameB = 'ama'
outgroup = 'num'
region = 'coding'
Lmin = 250
max_N_tolerated = 0.007
nMin = 12
Nref = 100000
mu = 0.000000003
rho_over_theta = 1

#ruleorder: fasta2ABC_2pops > RNAseqFGT > simulations_SI > simulations_AM
ruleorder: fasta2ABC_2pops > RNAseqFGT > simulations > model_comparison 

MODELS_2POPS = ['SC_1M_1N', 'SC_1M_2N', 'SC_2M_1N', 'SC_2M_2N', 'AM_1M_1N', 'AM_1M_2N', 'AM_2M_1N', 'AM_2M_2N', 'IM_1M_1N', 'IM_1M_2N', 'IM_2M_1N', 'IM_2M_2N', 'SI_1N', 'SI_2N']
#MODELS_2POPS = ['SC_2M_2N', 'AM_2M_2N', 'IM_2M_2N', 'SI_2N']
ITERATIONS_MODEL_COMP = range(nIterations_model_comp)

rule targets:
	input:
		expand("ABC_{nameA}_{nameB}/ABCstat_global.txt", nameA=nameA, nameB=nameB),
		expand("ABC_{nameA}_{nameB}/ABCstat_loci.txt", nameA=nameA, nameB=nameB),
		expand("ABC_{nameA}_{nameB}/bpfile", nameA=nameA, nameB=nameB),
		expand("ABC_{nameA}_{nameB}/nLoci.txt", nameA=nameA, nameB=nameB),
		expand("ABC_{nameA}_{nameB}/{nameA}_{nameB}_infos.txt", nameA=nameA, nameB=nameB),
		expand("ABC_{nameA}_{nameB}/{nameA}_{nameB}.ms", nameA=nameA, nameB=nameB),
		expand("ABC_{nameA}_{nameB}/results_recombination.txt", nameA=nameA, nameB=nameB),
		expand("ABC_{nameA}_{nameB}/{model}_{i}/ABCstat.txt", nameA=nameA, nameB=nameB, model=MODELS_2POPS, i=ITERATIONS_MODEL_COMP),
		expand("ABC_{nameA}_{nameB}/report_{nameA}_{nameB}.txt", nameA=nameA, nameB=nameB)


rule fasta2ABC_2pops:
	# fasta2ABC_2pops.py ../../all_loci.fasta txn mal num coding 250 0.007 12 100000 0.000000003 1
	params:
		nameA={nameA},
		nameB={nameB},
		outgroup={outgroup},
		region={region},
		Lmin={Lmin},
		max_N_tolerated={max_N_tolerated},
		nMin={nMin},
		Nref={Nref},
		mu={mu},
		rho_over_theta={rho_over_theta}
	input:
		expand("{infile}", infile=infile)
	output:
		expand("ABC_{nameA}_{nameB}/ABCstat_global.txt", nameA=nameA, nameB=nameB),
		expand("ABC_{nameA}_{nameB}/ABCstat_loci.txt", nameA=nameA, nameB=nameB),
		expand("ABC_{nameA}_{nameB}/bpfile", nameA=nameA, nameB=nameB),
		expand("ABC_{nameA}_{nameB}/nLoci.txt", nameA=nameA, nameB=nameB),
		expand("ABC_{nameA}_{nameB}/{nameA}_{nameB}_infos.txt", nameA=nameA, nameB=nameB),
		expand("ABC_{nameA}_{nameB}/{nameA}_{nameB}.ms", nameA=nameA, nameB=nameB)
	threads: 8
	shell:
		"fasta2ABC_2pops.py {infile} {params.nameA} {params.nameB} {params.outgroup} {params.region} {params.Lmin} {params.max_N_tolerated} {params.nMin} {params.Nref} {params.mu} {params.rho_over_theta}"


rule RNAseqFGT:
	params:
		nameA={nameA},
		nameB={nameB}
	input:
		expand("{infile}", infile=infile)
	output:
		expand("ABC_{nameA}_{nameB}/results_recombination.txt", nameA=nameA, nameB=nameB)
	threads: 8
	shell:	
		"RNAseqFGT {infile} ABC_{params.nameA}_{params.nameB}/results_recombination.txt"


rule simulations:
	#commande = 'submit_simulations_2pop.py {0} {1} {2} {3} {4}'.format(nMultilocus, nCPU, model, args['nameA'], args['nameB'])
	params:
		nmultilocus={nmultilocus}
	input:
		"ABC_{nameA}_{nameB}/bpfile",
		"ABC_{nameA}_{nameB}/nLoci.txt"
	output:
		"ABC_{nameA}_{nameB}/{model}_{i}/bpfile",
		"ABC_{nameA}_{nameB}/{model}_{i}/priorfile.txt",
		"ABC_{nameA}_{nameB}/{model}_{i}/ABCstat.txt"
	threads: 1
	shell:
		"""
		submit_simulations_2pop.py {params.nmultilocus} {wildcards.i} {wildcards.model} {nameA} {nameB}
		"""


rule model_comparison:
	input:
		expand("ABC_{nameA}_{nameB}/ABCstat_global.txt", nameA=nameA, nameB=nameB),
		expand("ABC_{nameA}_{nameB}/{model}_{i}/ABCstat.txt", nameA=nameA, nameB=nameB, model=MODELS_2POPS, i=ITERATIONS_MODEL_COMP)
	output:
		"ABC_{nameA}_{nameB}/report_{nameA}_{nameB}.txt"
	threads: nCPU_modelComp
	shell:
		"""
		model_comp_2pop.R nameA={nameA} nameB={nameB} nMin={nMin} nSubdir={nIterations_model_comp} ntree={ntree} ncores={nCPU_modelComp} outgroup=1
		#cat {input} >{output}
		"""

