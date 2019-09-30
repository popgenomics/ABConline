#!/usr/bin/python
# #!/home/roux/python/Python-2.7.14/python
# -*- coding: utf-8 -*-

import sys
from numpy.random import randint
from numpy.random import uniform
from numpy.random import binomial
from numpy.random import beta
from numpy import log
from numpy import median
from random import shuffle
help = "\t\033[1;31;40mTakes one model specifier, a number of multilocus simulations and a config.yaml file containing prior boundaries as arguments:\033[0m\n\t\t"
#help += "\n\t\t".join(["Constant_1N", "Constant_2N", "Discrete_1N", "Discrete_2N", "Expo_1N", "Expo_2N"])
help += "\n\t\t".join(["Constant_1N", "Constant_2N", "Expansion_1N", "Expansion_2N", "Contraction_1N", "Contraction_2N"])
help += "\n\n"
help += "\t\033[1;32;40m#Constant\033[0m\n\tmsnsam tbs 10000 -t tbs -r tbs tbs\n"
help += "\t\033[1;32;40m#Discrete\033[0m\n\tmsnsam tbs 10000 -t tbs -r tbs tbs -eN tbs tbs\n"
help += "\t\033[1;32;40m#Expo\033[0m\n\tmsnsam tbs 10000 -t tbs -r tbs tbs -G tbs -eG tbs 0.0 -eN tbs tbs\n"

help += "\t\033[1;32;40mExample: ./priorgen_gof_1pop.py Constant_1N 1000 posterior_file modePrior (joint; disjoint; randomBeta)\033[0m\n"

if len(sys.argv) != 5:
	print(help)
	sys.exit()

def randomBeta(posterior, nMultilocus):
	estimate = median(posterior)
	a = 50.0
	b = 50.0
	scalar_tmp = beta(a=a, b=b, size=nMultilocus)
	scalar = [ i/(a/(a+b)) for i in scalar_tmp ]
	
	res = [ estimate * i for i in scalar ]
	return(res)

# Configuration of the prior distribution
nMultilocus = int(sys.argv[2])
modePrior = sys.argv[4]

# read bpfile
infile = open("bpfile", "r")
tmp = infile.readline()
L = infile.readline().strip().split("\t") # in number of nucleotides
nsamA = infile.readline().strip().split("\t") # number of individuals within the population 
mu = infile.readline().strip().split("\t") # mutation rate per generation and per base pair
rec = infile.readline().strip().split("\t") # recombiation rate per generation and per base pair
infile.close()

# converts in integer/floating values
L = [ float(i) for i in L ]
nsamA = [ int(i) for i in nsamA ]
mu = [ float(i) for i in mu ]
rec = [ float(i) for i in rec ]

# number of loci
nLoci = len(L)

# get the posterior
infile = open(sys.argv[3], 'r')
posterior = {}
params_posterior = []
header = infile.readline().strip().split('\t')
for i in header:
	params_posterior.append(i)
	if i not in posterior:
		posterior[i] = []
for line in infile:
	line = line.strip().split('\t')
	cnt=0
	for i in line:
		cnt += 1
		posterior[ params_posterior[cnt-1] ].append(float(i))
infile.close()

# get the lines of the posterior used for the simulations: vector of length nMultilocus
# used_posterior = randint(cnt, size=nMultilocus)

model_changes = ['Contraction_1N', 'Contraction_2N', 'Expansion_1N', 'Expansion_2N']

if modePrior == "joint": # modePrior in {joint; disjoint; randomBeta}, where joint takes the exact joint values of the posterior as a prior; disjoint takes random associations of parameter values from posterior; randombeta simulates a beta distribution around the median of the posterior
	used_posterior = [ randint(0, cnt-1) for i in range(nMultilocus) ]
	N = [ posterior['N'][i] for i in used_posterior ]
	
	if sys.argv[1] in model_changes:
		Npast = [ posterior['Npast'][i] for i in used_posterior ]
		Tdem = [ posterior['Tdem'][i] for i in used_posterior ]
	if '2N' in sys.argv[1]:
		shape_N_a = [ posterior['shape_N_a'][i] for i in used_posterior ]
		shape_N_b = [ posterior['shape_N_b'][i] for i in used_posterior ]
else:
	if modePrior == "disjoint":
		used_posterior = [ randint(0, cnt-1) for i in range(nMultilocus) ]
		N = [ posterior['N'][i] for i in used_posterior ]
		
		if sys.argv[1] in model_changes:
			used_posterior = [ randint(0, cnt-1) for i in range(nMultilocus) ]
			Npast = [ posterior['Npast'][i] for i in used_posterior ]
			used_posterior = [ randint(0, cnt-1) for i in range(nMultilocus) ]
			Tdem = [ posterior['Tdem'][i] for i in used_posterior ]
		if '2N' in sys.argv[1]:
			used_posterior = [ randint(0, cnt-1) for i in range(nMultilocus) ]
			shape_N_a = [ posterior['shape_N_a'][i] for i in used_posterior ]
			used_posterior = [ randint(0, cnt-1) for i in range(nMultilocus) ]
			shape_N_b = [ posterior['shape_N_b'][i] for i in used_posterior ]

	else: # if modeprior == 'randombeta'
		N = randomBeta(posterior['N'], nMultilocus)
		
		if sys.argv[1] in model_changes:
			Npast = randomBeta(posterior['Npast'], nMultilocus)
			Tdem = randomBeta(posterior['Tdem'], nMultilocus)
		if '2N' in sys.argv[1]:
			shape_N_a = randomBeta(posterior['shape_N_a'], nMultilocus)
			shape_N_b = randomBeta(posterior['shape_N_b'], nMultilocus)


if sys.argv[1] == "Constant_1N":
	# msnsam tbs 10000 -t tbs
	# param multilocus: values that will be printed in priorfile.txt
	## N = N_pop_i / Nref
	theta = [ 4*L[i]*mu[i] for i in range(nLoci) ] # vector of theta/N
	rho = [ 4*L[i]*rec[i] for i in range(nLoci) ] # vector of rho/N

	# param monolocus: values that will be read by ms
	priorfile = "N\n"
	for sim in range(nMultilocus):
		priorfile += "{0:.5f}\n".format(N[sim])
		for locus in range(nLoci):
			theta_locus = theta[locus]*N[sim]
			rho_locus = rho[locus]*N[sim]
			if theta_locus < 0.00001:
				theta_locus = 0.00001
			print("{0}\t{1}\t{2}\t{3}".format(nsamA[locus], theta_locus, rho_locus, L[locus]))
	outfile = open("priorfile.txt", "w")
	outfile.write(priorfile)
	outfile.close()

if sys.argv[1] == "Constant_2N":
	# msnsam tbs 10000 -t tbs
	# param multilocus: values that will be printed in priorfile.txt
	## N = N_pop_i / Nref
	theta = [ 4*L[i]*mu[i] for i in range(nLoci) ] # vector of theta/N
	rho = [ 4*L[i]*rec[i] for i in range(nLoci) ] # vector of rho/N

	# param monolocus: values that will be read by ms
	priorfile = "N\tshape_N_a\tshape_N_b\n"
	for sim in range(nMultilocus):
		priorfile += "{0:.5f}\t{1:.5f}\t{2:.5f}\n".format(N[sim], shape_N_a[sim], shape_N_b[sim])
		# vectors of size 'nLoci' containing parameters
                scalar_N = beta(shape_N_a[sim], shape_N_b[sim], size=nLoci)
                rescale = shape_N_a[sim] / (shape_N_a[sim] + shape_N_b[sim])
		
		for locus in range(nLoci):
			theta_locus = theta[locus]*N[sim]*scalar_N[locus]/rescale # theta = (4.mu.N) * Beta / ( a / (a+b) )
			rho_locus = rho[locus]*N[sim]*scalar_N[locus]/rescale # theta = (4.mu.N) * Beta / ( a / (a+b) )
			if theta_locus < 0.00001:
				theta_locus = 0.00001
			print("{0}\t{1}\t{2}\t{3}".format(nsamA[locus], theta_locus, rho_locus, L[locus]))
	outfile = open("priorfile.txt", "w")
	outfile.write(priorfile)
	outfile.close()

if sys.argv[1] == "Expansion_1N":
	# msnsam tbs 10000 -t tbs -r tbs tbs -eN tbs tbs
	# nsamA theta 
	# param multilocus: values that will be printed in priorfile.txt
	## N = N_pop_i / Nref
	theta = [ 4*L[i]*mu[i] for i in range(nLoci) ] # vector of theta/N
	rho = [ 4*L[i]*rec[i] for i in range(nLoci) ] # vector of rho/N
	
	# param monolocus: values that will be read by ms
	priorfile = "N\tNpast\tTdem\n"
	for sim in range(nMultilocus):
		priorfile += "{0:.5f}\t{1:.5f}\t{2:.5f}\n".format(N[sim], Npast[sim], Tdem[sim])
		for locus in range(nLoci):
			theta_locus = theta[locus]*N[sim]
			rho_locus = rho[locus]*N[sim]
			if theta_locus < 0.00001:
				theta_locus = 0.00001
			print("{0}\t{1}\t{2}\t{3}\t{4}\t{5}".format(nsamA[locus], theta_locus, rho_locus, L[locus], Tdem[sim]/(4.0*N[sim]), Npast[sim]/(1.0*N[sim])) )
	outfile = open("priorfile.txt", "w")
	outfile.write(priorfile)
	outfile.close()

if sys.argv[1] == "Expansion_2N":
	# msnsam tbs 10000 -t tbs -r tbs tbs -eN tbs tbs
	# nsamA theta 
	# param multilocus: values that will be printed in priorfile.txt
	## N = N_pop_i / Nref
	theta = [ 4*L[i]*mu[i] for i in range(nLoci) ] # vector of theta/N
	rho = [ 4*L[i]*rec[i] for i in range(nLoci) ] # vector of rho/N
	
	# param monolocus: values that will be read by ms
	priorfile = "N\tNpast\tshape_N_a\tshape_N_b\tTdem\n"
	for sim in range(nMultilocus):
		priorfile += "{0:.5f}\t{1:.5f}\t{2:.5f}\t{3:.5f}\t{4:.5f}\n".format(N[sim], Npast[sim], shape_N_a[sim], shape_N_b[sim], Tdem[sim])
                scalar_N = beta(shape_N_a[sim], shape_N_b[sim], size=nLoci)
                rescale = shape_N_a[sim] / (shape_N_a[sim] + shape_N_b[sim])
		
		for locus in range(nLoci):
			theta_locus = theta[locus]*N[sim]*scalar_N[locus]/rescale
			rho_locus = rho[locus]*N[sim]*scalar_N[locus]/rescale
			if theta_locus < 0.00001:
				theta_locus = 0.00001
			print("{0}\t{1}\t{2}\t{3}\t{4}\t{5}".format(nsamA[locus], theta_locus, rho_locus, L[locus], Tdem[sim]/(4.0*N[sim]), Npast[sim]/(1.0*N[sim])))
	outfile = open("priorfile.txt", "w")
	outfile.write(priorfile)
	outfile.close()

if sys.argv[1] == "Contraction_1N":
	# msnsam tbs 10000 -t tbs -r tbs tbs -eN tbs tbs
	# nsamA theta 
	# param multilocus: values that will be printed in priorfile.txt
	## N = N_pop_i / Nref
	theta = [ 4*L[i]*mu[i] for i in range(nLoci) ] # vector of theta/N
	rho = [ 4*L[i]*rec[i] for i in range(nLoci) ] # vector of rho/N
	
	# param monolocus: values that will be read by ms
	priorfile = "N\tNpast\tTdem\n"
	for sim in range(nMultilocus):
		priorfile += "{0:.5f}\t{1:.5f}\t{2:.5f}\n".format(N[sim], Npast[sim], Tdem[sim])
		for locus in range(nLoci):
			theta_locus = theta[locus]*N[sim]
			rho_locus = rho[locus]*N[sim]
			if theta_locus < 0.00001:
				theta_locus = 0.00001
			#print("{0}\t{1}\t{2}\t{3}\t{4}\t{5}".format(nsamA[locus], theta[locus]*N[sim], rho[locus]*N[sim], L[locus], T[sim], Nanc[sim]))
			print("{0}\t{1}\t{2}\t{3}\t{4}\t{5}".format(nsamA[locus], theta_locus, rho_locus, L[locus], Tdem[sim]/(4.0*N[sim]), Npast[sim]/(1.0*N[sim])))
	outfile = open("priorfile.txt", "w")
	outfile.write(priorfile)
	outfile.close()

if sys.argv[1] == "Contraction_2N":
	# msnsam tbs 10000 -t tbs -r tbs tbs -eN tbs tbs
	# nsamA theta 
	# param multilocus: values that will be printed in priorfile.txt
	## N = N_pop_i / Nref
	theta = [ 4*L[i]*mu[i] for i in range(nLoci) ] # vector of theta/N
	rho = [ 4*L[i]*rec[i] for i in range(nLoci) ] # vector of rho/N
	
	# param monolocus: values that will be read by ms
	priorfile = "N\tNpast\tshape_N_a\tshape_N_b\tTdem\n"
	for sim in range(nMultilocus):
		priorfile += "{0:.5f}\t{1:.5f}\t{2:.5f}\t{3:.5f}\t{4:.5f}\n".format(N[sim], Npast[sim], shape_N_a[sim], shape_N_b[sim], Tdem[sim])
                scalar_N = beta(shape_N_a[sim], shape_N_b[sim], size=nLoci)
                rescale = shape_N_a[sim] / (shape_N_a[sim] + shape_N_b[sim])
		
		for locus in range(nLoci):
			theta_locus = theta[locus]*N[sim]*scalar_N[locus]/rescale
			rho_locus = rho[locus]*N[sim]*scalar_N[locus]/rescale
			if theta_locus < 0.00001:
				theta_locus = 0.00001
			print("{0}\t{1}\t{2}\t{3}\t{4}\t{5}".format(nsamA[locus], theta_locus, rho_locus, L[locus], Tdem[sim]/(4.0*N[sim]), Npast[sim]/(1.0*N[sim])))
	outfile = open("priorfile.txt", "w")
	outfile.write(priorfile)
	outfile.close()
