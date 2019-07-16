#!/usr/bin/python
# #!/home/roux/python/Python-2.7.14/python
# -*- coding: utf-8 -*-

import sys
from numpy.random import randint
from numpy.random import uniform
from numpy.random import binomial
from numpy.random import beta
from random import shuffle
help = "\t\033[1;31;40mTakes one model specifier, a number of multilocus simulations and a config.yaml file containing prior boundaries as arguments:\033[0m\n\t\t"
help += "\n\t\t".join(["SC_1M_1N", "SC_1M_2N", "SC_2M_1N", "SC_2M_2N", "AM_1M_1N", "AM_1M_2N", "AM_2M_1N", "AM_2M_2N", "IM_1M_1N", "IM_1M_2N", "IM_2M_1N", "IM_2M_2N", "SI_1N", "SI_2N"])
help += "\n\n"
help += "\t\033[1;32;40m#SI\033[0m\n\tmsnsam tbs 10000 -t tbs -r tbs tbs -I 2 tbs tbs 0 -ej tbs 2 1 -eN tbs tbs -g 1 tbs -g 2 tbs\n"
help += "\t\033[1;32;40m#AM\033[0m\n\tmsnsam tbs 10000 -t tbs -r tbs tbs -I 2 tbs tbs 0 -ema tbs 2 0 tbs tbs 0 -ej tbs 2 1 -eN tbs tbs -g 1 tbs -g 2 tbs\n"
help += "\t\033[1;32;40m#IM\033[0m\n\tmsnsam tbs 10000 -t tbs -r tbs tbs -I 2 tbs tbs 0 -m 1 2 tbs -m 2 1 tbs -ej tbs 2 1 -eN tbs tbs -g 1 tbs -g 2 tbs\n"
help += "\t\033[1;32;40m#SC\033[0m\n\tmsnsam tbs 10000 -t tbs -r tbs tbs -I 2 tbs tbs 0 -m 1 2 tbs -m 2 1 tbs -eM tbs 0 -ej tbs 2 1 -eN tbs tbs -g 1 tbs -g 2 tbs\n"
help += "\t\033[1;32;40mExample: ./priorgen_gof_2pop.py SC_2M_2N 1000 posterior_file config_yaml\033[0m\n"

if len(sys.argv) != 5:
	print(help)
	sys.exit()

# Configuration of the prior distribution
nMultilocus = int(sys.argv[2])


# read bpfile
infile = open("bpfile", "r")
tmp = infile.readline()
L = infile.readline().strip().split("\t")
nsamA = infile.readline().strip().split("\t")
nsamB = infile.readline().strip().split("\t")
theta = infile.readline().strip().split("\t")
rho = infile.readline().strip().split("\t")
infile.close()

# number of loci
nLoci = len(L)

# sum of nsamA + nsamB
nsam_tot = [ int(nsamA[i]) + int(nsamB[i]) for i in range(nLoci) ]

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
used_posterior = randint(cnt, size=nMultilocus)

# read the yaml
config_yaml = open(sys.argv[4], 'r')
for i in config_yaml:
	i = i.strip().split(':')
	if(i[0] == 'population_growth'): # =='constant' or 'variable'
		pop_growth = i[1]
config_yaml.close()


def alpha(Npresent, Nancestral, Tsplit, pop_growth):
	# return alpha, the population growth rate used in N(t) = N0.exp(-alpha x t)
	if pop_growth == 'variable'
		num = log(Npresent/(1.0*Nancestral))
		denom = -1.0*Tsplit
		return(num/denom)
	else:
		return(0)


if sys.argv[1] == "SC_1M_1N":
	# secondary contact
	# ms tbs 10000 -t tbs -r tbs tbs -I 2 tbs tbs 0 -m 1 2 tbs -m 2 1 tbs -n 1 tbs -n 2 tbs -eM tbs 0 -ej tbs 2 1 -eN tbs tbs
	# nsamtot theta rho L nsamA nsamB M12 M21 N1 N2 Tsc Tsplit Tsplit Na

	# param multilocus: values that will be printed in priorfile.txt
	N1 = [ posterior['N1'][i] for i in used_posterior ]
	N2 = [ posterior['N2'][i] for i in used_posterior ]
	Na = [ posterior['Na'][i] for i in used_posterior ]

	## Miration rates
	M12 = [ posterior['M12'][i] for i in used_posterior ]
	M21 = [ posterior['M21'][i] for i in used_posterior ]

	## times
	Tsplit = [ posterior['Tsplit'][i] for i in used_posterior ]
	Tsc = [ posterior['Tsc'][i] if posterior['Tsc'][i]<posterior['Tsplit'][i] else posterior['Tsplit'][i] for i in used_posterior ]

	## alpha
	alpha_1 = [ alpha(N1[i], Na[i], Tsplit[i], pop_growth) for i in range(nMultilocus) ]
	alpha_2 = [ alpha(N2[i], Na[i], Tsplit[i], pop_growth) for i in range(nMultilocus) ]

	# param monolocus: values that will be read by ms
	priorfile = "N1\tN2\tNa\tTsplit\tTsc\tM12\tM21\n"
	for sim in range(nMultilocus):
		priorfile += "{0:.5f}\t{1:.5f}\t{2:.5f}\t{3:.5f}\t{4:.5f}\t{5:.5f}\t{6:.5f}\n".format(N1[sim], N2[sim], Na[sim], Tsplit[sim], Tsc[sim], M12[sim], M21[sim])
		
		for locus in range(nLoci):
			print("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6:.5f}\t{7:.5f}\t{8:.5f}\t{9:.5f}\t{10:.5f}\t{11:.5f}\t{12:.5f}\t{13:.5f}".format(nsam_tot[locus], theta[locus], rho[locus], L[locus], nsamA[locus], nsamB[locus], M12[sim], M21[sim], Tsc[sim], Tsplit[sim], Tsplit[sim], Na[sim], alpha_1[sim], alpha_2[sim]))

	outfile = open("priorfile.txt", "w")
	outfile.write(priorfile)
	outfile.close()


if sys.argv[1] == "SC_1M_2N":
	# param multilocus: values that will be printed in priorfile.txt
	N1 = [ posterior['N1'][i] for i in used_posterior ]
	N2 = [ posterior['N2'][i] for i in used_posterior ]
	Na = [ posterior['Na'][i] for i in used_posterior ]

	## Miration rates
	M12 = [ posterior['M12'][i] for i in used_posterior ]
	M21 = [ posterior['M21'][i] for i in used_posterior ]

	## times
	Tsplit = [ posterior['Tsplit'][i] for i in used_posterior ]
	Tsc = [ posterior['Tsc'][i] if posterior['Tsc'][i]<posterior['Tsplit'][i] else posterior['Tsplit'][i] for i in used_posterior ]

	## alpha
	alpha_1 = [ alpha(N1[i], Na[i], Tsplit[i], pop_growth) for i in range(nMultilocus) ]
	alpha_2 = [ alpha(N2[i], Na[i], Tsplit[i], pop_growth) for i in range(nMultilocus) ]

	## factor of local reduction in Ne. Model of "background selection"
	shape_N_a = [ posterior['shape_N_a'][i] for i in used_posterior ]
	shape_N_b = [ posterior['shape_N_b'][i] for i in used_posterior ]

	# param monolocus: values that will be read by ms
	priorfile = "N1\tN2\tNa\tshape_N_a\tshape_N_b\tTsplit\tTsc\tM12\tM21\n"
	for sim in range(nMultilocus):
		priorfile += "{0:.5f}\t{1:.5f}\t{2:.5f}\t{3:.5f}\t{4:.5f}\t{5:.5f}\t{6:.5f}\t{7:.5f}\t{8:.5f}\n".format(N1[sim], N2[sim], Na[sim], shape_N_a[sim], shape_N_b[sim], Tsplit[sim], Tsc[sim], M12[sim], M21[sim])
		# vectors of size 'nLoci' containing parameters
		scalar_N = beta(shape_N_a[sim], shape_N_b[sim], size=nLoci)
		Na_vec = [ Na[sim]*i for i in scalar_N ]
		
		for locus in range(nLoci):
			print("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6:.5f}\t{7:.5f}\t{8:.5f}\t{9:.5f}\t{10:.5f}\t{11:.5f}\t{12:.5f}\t{13:.5f}".format(nsam_tot[locus], theta[locus], rho[locus], L[locus], nsamA[locus], nsamB[locus], M12[sim], M21[sim], Tsc[sim], Tsplit[sim], Tsplit[sim], Na_vec[locus], alpha_1[sim], alpha_2[sim]))

	outfile = open("priorfile.txt", "w")
	outfile.write(priorfile)
	outfile.close()


if sys.argv[1] == "SC_2M_1N":
	# param multilocus: values that will be printed in priorfile.txt
	N1 = [ posterior['N1'][i] for i in used_posterior ]
	N2 = [ posterior['N2'][i] for i in used_posterior ]
	Na = [ posterior['Na'][i] for i in used_posterior ]

	## Miration rates
	M12 = [ posterior['M12'][i] for i in used_posterior ]
	M21 = [ posterior['M21'][i] for i in used_posterior ]

	## times
	Tsplit = [ posterior['Tsplit'][i] for i in used_posterior ]
	Tsc = [ posterior['Tsc'][i] if posterior['Tsc'][i]<posterior['Tsplit'][i] else posterior['Tsplit'][i] for i in used_posterior ]

	## alpha
	alpha_1 = [ alpha(N1[i], Na[i], Tsplit[i], pop_growth) for i in range(nMultilocus) ]
	alpha_2 = [ alpha(N2[i], Na[i], Tsplit[i], pop_growth) for i in range(nMultilocus) ]

	## factor of local reduction in Ne. Model of "background selection"
	shape_M12_a = [ posterior['shape_M12_a'][i] for i in used_posterior ]
	shape_M12_b = [ posterior['shape_M12_b'][i] for i in used_posterior ]
	shape_M21_a = [ posterior['shape_M21_a'][i] for i in used_posterior ]
	shape_M21_b = [ posterior['shape_M21_b'][i] for i in used_posterior ]
	
	
	# param monolocus: values that will be read by ms
	priorfile = "N1\tN2\tNa\tTsplit\tTsc\tM12\tshape_M12_a\tshape_M12_b\tM21\tshape_M21_a\tshape_M21_b\n"
	for sim in range(nMultilocus):
		priorfile += "{0:.5f}\t{1:.5f}\t{2:.5f}\t{3:.5f}\t{4:.5f}\t{5:.5f}\t{6:.5f}\t{7:.5f}\t{8:.5f}\t{9:.5f}\t{10:.5f}\n".format(N1[sim], N2[sim], Na[sim], Tsplit[sim], Tsc[sim], M12[sim], shape_M12_a[sim], shape_M12_b[sim], M21[sim], shape_M21_a[sim], shape_M21_b[sim])
		
		# vectors of size 'nLoci' containing parameters
		scalar_M12 = beta(shape_M12_a[sim], shape_M12_b[sim], size = nLoci)
		scalar_M21 = beta(shape_M21_a[sim], shape_M21_b[sim], size = nLoci)
		M12_vec = [ M12[sim] * i for i in scalar_M12 ]
		M21_vec = [ M21[sim] * i for i in scalar_M21 ]
		
		for locus in range(nLoci):
			print("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6:.5f}\t{7:.5f}\t{8:.5f}\t{9:.5f}\t{10:.5f}\t{11:.5f}\t{12:.5f}\t{13:.5f}".format(nsam_tot[locus], theta[locus], rho[locus], L[locus], nsamA[locus], nsamB[locus], M12_vec[locus], M21_vec[locus], Tsc[sim], Tsplit[sim], Tsplit[sim], Na[sim], alpha_1[sim], alpha_2[sim]))
	
	outfile = open("priorfile.txt", "w")
	outfile.write(priorfile)
	outfile.close()



if sys.argv[1] == "SC_2M_2N":
	# param multilocus: values that will be printed in priorfile.txt
	N1 = [ posterior['N1'][i] for i in used_posterior ]
	N2 = [ posterior['N2'][i] for i in used_posterior ]
	Na = [ posterior['Na'][i] for i in used_posterior ]

	## Miration rates
	M12 = [ posterior['M12'][i] for i in used_posterior ]
	M21 = [ posterior['M21'][i] for i in used_posterior ]

	## times
	Tsplit = [ posterior['Tsplit'][i] for i in used_posterior ]
	Tsc = [ posterior['Tsc'][i] if posterior['Tsc'][i]<posterior['Tsplit'][i] else posterior['Tsplit'][i] for i in used_posterior ]

	## alpha
	alpha_1 = [ alpha(N1[i], Na[i], Tsplit[i], pop_growth) for i in range(nMultilocus) ]
	alpha_2 = [ alpha(N2[i], Na[i], Tsplit[i], pop_growth) for i in range(nMultilocus) ]

	## factor of local reduction in Ne. Model of "background selection"
	shape_M12_a = [ posterior['shape_M12_a'][i] for i in used_posterior ]
	shape_M12_b = [ posterior['shape_M12_b'][i] for i in used_posterior ]
	shape_M21_a = [ posterior['shape_M21_a'][i] for i in used_posterior ]
	shape_M21_b = [ posterior['shape_M21_b'][i] for i in used_posterior ]

	## factor of local reduction in Ne. Model of "background selection"
	shape_N_a = [ posterior['shape_N_a'][i] for i in used_posterior ]
	shape_N_b = [ posterior['shape_N_b'][i] for i in used_posterior ]
	
	# param monolocus: values that will be read by ms
	priorfile = "N1\tN2\tNa\tshape_N_a\tshape_N_b\tTsplit\tTsc\tM12\tshape_M12_a\tshape_M12_b\tM21\tshape_M21_a\tshape_M21_b\n"
	for sim in range(nMultilocus):
		priorfile += "{0:.5f}\t{1:.5f}\t{2:.5f}\t{3:.5f}\t{4:.5f}\t{5:.5f}\t{6:.5f}\t{7:.5f}\t{8:.5f}\t{9:.5f}\t{10:.5f}\t{11:.5f}\t{12:.5f}\n".format(N1[sim], N2[sim], Na[sim], shape_N_a[sim], shape_N_b[sim], Tsplit[sim], Tsc[sim], M12[sim], shape_M12_a[sim], shape_M12_b[sim], M21[sim], shape_M21_a[sim], shape_M21_b[sim])
		# vectors of size 'nLoci' containing parameters
                scalar_N = beta(shape_N_a[sim], shape_N_b[sim], size=nLoci)
                Na_vec = [ Na[sim]*i for i in scalar_N ]

                # vectors of size 'nLoci' containing parameters
                scalar_M12 = beta(shape_M12_a[sim], shape_M12_b[sim], size = nLoci)
                scalar_M21 = beta(shape_M21_a[sim], shape_M21_b[sim], size = nLoci)
                M12_vec = [ M12[sim] * i for i in scalar_M12 ]
                M21_vec = [ M21[sim] * i for i in scalar_M21 ]
		
		for locus in range(nLoci):
			print("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6:.5f}\t{7:.5f}\t{8:.5f}\t{9:.5f}\t{10:.5f}\t{11:.5f}\t{12:.5f}\t{13:.5f}".format(nsam_tot[locus], theta[locus], rho[locus], L[locus], nsamA[locus], nsamB[locus], M12_vec[locus], M21_vec[locus], Tsc[sim], Tsplit[sim], Tsplit[sim], Na_vec[locus], alpha_1[sim], alpha_2[sim]))

	outfile = open("priorfile.txt", "w")
	outfile.write(priorfile)
	outfile.close()



if sys.argv[1] == "AM_1M_1N":
	# param multilocus: values that will be printed in priorfile.txt
	N1 = [ posterior['N1'][i] for i in used_posterior ]
	N2 = [ posterior['N2'][i] for i in used_posterior ]
	Na = [ posterior['Na'][i] for i in used_posterior ]

	## Miration rates
	M12 = [ posterior['M12'][i] for i in used_posterior ]
	M21 = [ posterior['M21'][i] for i in used_posterior ]

	## times
	Tsplit = [ posterior['Tsplit'][i] for i in used_posterior ]
	Tam = [ posterior['Tam'][i] if posterior['Tam'][i]<posterior['Tsplit'][i] else posterior['Tsplit'][i] for i in used_posterior ]

	## alpha
	alpha_1 = [ alpha(N1[i], Na[i], Tsplit[i], pop_growth) for i in range(nMultilocus) ]
	alpha_2 = [ alpha(N2[i], Na[i], Tsplit[i], pop_growth) for i in range(nMultilocus) ]

	# param monolocus: values that will be read by ms
	priorfile = "N1\tN2\tNa\tTsplit\tTam\tM12\tM21\n"
	for sim in range(nMultilocus):
		priorfile += "{0:.5f}\t{1:.5f}\t{2:.5f}\t{3:.5f}\t{4:.5f}\t{5:.5f}\t{6:.5f}\n".format(N1[sim], N2[sim], Na[sim], Tsplit[sim], Tam[sim], M12[sim], M21[sim])
		
		for locus in range(nLoci):
			print("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6:.5f}\t{7:.5f}\t{8:.5f}\t{9:.5f}\t{10:.5f}\t{11:.5f}\t{12:.5f}\t{13:.5f}".format(nsam_tot[locus], theta[locus], rho[locus], L[locus], nsamA[locus], nsamB[locus], Tam[sim], M12[sim], M21[sim], Tsplit[sim], Tsplit[sim], Na[sim], alpha_1[sim], alpha_2[sim]))

	outfile = open("priorfile.txt", "w")
	outfile.write(priorfile)
	outfile.close()



if sys.argv[1] == "AM_1M_2N":
	# param multilocus: values that will be printed in priorfile.txt
	N1 = [ posterior['N1'][i] for i in used_posterior ]
	N2 = [ posterior['N2'][i] for i in used_posterior ]
	Na = [ posterior['Na'][i] for i in used_posterior ]

	## Miration rates
	M12 = [ posterior['M12'][i] for i in used_posterior ]
	M21 = [ posterior['M21'][i] for i in used_posterior ]

	## times
	Tsplit = [ posterior['Tsplit'][i] for i in used_posterior ]
	Tam = [ posterior['Tam'][i] if posterior['Tam'][i]<posterior['Tsplit'][i] else posterior['Tsplit'][i] for i in used_posterior ]

	## alpha
	alpha_1 = [ alpha(N1[i], Na[i], Tsplit[i], pop_growth) for i in range(nMultilocus) ]
	alpha_2 = [ alpha(N2[i], Na[i], Tsplit[i], pop_growth) for i in range(nMultilocus) ]

	## factor of local reduction in Ne. Model of "background selection"
	shape_N_a = [ posterior['shape_N_a'][i] for i in used_posterior ]
	shape_N_b = [ posterior['shape_N_b'][i] for i in used_posterior ]
	
	# param monolocus: values that will be read by ms
	priorfile = "N1\tN2\tNa\tshape_N_a\tshape_N_b\tTsplit\tTam\tM12\tM21\n"
	for sim in range(nMultilocus):
		priorfile += "{0:.5f}\t{1:.5f}\t{2:.5f}\t{3:.5f}\t{4:.5f}\t{5:.5f}\t{6:.5f}\t{7:.5f}\t{8:.5f}\n".format(N1[sim], N2[sim], Na[sim], shape_N_a[sim], shape_N_b[sim], Tsplit[sim], Tam[sim], M12[sim], M21[sim])
		# vectors of size 'nLoci' containing parameters
                scalar_N = beta(shape_N_a[sim], shape_N_b[sim], size=nLoci)
                Na_vec = [ Na[sim]*i for i in scalar_N ]
		
		for locus in range(nLoci):
			print("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6:.5f}\t{7:.5f}\t{8:.5f}\t{9:.5f}\t{10:.5f}\t{11:.5f}\t{12:.5f}\t{13:.5f}".format(nsam_tot[locus], theta[locus], rho[locus], L[locus], nsamA[locus], nsamB[locus], Tam[sim], M12[sim], M21[sim], Tsplit[sim], Tsplit[sim], Na_vec[locus], alpha_1[sim], alpha_2[sim]))

	outfile = open("priorfile.txt", "w")
	outfile.write(priorfile)
	outfile.close()



if sys.argv[1] == "AM_2M_1N":
	# param multilocus: values that will be printed in priorfile.txt
	N1 = [ posterior['N1'][i] for i in used_posterior ]
	N2 = [ posterior['N2'][i] for i in used_posterior ]
	Na = [ posterior['Na'][i] for i in used_posterior ]

	## Miration rates
	M12 = [ posterior['M12'][i] for i in used_posterior ]
	M21 = [ posterior['M21'][i] for i in used_posterior ]

	## times
	Tsplit = [ posterior['Tsplit'][i] for i in used_posterior ]
	Tam = [ posterior['Tam'][i] if posterior['Tam'][i]<posterior['Tsplit'][i] else posterior['Tsplit'][i] for i in used_posterior ]

	## alpha
	alpha_1 = [ alpha(N1[i], Na[i], Tsplit[i], pop_growth) for i in range(nMultilocus) ]
	alpha_2 = [ alpha(N2[i], Na[i], Tsplit[i], pop_growth) for i in range(nMultilocus) ]

	## factor of local reduction in Ne. Model of "background selection"
	shape_M12_a = [ posterior['shape_M12_a'][i] for i in used_posterior ]
	shape_M12_b = [ posterior['shape_M12_b'][i] for i in used_posterior ]
	shape_M21_a = [ posterior['shape_M21_a'][i] for i in used_posterior ]
	shape_M21_b = [ posterior['shape_M21_b'][i] for i in used_posterior ]

	# param monolocus: values that will be read by ms
	priorfile = "N1\tN2\tNa\tTsplit\tTam\tM12\tshape_M12_a\tshape_M12_b\tM21\tshape_M21_a\tshape_M21_b\n"
	for sim in range(nMultilocus):
		priorfile += "{0:.5f}\t{1:.5f}\t{2:.5f}\t{3:.5f}\t{4:.5f}\t{5:.5f}\t{6:.5f}\t{7:.5f}\t{8:.5f}\t{9:.5f}\t{10:.5f}\n".format(N1[sim], N2[sim], Na[sim], Tsplit[sim], Tam[sim], M12[sim], shape_M12_a[sim], shape_M12_b[sim], M21[sim], shape_M21_a[sim], shape_M21_b[sim])
		
		# vectors of size 'nLoci' containing parameters
                scalar_M12 = beta(shape_M12_a[sim], shape_M12_b[sim], size = nLoci)
                scalar_M21 = beta(shape_M21_a[sim], shape_M21_b[sim], size = nLoci)
                M12_vec = [ M12[sim] * i for i in scalar_M12 ]
                M21_vec = [ M21[sim] * i for i in scalar_M21 ]
	
		for locus in range(nLoci):
			print("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6:.5f}\t{7:.5f}\t{8:.5f}\t{9:.5f}\t{10:.5f}\t{11:.5f}\t{12:.5f}\t{13:.5f}".format(nsam_tot[locus], theta[locus], rho[locus], L[locus], nsamA[locus], nsamB[locus], Tam[sim], M12_vec[locus], M21_vec[locus], Tsplit[sim], Tsplit[sim], Na[sim], alpha_1[sim], alpha_2[sim]))

	outfile = open("priorfile.txt", "w")
	outfile.write(priorfile)
	outfile.close()



if sys.argv[1] == "AM_2M_2N":
	# param multilocus: values that will be printed in priorfile.txt
	N1 = [ posterior['N1'][i] for i in used_posterior ]
	N2 = [ posterior['N2'][i] for i in used_posterior ]
	Na = [ posterior['Na'][i] for i in used_posterior ]

	## Miration rates
	M12 = [ posterior['M12'][i] for i in used_posterior ]
	M21 = [ posterior['M21'][i] for i in used_posterior ]

	## times
	Tsplit = [ posterior['Tsplit'][i] for i in used_posterior ]
	Tam = [ posterior['Tam'][i] if posterior['Tam'][i]<posterior['Tsplit'][i] else posterior['Tsplit'][i] for i in used_posterior ]

	## alpha
	alpha_1 = [ alpha(N1[i], Na[i], Tsplit[i], pop_growth) for i in range(nMultilocus) ]
	alpha_2 = [ alpha(N2[i], Na[i], Tsplit[i], pop_growth) for i in range(nMultilocus) ]

	## factor of local reduction in Ne. Model of "background selection"
	shape_M12_a = [ posterior['shape_M12_a'][i] for i in used_posterior ]
	shape_M12_b = [ posterior['shape_M12_b'][i] for i in used_posterior ]
	shape_M21_a = [ posterior['shape_M21_a'][i] for i in used_posterior ]
	shape_M21_b = [ posterior['shape_M21_b'][i] for i in used_posterior ]

	## factor of local reduction in Ne. Model of "background selection"
	shape_N_a = [ posterior['shape_N_a'][i] for i in used_posterior ]
	shape_N_b = [ posterior['shape_N_b'][i] for i in used_posterior ]

	# param monolocus: values that will be read by ms
	priorfile = "N1\tN2\tNa\tshape_N_a\tshape_N_b\tTsplit\tTam\tM12\tshape_M12_a\tshape_M12_b\tM21\tshape_M21_a\tshape_M21_b\n"
	for sim in range(nMultilocus):
		priorfile += "{0:.5f}\t{1:.5f}\t{2:.5f}\t{3:.5f}\t{4:.5f}\t{5:.5f}\t{6:.5f}\t{7:.5f}\t{8:.5f}\t{9:.5f}\t{10:.5f}\t{11:.5f}\t{12:.5f}\n".format(N1[sim], N2[sim], Na[sim], shape_N_a[sim], shape_N_b[sim], Tsplit[sim], Tam[sim], M12[sim], shape_M12_a[sim], shape_M12_b[sim], M21[sim], shape_M21_a[sim], shape_M21_b[sim])
		# vectors of size 'nLoci' containing parameters
                scalar_N = beta(shape_N_a[sim], shape_N_b[sim], size=nLoci)
                Na_vec = [ Na[sim]*i for i in scalar_N ]
                scalar_M12 = beta(shape_M12_a[sim], shape_M12_b[sim], size = nLoci)
                scalar_M21 = beta(shape_M21_a[sim], shape_M21_b[sim], size = nLoci)
                M12_vec = [ M12[sim] * i for i in scalar_M12 ]
                M21_vec = [ M21[sim] * i for i in scalar_M21 ]
	
		for locus in range(nLoci):
			print("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6:.5f}\t{7:.5f}\t{8:.5f}\t{9:.5f}\t{10:.5f}\t{11:.5f}\t{12:.5f}\t{13:.5f}".format(nsam_tot[locus], theta[locus], rho[locus], L[locus], nsamA[locus], nsamB[locus], Tam[sim], M12_vec[locus], M21_vec[locus], Tsplit[sim], Tsplit[sim], Na_vec[locus], alpha_1[sim], alpha_2[sim]))

	outfile = open("priorfile.txt", "w")
	outfile.write(priorfile)
	outfile.close()



if sys.argv[1] == "IM_1M_1N":
	# param multilocus: values that will be printed in priorfile.txt
	N1 = [ posterior['N1'][i] for i in used_posterior ]
	N2 = [ posterior['N2'][i] for i in used_posterior ]
	Na = [ posterior['Na'][i] for i in used_posterior ]

	## Miration rates
	M12 = [ posterior['M12'][i] for i in used_posterior ]
	M21 = [ posterior['M21'][i] for i in used_posterior ]

	## times
	Tsplit = [ posterior['Tsplit'][i] for i in used_posterior ]

	## alpha
	alpha_1 = [ alpha(N1[i], Na[i], Tsplit[i], pop_growth) for i in range(nMultilocus) ]
	alpha_2 = [ alpha(N2[i], Na[i], Tsplit[i], pop_growth) for i in range(nMultilocus) ]

	# param monolocus: values that will be read by ms
	priorfile = "N1\tN2\tNa\tTsplit\tM12\tM21\n"
	for sim in range(nMultilocus):
		priorfile += "{0:.5f}\t{1:.5f}\t{2:.5f}\t{3:.5f}\t{4:.5f}\t{5:.5f}\n".format(N1[sim], N2[sim], Na[sim], Tsplit[sim], M12[sim], M21[sim])
		
		for locus in range(nLoci):
			# SC print("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6:.5f}\t{7:.5f}\t{8:.5f}\t{9:.5f}\t{10:.5f}\t{11:.5f}\t{12:.5f}\t{13:.5f}".format(nsam_tot[locus], theta[locus], rho[locus], L[locus], nsamA[locus], nsamB[locus], M12_vec[locus], M21_vec[locus], N1_vec[locus], N2_vec[locus], Tsc[sim], Tsplit[sim], Tsplit[sim], Na_vec[locus]))
			print("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6:.5f}\t{7:.5f}\t{8:.5f}\t{9:.5f}\t{10:.5f}\t{11:.5f}\t{12:.5f}".format(nsam_tot[locus], theta[locus], rho[locus], L[locus], nsamA[locus], nsamB[locus], M12[sim], M21[sim], Tsplit[sim], Tsplit[sim], Na[sim], alpha_1[sim], alpha_2[sim]))

	outfile = open("priorfile.txt", "w")
	outfile.write(priorfile)
	outfile.close()



if sys.argv[1] == "IM_1M_2N":
	# param multilocus: values that will be printed in priorfile.txt
	N1 = [ posterior['N1'][i] for i in used_posterior ]
	N2 = [ posterior['N2'][i] for i in used_posterior ]
	Na = [ posterior['Na'][i] for i in used_posterior ]

	## Miration rates
	M12 = [ posterior['M12'][i] for i in used_posterior ]
	M21 = [ posterior['M21'][i] for i in used_posterior ]

	## times
	Tsplit = [ posterior['Tsplit'][i] for i in used_posterior ]

	## alpha
	alpha_1 = [ alpha(N1[i], Na[i], Tsplit[i], pop_growth) for i in range(nMultilocus) ]
	alpha_2 = [ alpha(N2[i], Na[i], Tsplit[i], pop_growth) for i in range(nMultilocus) ]

	## factor of local reduction in Ne. Model of "background selection"
	shape_N_a = [ posterior['shape_N_a'][i] for i in used_posterior ]
	shape_N_b = [ posterior['shape_N_b'][i] for i in used_posterior ]
	
	# param monolocus: values that will be read by ms
	priorfile = "N1\tN2\tNa\tshape_N_a\tshape_N_b\tTsplit\tM12\tM21\n"
	for sim in range(nMultilocus):
		priorfile += "{0:.5f}\t{1:.5f}\t{2:.5f}\t{3:.5f}\t{4:.5f}\t{5:.5f}\t{6:.5f}\t{7:.5f}\n".format(N1[sim], N2[sim], Na[sim], shape_N_a[sim], shape_N_b[sim], Tsplit[sim], M12[sim], M21[sim])
		# vectors of size 'nLoci' containing parameters
                scalar_N = beta(shape_N_a[sim], shape_N_b[sim], size=nLoci)
                Na_vec = [ Na[sim]*i for i in scalar_N ]
		
		for locus in range(nLoci):
			# SC print("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6:.5f}\t{7:.5f}\t{8:.5f}\t{9:.5f}\t{10:.5f}\t{11:.5f}\t{12:.5f}\t{13:.5f}".format(nsam_tot[locus], theta[locus], rho[locus], L[locus], nsamA[locus], nsamB[locus], M12_vec[locus], M21_vec[locus], N1_vec[locus], N2_vec[locus], Tsc[sim], Tsplit[sim], Tsplit[sim], Na_vec[locus]))
			print("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6:.5f}\t{7:.5f}\t{8:.5f}\t{9:.5f}\t{10:.5f}\t{11:.5f}\t{12:.5f}".format(nsam_tot[locus], theta[locus], rho[locus], L[locus], nsamA[locus], nsamB[locus], M12[sim], M21[sim], Tsplit[sim], Tsplit[sim], Na_vec[locus], alpha_1[sim], alpha_2[sim]))

	outfile = open("priorfile.txt", "w")
	outfile.write(priorfile)
	outfile.close()



if sys.argv[1] == "IM_2M_1N":
	# param multilocus: values that will be printed in priorfile.txt
	N1 = [ posterior['N1'][i] for i in used_posterior ]
	N2 = [ posterior['N2'][i] for i in used_posterior ]
	Na = [ posterior['Na'][i] for i in used_posterior ]

	## Miration rates
	M12 = [ posterior['M12'][i] for i in used_posterior ]
	M21 = [ posterior['M21'][i] for i in used_posterior ]

	## times
	Tsplit = [ posterior['Tsplit'][i] for i in used_posterior ]

	## alpha
	alpha_1 = [ alpha(N1[i], Na[i], Tsplit[i], pop_growth) for i in range(nMultilocus) ]
	alpha_2 = [ alpha(N2[i], Na[i], Tsplit[i], pop_growth) for i in range(nMultilocus) ]

	## factor of local reduction in Ne. Model of "background selection"
	shape_M12_a = [ posterior['shape_M12_a'][i] for i in used_posterior ]
	shape_M12_b = [ posterior['shape_M12_b'][i] for i in used_posterior ]
	shape_M21_a = [ posterior['shape_M21_a'][i] for i in used_posterior ]
	shape_M21_b = [ posterior['shape_M21_b'][i] for i in used_posterior ]

	# param monolocus: values that will be read by ms
	priorfile = "N1\tN2\tNa\tTsplit\tM12\tshape_M12_a\tshape_M12_b\tM21\tshape_M21_a\tshape_M21_b\n"
	for sim in range(nMultilocus):
		priorfile += "{0:.5f}\t{1:.5f}\t{2:.5f}\t{3:.5f}\t{4:.5f}\t{5}\t{6:.5f}\t{7:.5f}\t{8:.5f}\t{9:.5f}\n".format(N1[sim], N2[sim], Na[sim], Tsplit[sim], M12[sim], shape_M12_a[sim], shape_M12_b[sim], M21[sim], shape_M21_a[sim], shape_M21_b[sim])
		
		# vectors of size 'nLoci' containing parameters
                scalar_M12 = beta(shape_M12_a[sim], shape_M12_b[sim], size = nLoci)
                scalar_M21 = beta(shape_M21_a[sim], shape_M21_b[sim], size = nLoci)
                M12_vec = [ M12[sim] * i for i in scalar_M12 ]
                M21_vec = [ M21[sim] * i for i in scalar_M21 ]

	
		for locus in range(nLoci):
			# SC print("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6:.5f}\t{7:.5f}\t{8:.5f}\t{9:.5f}\t{10:.5f}\t{11:.5f}\t{12:.5f}\t{13:.5f}".format(nsam_tot[locus], theta[locus], rho[locus], L[locus], nsamA[locus], nsamB[locus], M12_vec[locus], M21_vec[locus], N1_vec[locus], N2_vec[locus], Tsc[sim], Tsplit[sim], Tsplit[sim], Na_vec[locus]))
			print("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6:.5f}\t{7:.5f}\t{8:.5f}\t{9:.5f}\t{10:.5f}\t{11:.5f}\t{12:.5f}".format(nsam_tot[locus], theta[locus], rho[locus], L[locus], nsamA[locus], nsamB[locus], M12_vec[locus], M21_vec[locus], Tsplit[sim], Tsplit[sim], Na[sim], alpha_1[sim], alpha_2[sim]))

	outfile = open("priorfile.txt", "w")
	outfile.write(priorfile)
	outfile.close()



if sys.argv[1] == "IM_2M_2N":
	# param multilocus: values that will be printed in priorfile.txt
	N1 = [ posterior['N1'][i] for i in used_posterior ]
	N2 = [ posterior['N2'][i] for i in used_posterior ]
	Na = [ posterior['Na'][i] for i in used_posterior ]

	## Miration rates
	M12 = [ posterior['M12'][i] for i in used_posterior ]
	M21 = [ posterior['M21'][i] for i in used_posterior ]

	## times
	Tsplit = [ posterior['Tsplit'][i] for i in used_posterior ]

	## alpha
	alpha_1 = [ alpha(N1[i], Na[i], Tsplit[i], pop_growth) for i in range(nMultilocus) ]
	alpha_2 = [ alpha(N2[i], Na[i], Tsplit[i], pop_growth) for i in range(nMultilocus) ]

	## factor of local reduction in Ne. Model of "background selection"
	shape_M12_a = [ posterior['shape_M12_a'][i] for i in used_posterior ]
	shape_M12_b = [ posterior['shape_M12_b'][i] for i in used_posterior ]
	shape_M21_a = [ posterior['shape_M21_a'][i] for i in used_posterior ]
	shape_M21_b = [ posterior['shape_M21_b'][i] for i in used_posterior ]

	## factor of local reduction in Ne. Model of "background selection"
	shape_N_a = [ posterior['shape_N_a'][i] for i in used_posterior ]
	shape_N_b = [ posterior['shape_N_b'][i] for i in used_posterior ]

	# param monolocus: values that will be read by ms
	priorfile = "N1\tN2\tNa\tshape_N_a\tshape_N_b\tTsplit\tM12\tshape_M12_a\tshape_M12_b\tM21\tshape_M21_a\tshape_M21_b\n"
	for sim in range(nMultilocus):
		priorfile += "{0:.5f}\t{1:.5f}\t{2:.5f}\t{3:.5f}\t{4:.5f}\t{5:.5f}\t{6:.5f}\t{7:.5f}\t{8:.5f}\t{9:.5f}\t{10:.5f}\t{11:.5f}\n".format(N1[sim], N2[sim], Na[sim], shape_N_a[sim], shape_N_b[sim], Tsplit[sim], M12[sim], shape_M12_a[sim], shape_M12_b[sim], M21[sim], shape_M21_a[sim], shape_M21_b[sim])
		# vectors of size 'nLoci' containing parameters
                scalar_N = beta(shape_N_a[sim], shape_N_b[sim], size=nLoci)
                Na_vec = [ Na[sim]*i for i in scalar_N ]
	
                # vectors of size 'nLoci' containing parameters
                scalar_M12 = beta(shape_M12_a[sim], shape_M12_b[sim], size = nLoci)
                scalar_M21 = beta(shape_M21_a[sim], shape_M21_b[sim], size = nLoci)
                M12_vec = [ M12[sim] * i for i in scalar_M12 ]
                M21_vec = [ M21[sim] * i for i in scalar_M21 ]

		for locus in range(nLoci):
			# SC print("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6:.5f}\t{7:.5f}\t{8:.5f}\t{9:.5f}\t{10:.5f}\t{11:.5f}\t{12:.5f}\t{13:.5f}".format(nsam_tot[locus], theta[locus], rho[locus], L[locus], nsamA[locus], nsamB[locus], M12_vec[locus], M21_vec[locus], N1_vec[locus], N2_vec[locus], Tsc[sim], Tsplit[sim], Tsplit[sim], Na_vec[locus]))
			print("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6:.5f}\t{7:.5f}\t{8:.5f}\t{9:.5f}\t{10:.5f}\t{11:.5f}\t{12:.5f}".format(nsam_tot[locus], theta[locus], rho[locus], L[locus], nsamA[locus], nsamB[locus], M12_vec[locus], M21_vec[locus], Tsplit[sim], Tsplit[sim], Na_vec[locus], alpha_1[sim], alpha_2[sim]))

	outfile = open("priorfile.txt", "w")
	outfile.write(priorfile)
	outfile.close()



if sys.argv[1] == "SI_1N":
	# param multilocus: values that will be printed in priorfile.txt
	N1 = [ posterior['N1'][i] for i in used_posterior ]
	N2 = [ posterior['N2'][i] for i in used_posterior ]
	Na = [ posterior['Na'][i] for i in used_posterior ]

	## times
	Tsplit = [ posterior['Tsplit'][i] for i in used_posterior ]

	## alpha
	alpha_1 = [ alpha(N1[i], Na[i], Tsplit[i], pop_growth) for i in range(nMultilocus) ]
	alpha_2 = [ alpha(N2[i], Na[i], Tsplit[i], pop_growth) for i in range(nMultilocus) ]

	# param monolocus: values that will be read by ms
	priorfile = "N1\tN2\tNa\tTsplit\n"
	for sim in range(nMultilocus):
		priorfile += "{0:.5f}\t{1:.5f}\t{2:.5f}\t{3:.5f}\n".format(N1[sim], N2[sim], Na[sim], Tsplit[sim])
		
		for locus in range(nLoci):
			print("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6:.5f}\t{7:.5f}\t{8:.5f}\t{9:.5f}\t{10:.5f}".format(nsam_tot[locus], theta[locus], rho[locus], L[locus], nsamA[locus], nsamB[locus], Tsplit[sim], Tsplit[sim], Na[sim], alpha_1[sim], alpha_2[sim]))

	outfile = open("priorfile.txt", "w")
	outfile.write(priorfile)
	outfile.close()



if sys.argv[1] == "SI_2N":
	# param multilocus: values that will be printed in priorfile.txt
	N1 = [ posterior['N1'][i] for i in used_posterior ]
	N2 = [ posterior['N2'][i] for i in used_posterior ]
	Na = [ posterior['Na'][i] for i in used_posterior ]

	## times
	Tsplit = [ posterior['Tsplit'][i] for i in used_posterior ]

	## alpha
	alpha_1 = [ alpha(N1[i], Na[i], Tsplit[i], pop_growth) for i in range(nMultilocus) ]
	alpha_2 = [ alpha(N2[i], Na[i], Tsplit[i], pop_growth) for i in range(nMultilocus) ]

	## factor of local reduction in Ne. Model of "background selection"
	shape_N_a = [ posterior['shape_N_a'][i] for i in used_posterior ]
	shape_N_b = [ posterior['shape_N_b'][i] for i in used_posterior ]

	# param monolocus: values that will be read by ms
	priorfile = "N1\tN2\tNa\tshape_N_a\tshape_N_b\tTsplit\n"
	for sim in range(nMultilocus):
		priorfile += "{0:.5f}\t{1:.5f}\t{2:.5f}\t{3:.5f}\t{4:.5f}\t{5:.5f}\n".format(N1[sim], N2[sim], Na[sim], shape_N_a[sim], shape_N_b[sim], Tsplit[sim])
		# vectors of size 'nLoci' containing parameters
                scalar_N = beta(shape_N_a[sim], shape_N_b[sim], size=nLoci)
                Na_vec = [ Na[sim]*i for i in scalar_N ]

	
		for locus in range(nLoci):
			print("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6:.5f}\t{7:.5f}\t{8:.5f}\t{9:.5f}\t{10:.5f}".format(nsam_tot[locus], theta[locus], rho[locus], L[locus], nsamA[locus], nsamB[locus], Tsplit[sim], Tsplit[sim], Na_vec[locus], alpha_1[sim], alpha_2[sim]))

	outfile = open("priorfile.txt", "w")
	outfile.write(priorfile)
	outfile.close()

