#!/usr/bin/python
# #!/home/roux/python/Python-2.7.14/python
# -*- coding: utf-8 -*-

import sys
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
help += "\t\033[1;32;40mExample: ./priorgen_2pop_popGrowth.py SC_2M_2N 1000 config_yaml\033[0m\n"

if len(sys.argv) != 4:
	print(help)
	sys.exit()

# Configuration of the prior distribution
nMultilocus = int(sys.argv[2])

shape_bound = [0.01, 50]
N_bound = [0, 0] # number of diploid individuals in the population
T_bound = [0, 0] # number of generations
M_bound = [0, 0] # 4.N.m , so the number of diploid migrant copies is 2.N.m
config_yaml = open(sys.argv[3], 'r')
for i in config_yaml:
	i = i.strip().split(':')
	if(i[0] == 'N_min'):
		N_bound[0] = float(i[1])
	if(i[0] == 'N_max'):
		N_bound[1] = float(i[1])
	if(i[0] == 'Tsplit_min'):
		T_bound[0] = float(i[1])
	if(i[0] == 'Tsplit_max'):
		T_bound[1] = float(i[1])
	if(i[0] == 'M_min'):
		M_bound[0] = float(i[1])
	if(i[0] == 'M_max'):
		M_bound[1] = float(i[1])
	if(i[0] == 'population_growth'): # =='constant' or 'variable'
		pop_growth = i[1]
config_yaml.close()

# convert parameter values in coalescent units
Nref = (N_bound[1]+N_bound[0])/2.0
N_bound[0] /= Nref
N_bound[1] /= Nref
T_bound[0] /= (4*Nref)
T_bound[1] /= (4*Nref)

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
	# msnsam tbs 10000 -t tbs -r tbs tbs -I 2 tbs tbs 0 -m 1 2 tbs -m 2 1 tbs -eM tbs 0 -ej tbs 2 1 -eN tbs tbs -g 1 tbs -g 2 tbs
	# nsamtot theta rho L nsamA nsamB M12 M21 Tsc Tsplit Tsplit Na alpha_1 alpha_2

	# param multilocus: values that will be printed in priorfile.txt
	## N = N_pop_i / Nref
	N1 = uniform(low = N_bound[0], high = N_bound[1], size = nMultilocus)
	N2 = uniform(low = N_bound[0], high = N_bound[1], size = nMultilocus)
	Na = uniform(low = N_bound[0], high = N_bound[1], size = nMultilocus)
#	Na = [ (N1[i]+N2[i])/2.0 for i in range(len(N1)) ]

	## Miration rates
	M12 = uniform(low = M_bound[0], high = M_bound[1], size = nMultilocus)
	M21 = uniform(low = M_bound[0], high = M_bound[1], size = nMultilocus)

	## times
	Tsplit = uniform(low = T_bound[0], high = T_bound[1], size = nMultilocus)
	Tsc = [ uniform(low = 0, high = Tsplit[i], size = 1)[0] for i in range(nMultilocus) ]

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
	## N = N_pop_i / Nref
	N1 = uniform(low = N_bound[0], high = N_bound[1], size = nMultilocus)
	N2 = uniform(low = N_bound[0], high = N_bound[1], size = nMultilocus)
	Na = uniform(low = N_bound[0], high = N_bound[1], size = nMultilocus)
#       Na = [ (N1[i]+N2[i])/2.0 for i in range(len(N1)) ]

	## Miration rates
	M12 = uniform(low = M_bound[0], high = M_bound[1], size = nMultilocus)
	M21 = uniform(low = M_bound[0], high = M_bound[1], size = nMultilocus)

	## times
	Tsplit = uniform(low = T_bound[0], high = T_bound[1], size = nMultilocus)
	Tsc = [ uniform(low = 0, high = Tsplit[i], size = 1)[0] for i in range(nMultilocus) ]

	## alpha
	alpha_1 = [ alpha(N1[i], Na[i], Tsplit[i], pop_growth) for i in range(nMultilocus) ]
	alpha_2 = [ alpha(N2[i], Na[i], Tsplit[i], pop_growth) for i in range(nMultilocus) ]

	## factor of local reduction in Ne. Model of "background selection"
	shape_N_a = uniform(low = shape_bound[0], high=shape_bound[1], size = nMultilocus)
	shape_N_b = uniform(low = shape_bound[0], high=shape_bound[1], size = nMultilocus)

	# param monolocus: values that will be read by ms
	priorfile = "N1\tN2\tNa\tshape_N_a\tshape_N_b\tTsplit\tTsc\tM12\tM21\n"
	for sim in range(nMultilocus):
		priorfile += "{0:.5f}\t{1:.5f}\t{2:.5f}\t{3:.5f}\t{4:.5f}\t{5:.5f}\t{6:.5f}\t{7:.5f}\t{8:.5f}\n".format(N1[sim], N2[sim], Na[sim], shape_N_a[sim], shape_N_b[sim], Tsplit[sim], Tsc[sim], M12[sim], M21[sim])
		# vectors of size 'nLoci' containing parameters
		scalar_N = beta(shape_N_a[sim], shape_N_b[sim], size=nLoci)
#		N1_vec = [ N1[sim]*i for i in scalar_N ]
#		N2_vec = [ N2[sim]*i for i in scalar_N ]
		Na_vec = [ Na[sim]*i for i in scalar_N ]
		
		for locus in range(nLoci):
			print("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6:.5f}\t{7:.5f}\t{8:.5f}\t{9:.5f}\t{10:.5f}\t{11:.5f}\t{12:.5f}\t{13:.5f}".format(nsam_tot[locus], theta[locus], rho[locus], L[locus], nsamA[locus], nsamB[locus], M12[sim], M21[sim], Tsc[sim], Tsplit[sim], Tsplit[sim], Na_vec[locus], alpha_1[sim], alpha_2[sim]))

	outfile = open("priorfile.txt", "w")
	outfile.write(priorfile)
	outfile.close()


if sys.argv[1] == "SC_2M_1N":
	# param multilocus: values that will be printed in priorfile.txt
	## N = N_pop_i / Nref
	N1 = uniform(low = N_bound[0], high = N_bound[1], size = nMultilocus)
	N2 = uniform(low = N_bound[0], high = N_bound[1], size = nMultilocus)
	Na = uniform(low = N_bound[0], high = N_bound[1], size = nMultilocus)
#       Na = [ (N1[i]+N2[i])/2.0 for i in range(len(N1)) ]

	## Miration rates
	M12 = uniform(low = M_bound[0], high = M_bound[1], size = nMultilocus)
	M21 = uniform(low = M_bound[0], high = M_bound[1], size = nMultilocus)
	
	## factor of variation in M.
	shape_M12_a = uniform(low = shape_bound[0], high=shape_bound[1], size = nMultilocus)
	shape_M12_b = uniform(low = shape_bound[0], high=shape_bound[1], size = nMultilocus)
	shape_M21_a = uniform(low = shape_bound[0], high=shape_bound[1], size = nMultilocus)
	shape_M21_b = uniform(low = shape_bound[0], high=shape_bound[1], size = nMultilocus)

	## times
	Tsplit = uniform(low = T_bound[0], high = T_bound[1], size = nMultilocus)
	Tsc = [ uniform(low = 0, high = Tsplit[i], size = 1)[0] for i in range(nMultilocus) ]

	## alpha
	alpha_1 = [ alpha(N1[i], Na[i], Tsplit[i], pop_growth) for i in range(nMultilocus) ]
	alpha_2 = [ alpha(N2[i], Na[i], Tsplit[i], pop_growth) for i in range(nMultilocus) ]
	
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
	## N = N_pop_i / Nref
	N1 = uniform(low = N_bound[0], high = N_bound[1], size = nMultilocus)
	N2 = uniform(low = N_bound[0], high = N_bound[1], size = nMultilocus)
	Na = uniform(low = N_bound[0], high = N_bound[1], size = nMultilocus)
	#Na = [ (N1[i]+N2[i])/2.0 for i in range(len(N1)) ]

	## Miration rates
	M12 = uniform(low = M_bound[0], high = M_bound[1], size = nMultilocus)
	M21 = uniform(low = M_bound[0], high = M_bound[1], size = nMultilocus)
        ## factor of variation in M.
        shape_M12_a = uniform(low = shape_bound[0], high=shape_bound[1], size = nMultilocus)
        shape_M12_b = uniform(low = shape_bound[0], high=shape_bound[1], size = nMultilocus)
        shape_M21_a = uniform(low = shape_bound[0], high=shape_bound[1], size = nMultilocus)
        shape_M21_b = uniform(low = shape_bound[0], high=shape_bound[1], size = nMultilocus)

	## times
	Tsplit = uniform(low = T_bound[0], high = T_bound[1], size = nMultilocus)
	Tsc = [ uniform(low = 0, high = Tsplit[i], size = 1)[0] for i in range(nMultilocus) ]

	## alpha
	alpha_1 = [ alpha(N1[i], Na[i], Tsplit[i], pop_growth) for i in range(nMultilocus) ]
	alpha_2 = [ alpha(N2[i], Na[i], Tsplit[i], pop_growth) for i in range(nMultilocus) ]

	## bf = factor of local reduction in Ne. Model of "background selection"
        shape_N_a = uniform(low = shape_bound[0], high=shape_bound[1], size = nMultilocus)
        shape_N_b = uniform(low = shape_bound[0], high=shape_bound[1], size = nMultilocus)
	
	
	# param monolocus: values that will be read by ms
	priorfile = "N1\tN2\tNa\tshape_N_a\tshape_N_b\tTsplit\tTsc\tM12\tshape_M12_a\tshape_M12_b\tM21\tshape_M21_a\tshape_M21_b\n"
	for sim in range(nMultilocus):
		priorfile += "{0:.5f}\t{1:.5f}\t{2:.5f}\t{3:.5f}\t{4:.5f}\t{5:.5f}\t{6:.5f}\t{7:.5f}\t{8:.5f}\t{9:.5f}\t{10:.5f}\t{11:.5f}\t{12:.5f}\n".format(N1[sim], N2[sim], Na[sim], shape_N_a[sim], shape_N_b[sim], Tsplit[sim], Tsc[sim], M12[sim], shape_M12_a[sim], shape_M12_b[sim], M21[sim], shape_M21_a[sim], shape_M21_b[sim])
		# vectors of size 'nLoci' containing parameters
                scalar_N = beta(shape_N_a[sim], shape_N_b[sim], size=nLoci)
             #   N1_vec = [ N1[sim]*i for i in scalar_N ]
             #   N2_vec = [ N2[sim]*i for i in scalar_N ]
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
	## N = N_pop_i / Nref
	N1 = uniform(low = N_bound[0], high = N_bound[1], size = nMultilocus)
	N2 = uniform(low = N_bound[0], high = N_bound[1], size = nMultilocus)
	Na = uniform(low = N_bound[0], high = N_bound[1], size = nMultilocus)
        #Na = [ (N1[i]+N2[i])/2.0 for i in range(len(N1)) ]

	## Miration rates
	M12 = uniform(low = M_bound[0], high = M_bound[1], size = nMultilocus)
	M21 = uniform(low = M_bound[0], high = M_bound[1], size = nMultilocus)

	## times
	Tsplit = uniform(low = T_bound[0], high = T_bound[1], size = nMultilocus)
	Tam = [ uniform(low = 0, high = Tsplit[i], size = 1)[0] for i in range(nMultilocus) ]

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
	## N = N_pop_i / Nref
	N1 = uniform(low = N_bound[0], high = N_bound[1], size = nMultilocus)
	N2 = uniform(low = N_bound[0], high = N_bound[1], size = nMultilocus)
	Na = uniform(low = N_bound[0], high = N_bound[1], size = nMultilocus)
	#Na = [ (N1[i]+N2[i])/2.0 for i in range(len(N1)) ]

	## Miration rates
	M12 = uniform(low = M_bound[0], high = M_bound[1], size = nMultilocus)
	M21 = uniform(low = M_bound[0], high = M_bound[1], size = nMultilocus)

	## times
	Tsplit = uniform(low = T_bound[0], high = T_bound[1], size = nMultilocus)
	Tam = [ uniform(low = 0, high = Tsplit[i], size = 1)[0] for i in range(nMultilocus) ]

	## alpha
	alpha_1 = [ alpha(N1[i], Na[i], Tsplit[i], pop_growth) for i in range(nMultilocus) ]
	alpha_2 = [ alpha(N2[i], Na[i], Tsplit[i], pop_growth) for i in range(nMultilocus) ]

	## number of neutral loci
        shape_N_a = uniform(low = shape_bound[0], high=shape_bound[1], size = nMultilocus)
        shape_N_b = uniform(low = shape_bound[0], high=shape_bound[1], size = nMultilocus)

	# param monolocus: values that will be read by ms
	priorfile = "N1\tN2\tNa\tshape_N_a\tshape_N_b\tTsplit\tTam\tM12\tM21\n"
	for sim in range(nMultilocus):
		priorfile += "{0:.5f}\t{1:.5f}\t{2:.5f}\t{3:.5f}\t{4:.5f}\t{5:.5f}\t{6:.5f}\t{7:.5f}\t{8:.5f}\n".format(N1[sim], N2[sim], Na[sim], shape_N_a[sim], shape_N_b[sim], Tsplit[sim], Tam[sim], M12[sim], M21[sim])
		# vectors of size 'nLoci' containing parameters
                scalar_N = beta(shape_N_a[sim], shape_N_b[sim], size=nLoci)
                #N1_vec = [ N1[sim]*i for i in scalar_N ]
                #N2_vec = [ N2[sim]*i for i in scalar_N ]
                Na_vec = [ Na[sim]*i for i in scalar_N ]
		
		for locus in range(nLoci):
			print("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6:.5f}\t{7:.5f}\t{8:.5f}\t{9:.5f}\t{10:.5f}\t{11:.5f}\t{12:.5f}\t{13:.5f}".format(nsam_tot[locus], theta[locus], rho[locus], L[locus], nsamA[locus], nsamB[locus], N1_vec[locus], N2_vec[locus], Tam[sim], M12[sim], M21[sim], Tsplit[sim], Tsplit[sim], Na_vec[locus]))

	outfile = open("priorfile.txt", "w")
	outfile.write(priorfile)
	outfile.close()



if sys.argv[1] == "AM_2M_1N":
	# param multilocus: values that will be printed in priorfile.txt
	## N = N_pop_i / Nref
	N1 = uniform(low = N_bound[0], high = N_bound[1], size = nMultilocus)
	N2 = uniform(low = N_bound[0], high = N_bound[1], size = nMultilocus)
	Na = uniform(low = N_bound[0], high = N_bound[1], size = nMultilocus)
	#Na = [ (N1[i]+N2[i])/2.0 for i in range(len(N1)) ]

	## Miration rates
	M12 = uniform(low = M_bound[0], high = M_bound[1], size = nMultilocus)
	M21 = uniform(low = M_bound[0], high = M_bound[1], size = nMultilocus)

	## times
	Tsplit = uniform(low = T_bound[0], high = T_bound[1], size = nMultilocus)
	Tam = [ uniform(low = 0, high = Tsplit[i], size = 1)[0] for i in range(nMultilocus) ]

	## alpha
	alpha_1 = [ alpha(N1[i], Na[i], Tsplit[i], pop_growth) for i in range(nMultilocus) ]
	alpha_2 = [ alpha(N2[i], Na[i], Tsplit[i], pop_growth) for i in range(nMultilocus) ]

	## number of neutral loci
        shape_M12_a = uniform(low = shape_bound[0], high=shape_bound[1], size = nMultilocus)
        shape_M12_b = uniform(low = shape_bound[0], high=shape_bound[1], size = nMultilocus)
        shape_M21_a = uniform(low = shape_bound[0], high=shape_bound[1], size = nMultilocus)
        shape_M21_b = uniform(low = shape_bound[0], high=shape_bound[1], size = nMultilocus)

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
	## N = N_pop_i / Nref
	N1 = uniform(low = N_bound[0], high = N_bound[1], size = nMultilocus)
	N2 = uniform(low = N_bound[0], high = N_bound[1], size = nMultilocus)
	Na = uniform(low = N_bound[0], high = N_bound[1], size = nMultilocus)
	#Na = [ (N1[i]+N2[i])/2.0 for i in range(len(N1)) ]

	## Miration rates
	M12 = uniform(low = M_bound[0], high = M_bound[1], size = nMultilocus)
	M21 = uniform(low = M_bound[0], high = M_bound[1], size = nMultilocus)
        shape_M12_a = uniform(low = shape_bound[0], high=shape_bound[1], size = nMultilocus)
        shape_M12_b = uniform(low = shape_bound[0], high=shape_bound[1], size = nMultilocus)
        shape_M21_a = uniform(low = shape_bound[0], high=shape_bound[1], size = nMultilocus)
        shape_M21_b = uniform(low = shape_bound[0], high=shape_bound[1], size = nMultilocus)

	## times
	Tsplit = uniform(low = T_bound[0], high = T_bound[1], size = nMultilocus)
	Tam = [ uniform(low = 0, high = Tsplit[i], size = 1)[0] for i in range(nMultilocus) ]

	## alpha
	alpha_1 = [ alpha(N1[i], Na[i], Tsplit[i], pop_growth) for i in range(nMultilocus) ]
	alpha_2 = [ alpha(N2[i], Na[i], Tsplit[i], pop_growth) for i in range(nMultilocus) ]

	## bf = factor of local reduction in Ne. Model of "background selection"
        shape_N_a = uniform(low = shape_bound[0], high=shape_bound[1], size = nMultilocus)
        shape_N_b = uniform(low = shape_bound[0], high=shape_bound[1], size = nMultilocus)

	# param monolocus: values that will be read by ms
	priorfile = "N1\tN2\tNa\tshape_N_a\tshape_N_b\tTsplit\tTam\tM12\tshape_M12_a\tshape_M12_b\tM21\tshape_M21_a\tshape_M21_b\n"
	for sim in range(nMultilocus):
		priorfile += "{0:.5f}\t{1:.5f}\t{2:.5f}\t{3:.5f}\t{4:.5f}\t{5:.5f}\t{6:.5f}\t{7:.5f}\t{8:.5f}\t{9:.5f}\t{10:.5f}\t{11:.5f}\t{12:.5f}\n".format(N1[sim], N2[sim], Na[sim], shape_N_a[sim], shape_N_b[sim], Tsplit[sim], Tam[sim], M12[sim], shape_M12_a[sim], shape_M12_b[sim], M21[sim], shape_M21_a[sim], shape_M21_b[sim])
		# vectors of size 'nLoci' containing parameters
                scalar_N = beta(shape_N_a[sim], shape_N_b[sim], size=nLoci)
                #N1_vec = [ N1[sim]*i for i in scalar_N ]
                #N2_vec = [ N2[sim]*i for i in scalar_N ]
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
	## N = N_pop_i / Nref
	N1 = uniform(low = N_bound[0], high = N_bound[1], size = nMultilocus)
	N2 = uniform(low = N_bound[0], high = N_bound[1], size = nMultilocus)
	Na = uniform(low = N_bound[0], high = N_bound[1], size = nMultilocus)
	#Na = [ (N1[i]+N2[i])/2.0 for i in range(len(N1)) ]

	## Miration rates
	M12 = uniform(low = M_bound[0], high = M_bound[1], size = nMultilocus)
	M21 = uniform(low = M_bound[0], high = M_bound[1], size = nMultilocus)

	## times
	Tsplit = uniform(low = T_bound[0], high = T_bound[1], size = nMultilocus)

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
	## N = N_pop_i / Nref
	N1 = uniform(low = N_bound[0], high = N_bound[1], size = nMultilocus)
	N2 = uniform(low = N_bound[0], high = N_bound[1], size = nMultilocus)
	Na = uniform(low = N_bound[0], high = N_bound[1], size = nMultilocus)
	#Na = [ (N1[i]+N2[i])/2.0 for i in range(len(N1)) ]

	## Miration rates
	M12 = uniform(low = M_bound[0], high = M_bound[1], size = nMultilocus)
	M21 = uniform(low = M_bound[0], high = M_bound[1], size = nMultilocus)

	## times
	Tsplit = uniform(low = T_bound[0], high = T_bound[1], size = nMultilocus)

	## alpha
	alpha_1 = [ alpha(N1[i], Na[i], Tsplit[i], pop_growth) for i in range(nMultilocus) ]
	alpha_2 = [ alpha(N2[i], Na[i], Tsplit[i], pop_growth) for i in range(nMultilocus) ]

	## bf = factor of local reduction in Ne. Model of "background selection"
        shape_N_a = uniform(low = shape_bound[0], high=shape_bound[1], size = nMultilocus)
        shape_N_b = uniform(low = shape_bound[0], high=shape_bound[1], size = nMultilocus)
	
	
	# param monolocus: values that will be read by ms
	priorfile = "N1\tN2\tNa\tshape_N_a\tshape_N_b\tTsplit\tM12\tM21\n"
	for sim in range(nMultilocus):
		priorfile += "{0:.5f}\t{1:.5f}\t{2:.5f}\t{3:.5f}\t{4:.5f}\t{5:.5f}\t{6:.5f}\t{7:.5f}\n".format(N1[sim], N2[sim], Na[sim], shape_N_a[sim], shape_N_b[sim], Tsplit[sim], M12[sim], M21[sim])
		# vectors of size 'nLoci' containing parameters
                scalar_N = beta(shape_N_a[sim], shape_N_b[sim], size=nLoci)
                #N1_vec = [ N1[sim]*i for i in scalar_N ]
                #N2_vec = [ N2[sim]*i for i in scalar_N ]
                Na_vec = [ Na[sim]*i for i in scalar_N ]
		
		for locus in range(nLoci):
			# SC print("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6:.5f}\t{7:.5f}\t{8:.5f}\t{9:.5f}\t{10:.5f}\t{11:.5f}\t{12:.5f}\t{13:.5f}".format(nsam_tot[locus], theta[locus], rho[locus], L[locus], nsamA[locus], nsamB[locus], M12_vec[locus], M21_vec[locus], N1_vec[locus], N2_vec[locus], Tsc[sim], Tsplit[sim], Tsplit[sim], Na_vec[locus]))
			print("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6:.5f}\t{7:.5f}\t{8:.5f}\t{9:.5f}\t{10:.5f}\t{11:.5f}\t{12:.5f}".format(nsam_tot[locus], theta[locus], rho[locus], L[locus], nsamA[locus], nsamB[locus], M12[sim], M21[sim], Tsplit[sim], Tsplit[sim], Na_vec[locus], alpha_1[sim], alpha_2[sim]))

	outfile = open("priorfile.txt", "w")
	outfile.write(priorfile)
	outfile.close()



if sys.argv[1] == "IM_2M_1N":
	# param multilocus: values that will be printed in priorfile.txt
	## N = N_pop_i / Nref
	N1 = uniform(low = N_bound[0], high = N_bound[1], size = nMultilocus)
	N2 = uniform(low = N_bound[0], high = N_bound[1], size = nMultilocus)
	Na = uniform(low = N_bound[0], high = N_bound[1], size = nMultilocus)
	#Na = [ (N1[i]+N2[i])/2.0 for i in range(len(N1)) ]

	## Miration rates
	M12 = uniform(low = M_bound[0], high = M_bound[1], size = nMultilocus)
	M21 = uniform(low = M_bound[0], high = M_bound[1], size = nMultilocus)

	## times
	Tsplit = uniform(low = T_bound[0], high = T_bound[1], size = nMultilocus)

	## alpha
	alpha_1 = [ alpha(N1[i], Na[i], Tsplit[i], pop_growth) for i in range(nMultilocus) ]
	alpha_2 = [ alpha(N2[i], Na[i], Tsplit[i], pop_growth) for i in range(nMultilocus) ]

	## number of neutral loci
        shape_M12_a = uniform(low = shape_bound[0], high=shape_bound[1], size = nMultilocus)
        shape_M12_b = uniform(low = shape_bound[0], high=shape_bound[1], size = nMultilocus)
        shape_M21_a = uniform(low = shape_bound[0], high=shape_bound[1], size = nMultilocus)
        shape_M21_b = uniform(low = shape_bound[0], high=shape_bound[1], size = nMultilocus)

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
	## N = N_pop_i / Nref
	N1 = uniform(low = N_bound[0], high = N_bound[1], size = nMultilocus)
	N2 = uniform(low = N_bound[0], high = N_bound[1], size = nMultilocus)
	Na = uniform(low = N_bound[0], high = N_bound[1], size = nMultilocus)
	#Na = [ (N1[i]+N2[i])/2.0 for i in range(len(N1)) ]

	## Miration rates
	M12 = uniform(low = M_bound[0], high = M_bound[1], size = nMultilocus)
	M21 = uniform(low = M_bound[0], high = M_bound[1], size = nMultilocus)

	## times
	Tsplit = uniform(low = T_bound[0], high = T_bound[1], size = nMultilocus)

	## alpha
	alpha_1 = [ alpha(N1[i], Na[i], Tsplit[i], pop_growth) for i in range(nMultilocus) ]
	alpha_2 = [ alpha(N2[i], Na[i], Tsplit[i], pop_growth) for i in range(nMultilocus) ]

	## bf = factor of local reduction in Ne. Model of "background selection"
        shape_N_a = uniform(low = shape_bound[0], high=shape_bound[1], size = nMultilocus)
        shape_N_b = uniform(low = shape_bound[0], high=shape_bound[1], size = nMultilocus)

	## number of neutral loci
        shape_M12_a = uniform(low = shape_bound[0], high=shape_bound[1], size = nMultilocus)
        shape_M12_b = uniform(low = shape_bound[0], high=shape_bound[1], size = nMultilocus)
        shape_M21_a = uniform(low = shape_bound[0], high=shape_bound[1], size = nMultilocus)
        shape_M21_b = uniform(low = shape_bound[0], high=shape_bound[1], size = nMultilocus)

	# param monolocus: values that will be read by ms
	priorfile = "N1\tN2\tNa\tshape_N_a\tshape_N_b\tTsplit\tM12\tshape_M12_a\tshape_M12_b\tM21\tshape_M21_a\tshape_M21_b\n"
	for sim in range(nMultilocus):
		priorfile += "{0:.5f}\t{1:.5f}\t{2:.5f}\t{3:.5f}\t{4:.5f}\t{5:.5f}\t{6:.5f}\t{7:.5f}\t{8:.5f}\t{9:.5f}\t{10:.5f}\t{11:.5f}\n".format(N1[sim], N2[sim], Na[sim], shape_N_a[sim], shape_N_b[sim], Tsplit[sim], M12[sim], shape_M12_a[sim], shape_M12_b[sim], M21[sim], shape_M21_a[sim], shape_M21_b[sim])
		# vectors of size 'nLoci' containing parameters
                scalar_N = beta(shape_N_a[sim], shape_N_b[sim], size=nLoci)
                #N1_vec = [ N1[sim]*i for i in scalar_N ]
                #N2_vec = [ N2[sim]*i for i in scalar_N ]
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
	## N = N_pop_i / Nref
	N1 = uniform(low = N_bound[0], high = N_bound[1], size = nMultilocus)
	N2 = uniform(low = N_bound[0], high = N_bound[1], size = nMultilocus)
	Na = uniform(low = N_bound[0], high = N_bound[1], size = nMultilocus)
	#Na = [ (N1[i]+N2[i])/2.0 for i in range(len(N1)) ]

	## times
	Tsplit = uniform(low = T_bound[0], high = T_bound[1], size = nMultilocus)

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
	## N = N_pop_i / Nref
	N1 = uniform(low = N_bound[0], high = N_bound[1], size = nMultilocus)
	N2 = uniform(low = N_bound[0], high = N_bound[1], size = nMultilocus)
	Na = uniform(low = N_bound[0], high = N_bound[1], size = nMultilocus)
	#Na = [ (N1[i]+N2[i])/2.0 for i in range(len(N1)) ]

	## times
	Tsplit = uniform(low = T_bound[0], high = T_bound[1], size = nMultilocus)

	## alpha
	alpha_1 = [ alpha(N1[i], Na[i], Tsplit[i], pop_growth) for i in range(nMultilocus) ]
	alpha_2 = [ alpha(N2[i], Na[i], Tsplit[i], pop_growth) for i in range(nMultilocus) ]

	## bf = factor of local reduction in Ne. Model of "background selection"
        shape_N_a = uniform(low = shape_bound[0], high=shape_bound[1], size = nMultilocus)
        shape_N_b = uniform(low = shape_bound[0], high=shape_bound[1], size = nMultilocus)

	# param monolocus: values that will be read by ms
	priorfile = "N1\tN2\tNa\tshape_N_a\tshape_N_b\tTsplit\n"
	for sim in range(nMultilocus):
		priorfile += "{0:.5f}\t{1:.5f}\t{2:.5f}\t{3:.5f}\t{4:.5f}\t{5:.5f}\n".format(N1[sim], N2[sim], Na[sim], shape_N_a[sim], shape_N_b[sim], Tsplit[sim])
		# vectors of size 'nLoci' containing parameters
                scalar_N = beta(shape_N_a[sim], shape_N_b[sim], size=nLoci)
                #N1_vec = [ N1[sim]*i for i in scalar_N ]
                #N2_vec = [ N2[sim]*i for i in scalar_N ]
                Na_vec = [ Na[sim]*i for i in scalar_N ]

	
		for locus in range(nLoci):
			print("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6:.5f}\t{7:.5f}\t{8:.5f}\t{9:.5f}\t{10:.5f}".format(nsam_tot[locus], theta[locus], rho[locus], L[locus], nsamA[locus], nsamB[locus], Tsplit[sim], Tsplit[sim], Na_vec[locus], alpha_1[sim], alpha_2[sim]))

	outfile = open("priorfile.txt", "w")
	outfile.write(priorfile)
	outfile.close()

