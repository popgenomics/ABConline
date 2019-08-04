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
help += "\t\033[1;32;40mExample: ./priorgen_gof_2pop_popGrowth.py SC_2M_2N 1000 posterior_file\033[0m\n"

if len(sys.argv) != 4:
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


if sys.argv[1] == "SC_1M_1N":
	# secondary contact
	# ms tbs 10000 -t tbs -r tbs tbs -I 2 tbs tbs 0 -m 1 2 tbs -m 2 1 tbs -n 1 tbs -n 2 tbs -eM tbs 0 -ej tbs 2 1 -eN tbs tbs
	# nsamtot theta rho L nsamA nsamB M12 M21 N1 N2 Tsc Tsplit Tsplit Na

	# param multilocus: values that will be printed in priorfile.txt
	N1 = [ posterior['N1'][i] for i in used_posterior ]
	N2 = [ posterior['N2'][i] for i in used_posterior ]
	Na = [ posterior['Na'][i] for i in used_posterior ]
	founders1 = [ posterior['founders1'][i] for i in used_posterior ]
	founders2 = [ posterior['founders2'][i] for i in used_posterior ]

	## Miration rates
	M12 = [ posterior['M12'][i] for i in used_posterior ]
	M21 = [ posterior['M21'][i] for i in used_posterior ]

	## times
	Tsplit = [ posterior['Tsplit'][i] for i in used_posterior ]
	Tsc = [ posterior['Tsc'][i] if posterior['Tsc'][i]<posterior['Tsplit'][i] else posterior['Tsplit'][i] for i in used_posterior ]

	Tdem1 = [ posterior['Tdem1'][i] for i in used_posterior ]
	Tdem2 = [ posterior['Tdem2'][i] for i in used_posterior ]

	# param monolocus: values that will be read by ms
	priorfile = "N1\tN2\tNa\tfounders1\tfounders2\tTdem1\tTdem2\tTsplit\tTsc\tM12\tM21\n"
	for sim in range(nMultilocus):
		priorfile += "{0:.5f}\t{1:.5f}\t{2:.5f}\t{3:.5f}\t{4:.5f}\t{5:.5f}\t{6:.5f}\t{7:.5f}\t{8:.5f}\t{9:.5f}\t{10:.5f}\n".format(N1[sim], N2[sim], Na[sim], founders1[sim], founders2[sim], Tdem1[sim], Tdem2[sim], Tsplit[sim], Tsc[sim], M12[sim], M21[sim])
		
		for locus in range(nLoci):
			print("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6:.5f}\t{7:.5f}\t{8:.5f}\t{9:.5f}\t{10:.5f}\t{11:.5f}\t{12:.5f}\t{13:.5f}\t{14:.5f}\t{15:.5f}\t{16:.5f}\t{17:.5f}".format(nsam_tot[locus], theta[locus], rho[locus], L[locus], nsamA[locus], nsamB[locus], M12[sim], M21[sim], N1[sim], N2[sim], Tdem1[sim], founders1[sim]*Na[sim], Tdem2[sim], founders2[sim]*Na[sim], Tsc[sim], Tsplit[sim], Tsplit[sim], Na[sim]))

	outfile = open("priorfile.txt", "w")
	outfile.write(priorfile)
	outfile.close()


if sys.argv[1] == "SC_1M_2N":
	# param multilocus: values that will be printed in priorfile.txt
	N1 = [ posterior['N1'][i] for i in used_posterior ]
	N2 = [ posterior['N2'][i] for i in used_posterior ]
	Na = [ posterior['Na'][i] for i in used_posterior ]
	founders1 = [ posterior['founders1'][i] for i in used_posterior ]
	founders2 = [ posterior['founders2'][i] for i in used_posterior ]

	## Miration rates
	M12 = [ posterior['M12'][i] for i in used_posterior ]
	M21 = [ posterior['M21'][i] for i in used_posterior ]

	## times
	Tsplit = [ posterior['Tsplit'][i] for i in used_posterior ]
	Tsc = [ posterior['Tsc'][i] if posterior['Tsc'][i]<posterior['Tsplit'][i] else posterior['Tsplit'][i] for i in used_posterior ]

	Tdem1 = [ posterior['Tdem1'][i] for i in used_posterior ]
	Tdem2 = [ posterior['Tdem2'][i] for i in used_posterior ]

	## factor of local reduction in Ne. Model of "background selection"
	shape_N_a = [ posterior['shape_N_a'][i] for i in used_posterior ]
	shape_N_b = [ posterior['shape_N_b'][i] for i in used_posterior ]

	# param monolocus: values that will be read by ms
	priorfile = "N1\tN2\tNa\tfounders1\tfounders2\tTdem1\tTdem2\tshape_N_a\tshape_N_b\tTsplit\tTsc\tM12\tM21\n"
	for sim in range(nMultilocus):
		priorfile += "{0:.5f}\t{1:.5f}\t{2:.5f}\t{3:.5f}\t{4:.5f}\t{5:.5f}\t{6:.5f}\t{7:.5f}\t{8:.5f}\t{9:.5f}\t{10:.5f}\t{11:.5f}\t{12:.5f}\n".format(N1[sim], N2[sim], Na[sim], founders1[sim], founders2[sim], Tdem1[sim], Tdem2[sim], shape_N_a[sim], shape_N_b[sim], Tsplit[sim], Tsc[sim], M12[sim], M21[sim])
		# vectors of size 'nLoci' containing parameters
		scalar_N = beta(shape_N_a[sim], shape_N_b[sim], size=nLoci)
		N1_vec = [ N1[sim]*i for i in scalar_N ]
		N2_vec = [ N2[sim]*i for i in scalar_N ]
		Na_vec = [ Na[sim]*i for i in scalar_N ]
		
		for locus in range(nLoci):
			print("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6:.5f}\t{7:.5f}\t{8:.5f}\t{9:.5f}\t{10:.5f}\t{11:.5f}\t{12:.5f}\t{13:.5f}\t{14:.5f}\t{15:.5f}\t{16:.5f}\t{17:.5f}".format(nsam_tot[locus], theta[locus], rho[locus], L[locus], nsamA[locus], nsamB[locus], M12[sim], M21[sim], N1_vec[locus], N2_vec[locus], Tdem1[sim], founders1[sim]*Na_vec[locus], Tdem2[sim], founders2[sim]*Na_vec[locus], Tsc[sim], Tsplit[sim], Tsplit[sim], Na_vec[locus]))

	outfile = open("priorfile.txt", "w")
	outfile.write(priorfile)
	outfile.close()


if sys.argv[1] == "SC_2M_1N":
	# param multilocus: values that will be printed in priorfile.txt
	N1 = [ posterior['N1'][i] for i in used_posterior ]
	N2 = [ posterior['N2'][i] for i in used_posterior ]
	Na = [ posterior['Na'][i] for i in used_posterior ]
	founders1 = [ posterior['founders1'][i] for i in used_posterior ]
	founders2 = [ posterior['founders2'][i] for i in used_posterior ]

	## Miration rates
	M12 = [ posterior['M12'][i] for i in used_posterior ]
	M21 = [ posterior['M21'][i] for i in used_posterior ]

	## times
	Tsplit = [ posterior['Tsplit'][i] for i in used_posterior ]
	Tsc = [ posterior['Tsc'][i] if posterior['Tsc'][i]<posterior['Tsplit'][i] else posterior['Tsplit'][i] for i in used_posterior ]

	Tdem1 = [ posterior['Tdem1'][i] for i in used_posterior ]
	Tdem2 = [ posterior['Tdem2'][i] for i in used_posterior ]

	## factor of local reduction in Ne. Model of "background selection"
	shape_M12_a = [ posterior['shape_M12_a'][i] for i in used_posterior ]
	shape_M12_b = [ posterior['shape_M12_b'][i] for i in used_posterior ]
	shape_M21_a = [ posterior['shape_M21_a'][i] for i in used_posterior ]
	shape_M21_b = [ posterior['shape_M21_b'][i] for i in used_posterior ]
	
	
	# param monolocus: values that will be read by ms
	priorfile = "N1\tN2\tNa\tfounders1\tfounders2\tTdem1\tTdem2\tTsplit\tTsc\tM12\tshape_M12_a\tshape_M12_b\tM21\tshape_M21_a\tshape_M21_b\n"
	for sim in range(nMultilocus):
		priorfile += "{0:.5f}\t{1:.5f}\t{2:.5f}\t{3:.5f}\t{4:.5f}\t{5:.5f}\t{6:.5f}\t{7:.5f}\t{8:.5f}\t{9:.5f}\t{10:.5f}\t{11:.5f}\t{12:.5f}\t{13:.5f}\t{14:.5f}\n".format(N1[sim], N2[sim], Na[sim], founders1[sim], founders2[sim], Tdem1[sim], Tdem2[sim], Tsplit[sim], Tsc[sim], M12[sim], shape_M12_a[sim], shape_M12_b[sim], M21[sim], shape_M21_a[sim], shape_M21_b[sim])
		
		# vectors of size 'nLoci' containing parameters
		scalar_M12 = beta(shape_M12_a[sim], shape_M12_b[sim], size = nLoci)
		scalar_M21 = beta(shape_M21_a[sim], shape_M21_b[sim], size = nLoci)
		M12_vec = [ M12[sim] * i for i in scalar_M12 ]
		M21_vec = [ M21[sim] * i for i in scalar_M21 ]
		
		for locus in range(nLoci):
			print("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6:.5f}\t{7:.5f}\t{8:.5f}\t{9:.5f}\t{10:.5f}\t{11:.5f}\t{12:.5f}\t{13:.5f}\t{14:.5f}\t{15:.5f}\t{16:.5f}\t{17:.5f}".format(nsam_tot[locus], theta[locus], rho[locus], L[locus], nsamA[locus], nsamB[locus], M12_vec[locus], M21_vec[locus], N1[sim], N2[sim], Tdem1[sim], founders1[sim]*Na[sim], Tdem2[sim], founders2[sim]*Na[sim], Tsc[sim], Tsplit[sim], Tsplit[sim], Na[sim]))
	
	outfile = open("priorfile.txt", "w")
	outfile.write(priorfile)
	outfile.close()



if sys.argv[1] == "SC_2M_2N":
	# param multilocus: values that will be printed in priorfile.txt
	N1 = [ posterior['N1'][i] for i in used_posterior ]
	N2 = [ posterior['N2'][i] for i in used_posterior ]
	Na = [ posterior['Na'][i] for i in used_posterior ]
	founders1 = [ posterior['founders1'][i] for i in used_posterior ]
	founders2 = [ posterior['founders2'][i] for i in used_posterior ]

	## Miration rates
	M12 = [ posterior['M12'][i] for i in used_posterior ]
	M21 = [ posterior['M21'][i] for i in used_posterior ]

	## times
	Tsplit = [ posterior['Tsplit'][i] for i in used_posterior ]
	Tsc = [ posterior['Tsc'][i] if posterior['Tsc'][i]<posterior['Tsplit'][i] else posterior['Tsplit'][i] for i in used_posterior ]

	Tdem1 = [ posterior['Tdem1'][i] for i in used_posterior ]
	Tdem2 = [ posterior['Tdem2'][i] for i in used_posterior ]

	## factor of local reduction in Ne. Model of "background selection"
	shape_M12_a = [ posterior['shape_M12_a'][i] for i in used_posterior ]
	shape_M12_b = [ posterior['shape_M12_b'][i] for i in used_posterior ]
	shape_M21_a = [ posterior['shape_M21_a'][i] for i in used_posterior ]
	shape_M21_b = [ posterior['shape_M21_b'][i] for i in used_posterior ]

	## factor of local reduction in Ne. Model of "background selection"
	shape_N_a = [ posterior['shape_N_a'][i] for i in used_posterior ]
	shape_N_b = [ posterior['shape_N_b'][i] for i in used_posterior ]
	
	# param monolocus: values that will be read by ms
	priorfile = "N1\tN2\tNa\tfounders1\tfounders2\tTdem1\tTdem2\tshape_N_a\tshape_N_b\tTsplit\tTsc\tM12\tshape_M12_a\tshape_M12_b\tM21\tshape_M21_a\tshape_M21_b\n"
	for sim in range(nMultilocus):
		priorfile += "{0:.5f}\t{1:.5f}\t{2:.5f}\t{3:.5f}\t{4:.5f}\t{5:.5f}\t{6:.5f}\t{7:.5f}\t{8:.5f}\t{9:.5f}\t{10:.5f}\t{11:.5f}\t{12:.5f}\t{12:.5f}\t{13:.5f}\t{14:.5f}\t{15:.5f}\n".format(N1[sim], N2[sim], Na[sim], founders1[sim], founders2[sim], Tdem1[sim], Tdem2[sim], shape_N_a[sim], shape_N_b[sim], Tsplit[sim], Tsc[sim], M12[sim], shape_M12_a[sim], shape_M12_b[sim], M21[sim], shape_M21_a[sim], shape_M21_b[sim])
		# vectors of size 'nLoci' containing parameters
                scalar_N = beta(shape_N_a[sim], shape_N_b[sim], size=nLoci)
                N1_vec = [ N1[sim]*i for i in scalar_N ]
                N2_vec = [ N2[sim]*i for i in scalar_N ]
                Na_vec = [ Na[sim]*i for i in scalar_N ]

                # vectors of size 'nLoci' containing parameters
                scalar_M12 = beta(shape_M12_a[sim], shape_M12_b[sim], size = nLoci)
                scalar_M21 = beta(shape_M21_a[sim], shape_M21_b[sim], size = nLoci)
                M12_vec = [ M12[sim] * i for i in scalar_M12 ]
                M21_vec = [ M21[sim] * i for i in scalar_M21 ]
		
		for locus in range(nLoci):
			print("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6:.5f}\t{7:.5f}\t{8:.5f}\t{9:.5f}\t{10:.5f}\t{11:.5f}\t{12:.5f}\t{13:.5f}\t{14:.5f}\t{15:.5f}\t{16:.5f}\t{17:.5f}".format(nsam_tot[locus], theta[locus], rho[locus], L[locus], nsamA[locus], nsamB[locus], M12_vec[locus], M21_vec[locus], N1_vec[locus], N2_vec[locus], Tdem1[sim], founders1[sim]*Na_vec[locus], Tdem2[sim], founders2[sim]*Na_vec[locus], Tsc[sim], Tsplit[sim], Tsplit[sim], Na_vec[locus]))

	outfile = open("priorfile.txt", "w")
	outfile.write(priorfile)
	outfile.close()



if sys.argv[1] == "AM_1M_1N":
	# param multilocus: values that will be printed in priorfile.txt
	N1 = [ posterior['N1'][i] for i in used_posterior ]
	N2 = [ posterior['N2'][i] for i in used_posterior ]
	Na = [ posterior['Na'][i] for i in used_posterior ]
	founders1 = [ posterior['founders1'][i] for i in used_posterior ]
	founders2 = [ posterior['founders2'][i] for i in used_posterior ]

	## Miration rates
	M12 = [ posterior['M12'][i] for i in used_posterior ]
	M21 = [ posterior['M21'][i] for i in used_posterior ]

	## times
	Tsplit = [ posterior['Tsplit'][i] for i in used_posterior ]
	Tam = [ posterior['Tam'][i] if posterior['Tam'][i]<posterior['Tsplit'][i] else posterior['Tsplit'][i] for i in used_posterior ]

	Tdem1 = [ posterior['Tdem1'][i] for i in used_posterior ]
	Tdem2 = [ posterior['Tdem2'][i] for i in used_posterior ]

	# param monolocus: values that will be read by ms
	priorfile = "N1\tN2\tNa\tfounders1\tfounders2\tTdem1\tTdem2\tTsplit\tTam\tM12\tM21\n"
	for sim in range(nMultilocus):
		priorfile += "{0:.5f}\t{1:.5f}\t{2:.5f}\t{3:.5f}\t{4:.5f}\t{5:.5f}\t{6:.5f}\t{7:.5f}\t{8:.5f}\t{9:.5f}\t{10:.5f}\n".format(N1[sim], N2[sim], Na[sim], founders1[sim], founders2[sim], Tdem1[sim], Tdem2[sim], Tsplit[sim], Tam[sim], M12[sim], M21[sim])
		
		for locus in range(nLoci):
			print("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6:.5f}\t{7:.5f}\t{8:.5f}\t{9:.5f}\t{10:.5f}\t{11:.5f}\t{12:.5f}\t{13:.5f}\t{14:.5f}\t{15:.5f}\t{16:.5f}\t{17:.5f}".format(nsam_tot[locus], theta[locus], rho[locus], L[locus], nsamA[locus], nsamB[locus], N1[sim], N2[sim], Tdem1[sim], founders1[sim]*Na[sim], Tdem2[sim], founders2[sim]*Na[sim], Tam[sim], M12[sim], M21[sim], Tsplit[sim], Tsplit[sim], Na[sim]))

	outfile = open("priorfile.txt", "w")
	outfile.write(priorfile)
	outfile.close()



if sys.argv[1] == "AM_1M_2N":
	# param multilocus: values that will be printed in priorfile.txt
	N1 = [ posterior['N1'][i] for i in used_posterior ]
	N2 = [ posterior['N2'][i] for i in used_posterior ]
	Na = [ posterior['Na'][i] for i in used_posterior ]
	founders1 = [ posterior['founders1'][i] for i in used_posterior ]
	founders2 = [ posterior['founders2'][i] for i in used_posterior ]

	## Miration rates
	M12 = [ posterior['M12'][i] for i in used_posterior ]
	M21 = [ posterior['M21'][i] for i in used_posterior ]

	## times
	Tsplit = [ posterior['Tsplit'][i] for i in used_posterior ]
	Tam = [ posterior['Tam'][i] if posterior['Tam'][i]<posterior['Tsplit'][i] else posterior['Tsplit'][i] for i in used_posterior ]

	Tdem1 = [ posterior['Tdem1'][i] for i in used_posterior ]
	Tdem2 = [ posterior['Tdem2'][i] for i in used_posterior ]

	## factor of local reduction in Ne. Model of "background selection"
	shape_N_a = [ posterior['shape_N_a'][i] for i in used_posterior ]
	shape_N_b = [ posterior['shape_N_b'][i] for i in used_posterior ]
	
	# param monolocus: values that will be read by ms
	priorfile = "N1\tN2\tNa\tfounders1\tfounders2\tTdem1\tTdem2\tshape_N_a\tshape_N_b\tTsplit\tTam\tM12\tM21\n"
	for sim in range(nMultilocus):
		priorfile += "{0:.5f}\t{1:.5f}\t{2:.5f}\t{3:.5f}\t{4:.5f}\t{5:.5f}\t{6:.5f}\t{7:.5f}\t{8:.5f}\t{9:.5f}\t{10:.5f}\t{11:.5f}\t{12:.5f}\n".format(N1[sim], N2[sim], Na[sim], founders1[sim], founders2[sim], Tdem1[sim], Tdem2[sim], shape_N_a[sim], shape_N_b[sim], Tsplit[sim], Tam[sim], M12[sim], M21[sim])
		# vectors of size 'nLoci' containing parameters
                scalar_N = beta(shape_N_a[sim], shape_N_b[sim], size=nLoci)
                N1_vec = [ N1[sim]*i for i in scalar_N ]
                N2_vec = [ N2[sim]*i for i in scalar_N ]
                Na_vec = [ Na[sim]*i for i in scalar_N ]
		
		for locus in range(nLoci):
			print("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6:.5f}\t{7:.5f}\t{8:.5f}\t{9:.5f}\t{10:.5f}\t{11:.5f}\t{12:.5f}\t{13:.5f}\t{14:.5f}\t{15:.5f}\t{16:.5f}\t{17:.5f}".format(nsam_tot[locus], theta[locus], rho[locus], L[locus], nsamA[locus], nsamB[locus], N1_vec[locus], N2_vec[locus], Tdem1[sim], founders1[sim]*Na_vec[locus], Tdem2[sim], founders2[sim]*Na_vec[locus], Tam[sim], M12[sim], M21[sim], Tsplit[sim], Tsplit[sim], Na_vec[locus]))

	outfile = open("priorfile.txt", "w")
	outfile.write(priorfile)
	outfile.close()



if sys.argv[1] == "AM_2M_1N":
	# param multilocus: values that will be printed in priorfile.txt
	N1 = [ posterior['N1'][i] for i in used_posterior ]
	N2 = [ posterior['N2'][i] for i in used_posterior ]
	Na = [ posterior['Na'][i] for i in used_posterior ]
	founders1 = [ posterior['founders1'][i] for i in used_posterior ]
	founders2 = [ posterior['founders2'][i] for i in used_posterior ]

	## Miration rates
	M12 = [ posterior['M12'][i] for i in used_posterior ]
	M21 = [ posterior['M21'][i] for i in used_posterior ]

	## times
	Tsplit = [ posterior['Tsplit'][i] for i in used_posterior ]
	Tam = [ posterior['Tam'][i] if posterior['Tam'][i]<posterior['Tsplit'][i] else posterior['Tsplit'][i] for i in used_posterior ]

	Tdem1 = [ posterior['Tdem1'][i] for i in used_posterior ]
	Tdem2 = [ posterior['Tdem2'][i] for i in used_posterior ]

	## factor of local reduction in Ne. Model of "background selection"
	shape_M12_a = [ posterior['shape_M12_a'][i] for i in used_posterior ]
	shape_M12_b = [ posterior['shape_M12_b'][i] for i in used_posterior ]
	shape_M21_a = [ posterior['shape_M21_a'][i] for i in used_posterior ]
	shape_M21_b = [ posterior['shape_M21_b'][i] for i in used_posterior ]

	# param monolocus: values that will be read by ms
	priorfile = "N1\tN2\tNa\tfounders1\tfounders2\tTdem1\tTdem2\tTsplit\tTam\tM12\tshape_M12_a\tshape_M12_b\tM21\tshape_M21_a\tshape_M21_b\n"
	for sim in range(nMultilocus):
		priorfile += "{0:.5f}\t{1:.5f}\t{2:.5f}\t{3:.5f}\t{4:.5f}\t{5:.5f}\t{6:.5f}\t{7:.5f}\t{8:.5f}\t{9:.5f}\t{10:.5f}\t{11:.5f}\t{12:.5f}\t{13:.5f}\t{14:.5f}\n".format(N1[sim], N2[sim], Na[sim], founders1[sim], founders2[sim], Tdem1[sim], Tdem2[sim], Tsplit[sim], Tam[sim], M12[sim], shape_M12_a[sim], shape_M12_b[sim], M21[sim], shape_M21_a[sim], shape_M21_b[sim])
		
		# vectors of size 'nLoci' containing parameters
                scalar_M12 = beta(shape_M12_a[sim], shape_M12_b[sim], size = nLoci)
                scalar_M21 = beta(shape_M21_a[sim], shape_M21_b[sim], size = nLoci)
                M12_vec = [ M12[sim] * i for i in scalar_M12 ]
                M21_vec = [ M21[sim] * i for i in scalar_M21 ]
	
		for locus in range(nLoci):
			print("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6:.5f}\t{7:.5f}\t{8:.5f}\t{9:.5f}\t{10:.5f}\t{11:.5f}\t{12:.5f}\t{13:.5f}\t{14:.5f}\t{15:.5f}\t{16:.5f}\t{17:.5f}".format(nsam_tot[locus], theta[locus], rho[locus], L[locus], nsamA[locus], nsamB[locus], N1[sim], N2[sim], Tdem1[sim], founders1[sim]*Na[sim], Tdem2[sim], founders2[sim]*Na[sim], Tam[sim], M12_vec[locus], M21_vec[locus], Tsplit[sim], Tsplit[sim], Na[sim]))

	outfile = open("priorfile.txt", "w")
	outfile.write(priorfile)
	outfile.close()



if sys.argv[1] == "AM_2M_2N":
	# param multilocus: values that will be printed in priorfile.txt
	N1 = [ posterior['N1'][i] for i in used_posterior ]
	N2 = [ posterior['N2'][i] for i in used_posterior ]
	Na = [ posterior['Na'][i] for i in used_posterior ]
	founders1 = [ posterior['founders1'][i] for i in used_posterior ]
	founders2 = [ posterior['founders2'][i] for i in used_posterior ]

	## Miration rates
	M12 = [ posterior['M12'][i] for i in used_posterior ]
	M21 = [ posterior['M21'][i] for i in used_posterior ]

	## times
	Tsplit = [ posterior['Tsplit'][i] for i in used_posterior ]
	Tam = [ posterior['Tam'][i] if posterior['Tam'][i]<posterior['Tsplit'][i] else posterior['Tsplit'][i] for i in used_posterior ]

	Tdem1 = [ posterior['Tdem1'][i] for i in used_posterior ]
	Tdem2 = [ posterior['Tdem2'][i] for i in used_posterior ]

	## factor of local reduction in Ne. Model of "background selection"
	shape_M12_a = [ posterior['shape_M12_a'][i] for i in used_posterior ]
	shape_M12_b = [ posterior['shape_M12_b'][i] for i in used_posterior ]
	shape_M21_a = [ posterior['shape_M21_a'][i] for i in used_posterior ]
	shape_M21_b = [ posterior['shape_M21_b'][i] for i in used_posterior ]

	## factor of local reduction in Ne. Model of "background selection"
	shape_N_a = [ posterior['shape_N_a'][i] for i in used_posterior ]
	shape_N_b = [ posterior['shape_N_b'][i] for i in used_posterior ]

	# param monolocus: values that will be read by ms
	priorfile = "N1\tN2\tNa\tfounders1\tfounders2\tTdem1\tTdem2\tshape_N_a\tshape_N_b\tTsplit\tTam\tM12\tshape_M12_a\tshape_M12_b\tM21\tshape_M21_a\tshape_M21_b\n"
	for sim in range(nMultilocus):
		priorfile += "{0:.5f}\t{1:.5f}\t{2:.5f}\t{3:.5f}\t{4:.5f}\t{5:.5f}\t{6:.5f}\t{7:.5f}\t{8:.5f}\t{9:.5f}\t{10:.5f}\t{11:.5f}\t{12:.5f}\t{13:.5f}\t{14:.5f}\t{15:.5f}\t{16:.5f}\n".format(N1[sim], N2[sim], Na[sim], founders1[sim], founders2[sim], Tdem1[sim], Tdem2[sim], shape_N_a[sim], shape_N_b[sim], Tsplit[sim], Tam[sim], M12[sim], shape_M12_a[sim], shape_M12_b[sim], M21[sim], shape_M21_a[sim], shape_M21_b[sim])
		# vectors of size 'nLoci' containing parameters
                scalar_N = beta(shape_N_a[sim], shape_N_b[sim], size=nLoci)
                N1_vec = [ N1[sim]*i for i in scalar_N ]
                N2_vec = [ N2[sim]*i for i in scalar_N ]
                Na_vec = [ Na[sim]*i for i in scalar_N ]
                scalar_M12 = beta(shape_M12_a[sim], shape_M12_b[sim], size = nLoci)
                scalar_M21 = beta(shape_M21_a[sim], shape_M21_b[sim], size = nLoci)
                M12_vec = [ M12[sim] * i for i in scalar_M12 ]
                M21_vec = [ M21[sim] * i for i in scalar_M21 ]
	
		for locus in range(nLoci):
			print("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6:.5f}\t{7:.5f}\t{8:.5f}\t{9:.5f}\t{10:.5f}\t{11:.5f}\t{12:.5f}\t{13:.5f}\t{14:.5f}\t{15:.5f}\t{16:.5f}\t{17:.5f}".format(nsam_tot[locus], theta[locus], rho[locus], L[locus], nsamA[locus], nsamB[locus], N1_vec[locus], N2_vec[locus], Tdem1[sim], founders1[sim]*Na_vec[locus], Tdem2[sim], founders2[sim]*Na_vec[locus], Tam[sim], M12_vec[locus], M21_vec[locus], Tsplit[sim], Tsplit[sim], Na_vec[locus]))

	outfile = open("priorfile.txt", "w")
	outfile.write(priorfile)
	outfile.close()



if sys.argv[1] == "IM_1M_1N":
	# param multilocus: values that will be printed in priorfile.txt
	N1 = [ posterior['N1'][i] for i in used_posterior ]
	N2 = [ posterior['N2'][i] for i in used_posterior ]
	Na = [ posterior['Na'][i] for i in used_posterior ]
	founders1 = [ posterior['founders1'][i] for i in used_posterior ]
	founders2 = [ posterior['founders2'][i] for i in used_posterior ]

	## Miration rates
	M12 = [ posterior['M12'][i] for i in used_posterior ]
	M21 = [ posterior['M21'][i] for i in used_posterior ]

	## times
	Tsplit = [ posterior['Tsplit'][i] for i in used_posterior ]

	Tdem1 = [ posterior['Tdem1'][i] for i in used_posterior ]
	Tdem2 = [ posterior['Tdem2'][i] for i in used_posterior ]

	# param monolocus: values that will be read by ms
	priorfile = "N1\tN2\tNa\tfounders1\tfounders2\tTdem1\tTdem2\tTsplit\tM12\tM21\n"
	for sim in range(nMultilocus):
		priorfile += "{0:.5f}\t{1:.5f}\t{2:.5f}\t{3:.5f}\t{4:.5f}\t{5:.5f}\t{6:.5f}\t{7:.5f}\t{8:.5f}\t{9:.5f}\n".format(N1[sim], N2[sim], Na[sim], founders1[sim], founders2[sim], Tdem1[sim], Tdem2[sim], Tsplit[sim], M12[sim], M21[sim])
		
		for locus in range(nLoci):
			# SC print("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6:.5f}\t{7:.5f}\t{8:.5f}\t{9:.5f}\t{10:.5f}\t{11:.5f}\t{12:.5f}\t{13:.5f}".format(nsam_tot[locus], theta[locus], rho[locus], L[locus], nsamA[locus], nsamB[locus], M12_vec[locus], M21_vec[locus], N1_vec[locus], N2_vec[locus], Tsc[sim], Tsplit[sim], Tsplit[sim], Na_vec[locus]))
			print("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6:.5f}\t{7:.5f}\t{8:.5f}\t{9:.5f}\t{10:.5f}\t{11:.5f}\t{12:.5f}\t{13:.5f}\t{14:.5f}\t{15:.5f}\t{16:.5f}".format(nsam_tot[locus], theta[locus], rho[locus], L[locus], nsamA[locus], nsamB[locus], N1[sim], N2[sim], Tdem1[sim], founders1[sim]*Na[sim], Tdem2[sim], founders2[sim]*Na[sim], M12[sim], M21[sim], Tsplit[sim], Tsplit[sim], Na[sim]))

	outfile = open("priorfile.txt", "w")
	outfile.write(priorfile)
	outfile.close()



if sys.argv[1] == "IM_1M_2N":
	# param multilocus: values that will be printed in priorfile.txt
	N1 = [ posterior['N1'][i] for i in used_posterior ]
	N2 = [ posterior['N2'][i] for i in used_posterior ]
	Na = [ posterior['Na'][i] for i in used_posterior ]
	founders1 = [ posterior['founders1'][i] for i in used_posterior ]
	founders2 = [ posterior['founders2'][i] for i in used_posterior ]

	## Miration rates
	M12 = [ posterior['M12'][i] for i in used_posterior ]
	M21 = [ posterior['M21'][i] for i in used_posterior ]

	## times
	Tsplit = [ posterior['Tsplit'][i] for i in used_posterior ]

	Tdem1 = [ posterior['Tdem1'][i] for i in used_posterior ]
	Tdem2 = [ posterior['Tdem2'][i] for i in used_posterior ]

	## factor of local reduction in Ne. Model of "background selection"
	shape_N_a = [ posterior['shape_N_a'][i] for i in used_posterior ]
	shape_N_b = [ posterior['shape_N_b'][i] for i in used_posterior ]
	
	# param monolocus: values that will be read by ms
	priorfile = "N1\tN2\tNa\tfounders1\tfounders2\tTdem1\tTdem2\tshape_N_a\tshape_N_b\tTsplit\tM12\tM21\n"
	for sim in range(nMultilocus):
		priorfile += "{0:.5f}\t{1:.5f}\t{2:.5f}\t{3:.5f}\t{4:.5f}\t{5:.5f}\t{6:.5f}\t{7:.5f}\t{8:.5f}\t{9:.5f}\t{10:.5f}\t{11:.5f}\n".format(N1[sim], N2[sim], Na[sim], founders1[sim], founders2[sim], Tdem1[sim], Tdem2[sim], shape_N_a[sim], shape_N_b[sim], Tsplit[sim], M12[sim], M21[sim])
		# vectors of size 'nLoci' containing parameters
                scalar_N = beta(shape_N_a[sim], shape_N_b[sim], size=nLoci)
                N1_vec = [ N1[sim]*i for i in scalar_N ]
                N2_vec = [ N2[sim]*i for i in scalar_N ]
                Na_vec = [ Na[sim]*i for i in scalar_N ]
		
		for locus in range(nLoci):
			# SC print("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6:.5f}\t{7:.5f}\t{8:.5f}\t{9:.5f}\t{10:.5f}\t{11:.5f}\t{12:.5f}\t{13:.5f}".format(nsam_tot[locus], theta[locus], rho[locus], L[locus], nsamA[locus], nsamB[locus], M12_vec[locus], M21_vec[locus], N1_vec[locus], N2_vec[locus], Tsc[sim], Tsplit[sim], Tsplit[sim], Na_vec[locus]))
			print("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6:.5f}\t{7:.5f}\t{8:.5f}\t{9:.5f}\t{10:.5f}\t{11:.5f}\t{12:.5f}\t{13:.5f}\t{14:.5f}\t{15:.5f}\t{16:.5f}".format(nsam_tot[locus], theta[locus], rho[locus], L[locus], nsamA[locus], nsamB[locus], N1_vec[locus], N2_vec[locus], Tdem1[sim], founders1[sim]*Na_vec[locus], Tdem2[sim], founders2[sim]*Na_vec[locus], M12[sim], M21[sim], Tsplit[sim], Tsplit[sim], Na_vec[locus]))

	outfile = open("priorfile.txt", "w")
	outfile.write(priorfile)
	outfile.close()



if sys.argv[1] == "IM_2M_1N":
	# param multilocus: values that will be printed in priorfile.txt
	N1 = [ posterior['N1'][i] for i in used_posterior ]
	N2 = [ posterior['N2'][i] for i in used_posterior ]
	Na = [ posterior['Na'][i] for i in used_posterior ]
	founders1 = [ posterior['founders1'][i] for i in used_posterior ]
	founders2 = [ posterior['founders2'][i] for i in used_posterior ]

	## Miration rates
	M12 = [ posterior['M12'][i] for i in used_posterior ]
	M21 = [ posterior['M21'][i] for i in used_posterior ]

	## times
	Tsplit = [ posterior['Tsplit'][i] for i in used_posterior ]

	Tdem1 = [ posterior['Tdem1'][i] for i in used_posterior ]
	Tdem2 = [ posterior['Tdem2'][i] for i in used_posterior ]

	## factor of local reduction in Ne. Model of "background selection"
	shape_M12_a = [ posterior['shape_M12_a'][i] for i in used_posterior ]
	shape_M12_b = [ posterior['shape_M12_b'][i] for i in used_posterior ]
	shape_M21_a = [ posterior['shape_M21_a'][i] for i in used_posterior ]
	shape_M21_b = [ posterior['shape_M21_b'][i] for i in used_posterior ]

	# param monolocus: values that will be read by ms
	priorfile = "N1\tN2\tNa\tfounders1\tfounders2\tTdem1\tTdem2\tTsplit\tM12\tshape_M12_a\tshape_M12_b\tM21\tshape_M21_a\tshape_M21_b\n"
	for sim in range(nMultilocus):
		priorfile += "{0:.5f}\t{1:.5f}\t{2:.5f}\t{3:.5f}\t{4:.5f}\t{5}\t{6:.5f}\t{7:.5f}\t{8:.5f}\t{9:.5f}\t{10:.5f}\t{11:.5f}\t{12:.5f}\t{13:.5f}\n".format(N1[sim], N2[sim], Na[sim], founders1[sim], founders2[sim], Tdem1[sim], Tdem2[sim], Tsplit[sim], M12[sim], shape_M12_a[sim], shape_M12_b[sim], M21[sim], shape_M21_a[sim], shape_M21_b[sim])
		
		# vectors of size 'nLoci' containing parameters
                scalar_M12 = beta(shape_M12_a[sim], shape_M12_b[sim], size = nLoci)
                scalar_M21 = beta(shape_M21_a[sim], shape_M21_b[sim], size = nLoci)
                M12_vec = [ M12[sim] * i for i in scalar_M12 ]
                M21_vec = [ M21[sim] * i for i in scalar_M21 ]

	
		for locus in range(nLoci):
			# SC print("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6:.5f}\t{7:.5f}\t{8:.5f}\t{9:.5f}\t{10:.5f}\t{11:.5f}\t{12:.5f}\t{13:.5f}".format(nsam_tot[locus], theta[locus], rho[locus], L[locus], nsamA[locus], nsamB[locus], M12_vec[locus], M21_vec[locus], N1_vec[locus], N2_vec[locus], Tsc[sim], Tsplit[sim], Tsplit[sim], Na_vec[locus]))
			print("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6:.5f}\t{7:.5f}\t{8:.5f}\t{9:.5f}\t{10:.5f}\t{11:.5f}\t{12:.5f}\t{13:.5f}\t{14:.5f}\t{15:.5f}\t{16:.5f}".format(nsam_tot[locus], theta[locus], rho[locus], L[locus], nsamA[locus], nsamB[locus], N1[sim], N2[sim], Tdem1[sim], founders1[sim]*Na[sim], Tdem2[sim], founders2[sim]*Na[sim], M12_vec[locus], M21_vec[locus], Tsplit[sim], Tsplit[sim], Na[sim]))

	outfile = open("priorfile.txt", "w")
	outfile.write(priorfile)
	outfile.close()



if sys.argv[1] == "IM_2M_2N":
	# param multilocus: values that will be printed in priorfile.txt
	N1 = [ posterior['N1'][i] for i in used_posterior ]
	N2 = [ posterior['N2'][i] for i in used_posterior ]
	Na = [ posterior['Na'][i] for i in used_posterior ]
	founders1 = [ posterior['founders1'][i] for i in used_posterior ]
	founders2 = [ posterior['founders2'][i] for i in used_posterior ]

	## Miration rates
	M12 = [ posterior['M12'][i] for i in used_posterior ]
	M21 = [ posterior['M21'][i] for i in used_posterior ]

	## times
	Tsplit = [ posterior['Tsplit'][i] for i in used_posterior ]

	Tdem1 = [ posterior['Tdem1'][i] for i in used_posterior ]
	Tdem2 = [ posterior['Tdem2'][i] for i in used_posterior ]

	## factor of local reduction in Ne. Model of "background selection"
	shape_M12_a = [ posterior['shape_M12_a'][i] for i in used_posterior ]
	shape_M12_b = [ posterior['shape_M12_b'][i] for i in used_posterior ]
	shape_M21_a = [ posterior['shape_M21_a'][i] for i in used_posterior ]
	shape_M21_b = [ posterior['shape_M21_b'][i] for i in used_posterior ]

	## factor of local reduction in Ne. Model of "background selection"
	shape_N_a = [ posterior['shape_N_a'][i] for i in used_posterior ]
	shape_N_b = [ posterior['shape_N_b'][i] for i in used_posterior ]

	# param monolocus: values that will be read by ms
	priorfile = "N1\tN2\tNa\tfounders1\tfounders2\tTdem1\tTdem2\tshape_N_a\tshape_N_b\tTsplit\tM12\tshape_M12_a\tshape_M12_b\tM21\tshape_M21_a\tshape_M21_b\n"
	for sim in range(nMultilocus):
		priorfile += "{0:.5f}\t{1:.5f}\t{2:.5f}\t{3:.5f}\t{4:.5f}\t{5:.5f}\t{6:.5f}\t{7:.5f}\t{8:.5f}\t{9:.5f}\t{10:.5f}\t{11:.5f}\t{12:.5f}\t{13:.5f}\t{14:.5f}\t{15:.5f}\n".format(N1[sim], N2[sim], Na[sim], founders1[sim], founders2[sim], Tdem1[sim], Tdem2[sim], shape_N_a[sim], shape_N_b[sim], Tsplit[sim], M12[sim], shape_M12_a[sim], shape_M12_b[sim], M21[sim], shape_M21_a[sim], shape_M21_b[sim])
		# vectors of size 'nLoci' containing parameters
                scalar_N = beta(shape_N_a[sim], shape_N_b[sim], size=nLoci)
                N1_vec = [ N1[sim]*i for i in scalar_N ]
                N2_vec = [ N2[sim]*i for i in scalar_N ]
                Na_vec = [ Na[sim]*i for i in scalar_N ]
	
                # vectors of size 'nLoci' containing parameters
                scalar_M12 = beta(shape_M12_a[sim], shape_M12_b[sim], size = nLoci)
                scalar_M21 = beta(shape_M21_a[sim], shape_M21_b[sim], size = nLoci)
                M12_vec = [ M12[sim] * i for i in scalar_M12 ]
                M21_vec = [ M21[sim] * i for i in scalar_M21 ]

		for locus in range(nLoci):
			# SC print("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6:.5f}\t{7:.5f}\t{8:.5f}\t{9:.5f}\t{10:.5f}\t{11:.5f}\t{12:.5f}\t{13:.5f}".format(nsam_tot[locus], theta[locus], rho[locus], L[locus], nsamA[locus], nsamB[locus], M12_vec[locus], M21_vec[locus], N1_vec[locus], N2_vec[locus], Tsc[sim], Tsplit[sim], Tsplit[sim], Na_vec[locus]))
			print("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6:.5f}\t{7:.5f}\t{8:.5f}\t{9:.5f}\t{10:.5f}\t{11:.5f}\t{12:.5f}\t{13:.5f}\t{14:.5f}\t{15:.5f}\t{16:.5f}".format(nsam_tot[locus], theta[locus], rho[locus], L[locus], nsamA[locus], nsamB[locus], N1_vec[locus], N2_vec[locus], Tdem1[sim], founders1[sim]*Na_vec[locus], Tdem2[sim], founders2[sim]*Na_vec[locus], M12_vec[locus], M21_vec[locus], Tsplit[sim], Tsplit[sim], Na_vec[locus]))

	outfile = open("priorfile.txt", "w")
	outfile.write(priorfile)
	outfile.close()



if sys.argv[1] == "SI_1N":
	# param multilocus: values that will be printed in priorfile.txt
	N1 = [ posterior['N1'][i] for i in used_posterior ]
	N2 = [ posterior['N2'][i] for i in used_posterior ]
	Na = [ posterior['Na'][i] for i in used_posterior ]
	founders1 = [ posterior['founders1'][i] for i in used_posterior ]
	founders2 = [ posterior['founders2'][i] for i in used_posterior ]

	## times
	Tsplit = [ posterior['Tsplit'][i] for i in used_posterior ]

	Tdem1 = [ posterior['Tdem1'][i] for i in used_posterior ]
	Tdem2 = [ posterior['Tdem2'][i] for i in used_posterior ]

	# param monolocus: values that will be read by ms
	priorfile = "N1\tN2\tNa\tfounders1\tfounders2\tTdem1\tTdem2\tTsplit\n"
	for sim in range(nMultilocus):
		priorfile += "{0:.5f}\t{1:.5f}\t{2:.5f}\t{3:.5f}\t{4:.5f}\t{5:.5f}\t{6:.5f}\t{7:.5f}\n".format(N1[sim], N2[sim], Na[sim], founders1[sim], founders2[sim], Tdem1[sim], Tdem2[sim], Tsplit[sim])
		
		for locus in range(nLoci):
			print("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6:.5f}\t{7:.5f}\t{8:.5f}\t{9:.5f}\t{10:.5f}\t{11:.5f}\t{12:.5f}\t{13:.5f}\t{14:.5f}".format(nsam_tot[locus], theta[locus], rho[locus], L[locus], nsamA[locus], nsamB[locus], N1[sim], N2[sim], Tdem1[sim], founders1[sim]*Na[sim], Tdem2[sim], founders2[sim]*Na[sim], Tsplit[sim], Tsplit[sim], Na[sim]))

	outfile = open("priorfile.txt", "w")
	outfile.write(priorfile)
	outfile.close()

	outfile = open("priorfile.txt", "w")
	outfile.write(priorfile)
	outfile.close()



if sys.argv[1] == "SI_2N":
	# param multilocus: values that will be printed in priorfile.txt
	N1 = [ posterior['N1'][i] for i in used_posterior ]
	N2 = [ posterior['N2'][i] for i in used_posterior ]
	Na = [ posterior['Na'][i] for i in used_posterior ]
	founders1 = [ posterior['founders1'][i] for i in used_posterior ]
	founders2 = [ posterior['founders2'][i] for i in used_posterior ]

	## times
	Tsplit = [ posterior['Tsplit'][i] for i in used_posterior ]

	Tdem1 = [ posterior['Tdem1'][i] for i in used_posterior ]
	Tdem2 = [ posterior['Tdem2'][i] for i in used_posterior ]

	## factor of local reduction in Ne. Model of "background selection"
	shape_N_a = [ posterior['shape_N_a'][i] for i in used_posterior ]
	shape_N_b = [ posterior['shape_N_b'][i] for i in used_posterior ]

	# param monolocus: values that will be read by ms
	priorfile = "N1\tN2\tNa\tfounders1\tfounders2\tTdem1\tTdem2\tshape_N_a\tshape_N_b\tTsplit\n"
	for sim in range(nMultilocus):
		priorfile += "{0:.5f}\t{1:.5f}\t{2:.5f}\t{3:.5f}\t{4:.5f}\t{5:.5f}\t{6:.5f}\t{7:.5f}\t{8:.5f}\t{9:.5f}\n".format(N1[sim], N2[sim], Na[sim], founders1[sim], founders2[sim], Tdem1[sim], Tdem2[sim], shape_N_a[sim], shape_N_b[sim], Tsplit[sim])
		# vectors of size 'nLoci' containing parameters
                scalar_N = beta(shape_N_a[sim], shape_N_b[sim], size=nLoci)
                N1_vec = [ N1[sim]*i for i in scalar_N ]
                N2_vec = [ N2[sim]*i for i in scalar_N ]
                Na_vec = [ Na[sim]*i for i in scalar_N ]
	
		for locus in range(nLoci):
			print("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6:.5f}\t{7:.5f}\t{8:.5f}\t{9:.5f}\t{10:.5f}\t{11:.5f}\t{12:.5f}\t{13:.5f}\t{14:.5f}".format(nsam_tot[locus], theta[locus], rho[locus], L[locus], nsamA[locus], nsamB[locus], N1_vec[locus], N2_vec[locus], Tdem1[sim], founders1[sim]*Na_vec[locus], Tdem2[sim], founders2[sim]*Na_vec[locus], Tsplit[sim], Tsplit[sim], Na_vec[locus]))

	outfile = open("priorfile.txt", "w")
	outfile.write(priorfile)
	outfile.close()

