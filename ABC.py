#!/usr/bin/python
import sys
import os

# Hard-coded parameters
# number of CPU used to simulate a given model:
nCPU = 20
# number of CPU used to run the random forest:
ncores = 6
# number of tress for the random forest
ntree = 1000

# Parameters provided by users
args = {}
# for fasta2ABC.py
args['infile'] = '' # name of inputfile containing the sequences in a fasta format
args['nspecies'] = '' # specifies the number of species to consider, equal to 1, 2 or 4
args['nameA'] = '' # name of species A. example: flo
args['nameB'] = '' # name of species B. example: mal
args['nameC'] = '' # name of species C. example: txn
args['nameD'] = '' # name of species D. example: ama
args['outgroup'] = 'NA' # name of the outgroup. If no outgroup, then set the value to 'NA'. Obviously, an outgroup exactly named 'NA' will be discarded.
args['region'] = '' # coding (will only consider synonymous mutations); noncoding (will consider all SNPs)
args['Lmin'] = '' # minimum number of individuals within a species. example: 10
args['max_N_tolerated'] = '' # if an allele has %N > threshold_N --> sequence is rejected
args['nMin'] = '' # minimum number of individuals within a species. example: 10
args['Nref'] = '' # size of the reference population, arbitrary fixed. i.e: Nref=100000
args['mu'] = '' # mutation rate by bp and by generation. example: 0.00000002
args['rho_over_theta'] = '' # ratio of the recombination rate over mutation. example: 1

# for submit_simulations.py
args['nSimulations'] = ''

help = '\n\tGeneral command line:\n\t\tABC.py infile=data.fasta nspecies=4 nameA=arabidopsis_A nameB=arabidopsis_B nameC=arabidopsis_C nameD=arabidopsis_D region=coding Lmin=30 max_N_tolerated=0.1 nMin=10 Nref=100000 mu=0.0000002 rho_over_theta=1 nSimulations=3000\n\n'
help += '\tExample for 2 populations without outgroup:\n\t\tABC.py infile=all_loci.fasta nspecies=2 nameA=flo nameB=mal region=coding Lmin=30 max_N_tolerated=0.5 nMin=15 Nref=100000 mu=0.00000002 rho_over_theta=1 nSimulations=3000\n'
help += '\tExample for 2 populations with outgroup:\n\t\tABC.py infile=all_loci.fasta nspecies=2 nameA=flo nameB=mal outgroup=num region=coding Lmin=30 max_N_tolerated=0.5 nMin=15 Nref=100000 mu=0.00000002 rho_over_theta=1 nSimulations=3000\n'

nArgs = 0
for i in sys.argv[1::]:
	i = i.split('=')
	if i[0] in args:
		nArgs += 1
		args[i[0]] = i[1]
	else:
		print('\n\tERROR: {0} is not a valid argument\n'.format(i[0]))
		sys.exit(help)

if len(sys.argv) <13 or len(sys.argv)>16:
	sys.exit(help)

if args['nspecies'] not in ['1', '2', '4']:
	print('\n\tERROR: The number of studied species can only be equal to 1, 2 or 4. Here, a value of {0} is provided\n'.format(args['nspecies']))
	sys.exit(help)

if args['nspecies'] == '1':
	if nArgs != 12:
		print('\n\tERROR: The number of arguments is {0} while {1} are expected\n'.format(nArgs, 12))
		sys.exit(help)

if args['nspecies'] == '2':
	if nArgs != 13:
		print('\n\tERROR: The number of arguments is {0} while {1} are expected\n'.format(nArgs, 13))
		sys.exit(help)
	models = ['SC_1M_1N', 'SC_1M_2N', 'SC_2M_1N', 'SC_2M_2N', 'AM_1M_1N', 'AM_1M_2N', 'AM_2M_1N', 'AM_2M_2N', 'IM_1M_1N', 'IM_1M_2N', 'IM_2M_1N', 'IM_2M_2N', 'SI_1N', 'SI_2N']

if args['nspecies'] == '4':
	if nArgs != 15:
		print('\n\tERROR: The number of arguments is {0} while {1} are expected\n'.format(nArgs, 15))
		sys.exit(help)

if args['outgroup'] != 'NA':
	use_sfs = 1
else:
	use_sfs = 0
### END OF DATA CURATION ###


if args['nspecies'] == '2':
	# producing input files for the ABC pipeline from the provided sequences data
	# fasta2ABC_2pops.py
	print('\tPreparation of input files')
	commande = 'fasta2ABC_2pops.py {0} {1} {2} {3} {4} {5} {6} {7} {8} {9} {10}'.format(args['infile'], args['nameA'], args['nameB'], args['outgroup'], args['region'], args['Lmin'], args['max_N_tolerated'], args['nMin'], args['Nref'], args['mu'], args['rho_over_theta'])
	print(commande)
	os.system(commande)
	
	# recombination rates	
	print('\n\tRecombination rates')
	commande = 'RNAseqFGT {0} ABC_{1}_{2}/results_recombination.txt'.format(args['infile'], args['nameA'], args['nameB'])
	print(commande)
	os.system(commande)
		
	# submit simulations
	# submit_simulations_2pop.py
	# print("\n\tsubmit.py [nmultilocus] [niterations] [model: SI_x AM_x IM_x SC_x PSC_x PAM_x] [nameA] [nameB]")
	# print("\n\tex: submit_simulations.py 10000 150 SI_1N flo mal\n")
	nSimulations = int(args['nSimulations'])
	nMultilocus = nSimulations / nCPU
	print('\n\tSimulations')
	for model in models:
		commande = 'submit_simulations_2pop.py {0} {1} {2} {3} {4}'.format(nMultilocus, nCPU, model, args['nameA'], args['nameB'])
		print(commande)
		os.system(commande)

	# multilocus demographic inferences
	# model_comp_2pop.R
	print('\n\tDemographic inferences')
	commande = 'model_comp_2pop.R nameA={0} nameB={1} nreps={2} Nref={3} ntree={4} ncores={5} outgroup={6}'.format(args['nameA'], args['nameB'], nCPU, args['Nref'], ntree, ncores, use_sfs)
	print(commande)
	os.system(commande)

        # gather the results
	print('\n\tArchive the results')
        commande = 'tar -czvf results_{0}_{1}.tar.gz ABC_{0}_{1}/ABC* ABC_{0}_{1}/*{0}_{1}* ABC_{0}_{1}/results*'.format(args['nameA'], args['nameB'])
	print(commande)
        os.system(commande)

        # clean the space
	print('\n\tClean the space')
        commande = 'rm -rf ABC_{0}_{1}'.format(args['nameA'], args['nameB'])
	print(commande)
        os.system(commande)
print('\n\tEND OF THE ABC ANALYSIS\n')


