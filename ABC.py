#!/usr/bin/python
import sys
import os

args = {}
args['infile'] = '' # name of inputfile containing the sequences in a fasta format
args['nspecies'] = '' # specifies the number of species to consider, equal to 1, 2 or 4
args['nameA'] = '' # name of species A. example: flo
args['nameB'] = '' # name of species B. example: mal
args['nameC'] = '' # name of species B. example: txn
args['nameD'] = '' # name of species B. example: ama
args['region'] = '' # coding (will only consider synonymous mutations); noncoding (will consider all SNPs)
args['Lmin'] = '' # minimum number of individuals within a species. example: 10
args['max_N_tolerated'] = '' # if an allele has %N > threshold_N --> sequence is rejected
args['nMin'] = '' # minimum number of individuals within a species. example: 10
args['Nref'] = '' # size of the reference population, arbitrary fixed. i.e: Nref=100000
args['mu'] = '' # mutation rate by bp and by generation. example: 0.00000002
args['rho_over_theta'] = '' # ratio of the recombination rate over mutation. example: 1

help = '\n\tGeneral command line:\n\t\tABC.py infile=data.fasta nspecies=4 nameA=arabidopsis_A nameB=arabidopsis_B nameC=arabidopsis_C nameD=arabidopsis_D region=coding Lmin=30 max_N_tolerated=0.1 nMin=10 Nref=100000 mu=0.0000002 rho_over_theta=1\n\n'
help += '\tExample for 2 populations:\n\t\tABC.py infile=all_loci.fasta nspecies=2 nameA=flo nameB=mal region=coding Lmin=30 max_N_tolerated=0.5 nMin=15 Nref=100000 mu=0.00000002 rho_over_theta=1\n'

nArgs = 0
for i in sys.argv[1::]:
	i = i.split('=')
	if i[0] in args:
		nArgs += 1
		args[i[0]] = i[1]
	else:
		print('\n\tERROR: {0} is not a valid argument\n'.format(i[0]))
		sys.exit(help)

if len(sys.argv) <11 or len(sys.argv)>14:
	sys.exit(help)

if args['nspecies'] not in ['1', '2', '4']:
	print('\n\tERROR: The number of studied species can only be equal to 1, 2 or 4. Here, a value of {0} is provided\n'.format(args['nspecies']))
	sys.exit(help)

if args['nspecies'] == '1':
	if nArgs != 10:
		print('\n\tERROR: The number of arguments is {0} while {1} are expected\n'.format(nArgs, 10))
		sys.exit(help)

if args['nspecies'] == '2':
	if nArgs != 11:
		print('\n\tERROR: The number of arguments is {0} while {1} are expected\n'.format(nArgs, 11))
		sys.exit(help)

if args['nspecies'] == '4':
	if nArgs != 13:
		print('\n\tERROR: The number of arguments is {0} while {1} are expected\n'.format(nArgs, 13))
		sys.exit(help)

if args['nspecies'] == '2':
	# producing input files for the ABC pipeline from the provided sequences data
	commande = 'fasta2ABC_2pops.py {0} {1} {2} {3} {4} {5} {6} {7} {8} {9}'.format(args['infile'], args['nameA'], args['nameB'], args['region'], args['Lmin'], args['max_N_tolerated'], args['nMin'], args['Nref'], args['mu'], args['rho_over_theta'])
	os.system(commande)

print('\n\tEND OF THE ABC ANALYSIS\n')


