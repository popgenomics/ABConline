#!/usr/bin/pypy
import sys
import os
from math import ceil

# check the arguments
if len(sys.argv) != 11:
	print("\n\tfasta2ABC_2pops.py produces: bpfile (for simulations) and summary statistics (for inferences)")
	print("\n\033[1;33m\tExample: ./fasta2ABC_2pops.py all_loci.fasta flo mal coding 30 0.1 10 100000 0.00000002 1\033[0m\n")
	print("\t\targ1 =\tname of the fasta file containing all of the sequences")
	print("\t\targ2 =\tID of species A (example: flo)")
	print("\t\targ3 =\tID of species B (example: mal)")
	print("\t\targ4 =\t'coding' or 'noncoding', to study only synonymous polymorphisms (if coding) or all SNPs (if noncoding)")
	print("\t\targ5 =\tminimum length of a locus to be considered, i.e, number of total positions minus the positions containing a N")
	print("\t\targ6 =\tvalue in [0-1]. Corresponds to a threshold of %N above which a sequence is rejected")
	print("\t\targ7 =\tinteger, corresponding to the minimum number of retained sequences (after rejection).\n\t\t\tif not enough sequences are retained, the loci is excluded from the analysis")
	print("\t\targ8 =\tsize of the reference population, arbitrary fixed. i.e: Nref=100000")
	print("\t\targ9 =\tmutation rate by bp and by generation. example: 0.00000002")
	print("\t\targ10 =\tratio of the recombination rate over mutation. example: 1")
	if(len(sys.argv)<11):
		sys.exit("\n\033[1;31m ERROR: 10 arguments are required: {0} missing\033[0m\n".format(11-len(sys.argv)))
	if(len(sys.argv)>11):
		sys.exit("\n\033[1;31m ERROR: 10 arguments are required: {0} too much\033[0m\n".format(len(sys.argv)-11))

fileName = sys.argv[1] # example: all_loci.fasta
nameA = sys.argv[2] # name of species A. example: flo
nameB = sys.argv[3] # name of species B. example: mal
region = sys.argv[4] # if == coding: will only deal with synonymous codons; if == noncoding: will deal with all positions
Lmin = int(sys.argv[5]) # minimum length for a locus to be retained. example: 30
max_N_tolerated = float(sys.argv[6]) # if an allele has %N > threshold_N --> sequence is rejected
nMin = int(sys.argv[7]) # minimum number of individuals within a species. example: 10
Nref = int(sys.argv[8]) # size of the reference population, arbitrary fixed. i.e: Nref=100000
mu = float(sys.argv[9]) # mutation rate by bp and by generation. example: 0.00000002
rho_over_theta = float(sys.argv[10]) # ratio of the recombination rate over mutation. example: 1

test = os.path.isfile(fileName)
if test == False:
	sys.exit("\n\033[1;31m ERROR: alignement '{0}' is not found\033[0m\n".format(fileName))

if region not in ['coding', 'noncoding']:
	sys.exit("\n\033[1;31m ERROR: '{0}' is not an expected argument. 'coding' or 'noncoding' are expected here\033[0m\n".format(region))

if os.path.isdir('ABC_{0}_{1}'.format(nameA, nameB)) == True:
	commande = 'rm ABC_{0}_{1}'.format(nameA, nameB)
	os.system(commande)
commande = 'mkdir ABC_{0}_{1}'.format(nameA, nameB)
os.system(commande)

def coloredSeq(seq):
	# print sequences with the standard color code
	seq = seq.replace("A", '\x1b[5;30;41m' + 'A' + '\x1b[0m')
	seq = seq.replace("T", '\x1b[5;30;44m' + 'T' + '\x1b[0m')
	seq = seq.replace("G", '\x1b[6;30;43m' + 'G' + '\x1b[0m')
	seq = seq.replace("C", '\x1b[5;30;42m' + 'C' + '\x1b[0m')
	return(seq)


def trunc2triplets(size):
	# trunc a value to its closest and smaller multiple of 3
	size = int(size)
	for i in range(2):
		if size%3 != 0:
			size -= 1
	return(size)


# nN = number of non-synonymous sites in the codon i: example for CGG -> nN = 2/3 + 3/3 + 0/3
# nS = number of synonymous sites in the codon i: example for CGG -> n> = 1/3 + 0/3 + 3/3
codonTable = {'AAA': {'aa': 'K', 'nN': 2.6666666666666665, 'nS': 0.3333333333333333}, 'AAC': {'aa': 'N', 'nN': 2.6666666666666665, 'nS': 0.3333333333333333}, 'AAG': {'aa': 'K', 'nN': 2.6666666666666665, 'nS': 0.3333333333333333}, 'AAT': {'aa': 'N', 'nN': 2.6666666666666665, 'nS': 0.3333333333333333}, 'ACA': {'aa': 'T', 'nN': 2.0, 'nS': 1.0}, 'ACC': {'aa': 'T', 'nN': 2.0, 'nS': 1.0}, 'ACG': {'aa': 'T', 'nN': 2.0, 'nS': 1.0}, 'ACT': {'aa': 'T', 'nN': 2.0, 'nS': 1.0}, 'AGA': {'aa': 'R', 'nN': 2.3333333333333335, 'nS': 0.6666666666666666}, 'AGC': {'aa': 'S', 'nN': 2.6666666666666665, 'nS': 0.3333333333333333}, 'AGG': {'aa': 'R', 'nN': 2.3333333333333335, 'nS': 0.6666666666666666}, 'AGT': {'aa': 'S', 'nN': 2.6666666666666665, 'nS': 0.3333333333333333}, 'ATA': {'aa': 'I', 'nN': 2.3333333333333335, 'nS': 0.6666666666666666}, 'ATC': {'aa': 'I', 'nN': 2.3333333333333335, 'nS': 0.6666666666666666}, 'ATG': {'aa': 'M', 'nN': 3.0, 'nS': 0.0}, 'ATT': {'aa': 'I', 'nN': 2.3333333333333335, 'nS': 0.6666666666666666}, 'CAA': {'aa': 'Q', 'nN': 2.6666666666666665, 'nS': 0.3333333333333333}, 'CAC': {'aa': 'H', 'nN': 2.6666666666666665, 'nS': 0.3333333333333333}, 'CAG': {'aa': 'Q', 'nN': 2.6666666666666665, 'nS': 0.3333333333333333}, 'CAT': {'aa': 'H', 'nN': 2.6666666666666665, 'nS': 0.3333333333333333}, 'CCA': {'aa': 'P', 'nN': 2.0, 'nS': 1.0}, 'CCC': {'aa': 'P', 'nN': 2.0, 'nS': 1.0}, 'CCG': {'aa': 'P', 'nN': 2.0, 'nS': 1.0}, 'CCT': {'aa': 'P', 'nN': 2.0, 'nS': 1.0}, 'CGA': {'aa': 'R', 'nN': 1.6666666666666667, 'nS': 1.3333333333333333}, 'CGC': {'aa': 'R', 'nN': 2.0, 'nS': 1.0}, 'CGG': {'aa': 'R', 'nN': 1.6666666666666667, 'nS': 1.3333333333333333}, 'CGT': {'aa': 'R', 'nN': 2.0, 'nS': 1.0}, 'CTA': {'aa': 'L', 'nN': 1.6666666666666667, 'nS': 1.3333333333333333}, 'CTC': {'aa': 'L', 'nN': 2.0, 'nS': 1.0}, 'CTG': {'aa': 'L', 'nN': 1.6666666666666667, 'nS': 1.3333333333333333}, 'CTT': {'aa': 'L', 'nN': 2.0, 'nS': 1.0}, 'GAA': {'aa': 'E', 'nN': 2.6666666666666665, 'nS': 0.3333333333333333}, 'GAC': {'aa': 'D', 'nN': 2.6666666666666665, 'nS': 0.3333333333333333}, 'GAG': {'aa': 'E', 'nN': 2.6666666666666665, 'nS': 0.3333333333333333}, 'GAT': {'aa': 'D', 'nN': 2.6666666666666665, 'nS': 0.3333333333333333}, 'GCA': {'aa': 'A', 'nN': 2.0, 'nS': 1.0}, 'GCC': {'aa': 'A', 'nN': 2.0, 'nS': 1.0}, 'GCG': {'aa': 'A', 'nN': 2.0, 'nS': 1.0}, 'GCT': {'aa': 'A', 'nN': 2.0, 'nS': 1.0}, 'GGA': {'aa': 'G', 'nN': 2.0, 'nS': 1.0}, 'GGC': {'aa': 'G', 'nN': 2.0, 'nS': 1.0}, 'GGG': {'aa': 'G', 'nN': 2.0, 'nS': 1.0}, 'GGT': {'aa': 'G', 'nN': 2.0, 'nS': 1.0}, 'GTA': {'aa': 'V', 'nN': 2.0, 'nS': 1.0}, 'GTC': {'aa': 'V', 'nN': 2.0, 'nS': 1.0}, 'GTG': {'aa': 'V', 'nN': 2.0, 'nS': 1.0}, 'GTT': {'aa': 'V', 'nN': 2.0, 'nS': 1.0}, 'TAC': {'aa': 'Y', 'nN': 2.6666666666666665, 'nS': 0.3333333333333333}, 'TAT': {'aa': 'Y', 'nN': 2.6666666666666665, 'nS': 0.3333333333333333}, 'TCA': {'aa': 'S', 'nN': 2.0, 'nS': 1.0}, 'TCC': {'aa': 'S', 'nN': 2.0, 'nS': 1.0}, 'TCG': {'aa': 'S', 'nN': 2.0, 'nS': 1.0}, 'TCT': {'aa': 'S', 'nN': 2.0, 'nS': 1.0}, 'TGC': {'aa': 'C', 'nN': 2.6666666666666665, 'nS': 0.3333333333333333}, 'TGG': {'aa': 'W', 'nN': 3.0, 'nS': 0.0}, 'TGT': {'aa': 'C', 'nN': 2.6666666666666665, 'nS': 0.3333333333333333}, 'TTA': {'aa': 'L', 'nN': 2.3333333333333335, 'nS': 0.6666666666666666}, 'TTC': {'aa': 'F', 'nN': 2.6666666666666665, 'nS': 0.3333333333333333}, 'TTG': {'aa': 'L', 'nN': 2.3333333333333335, 'nS': 0.6666666666666666}, 'TTT': {'aa': 'F', 'nN': 2.6666666666666665, 'nS': 0.3333333333333333}}


#def fasta2dic(fastaFile):
#	fasta = open(fastaFile).readlines()
#	seqName = [x.split(" ")[0].rstrip().replace('>','') for x in fasta if x[0] == '>']
#	seq = ''.join([x.rstrip() if x[0]!='>' else '@' for x in fasta])[1:].split('@')
#	res = {}
#	for i in range(len(seq)):
#		res[seqName[i]] = seq[i]
#	return (res)

#alignA = fasta2dic(seqA)
#alignB = fasta2dic(seqB)


def fasta2list(fastaFile, nameA, nameB, nMin, max_N_tolerated):
	L = {}
	res = {}
	res[nameA] = {}
	res[nameB] = {}
	fasta = open(fastaFile).readlines()
	seqName = [x.split(" ")[0].rstrip().replace('>','') for x in fasta if x[0] == '>']
	seq = ''.join([x.rstrip() if x[0]!='>' else '@' for x in fasta])[1:].split('@')
	
	nsam = {} # number of individuals in both species
	for i in range(len(seqName)):
		tmp = seqName[i].split('|') # split a string similar to : Hmel219002_6|flo|flo.CS12|allele1
		locus = tmp[0]
		species = tmp[1]
		if locus not in nsam:
			nsam[locus] = {}
			nsam[locus][nameA] = 0
			nsam[locus][nameB] = 0

		if species == nameA:
			if locus not in res[species]:
				res[species][locus] = {}
				res[species][locus]['seq'] = []
				res[species][locus]['id'] = []
		if species == nameB:
			if locus not in res[species]:
				res[species][locus] = {}
				res[species][locus]['seq'] = []
				res[species][locus]['id'] = []
	
		# remove sequences with too many N
		propN = seq[i].count("N")/(1.0 * len(seq[i]))
		if propN <= max_N_tolerated:
			res[species][locus]['seq'].append(seq[i])
			res[species][locus]['id'].append(seqName[i])
			nsam[locus][species] += 1
	
	# remove loci not found in sufficient individuals in both species
	for i in nsam: # loop along loci
		test = 0
		if len(res[nameA][i]['seq'])==0 or len(res[nameB][i]['seq'])==0:
			test = 1
			del res[nameA][i] # delete locus i
			del res[nameB][i]
		else:
			if nsam[i][nameA] < nMin or nsam[i][nameB] < nMin:
				test = 1
				del res[nameA][i] # delete locus i
				del res[nameB][i]
		if test == 0:
			L[i] = len(res[nameA][i]['seq'][0]) - 3 # remove the last 3 bases to excluse final stop codon
			L[i] = trunc2triplets(L[i]) # convert the remaining length into a multiple of 3
	return ({'align': res, 'L': L})

# read the input file
align = fasta2list(fileName, nameA, nameB, nMin,max_N_tolerated)  # align[species][locus]['id', 'seq']

# treat the input file
bpfile_L1 = '# {0} {1} {2}'.format(nameA, nameB, mu)
bpfile_L2 = [] # L
bpfile_L3 = [] # nA
bpfile_L4 = [] # nB
bpfile_L5 = [] # theta
bpfile_L6 = [] # rho

output_ms = "./msnsam tbs 20 -t tbs -r tbs tbs -I 2 tbs tbs 0 -m 1 2 tbs -m 2 1 tbs -n 1 tbs -n 2 tbs -ej tbs 1 2 -eN tbs tbs\n3579 27011 59243\n\n"
outfile_ms = open('ABC_{0}_{1}/{0}_{1}.ms'.format(nameA, nameB), 'w')
outfile_ms.write(output_ms)


# For coding loci
if region == 'coding':
	output_info = "locusName\tL_including_N\tLsyno\tnSynSegSite\tnsamA\tnsamB\n"
	outfile_info = open('ABC_{0}_{1}/{0}_{1}_infos.txt'.format(nameA, nameB), 'w')
	outfile_info.write(output_info)
	for locus_i in align['L'].keys():
		geneName = locus_i
		
		nA = len(align['align'][nameA][locus_i]['id'])
		nB = len(align['align'][nameB][locus_i]['id'])
		
		L = align['L'][locus_i] 
		interspe = [] # contains the interspecific alignment
		interspeName = [] # contains the id of sequences in the interspecific alignment
		for i in range(nA):
			interspe.append(align['align'][nameA][locus_i]['seq'][i])
			interspeName.append(align['align'][nameA][locus_i]['id'][i])
		for i in range(nB):
			interspe.append(align['align'][nameB][locus_i]['seq'][i])
			interspeName.append(align['align'][nameB][locus_i]['id'][i])

		nSites = 0 # total number of synonymous sites within the sequence, computed using codonTable
		nSynSegSite = 0 # number of synonymous segregating sites among the nSites
		positions = [] # list of synonymous polymorphic positions: doesn't correspond to the SNP position, but to the first codon position
		msStyle = [] # contains the msStyle format
		for ind in range(nA):
			msStyle.append([])
		for ind in range(nB):
			msStyle.append([])

		# loop over codons:
		for pos in range(L)[::3]:
			alignmentOfCodons = [] # set of codons in the alignment, starting at the position 'pos1'
			# loop over individuals:
			# get all codons in the alignment
			for ind in range(nA + nB):
				pos1 = interspe[ind][pos]
				pos2 = interspe[ind][pos + 1]
				pos3 = interspe[ind][pos + 2]
				base = pos1 + pos2 + pos3 
				alignmentOfCodons.append(base)
			
			polyMcodons = list(set(alignmentOfCodons)) # list of codons found in the alignment
			nCodons = 0
			nCodons = len(polyMcodons)
			testN = False # False if no codon with 'N'; True if a 'N' is found in at least one codon for one individual
			testStopCodon = False # False if no stop codon was found; True if a stop codon was found
			for i in polyMcodons: # loop to test for some 'N'
				if 'N' in i:
					testN = True
				if i not in codonTable:
					testStopCodon = True
			
			# if: 1) a maximum of 2 polymorphic codons, and, 2) no codon with 'N', and, 3) all codons effectively code for an amino acid
			if nCodons <= 2 and testN==False and testStopCodon==False: 
				nSites_pos = 0.0
				for i in alignmentOfCodons:
					nSites_pos += codonTable[i]['nS']
				nSites += nSites_pos/len(alignmentOfCodons)
				
				# if two codons --> there is a polymorphism
				if nCodons == 2:
					alignmentOfAminoAcids = []
					for i in alignmentOfCodons:
						alignmentOfAminoAcids.append(codonTable[i]['aa'])
					setOfAminoAcids = list(set(alignmentOfAminoAcids))
					
					# if two codons but one amino acids --> synonymous polymorphism
					if len(setOfAminoAcids) == 1:
						nSynSegSite += 1
						positions.append(pos) # positions: list of first codon position of polymorphic synonymous codons
						ancestralAllele = polyMcodons[0] # in absence of outgroup --> the ancestral allele is the first in the alignement
						derivedAllele = polyMcodons[1] # without outgroup --> the derived allele is the one who is not the first...
						for i in range(nA + nB):
							if alignmentOfCodons[i] == ancestralAllele:
								msStyle[i].append('0')
							if alignmentOfCodons[i] == derivedAllele:
								msStyle[i].append('1')

		if nSites >= Lmin: # if the locus is big enough to be considered
			# ms_like output files
			locus_ms = ''
			locus_ms = locus_ms + "//{0}\n".format(geneName)
			locus_ms = locus_ms + "segsites: {0}\n".format(int(nSynSegSite))
			if nSynSegSite != 0:
				locus_ms += "positions: {0}\n".format( " ".join([ str(round((1.0*i)/L, 4)) for i in positions ]))
				for i in msStyle:
					locus_ms = locus_ms + "".join( [ str(j) for j in i ] ) + "\n"
			
			outfile_ms.write(locus_ms + '\n')
			
			bpfile_L2.append(int(ceil(nSites)))
			bpfile_L3.append(nA)
			bpfile_L4.append(nB)
			bpfile_L5.append(4*Nref*mu*nSites)
			bpfile_L6.append(4*Nref*mu*nSites*rho_over_theta)
				
			# informations about locus
			res = "{0}\t{1}\t{2}\t{3}\t{4}\t{5}\n".format(geneName, L, nSites, nSynSegSite, nA, nB)
			outfile_info.write(res)
			
		#	res = ""
		#	for i in range(len(interspe)):
		#		res += ">{0}\n{1}\n".format(interspeName[i], interspe[i])
		#	outfile = open('{0}.fas'.format(geneName), "w")
		#	outfile.write(res)
		#	outfile.close()
	outfile_ms.close()
	outfile_info.close()

	bpfile = bpfile_L1 + '\n'
	bpfile += '\t'.join([ str(i) for i in bpfile_L2 ]) + '\n'
	bpfile += '\t'.join([ str(i) for i in bpfile_L3 ]) + '\n'
	bpfile += '\t'.join([ str(i) for i in bpfile_L4 ]) + '\n'
	bpfile += '\t'.join([ str(i) for i in bpfile_L5 ]) + '\n'
	bpfile += '\t'.join([ str(i) for i in bpfile_L6 ]) + '\n'

	outfile = open('ABC_{0}_{1}/bpfile'.format(nameA, nameB), 'w')
	outfile.write(bpfile)
	outfile.close()

# For non coding loci
if region == 'noncoding':
	output_info = "locusName\tL_including_N\tL\tnSegSite\tnsamA\tnsamB\n"
	outfile_info = open('ABC_{0}_{1}_infos.txt'.format(nameA, nameB), 'w')
	outfile_info.write(output_info)
	for locus_i in align['L'].keys():
		geneName = locus_i
		
		nA = len(align['align'][nameA][locus_i]['id'])
		nB = len(align['align'][nameB][locus_i]['id'])
		
		L = align['L'][locus_i] 
		interspe = [] # contains the interspecific alignment
		interspeName = [] # contains the id of sequences in the interspecific alignment
		for i in range(nA):
			interspe.append(align['align'][nameA][locus_i]['seq'][i])
			interspeName.append(align['align'][nameA][locus_i]['id'][i])
		for i in range(nB):
			interspe.append(align['align'][nameB][locus_i]['seq'][i])
			interspeName.append(align['align'][nameB][locus_i]['id'][i])

		nSites = 0 # total number of sites within the sequence
		nSegSite = 0 # number of segregating sites among the nSites
		positions = [] # list of polymorphic positions: correspond to the SNP position
		msStyle = [] # contains the msStyle format
		for ind in range(nA):
			msStyle.append([])
		for ind in range(nB):
			msStyle.append([])

		# loop over pos:
		for pos in range(L):
			alignmentOfPos = [] # set of pos in the alignment, starting at the position 'pos1'
			# loop over individuals:
			# get all pos in the alignment
			for ind in range(nA + nB):
				pos1 = interspe[ind][pos]
				base = pos1 
				alignmentOfPos.append(base)
			
			polyMpos = list(set(alignmentOfPos)) # list of pos found in the alignment
			nPos = 0
			nPos = len(polyMpos)
			testN = False # False if no codon with 'N'; True if a 'N' is found in at least one codon for one individual
			for i in polyMpos: # loop to test for some 'N'
				if 'N' in i:
					testN = True
			
			# if: 1) a maximum of 2 polymorphic pos, and, 2) no codon with 'N'
			if nPos <= 2 and testN==False: 
				nSites += 1
				
				# if two pos --> there is a polymorphism
				if nPos == 2:
					nSegSite += 1
					positions.append(pos) # positions: list of first codon position of polymorphic synonymous pos
					ancestralAllele = polyMpos[0] # in absence of outgroup --> the ancestral allele is the first in the alignement
					derivedAllele = polyMpos[1] # without outgroup --> the derived allele is the one who is not the first...
					for i in range(nA + nB):
						if alignmentOfPos[i] == ancestralAllele:
							msStyle[i].append('0')
						if alignmentOfPos[i] == derivedAllele:
							msStyle[i].append('1')

		if nSites >= Lmin: # if the locus is big enough to be considered
			# ms_like output files
			locus_ms = ''
			locus_ms = locus_ms + "//{0}\n".format(geneName)
			locus_ms = locus_ms + "segsites: {0}\n".format(int(nSegSite))
			if nSegSite != 0:
				locus_ms += "positions: {0}\n".format( " ".join([ str(round((1.0*i)/L, 4)) for i in positions ]))
				for i in msStyle:
					locus_ms = locus_ms + "".join( [ str(j) for j in i ] ) + "\n"
			
			outfile_ms.write(locus_ms + '\n')
			
			bpfile_L2.append(int(ceil(nSites)))
			bpfile_L3.append(nA)
			bpfile_L4.append(nB)
			bpfile_L5.append(4*Nref*mu*nSites)
			bpfile_L6.append(4*Nref*mu*nSites*rho_over_theta)
				
			# informations about locus
			res = "{0}\t{1}\t{2}\t{3}\t{4}\t{5}\n".format(geneName, L, nSites, nSegSite, nA, nB)
			outfile_info.write(res)
			
		#	res = ""
		#	for i in range(len(interspe)):
		#		res += ">{0}\n{1}\n".format(interspeName[i], interspe[i])
		#	outfile = open('{0}.fas'.format(geneName), "w")
		#	outfile.write(res)
		#	outfile.close()
	outfile_ms.close()
	outfile_info.close()

	bpfile = bpfile_L1 + '\n'
	bpfile += '\t'.join([ str(i) for i in bpfile_L2 ]) + '\n'
	bpfile += '\t'.join([ str(i) for i in bpfile_L3 ]) + '\n'
	bpfile += '\t'.join([ str(i) for i in bpfile_L4 ]) + '\n'
	bpfile += '\t'.join([ str(i) for i in bpfile_L5 ]) + '\n'
	bpfile += '\t'.join([ str(i) for i in bpfile_L6 ]) + '\n'

	outfile = open('ABC_{0}_{1}/bpfile'.format(nameA, nameB), 'w')
	outfile.write(bpfile)
	outfile.close()


# compute the summary statistics for ABC
commande = 'cat ABC_{0}_{1}/{0}_{1}.ms | mscalc_2pop_observedDataset.py ABC_{0}_{1}'.format(nameA, nameB)
os.system(commande)

# remove the useless ms file
commande = 'rm ABC_{0}_{1}/{0}_{1}.ms'.format(nameA, nameB)
os.system(commande)

