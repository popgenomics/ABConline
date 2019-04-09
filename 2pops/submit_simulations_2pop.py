#!/usr/bin/python
# #!/home/roux/python/Python-2.7.14/python
import sys
import os
import time

# Example to submit jobs using slurm:
#    for model in AM IM SC; do for M in 1M 2M; do for N in 1N 2N; do ./submit.py 100000 10 ${model}_${M}_${N}; done; done; done
#    for model in SI; do for N in 1N 2N; do ./submit.py 100000 10 ${model}_${N}; done; done


if len(sys.argv) != 6:
	print("\n\tsubmit.py [nmultilocus] [niterations] [model: SI_x AM_x IM_x SC_x PSC_x PAM_x] [nameA] [nameB]")
	print("\n\tex: submit_simulations.py 10000 150 SI_1N flo mal\n")
	sys.exit(0)

nmultilocus = int(sys.argv[1]) # 10000
niterations = int(sys.argv[2]) # 150
model = sys.argv[3]
nameA = sys.argv[4] # name of the species A
nameB = sys.argv[5] # name of the species B

path = os.getcwd() + '/ABC_{0}_{1}'.format(nameA, nameB)

test_bpfile = os.path.isfile('{0}/bpfile'.format(path))
if test_bpfile == False:
	sys.exit('\n\tERROR in submit_simulations_2pop.py : the file {0}/bpfile is not found\n'.format(path))
else:
	infile = open('{0}/bpfile'.format(path), 'r')
	tmp = infile.readline()
	tmp = infile.readline().strip().split('\t')
	nlocus = len(tmp)
	infile.close()

mscommand = ""
if "SI" in model:
	mscommand = "-t tbs -r tbs tbs -I 2 tbs tbs 0 -n 1 tbs -n 2 tbs -ej tbs 2 1 -eN tbs tbs"
if "AM" in model:
	mscommand = "-t tbs -r tbs tbs -I 2 tbs tbs 0 -n 1 tbs -n 2 tbs -ema tbs 2 0 tbs tbs 0 -ej tbs 2 1 -eN tbs tbs"
if "SC" in model:
	mscommand = "-t tbs -r tbs tbs -I 2 tbs tbs 0 -m 1 2 tbs -m 2 1 tbs -n 1 tbs -n 2 tbs -eM tbs 0 -ej tbs 2 1 -eN tbs tbs"
if "IM" in model:
	mscommand = "-t tbs -r tbs tbs -I 2 tbs tbs 0 -n 1 tbs -n 2 tbs -m 1 2 tbs -m 2 1 tbs -ej tbs 2 1 -eN tbs tbs"
if "PSC" in model:
	mscommand = "-t tbs -r tbs tbs -I 2 tbs tbs 0 -n 1 tbs -n 2 tbs -m 1 2 tbs -m 2 1 tbs -ema tbs 2 0 0 0 0 -ema tbs 2 0 tbs tbs 0 -ema tbs 2 0 0 0 0 -ej tbs 2 1 -eN tbs tbs"
if "PAM" in model:
	mscommand = "-t tbs -r tbs tbs -I 2 tbs tbs 0 -n 1 tbs -n 2 tbs -m 1 2 0 -m 2 1 0 -ema tbs 2 0 tbs tbs 0 -ema tbs 2 0 0 0 0 -ema tbs 2 0 tbs tbs 0 -ej tbs 2 1 -eN tbs tbs"

if mscommand == "":
	print("You specified a wrong model: SI_x, AM_x, AM_x or SC_x\n")
	sys.exit()

for i in range(niterations):
#for i in [66,68,69,70,76,77,83,84,90,91,96,97]:
	tmp = "mkdir {0}/{1}_{2}_beta; ".format(path, model, i)
	tmp += "cp {0}/bpfile {0}/{1}_{2}_beta; ".format(path, model, i)
	tmp += "cd {0}/{1}_{2}_beta; ".format(path, model, i)
	tmp += "module load python/2.7.12; "
	tmp += "module load java; "
	tmp += "priorgen_2pop.py {0} {1} | msnsam tbs {2} {3} | mscalc_2pop.py".format(model, nmultilocus, nmultilocus*nlocus, mscommand)
	tmp2 = 'sbatch --nodes=1 --ntasks-per-node=1 --time=24:00:00 -J {0}_{1} --wrap="{2}"'.format(model, i, tmp)
#	os.system(tmp2)
	print(tmp2)

