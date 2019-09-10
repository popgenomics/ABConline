#!/usr/bin/python
import sys
import os
import time

# Example to submit jobs using slurm:
#    for model in AM IM SC; do for M in 1M 2M; do for N in 1N 2N; do ./submit.py 100000 10 ${model}_${M}_${N}; done; done; done
#    for model in SI; do for N in 1N 2N; do ./submit.py 100000 10 ${model}_${N}; done; done


if len(sys.argv) != 12:
	print("\n\tsubmit_simulations_2pop.py [outgroup] [nmultilocus] [iteration] [model: SI_x AM_x IM_x SC_x PSC_x PAM_x] [nameA] [nameB] [sub_dir_sim] [sub_dir_model] [config_yaml] [project's directory name, i.e, timeStamp] [beta or bimodal]")
	print("\n\tex: submit_simulations_2pop.py 1 1000 2 SI_1N flo mal sim_SI_1N SI_1N config.yaml Ng4PymB1dy beta\n\tto simulate 1000 multilocus simulations at the second iteration, in the folder sim_SI_1N, with outgroup") 
	sys.exit(0)

outgroup = int(sys.argv[1])
nmultilocus = int(sys.argv[2]) # 10000
iteration = int(sys.argv[3]) # 2
model = sys.argv[4]
nameA = sys.argv[5] # name of the species A
nameB = sys.argv[6] # name of the species B
sub_dir_sim = sys.argv[7] # name of the subdir where the simulations will be run
sub_dir_model = sys.argv[8] # name of the sub_sub_dir containing ABCstat.txt
config_yaml = sys.argv[9]
timeStamp = sys.argv[10]
modeBarrier = sys.argv[11]

path = os.getcwd() + '/{0}'.format(timeStamp)

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


#tmp = "mkdir {0}/{1}; ".format(path, sub_dir_sim)
#tmp += "mkdir {0}/{1}/{2}_{3}; ".format(path, sub_dir_sim, sub_dir_model, iteration)
tmp = "cp {0}/bpfile {0}/{1}/{2}_{3}; ".format(path, sub_dir_sim, sub_dir_model, iteration)
tmp += "cd {0}/{1}/{2}_{3}; ".format(path, sub_dir_sim, sub_dir_model, iteration)


#tmp += "module load python/2.7.12; "
#tmp += "module load java; "
#tmp += "priorgen_2pop.py {0} {1} {2} | msnsam tbs {3} {4} | mscalc_2pop.py".format(model, nmultilocus, config_yaml, nmultilocus*nlocus, mscommand)
tmp += "priorgen_2pop.py {0} {1} {2} | msnsam tbs {3} {4} | mscalc_2pop_SFS.py {5}".format(model, nmultilocus, config_yaml, nmultilocus*nlocus, mscommand, outgroup)
tmp2 = 'sbatch --nodes=1 --ntasks-per-node=1 --time=02:00:00 -J {0}_{1} --wrap="{2}"\n'.format(model, iteration, tmp)

print(tmp)
os.system(tmp) # to submit the job using slurm

