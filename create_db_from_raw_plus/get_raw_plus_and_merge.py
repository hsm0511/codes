import sys
import os, time
from shutil import copyfile
from multiprocessing import Process, Manager, Pool

start_time=time.time()
raw_plus_path_file=sys.argv[1]
raw_plus_paths_l=[]
task_path=sys.argv[2]

# get kor20k samples' sampleID_Raw_Plus.txt file paths
RP_path_f=open(raw_plus_path_file,'r')
for line in RP_path_f:
    if line.rstrip()!='':
        raw_plus_paths_l.append(line.rstrip())
RP_path_f.close()

sh_cmd_l=[]
number_of_cores=48
script_path='/clinix/Analysis/Projects/PersonalGenomics_RnD/add_DB_from_RawPlus/Scripts/'
merge_script=script_path+'merge_genotype.py'

# get number of files for each process
if len(raw_plus_paths_l)<number_of_cores:
    number_of_cores=len(raw_plus_paths_l)

rp_number_per_process=len(raw_plus_paths_l)/number_of_cores
rp_number_per_process_list=[rp_number_per_process]*number_of_cores
rp_remainder=len(raw_plus_paths_l)%number_of_cores
for i in range(rp_remainder):
    rp_number_per_process_list[i]+=1

# split raw plus paths for multiprocess
rp_for_merge_split=[]
start_idx_rp=0
for i in range(number_of_cores):
    last_idx_rp=start_idx_rp+rp_number_per_process_list[i]
    lb=raw_plus_paths_l[start_idx_rp:last_idx_rp]
    rp_for_merge_split.append(lb)
    start_idx_rp=last_idx_rp

# get shell command list
for i in range(number_of_cores):
    rp_l=rp_for_merge_split[i]
    sh_cmd='/mnt/clinix/Tools/Dev/python/bin/python '+merge_script+' '+','.join(rp_l)+' '+str(i)+' '+task_path
    sh_cmd_l.append(sh_cmd)

# execute multiprocess
pool=Pool(number_of_cores)
pool.map(os.system,sh_cmd_l)
pool.close()
pool.join()






