import os, time
from shutil import copyfile
from multiprocessing import Process, Manager, Pool

start_time=time.time()
cp_path='/mnt/clinix/Analysis/Projects/PersonalGenomics_RnD/merge_raw_plus_20k/'
raw_plus_paths_l=[]

# get kor20k samples' sampleID_Raw_Plus.txt file paths
FR_path_f=open('/mnt/clinix/Analysis/Projects/PersonalGenomics_RnD/merge_test_416/20k_samples_FR_paths.txt','r')
for line in FR_path_f:
    if line.rstrip()!='':
        tmp=line.rstrip().split('/')
        sample_name=tmp[-1].split('_Final')[0]
        path_front='/'.join(tmp[:8])
        raw_plus_path=path_front+'/Plus_Genotype/'+sample_name+'_Raw_Plus.txt'
        raw_plus_paths_l.append(raw_plus_path)
FR_path_f.close()

sh_cmd_l=[]
number_of_cores=48

# get number of files for each process
rp_number_per_process=len(raw_plus_paths_l)/number_of_cores
rp_number_per_process_list=[rp_number_per_process]*number_of_cores
rp_remainder=len(raw_plus_paths_l)%number_of_cores
for i in range(rp_remainder):
    rp_number_per_process_list[i]+=1

# split raw plus paths for multiprocess
rp_for_merge_split=[]
for i in range(number_of_cores):
    start_idx_rp=sum(rp_number_per_process_list[:i])
    last_idx_rp=start_idx_rp+rp_number_per_process_list[i]
    lb=raw_plus_paths_l[start_idx_rp:last_idx_rp]
    rp_for_merge_split.append(lb)

# get shell command list
for i in range(number_of_cores):
    rp_l=rp_for_merge_split[i]
    sh_cmd='python merge_genotype.py '+','.join(rp_l)+' '+str(i)
    sh_cmd_l.append(sh_cmd)

# execute multiprocess
pool=Pool(number_of_cores)
pool.map(os.system,sh_cmd_l)
pool.close()
pool.join()

# calculate time to complete this script
print('It took '+str(time.time()-start_time)+' seconds to complete the task.')




