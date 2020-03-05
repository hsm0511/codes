import os, sys, re, getopt, datetime
from pymongo import MongoClient
import subprocess
import logging, logging.handlers
import copy
import time
import pandas as pd
from multiprocessing import Process, Manager, Pool
import glob

start_time=time.time()

mongo = MongoClient (host='172.19.87.50', port=27017)
db = mongo ['Repository']
db_split = mongo ['kor20k_split_test']
probe_coll=db ['MGv2_probe_QC']

def split_tasks(number_of_cores,task_l):
    task_split_l=[]
    q=len(task_l)/number_of_cores
    r=len(task_l)%number_of_cores
    num_task_l=[q]*number_of_cores
    for i in range(r):
        num_task_l[i]+=1
    start_idx=0
    for i in range(number_of_cores):
        last_idx=start_idx+num_task_l[i]
        tasks_per_num=task_l[start_idx:last_idx]
        task_split_l.append(tasks_per_num)
        start_idx=last_idx
    return task_split_l

def get_probes_encoding_gene(gene):
    col_p={}
    p_find=probe_coll.find({'GeneSymbol':gene,'probe_col':{'$exists':True}},{'Name':1,'probe_col':1})
    for p_obj in p_find:
        probe_name=p_obj['Name']
        probe_name=probe_name.replace('.','_')
        col_name=p_obj['probe_col']
        if col_name not in col_p:
            col_p[col_name]=[]
        col_p[col_name].append(probe_name)
    return col_p

def pon_dict_to_string(pon_dict):
    pon_l=[]
    for gt in pon_dict:
        s=gt+':'+str(pon_dict[gt])
        pon_l.append(s)
    return '|'.join(pon_l)

gene_name=sys.argv[1]
probes_for_gene=get_probes_encoding_gene(gene_name)
task_script='get_stat_from_probes.py'
print(probes_for_gene)
num_cores=48
colls_for_probe=probes_for_gene.keys()
if len(colls_for_probe) >= num_cores:
    task_colls_splitted=split_tasks(num_cores,colls_for_probe)
else:
    num_cores=len(colls_for_probe)
    task_colls_splitted=split_tasks(num_cores,colls_for_probe)

task_commands_l=[]
file_num=1
for col_l in task_colls_splitted:
    cmd='python '+task_script+' '+gene_name+' '+','.join(col_l)+' '+str(file_num)
    task_commands_l.append(cmd)
    file_num+=1

# execute multiprocess
pool=Pool(num_cores)
pool.map(os.system,task_commands_l)
pool.close()
pool.join()

result_file=open(gene_name+'_kor20k_stat.txt','w')
output_f_lines=[]
# merge files row wise
split_files = glob.glob("/home/hsm0511/"+gene_name+'_probe_stat*.txt')
idx=0
for f in split_files:
    ff=open(f,'r')
    lines=ff.readlines()
    if idx==0:
        header='Probe\tGT_stat\n'
        result_file.write(header)
    output_f_lines.extend(lines[1:])
    ff.close()
    idx+=1
result_file.write(''.join(output_f_lines))
result_file.close()
os.system('rm '+gene_name+'_probe_stat_*')

print('It took '+str(time.time()-start_time)+' seconds to complete task.')