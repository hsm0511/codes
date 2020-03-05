import os, sys, re, getopt, datetime, time
import pymongo
import subprocess
import logging, logging.handlers

start_time=time.time()

# DB info
connection = pymongo.MongoClient('172.19.87.50', 27017)
db = connection.Repository
gt_coll = db['kor20k_split_8']
db_cathy=connection.kor20k_split_test
num_colls=int(sys.argv[1])
avoid_time_start=int(sys.argv[2])
avoid_time_end=int(sys.argv[3])

col_names=[]
for i in range(num_colls):
    ii=i+1
    col_name='kor20k_split_8_'+'num_'+str(ii)
    col_names.append(col_name)

def split_colls_list(gt_len,num_colls):
    q=gt_len/num_colls
    rm=gt_len%num_colls
    split_probes_num_l=[q]*num_colls
    for i in range(rm):
        split_probes_num_l[i]+=1
    return split_probes_num_l

gt_find=gt_coll.find({},no_cursor_timeout=True)
pass_idx=1
pass_limit=60645
for g_obj in gt_find:
    # pass for already inserted docs
    if pass_idx <= 60645:
        print('insert pass : '+str(pass_idx))
        pass_idx+=1
        continue
    # avoid work time for insertion
    avoid_start=time.time()
    now_time=int(str(datetime.datetime.now().time()).replace(':','')[:4])
    while now_time >= avoid_time_start and now_time <= avoid_time_end:
        now_time=int(str(datetime.datetime.now().time()).replace(':','')[:4])
        avoiding_t=time.time()
        diff_t=avoiding_t-avoid_start
        if diff_t >=3600:
            print('work hour avoiding : '+str(now_time))
            avoid_start=time.time()
    sample_name=str(g_obj['name'])
    gt_idx_f=g_obj['genotypes_indexed']
    gt_array_len=len(gt_idx_f)
    probes_num_split_l=split_colls_list(gt_array_len,num_colls)
    start_idx=0
    for i in range(num_colls):
        insert_col=db_cathy[col_names[i]]
        p_num=probes_num_split_l[i]
        last_idx=start_idx+p_num
        insert_l=gt_idx_f[start_idx:last_idx]
        
        insert_col.insert({'name':sample_name,'genotypes_indexed':insert_l})
        print('probe gt inserted '+str(p_num)+'\t'+sample_name)
        start_idx=last_idx
        print('col inserted time : '+str(time.time()))


print('It takes '+time.time()-start_time+' seconds to complete the task.')
