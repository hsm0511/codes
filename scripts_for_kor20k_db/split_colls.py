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
create_idx_flag=[False]*num_colls
for g_obj in gt_find:
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
        if create_idx_flag[i]==False:
            insert_col.create_index([("genotypes_indexed.probe_name", pymongo.ASCENDING), ("genotypes_indexed.gt", pymongo.ASCENDING)])
            create_idx_flag[i]=True
            print('create compound index : '+col_names[i])
        print('col inserted time : '+str(time.time()))


print('It takes '+time.time()-start_time+' seconds to complete the task.')
