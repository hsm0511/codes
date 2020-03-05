import os, sys, re, getopt, datetime
from pymongo import MongoClient
import subprocess
import logging, logging.handlers
import copy
import time

start_time=time.time()

mongo = MongoClient (host='172.19.87.50', port=27017)
dbname = mongo ['Repository']
gt_coll = dbname ['kor20k_split_8_probe_num']
probe_coll=dbname ['MGv2_probe_QC']

probe_name=sys.argv[1]
probe_name_rp=probe_name.replace('.','_')
header='Genotype\tSample Name\tProbe Name\n'

result_f=open('gt_sample_'+probe_name+'.txt','w')
result_f.write(header)

num_find=probe_coll.find_one({'Name':probe_name_rp},{'probe_num_DB':1})
probe_num=num_find['probe_num_DB']
gt_find=gt_coll.find({'genotypes_indexed_num':{'$elemMatch':{'probe_num': probe_num }}},{'genotypes_indexed_num.$':1,'name':1})
for g_obj in gt_find:
    sample_name=str(g_obj['name'])
    gt=g_obj['genotypes_indexed_num'][0]['gt']
    result_f.write(gt+'\t'+sample_name+'\t'+probe_name+'\n')
    print(gt+'\t'+sample_name+'\t'+probe_name)

result_f.close()

print('It took '+str(time.time()-start_time)+' seconds to complete task.')