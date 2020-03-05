import os, sys, re, getopt, datetime
from pymongo import MongoClient
import subprocess
import logging, logging.handlers
import copy
import time

start_time=time.time()

mongo = MongoClient (host='172.19.87.50', port=27017)
dbname = mongo ['Repository']
gt_coll = dbname ['kor20k_split_8']
probe_coll=dbname ['MGv2_probe_QC']

probe_name=sys.argv[1]
probe_name_rp=probe_name.replace('.','_')
header='Genotype\tSample Name\tProbe Name\n'

gt_ct={}
gt_find=gt_coll.find({'genotypes_indexed':{'$elemMatch':{'probe_name': probe_name_rp }}},{'genotypes_indexed.$':1})
for g_obj in gt_find:
    gt=g_obj['genotypes_indexed'][0]['gt']
    if gt not in gt_ct:
        gt_ct[gt]=0
    gt_ct[gt]+=1

print(gt_ct)

print('It took '+str(time.time()-start_time)+' seconds to complete task.')