import os, sys, re, getopt, datetime
from pymongo import MongoClient
import subprocess
import logging, logging.handlers
import copy
import time
from multiprocessing import Process, Manager, Pool

start_time=time.time()

mongo = MongoClient (host='172.19.87.50', port=27017)
db = mongo ['Repository']
db_split = mongo ['kor20k_split_test']
probe_coll=db ['MGv2_probe_QC']

gene_name=sys.argv[1]
col_name_l=sys.argv[2].split(',')
file_num=sys.argv[3]

gt_ct={}
for col_name in col_name_l:
    p_find=probe_coll.find({'GeneSymbol':gene_name,'probe_col':col_name},{'Name':1})
    c=col_name
    probes_tmp=[p_obj['Name'] for p_obj in p_find]
    probes=list(map(lambda x:x.replace('.','_'),probes_tmp))
    gt_coll=db_split[c]
    aggre_result=gt_coll.aggregate([
        {'$match':{'genotypes_indexed':{'$elemMatch':{'probe_name':{'$in':probes}}}}},
        {'$unwind':'$genotypes_indexed'},
        {'$match':{'genotypes_indexed.probe_name':{'$in':probes}}},
        {'$group':{'_id':{'p':'$genotypes_indexed.probe_name','gt':'$genotypes_indexed.gt'},'count':{'$sum':1}}}
    ])
    print(c)
    for ag_obj in aggre_result:
        probe_name=ag_obj['_id']['p']
        gt=ag_obj['_id']['gt']
        ct=ag_obj['count']
        print(probe_name+'\t'+gt+'\t'+str(ct))
        if probe_name not in gt_ct:
            gt_ct[probe_name]={}
        gt_ct[probe_name][gt]=ct

gene_probe_stat_f=open(gene_name+'_probe_stat_'+file_num+'.txt','w')
header='Probe_name\tGenotype_Count\n'
for probe_name in gt_ct:
    pon_info=[]
    for gt in gt_ct[probe_name]:
        ct=gt_ct[probe_name][gt]
        pon_info.append(gt+':'+str(ct))
    pon_string='|'.join(pon_info)
    gene_probe_stat_f.write(probe_name+'\t'+pon_string+'\n')
gene_probe_stat_f.close()