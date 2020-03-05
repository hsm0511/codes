import os, sys, re, getopt, datetime, json, time
from pymongo import MongoClient
import pymongo
import logging, logging.handlers
import copy

start_time=time.time()
connection = MongoClient('172.19.87.50', 27017)
db = connection.Repository
sCollection = db['cathy_test']
copy_col=db['cathy_test_merge_1']
merge_col=db['cathy_test_merge_2']

# get number of group which is inserted eight times by each sample
db_insert_log=open('insert_document_numbers.txt','r')
repeated_nums=[]
for line in db_insert_log:
    if line.rstrip()!='':
        repeated_num=line.rstrip().split('imported ')[1][:3]
        repeated_nums.append(repeated_num)
db_insert_log.close()

# get repeated group number
repeated_nums_set=[]
for i in range(len(repeated_nums)/8):
    k=8*i
    repeated_nums_set.append(repeated_nums[k])


def merge_aggre_obj(inserted_names):
    for name in inserted_names:
        # collecting genotypes for same sample name
        aggre_obj=db.cathy_test_merge_1.aggregate([
            { '$match' : {'Name':name}} ,
            { '$group' : {'_id':'$Name', 'genotypes': {'$mergeObjects':"$genotypes"}}}
        ])
        # insert aggregated genotypes object to new collection
        merge_col.insert(aggre_obj)

# get sample name and genotypes from DB
gt_find=sCollection.find({'genotypes':{'$exists':True}},{'Name':1,'genotypes':1}, no_cursor_timeout=True)
ct=0
repeat_idx=0
inserted_Name_l=[]
for g_obj in gt_find:
    repeat_num=repeated_nums_set[repeat_idx]
    # copy documents to temporary collection
    copy_col.insert(g_obj)
    sample_name=g_obj['Name']
    inserted_Name_l.append(sample_name)
    print('inserted '+str(sample_name))
    ct+=1
    # insertion complete for a group which includes eight documents for each sample
    if ct==int(repeat_num)*8:
        # inserted samples unique
        inserted_Name_l=list(set(inserted_Name_l))
        print(inserted_Name_l)
        # aggregate genotypes and insert to new coll
        merge_aggre_obj(inserted_Name_l)
        # remove documents in temporary coll
        copy_col.remove({'Name':{'$exists':True}})
        # initialize inserted sample list and insert count, but add count to group index
        inserted_Name_l=[]
        repeat_idx+=1
        ct=0

gt_find.close()
    
print('It took '+str(time.time()-start_time)+' seconds to complete task.')
