import os, sys, re, getopt, datetime, json
from pymongo import MongoClient
import logging, logging.handlers
import copy

connection = MongoClient('172.19.87.50', 27017)
db = connection.Repository
sCollection = db['cathy_test']

gt_find=sCollection.find({})

# make probe name field nested to genotypes field
for gt_obj in gt_find:
    gt_nested_dict={}
    print(gt_obj['Name'])
    for field in gt_obj:
        if field=='_id':
            id_val=gt_obj[field]
            print(str(id_val))
        elif field=='Name':
            pass
        else:
            gt_nested_dict[field]=gt_obj[field]
    sCollection.update({'_id':id_val},{'$set':{'genotypes':gt_nested_dict}})





