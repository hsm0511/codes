import os, sys, re, getopt, datetime
import pandas as pd
from pymongo import MongoClient
import json

# get input file name and output file name
tsv_file = sys.argv[1]
ref_num=tsv_file.split('_')[1]
raw_file=sys.argv[2].split('.')[0]

# transpose content of file and write result to tsv
pd.read_csv(tsv_file, sep='\t',header=None).T.to_csv(raw_file+'_'+ref_num+'_T.tsv', sep='\t',header=False, index=False)
