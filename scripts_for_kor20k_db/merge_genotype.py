import os, sys, re, getopt, datetime
import pandas as pd
from pymongo import MongoClient
import json

# get Raw_Plus file names to merge 
gt_files = sys.argv[1]
gt_files_l=gt_files.rstrip().split(',')

# final merge file's index number
file_num=sys.argv[2]

# insert first Raw_Plus file and get header
merged_data = pd.read_table(gt_files_l[0], sep="\t")
header_column=merged_data ['Name'].str.replace('.','_')
print(header_column)

# merge Raw_Plus to first file except first column
for gt_f in gt_files_l[1:]:
    gsa_data = pd.read_table(gt_f, sep="\t")
    sample_name=gsa_data.columns.tolist()[1]
    gts=gsa_data[sample_name]
    merged_data[sample_name]=gts.tolist()

# read csv file for merged data
merged_data.to_csv('output_gt_'+str(file_num)+'.tsv',sep='\t',header=True,index=None)

