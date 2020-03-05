import sys
import os

script_path='/clinix/Analysis/Projects/PersonalGenomics_RnD/add_DB_from_RawPlus/'

# scripts
merge_script=script_path+'get_raw_plus_and_merge.py'
modify_gt_data_script=script_path+'split_and_transpose.sh'
mongo_insert_script=script_path+'mongo_insert_GT.sh'
split_genotype_350_collections_script=script_path+'split_colls_time_avoiding.py'

raw_plus_file_path_f=sys.argv[1]
task_path=sys.argv[2]   # the directory to execute task must be empty

# merge Raw_Plus.txt
os.system('python '+merge_script+' '+raw_plus_file_path_f+' '+task_path)

# split rows and transpose the merged files
os.system('python '+modify_gt_data_script+' '+task_path)

# create raw DB by mongoimport after changing dot in probe name to under bar
os.system('sh '+mongo_insert_script+' '+task_path)

# split data to 350 collections
os.system('python '+split_genotype_350_collections_script+' 350 0800 1900')



