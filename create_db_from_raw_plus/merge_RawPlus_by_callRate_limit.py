import sys
import os
import subprocess
import datetime
import time

start_time=time.time()
script_path='/clinix/Analysis/Projects/PersonalGenomics_RnD/add_DB_from_RawPlus/Scripts/'
check_sample_f='/clinix/Analysis/Projects/PersonalGenomics_RnD/add_DB_from_RawPlus/merged_samples.txt'
task_path='/clinix/Analysis/Projects/PersonalGenomics_RnD/add_DB_from_RawPlus/task/'
merged_path_written_dir='/clinix/Analysis/Projects/PersonalGenomics_RnD/add_DB_from_RawPlus/merged_Raw_Plus_path/'

already_inserted_samples=[]
try:
    check_f=open(check_sample_f,'r')
    for line in check_f:
        if line.rstrip()!='':
            already_inserted_samples.append(line.rstrip())
    check_f.close()
except:
    pass

today_insert_samples=[]
def select_RawPlus_by_callRate(CHQ_path,call_rate_limit):
    selected_rawplus_paths=[]
    selected_sample_l=[]
    chq_num=CHQ_path.split('/')[-2]
    sample_table_f=open(CHQ_path+'Output/'+chq_num+'_Samples_Table_Info.txt','r')
    for line in sample_table_f:
        if line.rstrip()!='':
            l=line.rstrip().split('\t')
            if 'Call Rate' in line:
                call_rate_idx=l.index('Call Rate')
                manifest_idx=l.index('Manifest file')
                sample_name_idx=l.index('Sample ID')
            else:
                call_rate=float(l[call_rate_idx])
                manifest=l[manifest_idx]
                if manifest=='Macrogen_chip_L1_20039998_B1.bpm':
                    return []
                sample_name=l[sample_name_idx]
                if call_rate >= call_rate_limit:
                    selected_sample_l.append(sample_name)
                    print(sample_name+'\t'+str(call_rate))
    sample_table_f.close()
    raw_plus_upper_path=CHQ_path+'Output/Plus_Genotype/'
    ls_result=subprocess.check_output('ls '+raw_plus_upper_path+' | grep Raw_Plus',shell=True)
    for rp_file in ls_result.rstrip().split('\n'):
        sample_name=rp_file.split('_Raw')[0]
        if sample_name in selected_sample_l and sample_name not in today_insert_samples and sample_name not in already_inserted_samples:
            selected_rawplus_paths.append(raw_plus_upper_path+rp_file)
            today_insert_samples.append(sample_name)
    return selected_rawplus_paths

def get_chq_list(former_days):
    today=datetime.date.today()
    CHQs_cmd='ls -d /mnt/clinix/Analysis/Projects/PersonalGenomics/*'
    CHQs_cmd_output=subprocess.check_output(CHQs_cmd,shell=True)
    CHQs_list=CHQs_cmd_output.split('\n')[:-1]
    chq_list_to_get_RawPlus=[]
    for chq_path in CHQs_list:
        chq_num=chq_path.split('/')[-1]
        if 'CHQ' in chq_num:
            chq_date=chq_num[3:9]
            chq_year=int('20'+chq_date[:2])
            chq_month=int(chq_date[2:4])
            chq_day=int(chq_date[4:])
            try:
                day_to_calcul=datetime.date(chq_year,chq_month,chq_day)
            except:
                continue
            diff=today-day_to_calcul
            diff_num=diff.days
            if diff_num > 0 and diff_num <= former_days:
                chq_list_to_get_RawPlus.append(chq_path)
    return chq_list_to_get_RawPlus



today='%s'%(datetime.datetime.today().strftime("%Y%m%d"))
fw=open(merged_path_written_dir+today+'_Raw_Plus_paths.txt','w')
chq_list_to_get_sample=get_chq_list(7)
raw_plus_paths_to_merge=[]
for chq_path in chq_list_to_get_sample:
    if chq_path[-1]!='/':
        chq_path=chq_path+'/'
    raw_plus_paths_to_merge.extend(select_RawPlus_by_callRate(chq_path,0.99))
if len(raw_plus_paths_to_merge)==0:
    print('no sample which has call rate greater than or equal to 0.99')
    fw.close()
    sys.exit(0)
else:
    for rp in raw_plus_paths_to_merge:
        fw.write(rp+'\n')
    fw.close()

total_insert_samples=already_inserted_samples+today_insert_samples
check_f=open(check_sample_f,'w')
for s in total_insert_samples:
    check_f.write(s+'\n')
check_f.close()

# scripts
merge_script=script_path+'get_raw_plus_and_merge.py'
modify_gt_data_script=script_path+'split_and_transpose.sh'
extend_tmp_to_tsv=script_path+'extend_tsv.py'
raw_plus_file_path_f=merged_path_written_dir+today+'_Raw_Plus_paths.txt'

# merge Raw_Plus.txt
os.system('/mnt/clinix/Tools/Dev/python/bin/python '+merge_script+' '+raw_plus_file_path_f+' '+task_path)

# split rows and transpose the merged files
os.system('sh '+modify_gt_data_script+' '+task_path)

# merge *_T_tmp.tsv to *_T.tsv
os.system('/mnt/clinix/Tools/Dev/python/bin/python '+extend_tmp_to_tsv+' '+task_path)

# remove *_T_tmp.tsv files
os.system('rm '+task_path+'*_T_tmp.tsv')

# calculate time to complete this script
print('It took '+str(time.time()-start_time)+' seconds to complete the task.')
