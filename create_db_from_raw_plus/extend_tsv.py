import sys
import subprocess
import os.path

def remove_blank_in_list(l):
    while '' in l:
        l.remove('')
    return l

def merge_tsv_by_row(tsv_tmp,tsv_final):
    merge_lines_l=[]
    tmp_f=open(tsv_tmp,'r')
    merge_lines_2_l=tmp_f.readlines()[1:]
    merge_lines_2_l=list(map(lambda x:x.strip('\n'),merge_lines_2_l))
    merge_lines_2_l=remove_blank_in_list(merge_lines_2_l)
    tmp_f.close()
    final_f=open(tsv_final,'r')
    merge_lines_1_l=final_f.readlines()
    merge_lines_1_l=list(map(lambda x:x.strip('\n'),merge_lines_1_l))
    merge_lines_1_l=remove_blank_in_list(merge_lines_1_l)
    final_f.close()
    merge_lines_l=merge_lines_1_l+merge_lines_2_l
    final_f=open(tsv_final,'w')
    final_f.write('\n'.join(merge_lines_l))
    final_f.close()

task_path=sys.argv[1]
ls_cmd='ls '+task_path+' | grep _T_tmp'
tmp_tsv_l=[]
try:
    ls_result=subprocess.check_output(ls_cmd,shell=True)
    tmp_tsv_l=ls_result.rstrip().split('\n')
except:
    pass

if len(tmp_tsv_l)==0:
    print('no extending')
    sys.exit(0)

for tmp in tmp_tsv_l:
    tmp_with_path=task_path+tmp
    final_tsv_name=tmp.split('_tmp')[0]+'.tsv'
    final_with_path=task_path+final_tsv_name
    if os.path.isfile(final_with_path):
        merge_tsv_by_row(tmp_with_path,final_with_path)
        print('----------> merging : '+tmp_with_path+'\t'+final_with_path)
    else:
        print('there is no '+final_with_path)
        sys.exit(0)

