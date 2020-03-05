import pymongo; import sys; import time; import os; import re


connection = pymongo.MongoClient('172.19.87.50', 27017)
db_R = connection.Repository
k_col = db_R['kor20k_split_8']
result_col=db_R['kor20k_split_8_probe_num']
num_col=db_R['MGv2_probe_QC']


start_time=time.time()

# get probe number
qc_find=num_col.find({},{'Name':1,'probe_num_DB':1})
qc_num={}
for q_obj in qc_find:
    probe_name=q_obj['Name']
    probe_name=probe_name.replace('.','_')
    probe_num=q_obj['probe_num_DB']
    qc_num[probe_name]=probe_num

gt_find=k_col.find({},no_cursor_timeout=True)

create_idx_flag=False
for g_obj in gt_find:
    sample_name=g_obj['name']
    gt_obj=g_obj['genotypes_indexed']
    new_gt_field='genotypes_indexed_num'
    push_gt_list=[]
    for p_gt in gt_obj:
        probe_name=p_gt['probe_name']
        gt=p_gt['gt']
        probe_num=qc_num[probe_name]
        gt_dict={'probe_num':probe_num,'gt':gt}
        push_gt_list.append(gt_dict)
    result_col.insert({'name':sample_name,new_gt_field:push_gt_list})
    print(str(sample_name)+'\t'+str(probe_num)+'\tgenotypes_indexed_num field inserted'+' ----------- '+str(time.time()))
    # create index
    if create_idx_flag==False:
        result_col.create_index(
            [("genotypes_indexed_num.probe_num", pymongo.ASCENDING), ("genotypes_indexed_num.gt", pymongo.ASCENDING)]
        )
        create_idx_flag=True

print('It took '+str(time.time()-start_time)+' seconds to complete task.')

