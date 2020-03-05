import pymongo; import sys; import time; import os; import re


connection = pymongo.MongoClient('172.19.87.50', 27017)
db_R = connection.Repository
coll_gt = db_R['cathy_test']
copy_col=db_R['kor20k_split_8']

start_time=time.time()
gt_find=coll_gt.find({},no_cursor_timeout=True)

create_idx_flag=False
for g_obj in gt_find:
    sample_name=g_obj['Name']
    new_gt_field='genotypes_indexed'
    # make array form for each probe's genotype
    push_gt_list=[]
    for p in g_obj:
        if p=='_id' or p=='Name' or p=='genotypes':
            pass
        else:
            # change genotype field to probe_name, gt field and add to array
            gt_dict={'probe_name':p,'gt':g_obj[p]}
            push_gt_list.append(gt_dict)
    copy_col.insert({'name':sample_name,new_gt_field:push_gt_list})
    print(str(sample_name)+'genotypes_indexed field inserted'+' ----------- '+str(time.time()))
    # create index for probe_name, gt fields after inserting one document
    if create_idx_flag==False:
        copy_col.create_index(
            [("genotypes_indexed.probe_name", pymongo.ASCENDING), ("genotypes_indexed.gt", pymongo.ASCENDING)]
        )
        create_idx_flag=True

print('It took '+str(time.time()-start_time)+' seconds to complete task.')





