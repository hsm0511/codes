import pandas as pd
import sys
import pingouin as pg
from scipy import stats
import numpy as np
import itertools
from matplotlib import pyplot as plt
import matplotlib.cm as cm
import seaborn as sns
import os
import math
from scipy.stats import norm
import seaborn as sns

def check_err_delete(vaf,normal_vaf,normal_num):
    p_val='None'
    if normal_vaf==0 or normal_num<30:
        remove_flag='pass'
    elif vaf<normal_vaf:
        remove_flag='fail'
    else:
        normal_std=math.sqrt(normal_vaf*(1-normal_vaf)/float(normal_num))
        z_val=(vaf-normal_vaf)/float(normal_std)
        p_val=1-norm.cdf(z_val)
        if p_val<0.05:
            remove_flag='pass'
        else:
            remove_flag='fail'
    return remove_flag,p_val

tumor_var_info=pd.read_csv(sys.argv[1],sep='\t',index_col=None)
normal_var_info=pd.read_csv(sys.argv[2],sep='\t',index_col=None)

err_remove_info={}
p_info={}
for pat in tumor_var_info['patient_id'].unique():
    pat_var=tumor_var_info[tumor_var_info['patient_id']==pat]
    for cp in pat_var['locus'].unique():
        loc_var=pat_var[pat_var['locus']==cp]
        normal_var=normal_var_info[normal_var_info['locus']==cp]
        maf=loc_var['maf'].values[0]
        detect_rat=loc_var['detect_ratio'].values[0]
        
        if len(normal_var)==0:
            normal_vaf=0
            normal_tot_depth='None'
        else:
            normal_vaf=normal_var['vaf'].values[0]
            normal_tot_depth=normal_var['depth'].values[0]
        err_remove,p_val=check_err_delete(maf,normal_vaf,normal_tot_depth)
        
        if detect_rat<0.5:
            err_remove='fail'
        
        err_remove_info[pat+'\t'+cp]=err_remove
        p_info[pat+'\t'+cp]=p_val

tumor_var_info['error_flag']=tumor_var_info.apply(lambda x:err_remove_info[x['patient_id']+'\t'+x['locus']],axis=1)
tumor_var_info['p_value']=tumor_var_info.apply(lambda x:p_info[x['patient_id']+'\t'+x['locus']],axis=1)

tumor_var_info = tumor_var_info[~tumor_var_info.duplicated(subset=['patient_id','locus'], keep='first')]

hue_order=['pass','fail']
fg = sns.FacetGrid(data=tumor_var_info, hue='error_flag',hue_order=hue_order, aspect=1.61)
fg.map(plt.scatter, 'detect_ratio', 'maf',s=10).add_legend()
plt.savefig('error_suppress_prob_dist.png')
tumor_var_info.to_csv('tumor_somatic_var_info_add_error.txt',sep='\t',index=None)








