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

normal_var_info=pd.read_csv(sys.argv[1],sep='\t',index_col=None) #locus    maf    depth    sample_id
normal_var_info['mut_depth']=normal_var_info.apply(lambda x:x['depth']*x['maf'],axis=1)
normal_prob_d={}

for cp in normal_var_info['locus'].unique():
    loc_var=normal_var_info[normal_var_info['locus']==cp]
    tot_depth=loc_var['depth'].sum()
    loc_maf=loc_var['mut_depth'].sum()/float(tot_depth)
    normal_prob_d[cp]={'vaf':loc_maf,'depth':tot_depth}

normal_prob_dist=pd.DataFrame(normal_prob_d).T
normal_prob_dist['locus']=normal_prob_dist.index

normal_prob_dist.to_csv('normal_var_vaf_dp.txt',sep='\t',index=None)






