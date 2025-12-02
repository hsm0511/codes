import sys
from scipy.stats import norm
import math

# get model auc as input
auc=sys.argv[1]

# trueAUC is inside the interval sampleAUC-d ~ sample_AUC+d with probability of 90%
d=0.05
stat_power=0.9

# sample_control_ratio=treated_sample_size/control_sample_size
sample_control_ratio=1

alpha=1-stat_power
alpha_half=alpha/2
z_alpha_half=norm.ppf(1-alpha_half)

# binormal method for variance estimation by Obuchowski et al. (1998)
a = norm.ppf(auc)
v_auc=(0.0099*math.exp(-0.5*a*a))*((5*a*a+8)/sample_control_ratio+a*a+8)

control_sample_size=z_alpha_half*v_auc/sqrt(d)
treated_sample_size=control_sample_size*sample_control_ratio

print('control sample size : '+str(control_sample_size)+'\n'
    +'treated sample size : '+str(treated_sample_size))
