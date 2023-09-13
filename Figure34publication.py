# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.10.3
#   kernelspec:
#     display_name: Python 3.6.9 64-bit
#     language: python
#     name: python3
# ---

from scipy.io import loadmat
import scipy.io
#import seaborn as sns
from scipy import stats, optimize, interpolate
import pandas as pd
import numpy as np
from collections import defaultdict 
import math
import matplotlib.pyplot as plt
from scipy.stats import norm, lognorm
from scipy import stats
import matplotlib.cm as cm
import matplotlib.mlab as mlab
import os, fnmatch
# %matplotlib inline

# +
data=loadmat('rel_fit.mat')
low_error = loadmat('error_lower_bands.mat')
upper_error = loadmat('error_upper_bands.mat')

low_error=low_error['rel_fit_lb']
upper_error=upper_error['rel_fit_ub']
data=data['rel_fit']

# +
conditions=['0% Gal','0.01% Gal','0.015% Gal','0.02% Gal','0.03% Gal','0.04% Gal','0.05% Gal',
   '0.06% Gal','0.08% Gal','0.1% Gal','0.2% Gal','2% Gal']
   
data_new=data.copy()
low_error_new=low_error.copy()
upper_error_new=upper_error.copy()

np.nan_to_num(data_new,nan=0)
np.nan_to_num(low_error_new,nan=0)
np.nan_to_num(upper_error_new,nan=0)


# -

conditions=['0% Gal','0.01% Gal','0.015% Gal','0.02% Gal','0.03% Gal','0.04% Gal','0.05% Gal',
   '0.06% Gal','0.08% Gal','0.1% Gal','0.2% Gal','2% Gal']


# +
x_wt=np.arange(0,len(conditions))
y_wt= data_new[:,0]
yerr_wt=[upper_error_new[:,0],low_error_new[:,0]]
x=np.arange(0,len(conditions))
y=data_new[:,1]
yerror=[upper_error_new[:,1],low_error_new[:,1]]

#dbem1
x_dbem1=np.arange(0,len(conditions))
y_dbem1= data_new[:,2]
yerr_dbem1=[upper_error_new[:,2],low_error_new[:,2]]

#dbem1dbem3
x_dbem1dbem3=np.arange(0,len(conditions))
y_dbem1dbem3= data_new[:,3]
yerr_dbem1dbem3=[upper_error_new[:,3],low_error_new[:,3]]

#pos=[0,1,3,6,7,8,9,10,11]
pos=np.arange(0,len(conditions))

x_wt_pos=[]
x_pos=[]
x_dbem1_pos=[]
x_dbem1dbem3_pos=[]
y_wt_pos=[]
y_pos=[]
y_dbem1_pos=[]
y_dbem1dbem3_pos=[]

yerr_wt_pos=[]
yerror_pos=[]
yerr_dbem1_pos=[]
yerr_dbem1dbem3_pos=[]
conditions_pos=[]

for i in pos:
    x_wt_pos.append(x_wt[i])
    x_pos.append(x[i])
    x_dbem1_pos.append(x_dbem1[i])
    x_dbem1dbem3_pos.append(x_dbem1dbem3[i])
    y_wt_pos.append(y_wt[i])
    y_pos.append(y[i])
    y_dbem1_pos.append(y_dbem1[i])
    y_dbem1dbem3_pos.append(y_dbem1dbem3[i])
    
    yerr_wt_pos.append([yerr_wt[0][i],yerr_wt[1][i]])
    yerror_pos.append([yerror[0][i],yerror[1][i]])
    yerr_dbem1_pos.append([yerr_dbem1[0][i],yerr_dbem1[1][i]])
    yerr_dbem1dbem3_pos.append([yerr_dbem1dbem3[0][i],yerr_dbem1dbem3[1][i]])
    
    conditions_pos.append(conditions[i])

# +
#x0=np.arange(0,0.11,0.01)
x0=[0,0.01,0.015,0.02,0.03,0.04,0.05,0.06,0.08,0.1,0.2,0.25]
#y0=np.array([0.2,0.25])

d001=defaultdict(dict)
d065=defaultdict(dict)
d069=defaultdict(dict)
d070=defaultdict(dict)

for i in np.arange(0,len(x0)):
    
    d001[x0[i]]=y_wt_pos[i]

    d065[x0[i]]=y_pos[i]
    
    d069[x0[i]]=y_dbem1_pos[i]
    
    d070[x0[i]]=y_dbem1dbem3_pos[i]






# +
## replacing nan values with zeros , brut force

d065[0]=0
d069[0]=0
d069[0.01]=0
d069[0.015]=0
d069[0.02]=0
d069[0.04]=0
d069[0.05]=0
d069[0.08]=0
d070[0]=0
d070[0.01]=0
d070[0.015]=0
d070[0.02]=0
d070[0.03]=0
d070[0.04]=0
d070[0.05]=0



# +
fig, axes = plt.subplots(nrows=1, ncols=1, figsize=(9, 6))

plt.errorbar(*zip(*sorted(d001.items())),yerr=np.transpose(yerr_wt_pos),marker='o',color='#ff4444',alpha=0.8,capsize=10,linewidth=0.4,label='CDC42 BEM1 BEM3'),
plt.errorbar(*zip(*sorted(d065.items())),yerr=np.transpose(yerror_pos),marker='o',color="#009245",
alpha=0.8,capsize=10,linewidth=0.4,label='pGal1-CDC42-sfGFP$^{SW}$ BEM1 BEM3'),
plt.errorbar(*zip(*sorted(d069.items())),yerr=np.transpose(yerr_dbem1_pos),marker='o',color='#0400FF',alpha=0.8,capsize=10,linewidth=0.4,label='pGal1-CDC42-sfGFP$^{SW}$ $\Delta$bem1 BEM3'),
plt.errorbar(*zip(*sorted(d070.items())),yerr=np.transpose(yerr_dbem1dbem3_pos),marker='o',color='#600C3D',alpha=0.8,capsize=10,linewidth=0.4,label='pGal1-CDC42-sfGFP$^{SW}$ $\Delta$bem1$\Delta$bem3'),


plt.ylim([0,2])
plt.legend(bbox_to_anchor=(-0.08, 1.05, 1.05, 1), loc=3,ncol=2, mode="expand", borderaxespad=0.,
fontsize=14)
plt.ylabel('Relative fitness compared with WT',fontsize=14)
plt.xlabel('Galactose concentrations',fontsize=14)
#plt.xticks([x0[0],x0[1],x0[2],x0[5],x0[6],x0[8],x0[10],y0[0],y0[1]],rotation=50);
xt=[0,0.01,0.015,0.02,0.03,0.04,0.05,0.06,0.08,0.1,0.2,0.25]
#plt.xticks(xt,rotation=50)
plt.xticks(ticks=xt,labels=conditions,rotation=90,fontsize=14)
plt.yticks(fontsize=14)

plt.tight_layout()

fig.savefig("Fig-3A-python.svg",dpi=300)
