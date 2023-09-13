# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.10.3
#   kernelspec:
#     display_name: Python 3.9.7 ('transposonmapper')
#     language: python
#     name: python3
# ---

import scipy.io
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
import seaborn as sns
# %matplotlib inline
plt.rc('font', family='serif',size=14)
plt.rc('xtick',labelsize=14)
plt.rc('ytick',labelsize=14)

# +
## Reading the data 

data=pd.read_excel('DATA',sheet_name='Cell radii')
# -

data_bem1_all=data[data['strain']=='bem1d']
data_bem13_all=data[data['strain']=='bem1-3d']
data_wt_all=data[data['strain']=='WT']
data_wtg_all=data[data['strain']=='WT-pgal']


# +
# make boxplots of bem1d data across all three conditions by pulling all experiments together in one figure


data_bem1_all_0=data_bem1_all.query('condition==0')['cell radius'].tolist()
data_bem1_all_06=data_bem1_all.query('condition==0.06')['cell radius'].tolist()
data_bem1_all_1=data_bem1_all.query('condition==0.1')['cell radius'].tolist()

data_bem13_all_0=data_bem13_all.query('condition==0')['cell radius'].tolist()
data_bem13_all_06=data_bem13_all.query('condition==0.06')['cell radius'].tolist()
data_bem13_all_1=data_bem13_all.query('condition==0.1')['cell radius'].tolist()

data_wt_all_0=data_wt_all.query('condition==0')['cell radius'].tolist()
data_wt_all_06=data_wt_all.query('condition==0.06')['cell radius'].tolist()
data_wt_all_1=data_wt_all.query('condition==0.1')['cell radius'].tolist()

data_wtg_all_0=data_wtg_all.query('condition==0')['cell radius'].tolist()
data_wtg_all_06=data_wtg_all.query('condition==0.06')['cell radius'].tolist()
data_wtg_all_1=data_wtg_all.query('condition==0.1')['cell radius'].tolist()

# +
## Taken the sem of the cell radius for each condition and strain. There is no data for bem1-3d at 0.06% galactose and 0.1% galactose for 
# the first experiment.

data_bem1_all_0_mean_1=np.mean(data_bem1_all_0[0])
data_bem1_all_0_mean_2=np.mean(data_bem1_all_0[1])
data_bem1_all_0_mean_3=np.mean(data_bem1_all_0[2])


data_bem1_all_0_std_1=np.std(data_bem1_all_0[0])
data_bem1_all_0_std_2=np.std(data_bem1_all_0[1])
data_bem1_all_0_std_3=np.std(data_bem1_all_0[2])


data_bem1_all_06_mean_2=np.mean(data_bem1_all_06[1])
data_bem1_all_06_mean_3=np.mean(data_bem1_all_06[2])



data_bem1_all_06_std_2=np.std(data_bem1_all_06[1])
data_bem1_all_06_std_3=np.std(data_bem1_all_06[2])
data_bem1_all_1_mean_2=np.mean(data_bem1_all_1[1])
data_bem1_all_1_mean_3=np.mean(data_bem1_all_1[2])


data_bem1_all_1_std_2=np.std(data_bem1_all_1[1])
data_bem1_all_1_std_3=np.std(data_bem1_all_1[2])

data_bem13_all_0_mean_1=np.mean(data_bem13_all_0[0])
data_bem13_all_0_mean_2=np.mean(data_bem13_all_0[1])
data_bem13_all_0_mean_3=np.mean(data_bem13_all_0[2])



data_bem13_all_0_std_1=np.std(data_bem13_all_0[0])
data_bem13_all_0_std_2=np.std(data_bem13_all_0[1])
data_bem13_all_0_std_3=np.std(data_bem13_all_0[2])


data_bem13_all_06_mean_1=np.mean(data_bem13_all_06[0])
data_bem13_all_06_mean_2=np.mean(data_bem13_all_06[1])
data_bem13_all_06_mean_3=np.mean(data_bem13_all_06[2])



data_bem13_all_06_std_1=np.std(data_bem13_all_06[0])
data_bem13_all_06_std_2=np.std(data_bem13_all_06[1])
data_bem13_all_06_std_3=np.std(data_bem13_all_06[2])




data_bem13_all_1_mean_2=np.mean(data_bem13_all_1[1])
data_bem13_all_1_mean_3=np.mean(data_bem13_all_1[2])



data_bem13_all_1_std_2=np.std(data_bem13_all_1[1])
data_bem13_all_1_std_3=np.std(data_bem13_all_1[2])


data_wt_all_0_mean_1=np.mean(data_wt_all_0[0])
data_wt_all_0_mean_2=np.mean(data_wt_all_0[1])
data_wt_all_0_mean_3=np.mean(data_wt_all_0[2])



data_wt_all_0_std_1=np.std(data_wt_all_0[0])
data_wt_all_0_std_2=np.std(data_wt_all_0[1])
data_wt_all_0_std_3=np.std(data_wt_all_0[2])



data_wt_all_06_mean_2=np.mean(data_wt_all_06[1])
data_wt_all_06_mean_3=np.mean(data_wt_all_06[2])



data_wt_all_06_std_2=np.std(data_wt_all_06[1])
data_wt_all_06_std_3=np.std(data_wt_all_06[2])



data_wt_all_1_mean_2=np.mean(data_wt_all_1[1])
data_wt_all_1_mean_3=np.mean(data_wt_all_1[2])


data_wt_all_1_std_2=np.std(data_wt_all_1[1])
data_wt_all_1_std_3=np.std(data_wt_all_1[2])


data_wtg_all_0_mean_1=np.mean(data_wtg_all_0[0])
data_wtg_all_0_mean_2=np.mean(data_wtg_all_0[1])
data_wtg_all_0_mean_3=np.mean(data_wtg_all_0[2])


data_wtg_all_0_std_1=np.std(data_wtg_all_0[0])
data_wtg_all_0_std_2=np.std(data_wtg_all_0[1])
data_wtg_all_0_std_3=np.std(data_wtg_all_0[2])


data_wtg_all_06_mean_2=np.mean(data_wtg_all_06[1])
data_wtg_all_06_mean_3=np.mean(data_wtg_all_06[2])

data_wtg_all_06_std_2=np.std(data_wtg_all_06[1])
data_wtg_all_06_std_3=np.std(data_wtg_all_06[2])

data_wtg_all_1_mean_2=np.mean(data_wtg_all_1[1])
data_wtg_all_1_mean_3=np.mean(data_wtg_all_1[2])


data_wtg_all_1_std_2=np.std(data_wtg_all_1[1])
data_wtg_all_1_std_3=np.std(data_wtg_all_1[2])









# +
exp_1_mean_dbem1=[data_bem1_all_0_mean_1,np.nan,np.nan]
exp_1_dbem1_std=[data_bem1_all_0_std_1,np.nan,np.nan]

exp_2_mean_dbem1=[data_bem1_all_0_mean_2,data_bem1_all_06_mean_2,data_bem1_all_1_mean_2]
exp_2_dbem1_std=[data_bem1_all_0_std_2,data_bem1_all_06_std_2,data_bem1_all_1_std_2]

exp_3_mean_dbem1=[data_bem1_all_0_mean_3,data_bem1_all_06_mean_3,data_bem1_all_1_mean_3]
exp_3_dbem1_std=[data_bem1_all_0_std_3,data_bem1_all_06_std_3,data_bem1_all_1_std_3]

exp_1_mean_dbem1dbem3=[data_bem13_all_0_mean_1,np.nan,np.nan]
exp_1_dbem1dbem3_std=[data_bem13_all_0_std_1,np.nan,np.nan]

exp_2_mean_dbem1dbem3=[data_bem13_all_0_mean_2,data_bem13_all_06_mean_2,data_bem13_all_1_mean_2]
exp_2_dbem1dbem3_std=[data_bem13_all_0_std_2,data_bem13_all_06_std_2,data_bem13_all_1_std_2]

exp_3_mean_dbem1dbem3=[data_bem13_all_0_mean_3,data_bem13_all_06_mean_3,data_bem13_all_1_mean_3]
exp_3_dbem1dbem3_std=[data_bem13_all_0_std_3,data_bem13_all_06_std_3,data_bem13_all_1_std_3]

exp_1_mean_wtg=[data_wtg_all_0_mean_1,np.nan,np.nan]
exp_1_wtg_std=[data_wtg_all_0_std_1,np.nan,np.nan]

exp_2_mean_wtg=[data_wtg_all_0_mean_2,data_wtg_all_06_mean_2,data_wtg_all_1_mean_2]
exp_2_wtg_std=[data_wtg_all_0_std_2,data_wtg_all_06_std_2,data_wtg_all_1_std_2]

exp_3_mean_wtg=[data_wtg_all_0_mean_3,data_wtg_all_06_mean_3,data_wtg_all_1_mean_3]
exp_3_wtg_std=[data_wtg_all_0_std_3,data_wtg_all_06_std_3,data_wtg_all_1_std_3]

exp_1_mean_wt=[data_wt_all_0_mean_1,np.nan,np.nan]
exp_1_wt_std=[data_wt_all_0_std_1,np.nan,np.nan]

exp_2_mean_wt=[data_wt_all_0_mean_2,data_wt_all_06_mean_2,data_wt_all_1_mean_2]
exp_2_wt_std=[data_wt_all_0_std_2,data_wt_all_06_std_2,data_wt_all_1_std_2]

exp_3_mean_wt=[data_wt_all_0_mean_3,data_wt_all_06_mean_3,data_wt_all_1_mean_3]
exp_3_wt_std=[data_wt_all_0_std_3,data_wt_all_06_std_3,data_wt_all_1_std_3]

# +
plt.figure(figsize=(5,5))

plt.xticks([0,1,2],['0','0.06','0.1'])

plt.xlabel('Galactose (%)',fontsize=14)
plt.ylabel('Cell radii ($\mu$m)',fontsize=14)
plt.errorbar([0,1,2],exp_1_mean_dbem1,yerr=exp_1_dbem1_std,fmt='o',color='blue',capsize=5,alpha=0.6)
plt.errorbar([0,1,2],exp_2_mean_dbem1,yerr=exp_2_dbem1_std,fmt='o',color='blue',capsize=5,alpha=0.6)
plt.errorbar([0,1,2],exp_3_mean_dbem1,yerr=exp_3_dbem1_std,fmt='o',color='blue',capsize=5,alpha=0.6)

plt.errorbar([0,1,2],exp_1_mean_dbem1dbem3,yerr=exp_1_dbem1dbem3_std,fmt='o',color='#600C3D',capsize=5,alpha=0.6)
plt.errorbar([0,1,2],exp_2_mean_dbem1dbem3,yerr=exp_2_dbem1dbem3_std,fmt='o',color='#600C3D',capsize=5,alpha=0.6)
plt.errorbar([0,1,2],exp_3_mean_dbem1dbem3,yerr=exp_3_dbem1dbem3_std,fmt='o',color='#600C3D',capsize=5,alpha=0.6)

plt.errorbar([0,1,2],exp_1_mean_wtg,yerr=exp_1_wtg_std,fmt='o',color='green',capsize=5,alpha=0.6)
plt.errorbar([0,1,2],exp_2_mean_wtg,yerr=exp_2_wtg_std,fmt='o',color='green',capsize=5,alpha=0.6)
plt.errorbar([0,1,2],exp_3_mean_wtg,yerr=exp_3_wtg_std,fmt='o',color='green',capsize=5,alpha=0.6)

plt.errorbar([0,1,2],exp_1_mean_wt,yerr=exp_1_wt_std,fmt='o',color='red',capsize=5,alpha=0.6)
plt.errorbar([0,1,2],exp_2_mean_wt,yerr=exp_2_wt_std,fmt='o',color='red',capsize=5,alpha=0.6)
plt.errorbar([0,1,2],exp_3_mean_wt,yerr=exp_3_wt_std,fmt='o',color='red',capsize=5,alpha=0.6)


