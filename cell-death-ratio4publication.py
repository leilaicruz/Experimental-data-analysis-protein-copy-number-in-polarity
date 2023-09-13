# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.10.3
#   kernelspec:
#     display_name: Python 3.8.10 ('satay-dev')
#     language: python
#     name: python3
# ---

# +
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
import matplotlib

plt.rc("font", family="serif", size=14)
plt.rc("xtick", labelsize=14)
plt.rc("ytick", labelsize=14)
# -

## Reading the data
data = pd.read_excel("DATA", sheet_name="Cell death ratio")


# +
dbem1 = data[data["Strain"] == "dbem1"]
dbem1dbem3 = data[data["Strain"] == "bem1-3d"]
wg = data[data["Strain"] == "WT-pgal"]

dbem1_mean = dbem1.groupby(["Condition"]).mean().loc[:, "Ratio dead"]
dbem1dbem3_mean = dbem1dbem3.groupby(["Condition"]).mean().loc[:, "Ratio dead"]
wg_mean = wg.groupby(["Condition"]).mean().loc[:, "Ratio dead"]

dbem1_std = dbem1.groupby(["Condition"]).std().loc[:, "Ratio dead"]
dbem1dbem3_std = dbem1dbem3.groupby(["Condition"]).std().loc[:, "Ratio dead"]
wg_std = wg.groupby(["Condition"]).std().loc[:, "Ratio dead"]

# +
plt.figure(figsize=(5, 5))


sns.stripplot(
    x="Condition",
    y="Ratio dead",
    data=dbem1,
    jitter=True,
    size=10,
    color="blue",
    alpha=0.6,
)
sns.stripplot(
    x="Condition",
    y="Ratio dead",
    data=dbem1dbem3,
    jitter=True,
    size=10,
    color="#600C3D",
    alpha=0.6,
)
sns.stripplot(
    x="Condition",
    y="Ratio dead",
    data=wg,
    jitter=True,
    size=10,
    color="green",
    alpha=0.6,
)
plt.ylim(0, 1)
plt.xlabel("Galactose (%)")

plt.tight_layout()
