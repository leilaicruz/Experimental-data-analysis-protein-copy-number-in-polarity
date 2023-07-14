# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.10.3
#   kernelspec:
#     display_name: Python 3.9.7
#     language: python
#     name: python3
# ---

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

data = pd.read_excel("DATA", sheet_name="Multinucleated cell ratio")

# +
data.columns = [
    "Strain",
    "Condition",
    "total_number_of_cells",
    "cells_two_nuclei",
    "cells_three_nuclei",
    "cells_dead",
]

data.loc[:, "ratio_two_nuclei"] = data.loc[:, "cells_two_nuclei"] / (
    data.loc[:, "total_number_of_cells"] - data.loc[:, "cells_dead"]
)
data.loc[:, "ratio_three_nuclei"] = data.loc[:, "cells_three_nuclei"] / (
    data.loc[:, "total_number_of_cells"] - data.loc[:, "cells_dead"]
)
data.loc[:, "ratio_dead"] = (
    data.loc[:, "cells_dead"] / data.loc[:, "total_number_of_cells"]
)
data.loc[:, "alive_cells"] = (
    data.loc[:, "total_number_of_cells"] - data.loc[:, "cells_dead"]
)
# -

dbem1 = data[data["Strain"] == "ywkd069"]
dbem1dbem3 = data[data["Strain"] == "ywkd070"]
wt = data[data["Strain"] == "yll3a"]
wg = data[data["Strain"] == "ywkd065"]

# +
plt.figure(figsize=(5, 5))
colors = ["lightgray", "blue", "blue"]
sns.stripplot(
    x="Condition",
    y="ratio_two_nuclei",
    data=dbem1,
    jitter=True,
    size=8,
    palette=colors,
    alpha=0.6,
)

colors = ["lightgray", "#600C3D", "#600C3D"]
sns.stripplot(
    x="Condition",
    y="ratio_two_nuclei",
    data=dbem1dbem3,
    jitter=True,
    size=8,
    palette=colors,
    alpha=0.6,
)

colors = ["green", "green", "green"]
sns.stripplot(
    x="Condition",
    y="ratio_two_nuclei",
    data=wg,
    jitter=True,
    size=8,
    palette=colors,
    alpha=0.6,
)

colors = ["red", "red", "red"]
sns.stripplot(
    x="Condition",
    y="ratio_two_nuclei",
    data=wt,
    jitter=True,
    size=8,
    palette=colors,
    alpha=0.6,
)

plt.xlabel("Galactose (%)")
plt.ylabel("Ratio of multinucleated cells")

plt.savefig("ratio_two_nuclei_striplot.png", dpi=300, bbox_inches="tight")
