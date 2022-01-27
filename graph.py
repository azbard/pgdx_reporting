# %%
import datetime as dt
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import glob
import os
import re
from reading import parse_csv
from UploadToSql import addIdentifiers, caseInfo

# %%


textfile_dfs = []

files = glob.glob(
    "/mnt/pgdx_v1/ElioConnect_Output/*/TextFiles/*.rawfoldchange_ETC-RUO.csv",
    recursive=True,
)

for file in files:

    df = pd.read_csv(file)
    ids = caseInfo(file)
    addIdentifiers(df, *ids)
    textfile_dfs.append(df)

textfile_df = pd.concat(textfile_dfs, ignore_index=True)
textfile_df

# %%

report_dfs = []

files = glob.glob("/mnt/pgdx_v1/ElioConnect_Output/*/Reports/*.csv", recursive=True,)

for file in files:
    # print(file)
    report_dict = parse_csv(file)
    # print(report_dict.keys())

    df = report_dict.get("Sample Amplification Analysis")
    if df is not None:
        if not df.empty:
            ids = caseInfo(file)
            addIdentifiers(df, *ids)
            report_dfs.append(df)

report_df = pd.concat(report_dfs)
report_df

# %%

df = textfile_df.merge(report_df, how="left", on=["cid_number", "na", "batch", "Gene"])
repl = "0"  # "0","0.5", pd.NA
df.RawFoldChange = df.RawFoldChange.replace("<1", repl).astype(float, errors="ignore")
df.Status = df.Status.fillna("None")


# %%

# purity = pd.read_csv('purity.csv')
# purity

# # %%
# df.merge(purity, how='left', left_on='na', right_on='NA')


# %%

g = sns.FacetGrid(df, hue="Status", col="Gene", height=4, col_wrap=4)
g.map(
    sns.scatterplot, "RawFoldChange", "PercentRegions",
)
g.set(xlim=(None, 7))
g.add_legend()

# g.savefig('cnv - '+dt.datetime.now().strftime('%Y%m%d%H%M%s')+'.svg')
# %%
g = sns.scatterplot(data=df, x="RawFoldChange", y="PercentRegions", hue="Status")

# plt.savefig('cnv - '+dt.datetime.now().strftime('%Y%m%d%H%M%s')+'.svg')


# %% replace rawfoldchange = "<1" with what?


# %%
fig, ax = plt.subplots(figsize=(12, 3), dpi=300)
vp = sns.violinplot(x="Gene", y="RawFoldChange", data=df, ax=ax, orient="v")
ax.set_xticklabels(ax.get_xticklabels(), rotation=45)
ax.axhline(2.5)
ax.set_ylim(-1, 5)
plt.show()
# %%
# fig.savefig('cnv - '+dt.datetime.now().strftime('%Y%m%d%H%M%s')+'.svg')
# %%
