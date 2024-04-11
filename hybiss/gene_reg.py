import matplotlib as mpl
mpl.rcParams['pdf.fonttype'] = 42 # for editable text in AI

import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import seaborn as sns

import numpy as np
import pandas as pd
import scanpy as sc

import glob
import os, sys


from patsy import dmatrix
import statsmodels.api as sm


import pickle
with open('hybiss_h5ad.pkl', "rb") as fn:
    adata_list_all =  pickle.load(fn)


adata_list = []
for adata in adata_list_all:
    if adata.obs['loc'][0] != 'A':
        del adata.obsm
        adata_list.append(adata)
print(len(adata_list))
# 41


g2c = {'meso9':'#015492', 'mesoAER':'#8E0068','mesoSE':'#FF9300','':'lightgrey'}


genes = ['Msx1', 'Fgf10', 'Wnt5a', 'Axin2', 'Bmp4', 'Id2', 'Id3', 'Gli3']

fig, axes = plt.subplots(2, 4, figsize = (12, 5))
ax = axes.ravel()
for i, g in enumerate(genes):
    df_list = []

    for test in adata_list:
        if g in test.var_names:
            plotdf = sc.get.obs_df(test, keys=[g, 'major_coor_used','cond','sample'], use_raw = True)
            df_list.append(plotdf)

    df = pd.concat(df_list)
    df['cond'] = pd.Categorical(df['cond'], categories=['meso9', 'mesoSE', 'mesoAER'])


    for cond in df['cond'].unique().sort_values():
        df_sub = df[df['cond'] == cond]
        X = df_sub[['major_coor_used']]
        y = df_sub[[g]]

        X = sm.add_constant(X)
        model = sm.ZeroInflatedPoisson(y, X).fit()


        idx = X.duplicated()
        a = X[~idx]
        a['pred'] = model.predict(a)
        summary = model.get_prediction().summary_frame()
        summary.index = X.index
        plot_df = pd.concat([a, summary[~idx]], axis =1)
        plot_df = plot_df.sort_values('major_coor_used')


        ax[i].fill_between(plot_df['major_coor_used'], plot_df.ci_lower, plot_df.ci_upper, alpha=.3, color=g2c.get(cond), edgecolor = 'none')
        ax[i].plot(plot_df['major_coor_used'], plot_df['pred'], color=g2c.get(cond), label=i, lw = 3, solid_capstyle='round')

        ax[i].set_title(g)
        ax[i].set_xlabel('Morphological midline')
        ax[i].set_ylabel('')



legend_elements = [
    mlines.Line2D([0], [0],color='#015492', lw=4, label='Meso-meso'),
    mlines.Line2D([0], [0],color='#FF9300', lw=4, label='Meso-SE'),
    mlines.Line2D([0], [0],color='#8E0068', lw=4, label='Meso-AER'),
    ]

plt.legend(handles = legend_elements, bbox_to_anchor=(1.05, 1), loc='upper left', frameon = False)
plt.tight_layout()
plt.savefig('FigS14B_lineplot.pdf', transparent=True)
