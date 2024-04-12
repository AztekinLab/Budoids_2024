import matplotlib as mpl
mpl.rcParams['pdf.fonttype'] = 42

import matplotlib.pyplot as plt
import seaborn as sns

import numpy as np
import pandas as pd
import scanpy as sc

import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

import utilis


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

adata = adata_list[0].concatenate(adata_list[1:], join = 'outer', fill_value = 0)



# boxplot
genes = ['Msx1', 'Fgf10', 'Wnt5a', 'Axin2', 'Bmp4', 'Id2', 'Id3', 'Gli3']
gene_mat_PD = utilis.grouped_obs_single(adata, groupby = 'comp', use_raw = True, method = 'mean') # raw counts
gene_mat_PD = gene_mat_PD.loc[genes,:]

# average by sample, divide PD
gene_mat_long_PD = pd.melt(gene_mat_PD.reset_index(), id_vars = 'index')
feilds =  gene_mat_long_PD['variable'].str.split('_').str

gene_mat_long_PD['cond'] = feilds[0]
gene_mat_long_PD['cond'] = gene_mat_long_PD['cond'].astype("category")
gene_mat_long_PD['cond'] = gene_mat_long_PD['cond'].cat.reorder_categories([ 'meso9', 'mesoSE','mesoAER'])
gene_mat_long_PD['biorep'] = [''.join(x.split('_')[1:4]) for x in gene_mat_long_PD['variable']]

gene_mat_long_PD['comp'] =  feilds[0] + '_' + feilds[-1]
gene_mat_long_PD['comp'] = gene_mat_long_PD['comp'].astype("category")
gene_mat_long_PD['comp'] = gene_mat_long_PD['comp'].cat.reorder_categories([
   'meso9_P','meso9_D',
   'mesoSE_P','mesoSE_D',
   'mesoAER_P','mesoAER_D',
   ])

gene_mat_long_PD['loc'] = feilds[-1]
gene_mat_long_PD['loc'] = gene_mat_long_PD['loc'].astype("category")
gene_mat_long_PD['loc'] = gene_mat_long_PD['loc'].cat.reorder_categories(['P','D'
])



# plot averaging samples
fig, axes = plt.subplots(2, 4, figsize = (12, 5))
ax = axes.ravel()

for i, g in enumerate(genes):
    dat =  gene_mat_long_PD[gene_mat_long_PD['index'] == g]

    sns.boxplot(data = dat, x = 'loc', hue = 'cond', y = 'value',
    ax = ax[i],
    showfliers = False, showmeans=True,
    hue_order = dat['cond'].cat.categories,
    palette = {'meso9':'#015492', 'mesoAER':'#8E0068','mesoSE':'#FF9300'},
    boxprops={'alpha': 0.9},
    meanprops={'marker':'o','markerfacecolor':'white',  'markeredgecolor':'black','markersize':'8'}
    )

    sns.stripplot(data = dat, x = 'loc', hue = 'cond', y = 'value', dodge=True, ax=ax[i], color = 'black', edgecolor = 'none', alpha=.6, legend = False)
    if i < 7:
        ax[i].get_legend().remove()
    ax[i].set_title(g)
    ax[i].set_xlabel('')
    ax[i].set_ylabel('')
    ax[i].set_xticklabels(['Proximal','Distal'])


plt.legend(loc='upper left', bbox_to_anchor=(1, 1))
plt.tight_layout()
plt.savefig('FigS14B_boxplot.pdf')
