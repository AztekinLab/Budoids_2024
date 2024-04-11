import matplotlib as mpl
mpl.rcParams['pdf.fonttype'] = 42 # for editable text in AI

import scanpy as sc
import pandas as pd
import numpy as np
import glob, os

import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import matplotlib.lines as mlines
import seaborn as sns


import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

import utilis


# Paramters
d1 = 'Dim' + str(1)
d2 = 'Dim' + str(2)

color_dict = {'meso9':'#015492', 'mesoAER':'#8E0068','mesoSE':'#FF9300','':'grey'}
marker_dict = {'A':'o', 'P':'X', 'D':'s'}
line_dict = {'P': '--', 'D':'-'}



import pickle
with open('hybiss_h5ad.pkl', "rb") as fn:
    adata_list_all =  pickle.load(fn)

out_all = []
for adata in adata_list_all:
    sc.pp.filter_genes(adata, min_cells = 0.01 * adata.n_obs)
    # if adata.obs['loc'][0] != 'A':
    out = utilis.grouped_obs_single(adata, groupby = 'comp', method = 'mean')
    out = 100 * out / out.sum(axis = 0)
    out_all.append(out)
print(len(out_all))



## PCA on sample
mat = pd.concat(out_all, axis = 1)
X_reduced, pca = utilis.pca_budoids(mat.transpose())

label = [x.split('_')[4] for x in mat.columns]
col = [x.split('_')[0] for x in mat.columns]
shape = [x.split('_')[-1] for x in mat.columns]


X_reduced = pd.DataFrame(X_reduced)
X_reduced.columns = ['Dim' + str(i+1) for i in range(X_reduced.shape[1])]
# X_reduced['label'] = label
X_reduced['col'] = col
X_reduced['shape'] = shape


# plot PCA separately for each condition
idx1 = np.where([i.startswith('meso9_') for i in mat.columns])[0]
idx3 = np.where([i.startswith('mesoAER_') for i in mat.columns])[0]
idx4 = np.where([i.startswith('mesoSE_') for i in mat.columns])[0]


title = ['meso9','mesoSE', 'mesoAER']
fig, axes = plt.subplots(2, 3,figsize = (15,5), height_ratios=[0.5,3])
ax = axes.ravel()
main_ax = [3,4,5]

for i in [0,1,2]:
    ax[i].set_axis_off()

Centroids = {}
dist_all = []

for j, idx in enumerate([idx1, idx4, idx3]):
    labeli = [l if i in idx else '' for i, l in enumerate(label)]
    coli = [l if i in idx else '' for i, l in enumerate(col)]

    df = X_reduced
    df['label'] = labeli
    df['color'] = coli

    df.loc[df['label'] == '', 'Top'] = 1
    df.loc[df['label'] != '', 'Top'] = 2
    df = df.sort_values(by=['Top'], ascending=True)


    sns.scatterplot(x=d1, y=d2, data=df, hue="color", style = 'shape', ax = ax[main_ax[j]], palette=color_dict, markers = marker_dict)


    df_sub = df[df['color'] != '']
    col_ = color_dict.get(np.unique(coli)[1])
    for m in ['P','D']:
        ax[main_ax[j]-3].set_xlim(ax[main_ax[j]].get_xlim())
        ax[main_ax[j]-3].set_xticks([])
        ax[main_ax[j]-3].set_xlabel('')
        sns.kdeplot(df_sub.loc[df_sub['shape']==m, 'Dim1'], ax=ax[main_ax[j]-3], vertical=False, cut = 0, label = m, linestyle=line_dict.get(m), color = col_, linewidth = 4, solid_capstyle='round')
        # ax[main_ax[j]-3].set_solid_capstyle('round')


    for m in df_sub['shape'].unique():

        x = df_sub.loc[df_sub['shape'] == m, 'Dim1']
        y = df_sub.loc[df_sub['shape'] == m, 'Dim2']

        Centroids[np.unique(coli)[1] + '_' + m] = np.array((x.mean(), y.mean()))

        # Centroids
        ax[main_ax[j]].scatter(x.mean(), y.mean(), s = 100, marker = marker_dict.get(m) ,c = col_, edgecolor = 'black' )

    ax[main_ax[j]-3].set_title(title[j])
    ax[main_ax[j]].legend([],[], frameon=False)



dist = [np.linalg.norm(Centroids.get('meso9_D') - Centroids.get('meso9_P')),
np.linalg.norm(Centroids.get('mesoSE_D') - Centroids.get('mesoSE_P')),
np.linalg.norm(Centroids.get('mesoAER_D') - Centroids.get('mesoAER_P')),
]

for j, _ in enumerate([idx1,  idx4, idx3]):
    ax[main_ax[j]].text(-5, 12, 'P to D distance: ' + str(round(dist[j], 2)))



legend_elements = [
    Patch(facecolor='#015492', edgecolor=None, label='Meso only 9K'),
    Patch(facecolor='#8E0068', edgecolor=None, label='Meso+AER'),
    Patch(facecolor='#FF9300', edgecolor=None, label='Meso+SE'),
    mlines.Line2D([], [], marker='o',color='black',  linestyle='None', markersize=10, label='A'),
    mlines.Line2D([], [], marker='X', color='black',linestyle='None', markersize=10, label='P'),
    mlines.Line2D([], [], marker='s', color='black',linestyle='None', markersize=10, label='D'),

    mlines.Line2D([], [], color='black',linestyle='--', label='P'),
    mlines.Line2D([], [], color='black',linestyle='-', label='D')
    ]



plt.legend(handles = legend_elements, loc='upper left', bbox_to_anchor=(1, 1))
plt.tight_layout()
plt.subplots_adjust(hspace=0.01)
fn = 'Fig4K_pca.pdf'
plt.savefig(fn, transparent = True)
plt.clf()
