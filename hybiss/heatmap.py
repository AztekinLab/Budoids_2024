import matplotlib as mpl
mpl.rcParams['pdf.fonttype'] = 42 # for editable text in AI

import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt

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


fields = adata.obs['comp'].str.split('_').str
adata.obs['merge'] = fields[0] + '_' + fields[-1]
adata.obs['merge'] = adata.obs['merge'].astype('category')
adata.obs['merge'] = adata.obs['merge'].cat.reorder_categories([
    'meso9_P', 'mesoSE_P','mesoAER_P',
    'meso9_D', 'mesoSE_D','mesoAER_D',
    ])
adata_P = adata[adata.obs['loc'] == 'P']
adata_D = adata[adata.obs['loc'] == 'D']


# get gene categories
genelist = pd.read_table('hybiss_genes.txt')
genelist['Genes'] = [genelist['Genes'][i].split(',') for i in genelist.index]
genelist = genelist.explode(['Genes'], ignore_index = True)

# main: CT, Fibro, Cartilage
sub_list = genelist[genelist['Category_2'].isin(['Limb_bud_mesoderm','Fibroblast','Chondrocytes'])]
gene_dict = {i:sub_list[sub_list['Category_2'] == i]['Genes'].tolist() for i in sub_list['Category_2'].unique()}


fig, ax = plt.subplots(2, 1, figsize = (10,5))
mpp = sc.pl.matrixplot(adata_P, gene_dict, groupby='merge', log = True, vmax = 0.4, dendrogram = False, var_group_rotation = 0, show=False, save = False, cmap='Blues', ax = ax[0], colorbar_title='Mean expression')
mpd = sc.pl.matrixplot(adata_D, gene_dict, groupby='merge', log = True, vmax = 0.4, dendrogram = False, var_group_rotation = 0, show=False, save = False, cmap='Blues', ax = ax[1], colorbar_title='Mean expression')

plt.tight_layout()
plt.subplots_adjust(hspace=0.01)
fn = 'Fig4J_ht_main_log_max0.4.pdf'
plt.savefig(fn)



# supp: all
for i in genelist['Category_1'].unique():
    print(i)
    sub_list = genelist[genelist['Category_1'] == i]
    gene_dict = {i:sub_list[sub_list['Category_2'] == i]['Genes'].tolist() for i in sub_list['Category_2'].unique()}

    fn = 'FigS13_ht_log_' + i  + '.pdf'
    fig, ax = plt.subplots(2, 1)

    if i == 'Celltype_marker': # modify vmax
        mpp = sc.pl.matrixplot(adata_P, gene_dict, groupby='merge', log = True, vmax = 0.4,  var_group_rotation = 0, show=False, save = False, cmap='Blues', ax = ax[0], colorbar_title='Mean expression')
        mpd = sc.pl.matrixplot(adata_D, gene_dict, groupby='merge', log = True, vmax = 0.4,  var_group_rotation = 0, show=False, save = False, cmap='Blues', ax = ax[1], colorbar_title='Mean expression')

    else:
        mpp = sc.pl.matrixplot(adata_P, gene_dict, groupby='merge', log = True, var_group_rotation = 0, show=False, save = False, cmap='Blues', ax = ax[0], colorbar_title='Mean expression')
        mpd = sc.pl.matrixplot(adata_D, gene_dict, groupby='merge', log = True, var_group_rotation = 0, show=False, save = False, cmap='Blues', ax = ax[1], colorbar_title='Mean expression')

    # for i, label in enumerate(mpp['mainplot_ax'].get_xticklabels()):
    #     if (label.get_text() in meso_p) and (label.get_text() in se_p):
    #         mpp['mainplot_ax'].get_xticklabels()[i].set_color("orange")
    #     elif label.get_text() in meso_p:
    #         mpp['mainplot_ax'].get_xticklabels()[i].set_color("green")
    #     elif label.get_text() in se_p:
    #         mpp['mainplot_ax'].get_xticklabels()[i].set_color("red")

    # for i, label in enumerate(mpd['mainplot_ax'].get_xticklabels()):
    #     if (label.get_text() in meso_d) and (label.get_text() in se_d):
    #         mpd['mainplot_ax'].get_xticklabels()[i].set_color("orange")
    #     elif label.get_text() in meso_d:
    #         mpd['mainplot_ax'].get_xticklabels()[i].set_color("green")
    #     elif label.get_text() in se_d:
    #         mpd['mainplot_ax'].get_xticklabels()[i].set_color("red")

    plt.tight_layout()
    plt.subplots_adjust(hspace=0.01)
    plt.savefig(fn)
