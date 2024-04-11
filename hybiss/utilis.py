import numpy as np
import pandas as pd

import seaborn as sns
from adjustText import adjust_text

from sklearn.impute import SimpleImputer
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from sklearn.manifold import MDS


def grouped_obs_single(adata, groupby, method = 'sum', layer=None, use_raw = False):
    if layer is not None:
        getX = lambda x: x.layers[layer]
    else:
        getX = lambda x: x.X

    new_idx = adata.var_names
    grouped = adata.obs.groupby(groupby)
    out = pd.DataFrame(
        np.zeros((len(new_idx), len(grouped)), dtype=np.float64),
        columns = list(grouped.groups.keys()),
        index = new_idx
    )
    if method == 'sum':
        if use_raw:
            for group, idx in grouped.indices.items():
                X = getX(adata.raw[idx, new_idx])
                out[group] = np.ravel(X.sum(axis=0, dtype=np.float64)).tolist()
        else:
            for group, idx in grouped.indices.items():
                X = getX(adata[idx, new_idx])
                out[group] = np.ravel(X.sum(axis=0, dtype=np.float64)).tolist()

    else:
        if use_raw:
            for group, idx in grouped.indices.items():
                X = getX(adata.raw[idx, new_idx])
                out[group] = np.ravel(X.mean(axis=0, dtype=np.float64)).tolist()
        else:
            for group, idx in grouped.indices.items():
                X = getX(adata[idx, new_idx])
                out[group] = np.ravel(X.mean(axis=0, dtype=np.float64)).tolist()
    return out


def pca_budoids(X, n_components = 2):
    # X: row as sample, col as features
    sim = SimpleImputer(missing_values=np.nan, strategy='mean')
    pca = PCA(n_components=n_components, random_state = 1234)

    X_imp = sim.fit_transform(X)
    scaled_data = StandardScaler().fit_transform(X_imp)
    X_reduced = pca.fit_transform(scaled_data)

    return X_reduced, pca
