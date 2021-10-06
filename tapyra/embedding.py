"""
Dimensionality reduction functions for
Transcriptional Analysis with Python Imported from R (TAPIR)
F. Comitani     @2021
"""

import os

import numpy as np
import pandas as pd

import umap

def remove_collinear(df, thresh=0.75):

    """ Remove collinearity. WARNING: slow!
    
        Args:
            df (pandas dataframe): expression counts matrix with samples as row and 
                features (genes) as columns. All features within this matrix
                are compared.
            thresh (float): correlation cutoff. Features whose correlation
                is above this value will be candidate for removal (default 0.75).
                
        Returns:
            (pandas dataframe): reduced expression counts matrix.		
    """

    crmat = pd.DataFrame(np.corrcoef(df.astype(float).T), columns=df.columns, index=df.columns)
    crmat.index.name = None
    crmat2 = crmat.where(np.triu(np.ones(crmat.shape), k=1).astype(np.bool)).stack()
    crmat2 = crmat2.reset_index().sort_values(0, ascending=False)
    crmat2 = crmat2[crmat2[crmat2.columns[2]]>thresh]

    toremove = []
    while len(crmat2[crmat2.columns[2]])>0:
        a = crmat2[crmat2.columns[0]].iloc[0]
        b = crmat2[crmat2.columns[1]].iloc[0]
        meana = crmat.loc[a,crmat.columns.difference([a,b])].mean()
        meanb = crmat.loc[b,crmat.columns.difference([a,b])].mean()

        toremove.append([a,b][np.argmax([meana,meanb])])

        crmat2 = crmat2[(crmat2[crmat2.columns[0]] != toremove[-1]) & 
                        (crmat2[crmat2.columns[1]] != toremove[-1])]
        
    return df.drop(toremove, axis=1)

def near_zero_var_drop(df, thresh=0.99):

    """ Remove features with near-zero variance.

        Args:
            df (pandas dataframe): expression counts matrix with samples as row and 
                features (genes) as columns. All features within this matrix
                are compared.
            thresh (float): cumulative variance cutoff. Features whose variance sum up
                to this percentage of the total variance will be kept (default 0.99).
            
        Returns:
            (pandas dataframe): reduced expression counts matrix.
    """

    vVal   = df.var(axis=0).values
    cs     = pd.Series(vVal).sort_values(ascending=False).cumsum()
    remove = cs[cs>cs.values[-1]*thresh].index.values

    return df.drop(df.columns[remove],axis=1)

def get_umap(vals, collinear_thresh=None, var_drop_thresh=None, n_neighbors='sqrt', **kwargs):

    """ Wrapper function for UMAP dimensionality reduction.

        Args:
            vals (np.array or pandas dataframe): expression counts matrix with samples as row and 
                features (genes) as columns. 
            var_drop_thresh (float): cumulative variance cutoff for near_zero_var_drop. If None, 
                skip this step (default None).
            var_drop_thresh (float): correlation cutoff for remove_collinear. If None, 
                skip this step (default None).
            n_neighbors (int): number of nearest neighbours for UMAP dimensionality reduction.
                If 'sqrt' take the square root of the total number of samples (default 'sqrt').
            kwargs (dict): dictionary containing further arguments for UMAP.
            
        Returns:
            proj (pandas dataframe): coordinates of the projected points.
            trained_map (UMAP): the trained UMAP object.
    """

    vals_cut=vals.copy(deep=True)

    if collinear_thresh is not None:
        vals_cut = remove_collinear(vals_cut, thresh=collinear_thresh)
    if var_drop_thresh is not None:
        vals_cut = near_zero_var_drop(vals_cut, thresh=var_drop_thresh)

    if n_neighbors=='sqrt':
        n_neighbors = int(np.sqrt(vals_cut.shape[0]))

    mapping       = umap.UMAP(n_neighbors=n_neighbors, **kwargs)
    trained_map   = mapping.fit(vals_cut)
    proj          = pd.DataFrame(mapping.transform(vals_cut))

    return proj - proj.mean(), trained_map

if __name__ == "__main__":
    
        pass