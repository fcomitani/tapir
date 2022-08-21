"""
Gene sets enrichment analysis functions for
Transcriptional Analysis with Python Imported from R (TAPIR)
F. Comitani     @2021
"""

import os

import numpy as np
import pandas as pd
import gseapy as gp

from tapir.auxiliary import invert_dict

ref_path = os.path.join(os.path.dirname(__file__),'ref/')

def gset_as_dict(subsel=None, ref=None):

    """ Read in gene sets reference file and transforms it in
        a dictionary.

        Args:
            subsel (list of strings): gene sets to subselect. All sets
                containing any of the strings in this list will be kept
                (default None).
            ref (str): path to the reference files containing the gene sets to be 
                included in the analysis, if None use the provided file (default None).
           

        Returns:
            sets_dict (dictionary): a dictionary with gene sets names
                as keys and the list of member genes as value.
    """

    if ref is None:
        ref = os.path.join(ref_path,'all_gsets.gmt')
        
    with open(ref, 'r') as gmt:
        gset_list = [d.split('\t') for d in gmt.read().split('\n')]
    
    if not isinstance(subsel, list) and subsel!=None:
        subsel = [subsel]
        
    sets_dict = {}
    for gl in gset_list:
        if subsel is None or any(x in gl[0] for x in subsel):
            sets_dict[gl[0]] = gl[2:]
    
    return sets_dict

def connection_matrix_gsets(gset, sets_dict):

    """ Builds a matrix counting the number of times genes are
        found in together in gene sets.

        Args:
            gset (list of string): the list of genes to search
                and compare.
            sets_dict (dictionary): a dictionary of the gene sets
                and to be considered for the counting. The format should
                follow the output of gset_as_dict.
                
        Returns:
            (pandas dataframe): a square dataframe with the number
                of common gene sets between each provided gene pair.
    """

    inv_dict=invert_dict(sets_dict)

    mat=[]
    for i in np.arange(len(gset)):
        mat.append([])
        mat[-1]+=[np.nan]*i
        for j in np.arange(i,len(gset)):
            mat[-1].append(len(set(inv_dict[gset[i]]).intersection(inv_dict[gset[j]])))
            
    A = np.nan_to_num(np.array(mat),0)
    A = np.tril(A.T,1)+A
    np.fill_diagonal(A,np.diagonal(A)/2)

    return pd.DataFrame(A, columns=gset, index=gset)

def run_gsea(df, subsel=None, type='gsea', ref=None, tmp_path=r'./tmp_gsea', **kwargs):
    
    """ Run gene sets enrichment analysis.

        Args:
            df (pandas dataframe): count matrix or preranked list of genes to be used
                for enrichment analysis
            subsel (list of strings): list of pathways to subselect. All pathways containing
                the provided strings will be selected (e.g. 'HALLMARK_' will select all
                hallmark of cancer pathways, default None) .
            type (str): chose the type of analysis to run, single sample 'ssgsea', preranked
                'prerank' or standard 'gsea' (default 'gsea').
            ref (str): path to the reference files containing the gene sets to be 
                included in the analysis, if None use the provided file (default None).
            tmp_path (str): path where temporary files will be stored. If it doesn't exist
                the function will try and build it (default ./tmp_gsea).
            **kwargs: keyword parameters for the enrichment analysis functions.
            
        Returns:
            (pandas dataframe): dataframe with enrichment scores and p-values, 
                where available.
            
    """

    """ Set up temporary files directory. """

    if not os.path.exists(tmp_path):
        os.makedirs(tmp_path)

    if ref is None:
        ref = os.path.join(ref_path,'all_gsets.gmt')

    """ Set up gene sets files, select relevant sets. """

    with open(ref, 'r') as gmt:
        data = gmt.read().split('\n')

    for d in range(len(data)):
        data[d]=data[d].split('\t')
            
    if not isinstance(subsel, list) and subsel!=None:
            subsel=[subsel]
    with open(os.path.join(tmp_path,'subsel_gsets.gmt'), 'w') as gmt:
        for d in data:
            if subsel==None or any(x in d[0] for x in subsel):
                if isinstance(df, pd.DataFrame):
                    setdf=set(df.columns)
                elif isinstance(df, pd.Series):
                    setdf=set(df.index.values)
                else:
                    print('Need a pandas Series or DataFrame!')
                    return -1

                if len(setdf.intersection(d))>0:
                    for el in d:
                        gmt.write(el+'\t')
                    gmt.seek(0, os.SEEK_END)
                    gmt.seek(gmt.tell() - 1, os.SEEK_SET)
                    gmt.truncate()
                    gmt.write('\n')
    gmt.close()

    df_new             = df
    df_new.index.names = ['NAME']

    params={'gene_sets':           tmp_path+r'subsel_gsets.gmt',
            'outdir':              tmp_path,
            'min_size':            5,
            'max_size':            1000, 
            'permutation_num':     0,
            'weighted_score_type': 1,
            'processes':           8,
            'verbose':             True,
            'no_plot':             True}

    """ Override with keyword arguments if provided. """
    for key, value in kwargs.items():
        params[key] = value

    """ Run gene sets analysis. """

    if type=='ssgsea':
        params['sample_norm_method']  = 'rank'
        params['weighted_score_type'] = 0
        gs_res = gp.ssgsea(df_new.T, **params)
    elif type=='prerank':
        params['weighted_score_type'] = 0
        gs_res = gp.preranked(df_new.T, **params)
    else:
        gs_res = gp.gsea(df_new.T, **params)

    return gs_res.res2d.sort_index().T

if __name__ == "__main__":
    
        pass