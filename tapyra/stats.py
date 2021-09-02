"""
Statistical analysis functions for
Transcriptional Analysis with Python Imported from R (TAPIR)
F. Comitani     @2021
"""

import os
import warnings
import itertools

import numpy as np
import pandas as pd

from scipy.stats import mannwhitneyu as mwu
from scipy.stats import kruskal as kw
from scipy.stats import chisquare as chi2
from scipy.stats import fisher_exact as fet

from statsmodels.stats.multitest import multipletests

from lifelines import KaplanMeierFitter
from lifelines.statistics import multivariate_logrank_test as mlt

from tapyra.auxiliary import smart_selection

def compare(gene, classes):

    """ Compare the distribution of two or more groups and automatically selects
        the proper statistical test

        Args:
            gene (string): feature to be compared.
            classes (list of pandas dataframe): list of groups (classes) to compare.
            
        Returns:
            (list): set of test outputs, including median values,
                median values ratio (maximum if >2 groups are compared),
                median values difference (mainimum if >2 groups are compared),
                statistical test score and p-value.
            dunn (pandas dataframe): dunn post-hoc test values, None if only two groups
                are compared.
    """
    
    csem = [c[gene] for c in classes]
    meds = [c.median() for c in csem]

    if len(classes) == 2:
            
        """ Only two groups, Mann-Whitney U test. """

        try:
            tmp = mwu(csem[0], csem[1])
        except:
            tmp = [np.nan, 1]

        if meds[0] == 0:
            ratio = 0
        elif meds[1] == 0:
            ratio = np.inf
        else:
            ratio = np.abs(meds[0]/meds[1])

        diff = meds[0] - meds[1]
        dunn = None
    
    else:
    
        """ More than two groups, Kruskal-Wallis test. """

        try:
            tmp = kw(*csem)
        except:
            tmp = [np.nan, 1]

        ratio = []
        diff  = []
        for i in range(len(meds)-1):
            for j in range(i+1,len(meds)):
                if meds[j] != 0:
                    ratio.append(np.min([np.abs(meds[i]/meds[j])]))
                    diff.append(np.abs(meds[i]-meds[j]))

        if len(minratio) == 0:
            minratio = np.nan
            diff     = np.nan
        else:
            ratio = np.min(ratio)
            diff  = np.max(diff)

        """ Dunn post-hoc significance test. """    
        
        dunn = sp.posthoc_dunn(csem, p_adjust='bonferroni')        
        dunn = pd.DataFrame(dunn) #columns=classes, index=classes
        
    return [*meds, ratio, diff,
            tmp[0], tmp[1]], dunn
            
def multicompare(classes, labs, df, cutoff=0.05, multi_method='fdr_tsbh', multi_alpha=0.05):
    
    """ Run the comparison test between groups for multiple features (genes).

        Args:
            classes (list int or strings): list of group (class) names to compare.
            labs (panda dataframe): one-hot-encoded classes membership dataframe 
                with samples as rows and classes as columns.
            df (pandas dataframe): expression counts matrix with samples as row and 
                features (genes) as columns. All features within this matrix
                are compared.
            cutoff (float): significance cutoff. Features with a p-value/fdr above this
                value are automatically discarded (default 0.05).
            multi_method (string): method for multipletesting correction with statsmodels.
                (default 'fdr_tsbh').
            multi_alpha (float): FWER, family-wise error rate for multipletesting correction 
                with statsmodels (default 0.05).
            
        Returns:
            (pandas dataframe): dataframe summarizing the statistical test outputs
                for each feature which reached significance.
            (list of pandas dataframe): list of dunn post-hoc test outputs, None if 
                only two groups are compared.
    """
    
    """ Subset expression counts matrix by comparable groups. """

    if len(classes)==2:

        """ When two classes are compared, make sure they don't intersect. 
            The first class gets priority. """

        sel = [df.loc[labs[smart_selection(labs,classes[0],how='any',val=1)].index],
             df.loc[labs[(smart_selection(labs,classes[1],how='any',val=1))&\
                         (smart_selection(labs,classes[0],how='all',val=0))].index]]
    else:
        sel=[]
        for c in classes:
            sel.append(df.loc[labs[smart_selection(labs,c,how='any',val=1)].index])

    """ Run comparison for each gene, build dataframe. """   

    allgs = [compare(gene,sel) for gene in df.columns]
    
    if len(classes)>2:
        dunns = [x[1] for x in allgs]
    else:
        dunns = None

    allgs = [x[0] for x in allgs]
    allgs = pd.DataFrame(allgs, index=df.columns, 
            columns=['med_'+str(x) for x in np.arange(len(sel))]+\
            ['ratio','diff','value','pval'])

    """ Apply multiple-tests correction. """
    allgs['padj'] = multipletests(allgs['pval'], 
                    alpha=multi_alpha, method=multi_method)[1]
    
    """ Apply p-value cutoff. """

    allgs = allgs[allgs['pval']<cutoff]

    return allgs.sort_values('pval', ascending=True), dunns

def get_contingency(series, classes, labs):
    
    """ Build contingency table from a series of variables and samples
        groups.
    
        Args:
            series (list or panda series): list of variables (e.g. sex, status, therapy). 
            classes (list int or strings): list of group (class) names to compare.
            labs (panda dataframe): one-hot-encoded classes membership dataframe 
                with samples as rows and classes as columns.
                
        Returns:
            contab (pandas dataframe): contingency table with groups as rows
                and input variable values as columns.
    """

    dummy = pd.get_dummies(series)
    contab = [dummy.loc[[x for x in labs[labs[x]==1].index if x in dummy.index]].sum()
            for x in classes]

    return pd.DataFrame(contab)

def test_contingency(contab, method='auto'):
    
    """ Build contingency table from a series of variables and samples
        groups.
    
        Args:
            contab (pandas dataframe): contingency table with groups as rows
                and input variable values as columns.
            method (string): statistical test to run, to be chosen among 'chi2' 
                for chi square, 'fet' for Fisher Exact, and 'auto' to automatically
                select the test based on the contingency table shape (default 'auto'). 
                
        Returns:
            (pandas dataframe): dataframe containing the statistic and pvalue
                outputs obtained with the selected test.
        
    """
    
    """ Automatic method selection. """

    if method=='auto':
        method = 'chi2'
        if contab.shape==(2,2):
            method = 'fet'

    """ Chi2. """
    
    if method=='chi2':
        if (contab.sum()<5).any():
            warnings.warn('Insufficient population for Chi2 test, the total population '+\
                  'for each cell should be >5. The test will be interrupted.')
            return None
        else:
            return pd.DataFrame(chi2(contab), columns=contab.columns, index=['statistic','pval']).T
    
    """ FET. """

    return pd.DataFrame(fet(contab), columns=['FET'], index=['statistic','pval']).T

def st_curves(st_stats, groups, labs, st_name='survival_time', obs_name='event_observed'):
    
    """ Calculate survival curves and runs multivariate logrank test.
        
        Args:
            st_stats (pandas dataframe): dataframe containing survival times
                information. It must contain two columns, survival_time with the age
                of the sample and event_observed, a binary-encoded report of the patient 
                deceased (1) or alive (0) status.
                Alternatively the column names must be provided.
            groups (list int or strings): list of group names to compare.
            labs (panda dataframe): one-hot-encoded classes membership dataframe 
                with samples as rows and classes as columns.
            st_name (string): name of the column containing age information
                (default 'survival_time').
            obs_name (string): name of the column containing patient status information
                (default 'event_observed').
                
        Returns:
            ml (pandas dataframe): logrank test statistic and p-value.
            curves (pandas dataframe): survival curves for each group.
            
    """
    
    kmf = KaplanMeierFitter()

    st, obs, gps, curves = [], [], [], []

    """ Subset by group and fit Kaplan-Meier Curve. """

    for group in groups:
        
        st_group = st_stats[labs[group]==1]
        st.append(st_group[st_name])
        obs.append(st_group[obs_name])
        gps.append([group]*st_group.shape[0])

        kmf.fit(st[-1], event_observed=obs[-1], 
                label=str(group)+' n='+str(st[-1].shape[0]))
        curves.append(pd.DataFrame(kmf.survival_function_at_times(kmf.timeline)))
    
    curves=pd.concat(curves,axis=1).fillna(method='bfill')
    
    """ Return logrank test and curves. """

    evts = pd.DataFrame({
        'durations': list(itertools.chain.from_iterable(st)),
        'events': list(itertools.chain.from_iterable(obs)),
        'groups': list(itertools.chain.from_iterable(gps)),
        })
    
    return mlt(evts['durations'], evts['groups'], evts['events']).summary, \
           curves

if __name__ == "__main__":
    
        pass