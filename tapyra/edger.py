"""
Differential expression functions from EdgeR for
Transcriptional Analysis with Python Imported from R (TAPIR)
F. Comitani     @2021
J. O. Nash      @2021
"""

import pandas as pd

import rpy2.robjects as ro
from rpy2.robjects.packages import importr
from rpy2.robjects import pandas2ri
from rpy2.robjects.conversion import localconverter

""" Import R libraries. """

base = importr('base')
edgeR = importr('edgeR')
stats = importr('stats')
dplyr = importr('dplyr')

def build_dgelist(tab): 

	""" Builds a DGE list and normalizes counts as TMM. 

        Args:
            tab (pandas dataframe): expression counts matrix with samples as row and 
				genes as columns.
		Returns:
			dgelist (edgeR DGElist): a composite dataframe with samples and expression
				counts information to be used by diff_exp.
			tmmlog (pandas dataframe): log2(TMM)-normalized expression counts dataset,
				with samples as rows and genes as columns.
	"""

	""" Convert input dataframe to R. """
	
	with localconverter(ro.default_converter + pandas2ri.converter):
		tabr = ro.conversion.py2rpy(tab.T)

	""" Build DGEList and normalize counts. """
	
	dgelist = edgeR.DGEList(counts=tabr, genes=base.rownames(tabr), samples = base.colnames(tabr))
	dgelist = edgeR.calcNormFactors(dgelist, method="TMM")
	tmmlog = edgeR.cpm(dgelist, log=True, **{'normalized.lib.sizes': True})

	""" Convert TMM counts back to pandas. """

	with localconverter(ro.default_converter + pandas2ri.converter):
		tmmlog = pd.DataFrame(ro.conversion.rpy2py(tmmlog), index=dgelist.rx2('genes')['genes'], columns=dgelist.rx2('samples')['samples'])

	return dgelist, tmmlog.T

def diff_exp(dgelist, groups, batches=None, filter=False, min_count=2, min_total_count=3):

	""" Runs a differential expression analysis with edgeR on a given DGElist.

        Args:
            dgelist (edgeR DGElist): a composite dataframe with samples and expression
				counts information obtained with build_dgelist.
            groups (pandas dataframe): a single-column dataframe containing information
				on the groups to compare. Samples are as rows, "group" is the only column
				and contains the groups membership for each sample encoded as integers.
            batches (pandas dataframe): a single-column dataframe containing information
				on the batches. Samples are as rows, "batch" is the only column
				and contains the batch membership for each sample encoded as integers
				(optional, defaul None).
            filter (bool): if True, filter genes with counts lower than a set threshold,
				(defualt False)
			min_count (int): minimum number of per-group counts when filtering genes.
			min_count (int): minimum number of total counts when filtering genes.

		Returns:
			de_sorted (pandas dataframe): a dataframe summirizing the results of 
			the differential expression analysis, including log-normalized fold-change, average
			log-normalized count per millions over all the libraries, F-statistics, p-values and
			flase discovery ratio.
	"""

	
	""" Setup input dataframe and constants. """

	groups_tmp = groups.reset_index()
	groups_tmp.columns = ['samples','group']

	numcols=2
	formula='~groups'

	if batches is not None:
		
		""" If provided, add batch information. """

		batches_tmp = batches.reset_index(drop=True)
		groups_tmp = pd.concat([groups_tmp,batches_tmp],axis=1)
		groups_tmp.columns = ['samples','group','batch']
		
		numcols+=1
		formula+='+batches'

	""" Convert input groups dataframe to R. """

	with localconverter(ro.default_converter + pandas2ri.converter):
		groupsr = ro.conversion.py2rpy(groups_tmp)


	""" Add groups and subset the DGElist if necessary. """
	
	dgelist[1] = base.merge(x = groupsr, y = dgelist.rx2('samples'), by = "samples", **{'all.x' : True})
	dgelist[0] = base.subset(dgelist[0], select=dgelist[1][0])

	""" Clean up columns name and select only relevant columns. """

	base.colnames(dgelist[1])[1] = 'group'
	dgelist[1]=base.subset(dgelist[1],select=base.colnames(dgelist[1])[:numcols]+
										base.colnames(dgelist[1])[numcols+1:])

	""" Filter low expression genes if requested. """

	if filter:
		keep = edgeR.filterByExpr(dgelist, **{'min.count' : min_count, 'min.total.count' : min_total_count})
		dgelist = dgelist.rx(keep, True)

	""" Build the design matrix. """

	fmla = ro.Formula(formula)
	env = fmla.environment
	env['groups'] = dgelist.rx2('samples').rx2('group')

	if batches is not None:
		env['batches'] = dgelist.rx2('samples').rx2('batch')

	design=stats.model_matrix(fmla)
	set_method = ro.r("`colnames<-`")
	design = set_method(design, base.make_names(base.levels(dgelist.rx2('samples').rx2('group'))))
	
	""" Fit linear model, calculate DE statistics. """

	dge_disp = edgeR.estimateGLMRobustDisp(dgelist, design = design)
	fit = edgeR.glmQLFit(dge_disp, design, robust = True)
	de_result = edgeR.glmQLFTest(fit, coef = len(set(env['groups'])))

	""" Sort genes by p-value and reformat. """
	
	de_sorted = edgeR.topTags(de_result, n=base.nrow(de_result), **{'sort.by' : "p.value"})
	de_sorted = pd.DataFrame(ro.conversion.rpy2py(de_sorted[0]), index=base.colnames(de_sorted[0])).T
	de_sorted.set_index('genes', inplace=True, drop=True)

	return de_sorted


if __name__ == "__main__":

	pass