"""
Immune deconvolution from MCPcounter for
Transcriptional Analysis with Python Imported from R (TAPIR)
F. Comitani     @2021
J. O. Nash      @2021
"""

import pandas as pd
import os

import rpy2.robjects as ro
from rpy2.robjects.packages import importr
from rpy2.robjects import pandas2ri
from rpy2.robjects.conversion import localconverter

ref_path = os.path.join(os.path.dirname(__file__),'ref/')

""" Import R libraries. """

base = importr('base')
MCPcounter = importr('MCPcounter')

def mcpc_estimate(tab): 

	""" Runs an immune cell deconvolution estimation with MCPcounter. 

        Args:
            tab (pandas dataframe): TMM-normalized expression counts matrix with samples as row and 
				genes as columns.
		Returns:
			imm_counts (pandas dataframe): immune population counts obtained with
				MCPcounter
	"""

	""" Convert input dataframe and reference probes/genes to R. """
	
	with localconverter(ro.default_converter + pandas2ri.converter):
		tabr    = ro.conversion.py2rpy(tab.T)
		probesr = ro.conversion.py2rpy(pd.read_hdf(os.path.join(ref_path,'immune_probes.h5')))
		genesr  = ro.conversion.py2rpy(pd.read_hdf(os.path.join(ref_path,'immune_genes.h5')))

	""" Run the immune deconvolution. """
    
	imm_counts = MCPcounter.MCPcounter_estimate(expression = tabr, featuresType = 'HUGO_symbols', 
                                         probesets= probesr,
                                         genes= genesr)

	""" Convert cell populations counts back to pandas. """

	with localconverter(ro.default_converter + pandas2ri.converter):
		imm_counts = pd.DataFrame(ro.conversion.rpy2py(imm_counts), 
			index=base.rownames(imm_counts), columns=base.colnames(imm_counts))

	return imm_counts.T


if __name__ == "__main__":

	pass