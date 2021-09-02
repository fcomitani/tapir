
================================
Differential Expression Analysis
================================

There are two main functions for differential expression converted from EdgeR.
:code:`build_dgelist` takes as input a pandas dataframe with expression counts,
with samples as rows and genes as columns. It returns the log2-normalized TMM matrix 
and, as the name implies, a DGElist. The latter is to be used as an input for 
:code:`diff_exp` which will fit a glmQL [Robinson2010]_ model and return the results of the
differential expression analysis. This function also needs a pandas dataframe
containing information on the samples membership to the groups to be compared,
with samples as rows, and a single column, :code:`group`, with the group
number.
If activated, the :code:`filter` option allows to remove genes that fall below 
an expression threshold set by :code:`min_count` and :code:`min_total_count`, 
equivalent to the flags in :code:`EdgeR.filterByExpr`.

.. code-block:: python

  from tapyra.edger import build_dgelist, diff_exp

  dgelist, tmmlog = build_dgelist(input_table)
  de              = diff_exp(dgelist, groups, filter=True)
  
The output includes log-normalized fold change, 
average log-normalized CPM across all samples, the
quasi-likelihood F-statistics, p-values and FDR.

.. 
  ===== ========= ========== ========= ========= =========
  gene     logFC     logCPM         F    PValue       FDR 
  ===== ========= ========== ========= ========= =========
  MYCN  -2.492382  17.356105  3.951764  0.096183  0.384733
  MYC   -0.195281  18.062056  0.139693  0.741225  0.937882
  ALK    0.069065  18.439629  0.757448  0.889513  0.937882
  ERBB2 -0.035652  18.616506  0.203695  0.937882  0.937882
  ===== ========= ========== ========= ========= =========

References
----------
        
.. [Robinson2010] Robinson M.D., McCarthy D.J., Smyth G.K. (2010). “edgeR: a Bioconductor package for differential expression analysis of digital gene expression data”, Bioinformatics, 26(1), 139-140.