<img src="docs/figs/logo_tp.png" width=400, padding=100>


## Transcriptome Analysis in Python Imported from R
### v 0.1

[![Licence](https://img.shields.io/github/license/fcomitani/tapir?style=flat-square)](https://github.com/fcomitani/tapir/blob/main/LICENSE)
[![GitHub top language](https://img.shields.io/github/languages/top/fcomitani/tapir?style=flat-square)](https://github.com/fcomitani/tapir/search?l=python)
[![Documentation Status](https://readthedocs.org/projects/tapir/badge/?version=latest&style=flat-square)](https://tapir.readthedocs.io/en/latest/?badge=latest)
<!--
[![Build Status](https://img.shields.io/travis/com/fcomitani/tapir/main?style=flat-square)](https://travis-ci.com/fcomitani/tapir)
-->

`tapir` is a python 3 package for the analysis of gene expression data.
It includes a number of functions for statistical analysis, differential expression
and gene sets enrichment analysis.

WARNING: This library is still in active development and we hope to add more options in the future. Please feel free 
to leave feedback, suggestions or to contribute to this repository.

This library includes

* TMM normalization with [EdgeR](https://bioconductor.org/packages/release/bioc/html/edgeR.html)
* differential expression analysis with [EdgeR](https://bioconductor.org/packages/release/bioc/html/edgeR.html)
* gene sets enrichment analysis with [gseapy](https://github.com/zqfang/GSEApy)
* survival analysis with [lifelines](https://github.com/CamDavidsonPilon/lifelines)
* immune deconvolution with [MCPcounter](https://github.com/ebecht/MCPcounter)
* dimensionality reduction with [UMAP](https://github.com/lmcinnes/umap)
* plotting functions for distribution comparisons, heatmaps and gene sets networks.

Detailed documentation, API references and tutorials can be found at this [link](https://tapir.readthedocs.io/en/latest/).

### Dependencies

Besides basic scientific and plotting libraries, the current version requires

```
- gseapy
- lifelines
- rpy2
- seaborn
- scikit-learn
- statsmodels
- umap-learn
```

** R, EdgeR and MCPcounter need to be installed independently. **

### Installation

tapir releases can be easily installed through the python standard package manager  
`pip install tapir-rna`

To install the latest (unreleased) version you can download it from this repository by running 
 
    git clone https://github.com/fcomitani/tapir
    cd tapir
    python setup.py install

### Basic usage

Given an `input` dataset in pandas-like format (samples X genes), the `build_dgelist`  and `diff_exp` functions will allow you to normalize 
the samples as TMM and fit a glmQL model for differential expression
significance.

    from tapir.edger import build_dgelist, diff_exp

    dgelist, tmmlog = build_dgelist(input_table)
    de              = diff_exp(dgelist, groups, filter=True)

### Contact us

- federico.comitani at sickkids.ca
- josh.nash at sickkids.ca

<!--
### Citation

When using this library, please cite

> F. Comitani, J. O. Nash 
-->

### Contributions

This library is still a work in progress and we are striving to improve it, by adding more flexibility and increase the memory and time efficiency of the code. If you would like to be part of this effort, please fork the master branch and work from there. 

<!-- Make sure your code passes the travis build tests. -->

Contributions are always welcome.
