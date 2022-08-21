
====================
Immune Deconvolution
====================

Immune deconvolution is available through MCP counter, starting from
an expression matrix :code:`tmm` with samples as rows and genes as columns.

.. code-block:: python

  from tapir.immune import mcpc_estimate

  estimates = mcpc_estimate(tmm)


For details on the deconvolution process, please see the MCPcounter
publication [Becht2016]_.

References
----------
        
.. [Becht2016] Becht E., Giraldo N. A., Lacroix L., Buttard B., Elarouci N., Petitprez F., Selves J., Laurent-Puig P., Sautès-Fridman C., Fridman W. H. and de Reyniès A., (2016). “Estimating the population abundance of tissue-infiltrating immune and stromal cell populations using gene expression”, Genome Biol. 20; 17(1) 218.
