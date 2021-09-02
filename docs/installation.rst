============
Installation
============

..
   tapir can be easily installed through python standard 
   package managers without worrying about dependencies. 

   With pip

   .. code-block:: bash

      pip install tapir

To install the latest (unreleased) version 
you can download it from our GitHub repository by running 

.. code-block:: bash

   git clone https://github.com/fcomitani/tapir
   cd tapir
   python setup.py install


The current version requires the following 
packages and their inherited dependencies:

   - gseapy >= 0.9.5
   - lifelines >= 0.21.0
   - rpy2 >= 3.4.5
   - seaborn >= 0.11.1
   - scikit-learn >= 0.4.6
   - statsmodels >= 0.11.1
   - umap-learn >= 0.3.9

The code has been tested with these version, it may 
work with newer releases although not guaranteed.
For a full list of python dependencies see :code:`requirements.txt`

**It also requires R and EdgeR to be installed independently.**
