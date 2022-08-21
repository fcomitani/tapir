from setuptools import setup, find_packages
from codecs import open
from os import path

here    = path.abspath(path.dirname(__file__))
version = open("tapir/_version.py").readlines()[-1].split()[-1].strip("\"'")

with open(path.join(here, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()


setup(

    name='tapir-rna',

    version=version,

    description='Transcriptional Analysis in Python Imported from R',
    long_description=long_description,

    url='https://github.com/fcomitani/tapir',
	download_url = 'https://github.com/fcomitani/tapir/archive/'+version+'.tar.gz', 
    author='Federico Comitani, Josh O. Nash',
    author_email='federico.comitani@gmail.com',

    license='MIT',

    classifiers=[
        'Development Status :: 4 - Beta',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python'
        'Programming Language :: Python :: 3.7'
        'Programming Language :: Python :: 3.8'
		],

    keywords='clustering recursion dimension-reduction k-NN hiearchical-clustering optimal-clusters differential-evolution',

    packages=find_packages(exclude=['contrib', 'docs', 'tests']),

    install_requires=['numpy>=1.20.2',
	'pandas>=1.1.3',   
	'gseapy>=0.9.5',
	'lifelines >= 0.21.0',
	'rpy2 >= 3.4.5',
	'seaborn >= 0.11.1',
	'scikit-learn >= 0.4.6',
	'statsmodels >= 0.11.1',
	'umap-learn >= 0.3.9'],

    zip_safe=False,

)

