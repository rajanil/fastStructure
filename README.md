# fastStructure

*fastStructure* is a fast algorithm for inferring population structure from large SNP genotype data. 
It is based on a variational Bayesian framework for posterior inference and is written in Python2.x. 
This repo contains a *vars* library with definitions and methods for
the relevant variables, and a set of scripts to load the data and run the algorithm.

This document summarizes how to setup this software package, compile the C and 
[Cython](http://cython.org/) scripts and run the algorithm on a test simulated genotype dataset.

## Parts 

This repo has two components: a library of C and Cython scripts *vars* and
a set of Cython and pure Python scripts to load the data and run the algorithm.

## Dependencies

*fastStructure* depends on 
+ [Numpy](http://www.numpy.org/)
+ [Scipy](http://www.scipy.org/)
+ [Cython](http://cython.org/)

## Getting the source code

To obtain the source code from github, let us assume you want to clone this repo into a
directory named proj:

    mkdir ~/proj
    cd ~/proj
    git clone https://github.com/rajanil/fastStructure

To retrieve the latest code updates, you can do the following:

    cd ~/proj/fastStructure
    git fetch
    git merge origin/master

## Building Python extensions

To build library extensions, you can do the following:

    cd ~/proj/fastStructure/vars
    python setup.py build_ext --inplace

To compile the main cython scripts, you can do the following:

    cd ~/proj/fastStructure
    python setup.py build_ext --inplace

Each setup will create some *.c and *.so (shared object) files.
This setup may give some warnings, which are OK. If you get errors that indicate the 
build failed, this might be because the LD_LIBRARY_PATH, CFLAGS, LDFLAGS environment 
variables are set incorrectly.

## Executing the code

The main script you will need to execute is `structure.py`. To see command-line 
options that need to be passed to the script, you can do the following:

    python structure.py

The current implementation can import data from [plink bed](http://pngu.mgh.harvard.edu/~purcell/plink/data.shtml#bed) 
format and the original Structure format. If the data are in plink format, ensure that
bed, bim and fam files for the dataset are all present in the same path.

## Running on test data
