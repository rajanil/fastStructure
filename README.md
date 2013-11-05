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

A number of python distributions already have these modules packaged in them. It is also
straightforward to install them from platform-specific binary packages OR from source.

## Getting the source code

To obtain the source code from github, let us assume you want to clone this repo into a
directory named `proj`:

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

    $ python structure.py

    Here is how you can use this script

    Usage: python structure.py
         -K <int>
         --input=<file>
         --output=<file>
         --tol=<float>   (default: 10e-6)
         --prior={simple,logistic}   (default: simple)
         --cv=<int>   (default: 0)
         --full   (to output all variational parameters)
         --seed=<int>

The current implementation can import data from [plink bed](http://pngu.mgh.harvard.edu/~purcell/plink/data.shtml#bed) 
format and the original Structure format. If the data are in plink format, ensure that
bed, bim and fam files for the dataset are all present in the same path.

### Main options

The key options to pass to the scripts are the input file, the output file and the number of populations.
Assuming the input file is named `genotypes.bed` (with corresponding `genotypes.fam` and `genotypes.bim`),
the output file is named `genotypes_output` and the number of populations you would like is 3, 
you can run the algorithm as follows:

    python structure.py -K 3 --input=genotypes --output=genotypes_output

This generates a `genotypes_output.3.log` file that tracks how the algorithm proceeds, and files
`genotypes_output.3.meanQ` and `genotypes_output.3.meanP` containing the posterior mean of
admixture proportions and allele frequencies, respectively. Note that input file names should
not include suffixes (e.g., .bed) and are relative to the main project directory (unless a full
path is provided).

## Running on test data

A test simulated dataset is provided in `test/testdata.bed` with genotypes sampled for
200 individuals at 500 SNP loci. The output files in `test/` were generated as follows:

    $ python structure.py -K 3 --input=test/testdata --output=testoutput_simple --full --seed=100
    $ ls test/testoutput_simple*
    test/testoutput_simple.3.log  test/testoutput_simple.3.meanP  test/testoutput_simple.3.meanQ  
    test/testoutput_simple.3.varP  test/testoutput_simple.3.varQ

    $ python structure.py -K 3 --input=test/testdata --output=testoutput_logistic --full --seed=100 --prior=logistic
    $ ls test/testoutput_logistic*
    test/testoutput_logistic.3.log    test/testoutput_logistic.3.meanQ  test/testoutput_logistic.3.varQ
    test/testoutput_logistic.3.meanP  test/testoutput_logistic.3.varP

Executing the code with the provided test data should generate a log file identical to the ones in `test/`, 
as a final check that the source code has been downloaded and compiled correctly.
