# fastStructure

*fastStructure* is a fast algorithm for inferring population structure from large SNP genotype data, 
implemented in Python2.x. This repo contains a *vars* library with definitions and methods for
the relevant variables, and a set of scripts to load the data and run the algorithm.

This document summarizes how to setup this software package, compile the C and 
[Cython](http://cython.org/) scripts and run the algorithm on a test simulated genotype dataset.

## Setup

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


## Executing the code

## Running on test data
