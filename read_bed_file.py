#! /usr/bin/env python2

import numpy as np
#import fastStructure
import parse_bed
import parse_str
import getopt
import sys
import pdb

def parseopts(opts):

    """
    parses the command-line flags and options passed to the script
    """

    params = { 'format': 'bed' }

    for opt, arg in opts:
        if opt in ["--input"]:
            params['inputfile'] = arg

        elif opt in ["--format"]:
            params['format'] = arg

        elif opt in ["--output"]:
            params['outputfile'] = arg

    return params

def checkopts(params):

    """
    checks if some of the command-line options passed are valid.
    In the case of invalid options, an exception is always thrown.
    """

    if params['format'] not in ['bed','str']:
        print "%s data format is not currently implemented"
        raise ValueError

    if not params.has_key('inputfile'):
        print "an input file needs to be provided"
        raise KeyError

    if not params.has_key('outputfile'):
        print "an output file needs to be provided"
        raise KeyError

def usage():

    """
    brief description of various flags and options for this script
    """

    print "\nHere is how you can use this script\n"
    print "Usage: python %s"%sys.argv[0]
    print "\t --input=<file>"
    print "\t --format={bed,str} (default: bed)"
    print "\t --output=<file>"


if __name__=="__main__":

    # parse command-line options
    argv = sys.argv[1:]
    smallflags = "K:"
    bigflags = ["input=", "format=", "output="]
    try:
        opts, args = getopt.getopt(argv, smallflags, bigflags)
        if not opts:
            usage()
            sys.exit(2)
    except getopt.GetoptError:
        print "Incorrect options passed"
        usage()
        sys.exit(2)

    params = parseopts(opts)

    # check if command-line options are valid
    try:
        checkopts(params)
    except (ValueError,KeyError):
        sys.exit(2)

    # load data
    if params['format']=='bed':
        G = parse_bed.load(params['inputfile'])
    elif params['format']=='str':
        G = parse_str.load(params['inputfile'])
    G = np.require(G, dtype=np.uint8, requirements='C')

    # Write the genome file.
    handle = open('%s.genome' % (params['outputfile']), 'w')
    handle.write('\n'.join(['  '.join(['%d' % i for i in g]) for g in G])+'\n')
    handle.close()
