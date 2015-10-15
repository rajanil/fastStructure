
import numpy as np
import fastStructure 
import parse_bed
import parse_str
import random
import getopt
import sys
import pdb
import warnings

# ignore warnings with these expressions
warnings.filterwarnings('ignore', '.*divide by zero.*',)
warnings.filterwarnings('ignore', '.*invalid value.*',)

def parseopts(opts):

    """
    parses the command-line flags and options passed to the script
    """

    params = {'mintol': 1e-6,
            'prior': "simple",
            'cv': 0,
            'full': False,
            'format': 'bed'
            }

    for opt, arg in opts:

        if opt in ["-K"]:
            params['K'] = int(arg)

        elif opt in ["--input"]:
            params['inputfile'] = arg

        elif opt in ["--output"]:
            params['outputfile'] = arg

        elif opt in ["--prior"]:
            params['prior'] = arg

            if params['prior'] not in ['simple','logistic']:
                print("%s prior is not currently implemented, defaulting to the simple prior")
                params['prior'] = 'simple'

        elif opt in ["--format"]:
            params['format'] = arg

        elif opt in ["--cv"]:
            params['cv'] = int(arg)
        
        elif opt in ["--tol"]:
            params['mintol'] = float(arg)

        elif opt in ["--full"]:
            params['full'] = True

        elif opt in ["--seed"]:
            np.random.seed(int(arg))
            random.seed(int(arg))

    return params

def checkopts(params):

    """
    checks if some of the command-line options passed are valid.
    In the case of invalid options, an exception is always thrown.
    """

    if params['mintol']<=0:
        print("a non-positive value was provided as convergence criterion")
        raise ValueError
    
    if params['cv']<0:
        print("a negative value was provided for the number of cross-validations folds")
        raise ValueError

    if 'K' not in params:
        print("a positive integer should be provided for number of populations")
        raise KeyError

    if params['format'] not in ['bed','str']:
        print("%s data format is not currently implemented")
        raise ValueError

    if params['K']<=0:
        print("a negative value was provided for the number of populations")
        raise ValueError
    
    if 'inputfile' not in params:
        print("an input file needs to be provided")
        raise KeyError 

    if 'outputfile' not in params:
        print("an output file needs to be provided")
        raise KeyError
    
def write_output(Q, P, other, params):

    """
    write the posterior means and variational parameters
    to separate output files.
    """

    handle = open('%s.%d.meanQ'%(params['outputfile'],params['K']),'w')
    handle.write('\n'.join(['  '.join(['%.6f'%i for i in q]) for q in Q])+'\n')
    handle.close()

    handle = open('%s.%d.meanP'%(params['outputfile'],params['K']),'w')
    handle.write('\n'.join(['  '.join(['%.6f'%i for i in p]) for p in P])+'\n')
    handle.close()

    if params['full']:
        handle = open('%s.%d.varQ'%(params['outputfile'],params['K']),'w')
        handle.write('\n'.join(['  '.join(['%.6f'%i for i in q]) for q in other['varQ']])+'\n')
        handle.close()

        handle = open('%s.%d.varP'%(params['outputfile'],params['K']),'w')
        handle.write('\n'.join(['  '.join(['%.6f'%i for i in np.hstack((pb,pg))]) \
            for pb,pg in zip(other['varPb'],other['varPg'])])+'\n')
        handle.close()

def usage():
    
    """
    brief description of various flags and options for this script
    """

    print("\nHere is how you can use this script\n")
    print("Usage: python %s"%sys.argv[0])
    print("\t -K <int> (number of populations)")
    print("\t --input=<file> (/path/to/input/file)")
    print("\t --output=<file> (/path/to/output/file)")
    print("\t --tol=<float> (convergence criterion; default: 10e-6)")
    print("\t --prior={simple,logistic} (choice of prior; default: simple)")
    print("\t --cv=<int> (number of test sets for cross-validation, 0 implies no CV step; default: 0)")
    print("\t --format={bed,str} (format of input file; default: bed)")
    print("\t --full (to output all variational parameters; optional)")
    print("\t --seed=<int> (manually specify seed for random number generator; optional)")


if __name__=="__main__":

    # parse command-line options
    argv = sys.argv[1:]
    smallflags = "K:"
    bigflags = ["prior=", "tol=", "input=", "output=", "cv=", "seed=", "format=", "full"] 
    try:
        opts, args = getopt.getopt(argv, smallflags, bigflags)
        if not opts:
            usage()
            sys.exit(2)
    except getopt.GetoptError:
        print("Incorrect options passed")
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

    # run the variational algorithm
    Q, P, other = fastStructure.infer_variational_parameters(G, params['K'], \
                    params['outputfile'], params['mintol'], \
                    params['prior'], params['cv'])

    # write out inferred parameters
    write_output(Q, P, other, params)
