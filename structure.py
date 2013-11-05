import numpy as np
#np.random.seed(10)
import fastStructure 
import parse_bed
import getopt
import sys
import os, pdb

def parseopts(opts):

    params = {'mintol': 1e-6,
            'prior': "simple",
            'cv': 0,
            'full': False
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
                print "%s prior is not currently implemented, defaulting to the simple prior"
                params['prior'] = 'simple'

        elif opt in ["--cv"]:
            params['cv'] = int(arg)
        
        elif opt in ["--tol"]:
            params['mintol'] = float(arg)

        elif opt in ["--full"]:
            params['full'] = True

        elif opt in ["--seed"]:
            np.random.seed(int(arg))

    return params

def checkopts(params):

    if params['mintol']<=0:
        print "a non-positive value was provided as convergence criterion"
        raise ValueError
    
    if params['cv']<0:
        print "a negative value was provided for the number of cross-validations folds"
        raise ValueError

    if params['K']<=0:
        print "a negative value was provided for the number of populations"
        raise ValueError
    
    if not params.has_key('inputfile'):
        print "an input file needs to be provided"
        raise ValueError

    if not params.has_key('outputfile'):
        print "an output file needs to be provided"
        raise ValueError
    
def write_output(Q, P, other, params):

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

if __name__=="__main__":
    project_path = os.path.split(os.getcwd())[0]

    # parse command-line options
    argv = sys.argv[1:]
    smallflags = "K:"
    bigflags = ["prior=", "tol=", "input=", "output=", "cv=", "seed=", "full"] 
    try:
        opts, args = getopt.getopt(argv, smallflags, bigflags)
    except getopt.GetoptError:
        print "Incorrect arguments passed"
        sys.exit(2)

    params = parseopts(opts)

    # check if command-line options are valid
    try:
        checkopts(params)
    except ValueError:
        sys.exit(2)

    # load data
    G = parse_bed.load(params['inputfile'])
    G = np.require(G, dtype=np.uint8, requirements='C')

    Q, P, other = fastStructure.infer_variational_parameters(G, params['K'], \
                    params['outputfile'], params['mintol'], \
                    params['prior'], params['cv'])

    # write out inferred parameters
    write_output(Q, P, other, params)
