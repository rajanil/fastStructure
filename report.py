import numpy as np
import colorsys
import getopt
import sys, pdb


def get_admixture_proportions(params):

    # load admixture proportions
    handle = open('%s.%d.meanQ'%(params['outputfile'],params['K']),'r')
    admixture = np.array([line.strip().split() for line in handle]).astype('float')
    handle.close()
    N,K = admixture.shape
    admixture = admixture/admixture.sum(1).reshape(N,1)


    # get population labels
    if params.has_key('popfile'):
        handle = open(params['popfile'],'r')
        populations = [line.strip() for line in handle]
        handle.close()
        population_labels = list(set(populations))

        # group populations by cluster similarity
        population_cluster = [np.mean(admixture[[i for i,p in enumerate(populations) if p==label],:],0).argmax() \
            for label in population_labels]
        population_labels = [population_labels[j] for j in np.argsort(population_cluster)]
        population_indices = np.array([population_labels.index(pop) for pop in populations])

        # re-order samples in admixture matrix
        order = np.argsort(population_indices)
        population_indices = population_indices[order]
        admixture = admixture[order,:]

    else:

        print "file with population labels is not provided or does not exist .... \ncreating population labels based on inferred admixture proportions"
        population_labels = ['population %d'%i for i in xrange(1,K+1)]
        population_indices = np.argmax(admixture,1)

        # re-order samples in admixture matrix
        order = np.argsort(population_indices)
        population_indices = population_indices[order]
        admixture = admixture[order,:]
        order = np.arange(N)
        for k in xrange(K):
            order[population_indices==k] = order[population_indices==k][np.argsort(admixture[population_indices==k,:][:,k])[::-1]]
        admixture = admixture[order,:]

    return admixture, population_indices, population_labels

def generate_report(params, population_indices, population_labels):

    """
    generate report with accessions, population_indices, population_labels in an array
    """

    report_array = []
    i = 0

    with open('%s.fam'%(params['inputfile']), 'r') as reader:
        for line in reader:
            line_array = line.strip("\n").strip("\r").strip("\r\n").strip(" ").split(" ")

            report_array.append([
                str(line_array[0]), 
                str(line_array[1]), 
                str(population_indices[i]), 
                str(population_labels[population_indices[i]])
            ])

            i += 1
    
    return report_array


def parseopts(opts):

    """
    parses the command-line flags and options passed to the script
    """

    params = {}

    for opt, arg in opts:

        if opt in ["-K"]:
            params['K'] = int(arg)

        elif opt in ["--input"]:
            params['inputfile'] = arg

        elif opt in ["--output"]:
            params['outputfile'] = arg

        elif opt in ["--popfile"]:
            params['popfile'] = arg

    return params

def usage():

    """
    brief description of various flags and options for this script
    """

    print "\nHere is how you can use this script\n"
    print "Usage: python %s"%sys.argv[0]
    print "\t -K <int>  (number of populations)"
    print "\t --input=<file>  (/path/to/input/file; same as input flag passed to structure.py)"
    print "\t --output=<file> (/path/to/output/file; same as output flag passed to structure.py)"
    print "\t --popfile=<file> (file with known categorical labels; optional)"


if __name__=="__main__":

    # parse command-line options
    argv = sys.argv[1:]
    smallflags = "K:"
    bigflags = ["input=", "output=", "popfile="]
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
    
    # get the data to be plotted
    admixture, population_indices, population_labels = get_admixture_proportions(params)

    # generate report
    report_array = generate_report(params, population_indices, population_labels)

    # output the report
    if params.has_key('popfile'):
        output_file_path = '%s.%d.report.pop.txt'%(params['outputfile'],params['K'])
    else:
        output_file_path = '%s.%d.report.txt'%(params['outputfile'],params['K'])
    try:
        with open(output_file_path,'w') as writer:
            writer.write(''.join([ str('\t'.join(single_row)) + '\n' for single_row in report_array ]))
    except:
        print "Unexpected error when writing the report"
