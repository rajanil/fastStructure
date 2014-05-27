import numpy as np
import matplotlib.pyplot as plot
import colorsys
import getopt
import sys, pdb

def plot_admixture(admixture, population_indices, population_labels, title):

    N,K = admixture.shape
    colors = [colorsys.hsv_to_rgb(h,0.9,0.7) for h in np.linspace(0,1,K+1)[:-1]]
    text_color = 'k'
    bg_color = 'w'
    fontsize = 12

    figure = plot.figure(figsize=(5,3))

    xmin = 0.13
    ymin = 0.2
    height = 0.6
    width = 0.74
    indiv_width = width/N
    subplot = figure.add_axes([xmin,ymin,width,height])
    [spine.set_linewidth(0.001) for spine in subplot.spines.values()]

    for k in xrange(K):
        if k:
            bottoms = admixture[:,:k].sum(1)
        else:
            bottoms = np.zeros((N,),dtype=float)

        lefts = np.arange(N)*indiv_width
        subplot.bar(lefts, admixture[:,k], width=indiv_width, bottom=bottoms, facecolor=colors[k], edgecolor=colors[k], linewidth=0.4)

        subplot.axis([0, N*indiv_width, 0, 1])
        subplot.tick_params(axis='both', top=False, right=False, left=False, bottom=False)
        xtick_labels = tuple(map(str,['']*N))
        subplot.set_xticklabels(xtick_labels)
        ytick_labels = tuple(map(str,['']*K))
        subplot.set_yticklabels(ytick_labels)

    position = subplot.get_position()
    title_position = (0.5, 0.9)
    figure.text(title_position[0], title_position[1], title, fontsize=fontsize, \
        color='k', horizontalalignment='center', verticalalignment='center')

    for p,popname in enumerate(population_labels):
        indices = np.where(population_indices==p)[0]
        vline_pos = (indices.max()+1)*indiv_width 
        subplot.axvline(vline_pos, linestyle='-', linewidth=0.2, c='#888888')
        label_position = (xmin+(2*indices.min()+indices.size)*0.5*indiv_width, ymin-0.01)
        figure.text(label_position[0], label_position[1], popname, fontsize=6, color='k', \
            horizontalalignment='right', verticalalignment='top', rotation=70)

    return figure

def get_admixture_proportions(params):

    # load admixture proportions
    handle = open('%s.%d.meanQ'%(params['inputfile'],params['K']),'r')
    admixture = np.array([line.strip().split() for line in handle]).astype('float')
    handle.close()
    N,K = admixture.shape
    admixture = admixture/admixture.sum(1).reshape(N,1)

    # get population labels
    if params.has_key('popfile'):
        handle = open(params['popfile'],'r')
        populations = [line.strip().split() for line in handle]
        handle.close()
        population_labels = list(np.unique(populations))
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

        elif opt in ["--title"]:
            params['title'] = arg

    return params

def usage():

    """
    brief description of various flags and options for this script
    """

    print "\nHere is how you can use this script\n"
    print "Usage: python %s"%sys.argv[0]
    print "\t -K <int>"
    print "\t --input=<file>"
    print "\t --output=<file>"
    print "\t --popfile=<file> (optional)"
    print "\t --title=<figure title> (optional)"


if __name__=="__main__":

    # parse command-line options
    argv = sys.argv[1:]
    smallflags = "K:"
    bigflags = ["input=", "output=", "title="]
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
    if params.has_key('title'):
        title = params['title']
    else:
        title = params['inputfile']

    # plot the data
    figfile = '%s.png'%params['outputfile']
    figure = plot_admixture(admixture, population_indices, population_labels, title)
    figure.savefig(figfile, dpi=300, format='png')
