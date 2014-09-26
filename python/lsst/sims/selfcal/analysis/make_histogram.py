######
# Lynne Jones, ljones@astro.washington.edu
#   8/2010
#
# $Id$
#
# tools for making histograms with big data files
#
####

import numpy 
import pylab
from matplotlib.ticker import NullFormatter

import sys
from optparse import OptionParser

figformat = "png"

def make_histogram(filename, columns, bin_numbers=None, file_lineskip=1, titles=None, xlabels=None, savefig=False):
    """Read data from 'filename' in columns 'columns' and create a histogram from the result.

    This method reads the data and stores it in memory, to allow pylab to create the histogram
    with the optimal ranges for the data, while the number of bins can be specified with 'bin_numbers'.
    Can optionally subsample the file, reading only every file_lineskip line.
    Returns the data read from 'filename' in a dictionary of numpy arrays so if this was called from an
    interactive python shell, the histogram could be redone to visually examine the results of
    different ranges and number of bins. """
    # filename = name of data file.
    # columns = column numbers to read from data file. Count starts at 1. List of ints.
    # bin_numbers = can be specified (otherwise defaults to 100).
    #               gives number of bins in histogram. List of ints, same length as columns.
    # file_lineskip = read only lines which are multiples of file_lineskip.
    for c in columns:
        if not (isinstance(c, int)):
            raise Exception('Error - columns should be list of ints.')
    if bin_numbers != None:
        if ((len(bin_numbers) != len(columns)) & (len(bin_numbers)>1)):
            raise Exception('Error - bin_numbers should be same length as columns or length 1.')
    # Set up dictionary and lists to save data while being read.
    data = {}
    for c in columns:
        data[c] = []
    # Read the data.
    # Open data file. 
    f = open(filename, 'r')    
    line_num = 0
    for line in f:
        line_num = line_num + 1
        # Skip comment lines.
        if (line.startswith("#") | line.startswith("!")):
            continue
        # Skip the lines which are not multiples of file_lineskip (to enable subsampling file). 
        #  i.e. if file_lineskip = 10, this will read only every tenth line.
        if ((line_num % file_lineskip) != 0):
            continue
        values = line.split()
        # If there are not enough values in the line, quit reading (assume end of file). 
        if len(values)<max(columns):
            break
        # Assign data to dictionary. 
        for c in columns:
            data[c].append(values[c-1])
    # Close file.
    f.close()
    # Convert to numpy arrays.
    for c in columns:
        data[c] = numpy.array(data[c], dtype='float')
    # Set up to create histograms.
    # Set bin_numbers to default value, if not set by user.
    if bin_numbers == None:
        bin_numbers = [100]
    if len(bin_numbers)==1:
        for c in columns:
            bin_numbers.append(bin_numbers[0])
    # Create histograms - a new figure for each column.
    i = 0
    n = {}
    b = {}
    p = {}
    for c in columns:
        pylab.figure()
        n[c], b[c], p[c] = pylab.hist(data[c], bins=bin_numbers[i])
        if titles == None:
            pylab.title("%s Column %d" %(filename, c))
        else:
            pylab.title(titles[i])
        if xlabels != None:
            pylab.xlabel(xlabels[i])
        if savefig:
            figname = "hist_%d" %(c)
            pylab.savefig(figname+"."+figformat, format=figformat)
        i = i + 1
    # Calculate some basic statistics for output.
    stats = {}
    statlist = ('data_min', 'data_max', 'data_ave', 'hist_min', 'hist_max', 'hist_ave')
    for c in columns:
        stats[c] = {}        
        stats[c]['data_min'] = data[c].min()
        stats[c]['data_max'] = data[c].max()
        stats[c]['data_ave'] = data[c].sum() / float(len(data[c]))
        stats[c]['hist_min'] = n[c].min()
        stats[c]['hist_max'] = n[c].max()
        stats[c]['hist_ave'] = n[c].sum() / float(len(n[c]))
    print ""
    writestring = "# column "
    for key in statlist:
        writestring = writestring + " %s " %(key)
    print writestring
    for c in columns:
        writestring = "c %d " %(c)
        for key in statlist:
            writestring = writestring + "%g " %(stats[c][key])
        print writestring
    return data, stats

def make_2d_scatterhist(filename, columns, bin_numbers=None, file_lineskip=1, titles=None, xlabels=None, savefig=False):
    """Read data from 'filename' in columns 'columns' and create a histogram from the result.

    This method reads the data and stores it in memory, to allow pylab to create the histogram
    with the optimal ranges for the data, while the number of bins can be specified with 'bin_numbers'.
    Can optionally subsample the file, reading only every file_lineskip line.
    Returns the data read from 'filename' in a dictionary of numpy arrays so if this was called from an
    interactive python shell, the histogram could be redone to visually examine the results of
    different ranges and number of bins. """
    # filename = name of data file.
    # columns = column numbers to read from data file. Count starts at 1. List of ints. Can only be 2 columns.
    # bin_numbers = can be specified (otherwise defaults to 100).
    #               gives number of bins in histogram. List of ints, same length as columns.
    # file_lineskip = read only lines which are multiples of file_lineskip.
    if len(columns) > 2:
        raise Exception('This routine can only handle 2 columns.')
    for c in columns:
        if not (isinstance(c, int)):
            raise Exception('Error - columns should be list of ints.')
    if bin_numbers != None:
        if ((len(bin_numbers) != len(columns)) & (len(bin_numbers)>1)):
            raise Exception('Error - bin_numbers should be same length as columns or length 1.')
    # Set up dictionary and lists to save data while being read.
    data = {}
    for c in columns:
        data[c] = []
    # Read the data.
    # Open data file. 
    f = open(filename, 'r')    
    line_num = 0
    for line in f:
        line_num = line_num + 1
        # Skip comment lines.
        if (line.startswith("#") | line.startswith("!")):
            continue
        # Skip the lines which are not multiples of file_lineskip (to enable subsampling file). 
        #  i.e. if file_lineskip = 10, this will read only every tenth line.
        if ((line_num % file_lineskip) != 0):
            continue
        values = line.split()
        # If there are not enough values in the line, quit reading (assume end of file). 
        if len(values)<max(columns):
            break
        # Assign data to dictionary. 
        for c in columns:
            data[c].append(values[c-1])
    # Close file.
    f.close()
    # Convert to numpy arrays.
    for c in columns:
        data[c] = numpy.array(data[c], dtype='float')
    # Set up to create histograms.
    # Set bin_numbers to default value, if not set by user.
    if bin_numbers == None:
        bin_numbers = [100]
    if len(bin_numbers)==1:
        for c in columns:
            bin_numbers.append(bin_numbers[0])
    # Create scatter plot / histograms. One figure only.
    nullfmt   = NullFormatter()         # no labels    
    # definitions for the axes - set any values.
    left, width = 0.1, 0.65
    bottom, height = 0.1, 0.65
    bottom_h = left_h = left+width+0.02
    rect_scatter = [left, bottom, width, height]
    rect_histx = [left, bottom_h, width, 0.2]
    rect_histy = [left_h, bottom, 0.2, height]
    # start with a rectangular Figure
    pylab.figure(1, figsize=(8,8))
    axScatter = pylab.axes(rect_scatter)
    axHistx = pylab.axes(rect_histx)
    axHisty = pylab.axes(rect_histy)
    # no labels
    axHistx.xaxis.set_major_formatter(nullfmt)
    axHisty.yaxis.set_major_formatter(nullfmt)
    # Create the scatter plot.
    axScatter.scatter(data[columns[0]], data[columns[1]])
    # Now determine nice limits by hand.    
    xymax = numpy.max( [numpy.max(data[columns[0]]), numpy.max(data[columns[1]])] )
    xymin = numpy.min( [numpy.min(data[columns[0]]), numpy.min(data[columns[1]])] )
    binwidth = (xymax - xymin) / bin_numbers[0]
    lim = ( int(xymax/binwidth)+0.5) * binwidth
    axScatter.set_xlim( (-lim, lim) )
    axScatter.set_ylim( (-lim, lim) )
    bins = numpy.arange(-lim, lim + binwidth, binwidth)
    n = {}
    n[columns[0]], b, p = axHistx.hist(data[columns[0]], bins=bins)
    n[columns[1]], b, p = axHisty.hist(data[columns[1]], bins=bins, orientation='horizontal')    
    axHistx.set_xlim( axScatter.get_xlim() )
    axHisty.set_ylim( axScatter.get_ylim() )
    if xlabels != None:
        pylab.xlabel(xlabels[0])
        pylab.ylabel(xlabels[1])
    if savefig:
        figname = "hist_%d" %(c)
        pylab.savefig(figname+"."+figformat, format=figformat)
    # Calculate some basic statistics for output.
    stats = {}
    statlist = ('data_min', 'data_max', 'data_ave', 'hist_min', 'hist_max', 'hist_ave')
    for c in columns:
        stats[c] = {}        
        stats[c]['data_min'] = data[c].min()
        stats[c]['data_max'] = data[c].max()
        stats[c]['data_ave'] = data[c].sum() / float(len(data[c]))
        stats[c]['hist_min'] = n[c].min()
        stats[c]['hist_max'] = n[c].max()
        stats[c]['hist_ave'] = n[c].sum() / float(len(n[c]))
    print ""
    writestring = "# column "
    for key in statlist:
        writestring = writestring + " %s " %(key)
    print writestring
    for c in columns:
        writestring = "c %d " %(c)
        for key in statlist:
            writestring = writestring + "%g " %(stats[c][key])
        print writestring
    return data, stats

#########copy above routine and but in 2d plot.
                   
def make_histogram_predef(filename, columns, bin_ranges, bin_numbers, file_lineskip=1,
                          titles=None, xlabels=None,  savefig=False):
    """Read data from 'filename' in columns 'columns' and create a histogram from the result.

    This method reads the data and immediately bins the data to produce a histogram, without saving
    the data in memory. Uses bin_ranges and bin_numbers to create histogram bins.
    Can optionally subsample the file, reading only every file_lineskip line.
    Returns the histogram binned data in a dictionary of numpy arrays, so new plots could be created."""
    # filename = name of data file.
    # columns = column numbers to read from data file. Count starts at 1. List of ints.
    # bin_ranges = must be specified. gives lower/upper ranges for the histogram.
    #              List of pairs of values, same length as columns.
    # bin_numbers = must be specified. gives number of bins in histogram.
    #               gives number of bins in histogram. List of ints, same length as columns.
    # file_lineskip = read only lines which are multiples of file_lineskip.
    # Check input types.
    for c in columns:
        if not (isinstance(c, int)):
            raise Exception('Error - columns should be list of ints.')
    if (bin_numbers ==None) | (bin_ranges == None):
        raise Exception('Error - bin_ranges and bin_numbers have to be predefined.')
    if ((len(bin_numbers) != len(columns)) & (len(bin_numbers)>1)):
        raise Exception('Error - bin_numbers should be same length as columns or length 1.')
    if len(bin_ranges) != len(columns):
        raise Exception('Error - bin_ranges should be same length as columns. (2 pairs per list item).')
    # Set up histogram bins and values to help assign values into histogram bins.
    histobins = {}
    histovals = {}
    minval = {}
    stepval = {}
    nbinmax = {}
    i = 0
    for c in columns:
        minval[c] = bin_ranges[i][0]
        hival = bin_ranges[i][1]
        stepval[c] = (hival - minval[c])/float(bin_numbers[i])
        histobins[c] = numpy.arange(minval[c], hival+stepval[c], stepval[c], dtype='float')
        histovals[c] = numpy.zeros(len(histobins[c]), dtype='int')
        nbinmax[c] = len(histobins[c])
        i = i + 1
    # Set up for calculating some basic statistics for output.
    stats = {}
    statlist = ('data_min', 'data_max', 'data_ave', 'hist_min', 'hist_max', 'hist_ave')
    for c in columns:
        stats[c] = {}
        stats[c]['data_min'] = 1e9
        stats[c]['data_max'] = -1e9
        stats[c]['data_ave' ] = 0
    # Read the data.
    # Open data file. 
    f = open(filename, 'r')    
    line_num = 0
    for line in f:
        line_num = line_num + 1
        # Skip comment lines.
        if (line.startswith("#") | line.startswith("!")):
            continue
        # Skip the lines which are not multiples of file_lineskip (to enable subsampling file). 
        #  i.e. if file_lineskip = 10, this will read only every tenth line.
        if (line_num % file_lineskip != 0):
            continue
        values = line.split()        
        # If there are not enough values in the line, quit reading (assume end of file). 
        if len(values)<max(columns):
            break
        # Assign data to histogram bins.        
        for c in columns:
            dataval = float(values[c-1])
            histidx = min(int((dataval - minval[c])/stepval[c]), nbinmax[c]-1)
            histovals[c][histidx] = histovals[c][histidx] +  1
            # And calculate the min/max/ave of the data.
            stats[c]['data_min'] = min(dataval, stats[c]['data_min'])
            stats[c]['data_max'] = max(dataval, stats[c]['data_max'])
            stats[c]['data_ave'] = stats[c]['data_ave'] + dataval
    # Close file.
    f.close()
    # Print some basic output.
    for c in columns:
        print "# For column ", c, " used ", len(histovals[c]), " bins."
        print "# And ", histovals[c].sum(), " values from the data file."
        stats[c]['data_ave'] = stats[c]['data_ave'] / histovals[c].sum()
        stats[c]['hist_min'] = histovals[c].min()
        stats[c]['hist_max'] = histovals[c].max()
        stats[c]['hist_ave'] = histovals[c].sum() / float(len(histovals[c]))
    print ""
    writestring = "# column "
    for key in statlist:
        writestring = writestring + " %s " %(key)
    print writestring
    for c in columns:
        writestring = "c %d " %(c)
        for key in statlist:
            writestring = writestring + "%g " %(stats[c][key])
        print writestring            
    # Make histogram plots.
    i = 0
    for c in columns:
        pylab.figure()
        pylab.bar(histobins[c], histovals[c], width=stepval[c], linewidth=0)
        if titles == None:
            pylab.title("%s Column %d" %(filename, c))
        else:
            pylab.title(titles[i])
        if xlabels != None:
            pylab.xlabel(xlabels[i])
        if savefig:
            figname = "hist_%d" %(c)
            pylab.savefig(figname+"."+figformat, format=figformat)
        i = i + 1
    return histobins, histovals


def parse_options(argv):
    usage = "usage: %prog [options] datafilename"
    epilog ="Create histogram from a datafile, either by using pylab's built-in histogram function (which requires saving the data in memory but allows more dynamic setting of the histogram limits) or by assigning the data directly into histogram bins (which does not save the data in memory but requires predefining the histogram bins)."
    parser = OptionParser(usage, epilog=epilog)
    parser.add_option("-l", "--linesample", action="store", type="int", dest="filelines",
                      help="subsample the data file by reading only multiples of linesample", default=1)
    parser.add_option("-c", "--columns", action="store", dest="columns", type="string",
                      help="list of the columns to read from the data file, starting at 1. Enclose in ' '.",
                      default="1")
    parser.add_option("-x", "--xlabel_hist", action="store", dest="hist_xlabels", type="string",
                      help="list of labels to apply to x axis for each histogram - separate with whitespace and enclose in ' '.", default=None)
    parser.add_option("-p", "--predefine", action="store_true", dest="use_predef",
                      help="use predefined bins for histogram and don't save all data in memory",
                      default=False)
    parser.add_option("-b", "--bins_hist", action="store", dest="num_histbins", type="string",
                      help="the number of bins to use in each histogram - (list) - enclose in ''", default=None)
    parser.add_option("-r", "--range_hist", action="store",  dest="range_histbins", type="string",
                      help="the ranges of the bins to use in each histogram - (list of pairs) - enclose in ''",
                      default=None)
    parser.add_option("-t", "--titles_hist", action="store", dest="hist_titles", type="string",
                      help="either single title to apply to each histogram or list of titles - separate with whitespace and enclose in ' '.",
                      default=None)
    parser.add_option("-s", "--save_hists", action="store_true", dest="savefig",
                      help="save histograms to files", default=False)
    parser.add_option("-2", "--2d_scatter/hist", action="store_true", dest="do2d",
                      help="Create a 2-d scatter plot and histogram", default=False)
    (options, argv) = parser.parse_args()
    if len(argv)!=1:
        parser.error("Incorrect number of arguments. Must supply filename.")
    filename = argv[0]
    return options, filename

if __name__ == "__main__":
    """Create histogram from a datafile, either by using pylab's built-in histogram function
    (which requires saving the data in memory but allows more dynamic setting of the histogram
    limits) or by assigning the data directly into histogram bins (which does not save the data
    in memory but requires predefining the histogram bins).  """
    # jobname, pbsfilename, walltime and job command can be input from command line
    options, filename = parse_options(sys.argv)
    filelines = options.filelines
    strcolumns = options.columns
    strxlabels = options.hist_xlabels
    strtitles = options.hist_titles
    use_predef = options.use_predef
    strnum_histbins = options.num_histbins
    strrange_histbins = options.range_histbins
    savefig = options.savefig
    do2d = options.do2d

    # Need some parsing here to accomodate input from command line.
    # The problem is the options only reads on big string from the command line -
    # I have to separate it into individual list items here. 

    # Convert columns to list for input to routines.
    columns = []
    tmp = ""
    for char in strcolumns:
        if char == "[" :
            continue
        if ((char==" ") | (char==",") | (char=="]")):
            if ((tmp!="") & (tmp!=" ") & (tmp!=",") & (tmp!="]")):
                columns.append(int(tmp))
            tmp = ""
        else:
            tmp = tmp + char
    if ((tmp!="") & (tmp!=" ") & (tmp!=",") & (tmp!="]")):
        columns.append(int(tmp))

    # Convert num_histbins to list for input to routines.
    if strnum_histbins == None:
        num_histbins=None
    else:
        num_histbins = []
        tmp = ""
        for char in strnum_histbins:
            if char == "[" :
                continue
            if ((char==" ") | (char==",") | (char=="]")):
                if ((tmp!="") & (tmp!=" ") & (tmp!=",") & (tmp!="]")):
                    num_histbins.append(int(tmp))
                tmp = ""
            else:
                tmp = tmp + char
        if ((tmp!="") & (tmp!=" ") & (tmp!=",") & (tmp!="]")):
            num_histbins.append(int(tmp))

    # Convert range_histbins to list for input to routines.
    if strrange_histbins == None:
        range_histbins=None
    else:
        range_histbins = [[],]
        i  = 0
        tmp = ""
        for char in strrange_histbins:
            if char == "[":
                if i > 0:
                    range_histbins.append([])
                tmp = ""
            elif (char == "]") & (tmp!=""):
                range_histbins[i].append(float(tmp))
                tmp = ""
                i = i + 1
            elif ((char==" ") | (char==",")):
                if ((tmp!="") & (tmp!=" ") & (tmp!=",") & (tmp!="]")):
                    range_histbins[i].append(float(tmp))
                tmp = ""
            else:
                tmp = tmp + char
        if ((tmp!="") & (tmp!=" ") & (tmp!=",") & (tmp!="]") & (tmp!="[")):
            range_histbins[i].append(float(tmp))

    # Parse the xlabels.
    if strxlabels != None:
        xlabels = []
        tmp = strxlabels.split(" ")
        for t in tmp:
            xlabels.append(t.strip("[").strip(",").strip("]"))
    else:
        xlabels = None
    # and the titles.
    if strtitles != None:
        titles = []
        tmp = strtitles.split(" ")
        for t in tmp:
            titles.append(t.strip("[").strip(",").strip("]"))
    else:
        titles = None

    # Check on xlabels and titles to see if they match the number of columns.
    if xlabels != None:
        if len(xlabels) != len(columns):
            if len(xlabels) == 1:
                # then just duplicate xlabels so that matches number of columns
                newxlabels = []
                for i in range(len(columns)):
                    newlabels.append(xlabels[0])
                xlabels = newlabels
            else:
                # then the number of columns and the number of labels do not match easily
                raise Exception("xlabels (length %d) can be either length 1 or the length of 'columns' (%d)" \
                                %(len(xlabels), len(columns)))
    # Also for the titles.
    if titles != None:
        if len(titles) != len(columns):
            if len(titles) == 1:
                newtitles = []
                for i in range(len(columns)):
                    newtitles.append(titles[0])
                titles = newtitles
            else:
                raise Exception("titles can be either length 1 or the length of 'columns'")

    # Send some feedback about usage.
    print "# Will reading from file ", filename
    if use_predef:
        print "# Using predefined histogram bins / predef routine for creating histogram."
    else:
        print "# Not using predefined histogram bins / pylab version for creating histogram."
    print "# Creating histograms from columns ", columns
    if num_histbins != None:
        print "# Using ", num_histbins, " numbers for the histogram bins."
    if range_histbins != None:
        print "# Using ", range_histbins, " for the histogram bin ranges."
    if xlabels != None:
        print "# Using ", xlabels, " for the histogram x labels."
    if titles != None:
        print "# Using ", titles, " for the histogram titles."
    if do2d:
        print "# Going to create a 2-d scatter plot with histograms."

    # Read data and make histograms.
    if not use_predef:
        if do2d:
            make_2d_scatterhist(filename, columns, num_histbins, filelines, titles, xlabels, savefig)
        else:
            make_histogram(filename, columns, num_histbins, filelines, titles, xlabels, savefig)                    
    if use_predef:
        make_histogram_predef(filename, columns, range_histbins, num_histbins, filelines, titles, xlabels, savefig)

    # Show the plots to the screen.
    pylab.show()
