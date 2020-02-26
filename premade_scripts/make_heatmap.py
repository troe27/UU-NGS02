### Import libraries ########################################################################
import matplotlib                        # plotting library
import matplotlib.pyplot as plt          # plotting library
plt.switch_backend('TkAgg')              # switch plotting backend to enable X11 forwarding
import matplotlib.gridspec as gridspec   # more plotting utilities
import cyvcf2                            # python/cython wrapper around htslib, for reading/writing VCF files
import numpy as np                       # math
import argparse                          # library to make the command line interface pretty


### function to create the CLI ##############################################################
def cli_parser():
    '''Parse command line input.'''
    parser_main = argparse.ArgumentParser(prog='plot_heatmap.py',        # define name of the program
                                          description="plot a heatmap")  # define description
    parser_main.add_argument("-i","--input_vcf",                         # define the long and short flag for input
                             help="path/to/input.vcf",                   # define the associated help msg.
                             required = True)                            # make it a required argument
    args = parser_main.parse_args()                                      # execute parsing
    return args                                                          # return resulting arguments

### function to parse the VCF into a genotype matrix ########################################
def load_vcf(input_vcf, threads=1, aaf_thresh=0.0):
    """ """
    vcf = cyvcf2.VCF(input_vcf, gts012=True,threads=threads) # load the vcf
    gts = []                                                 #init empty list for genotype_entries
    aaf = []                                                 #init empty list for alternative allele frequencies
    chr_pos = []                                             #init empty list for positions

    for variant in vcf:                                      # for each variant/position, do:
        if variant.aaf > aaf_thresh:                         # if the alt. allele freq. is above threshold:
            gts.append(variant.gt_types.astype(int))         # append genotype array to gts
            chr_pos.append(variant.POS)                      # append position to position-list
            aaf.append(variant.aaf)                          # append aaf to alt. allele. freq. list

    gt_array = np.array(gts)                                 # make list of per-position arrays into SAMPLE x POS rectangular matrix
    samples = vcf.samples                                    # extract list of sample names from vcf

    return aaf, [chr_pos, samples, gt_array]                 # return aaf, and genotype matrix with column and row names

### function to make a rolling mean for the allele frequencies ##############################
def make_sliding_aaf_mean(aaf,winsize=10):
    """takes array and windowsize, default 10 SNPs """
    rmean = np.convolve(aaf, np.ones((winsize,))/winsize, mode='valid')  # https://en.wikipedia.org/wiki/Convolution
    return rmean                                                         # return rolling mean



def heatmap(data, row_labels, col_labels, ax=None, **kwargs):
    """
    Create a heatmap from a numpy array and two lists of labels.
    Parameters
    ----------
    data
        A 2D numpy array of shape (N, M).
    row_labels
        A list or array of length N with the labels for the rows.
    col_labels
        A list or array of length M with the labels for the columns.
    ax
        A `matplotlib.axes.Axes` instance to which the heatmap is plotted.  If
        not provided, use current axes or create a new one.  Optional.
    **kwargs
        All other arguments are forwarded to `imshow`.
    """
    if not ax:                                              # if no plotting ax exists:
        ax = plt.gca()                                      # create it
    im = ax.imshow(data, aspect="auto", **kwargs)           # Plot the heatmap
    ax.set_xticks([0,11800])                                # set one tick at the beginning, one at the end
    ax.set_yticks(np.arange(data.shape[0]))                 # set as many y-ticks as there are samples
    ax.set_xticklabels(["\"Beginning\"", "\"End\""])        # ignoring col_labels, merely plotting "beginning" and "end"
    ax.set_yticklabels(row_labels)                          # row-labels are the samples

    return im                                               # return the plot

### plotting ############################################################################
def make_plots(gt_array, rmean):
    """ """
    figure = plt.figure(constrained_layout=True, figsize=(20,20)) # create the figure with a size of 20x20
    gs = figure.add_gridspec(2, 1, height_ratios=[10,1])          # make a grid of 2x1 subplots,
                                                                  # i.e. two plots on top of each other
                                                                  # with a height ratio of 10:1
    figure_ax1 = figure.add_subplot(gs[0])                        # ax1 is the top plot
    figure_ax2 = figure.add_subplot(gs[1])                        # ax2 is the smaller bottom plot


    heatmap(gt_array[2].transpose(),                              # plot the heatmap on the big-plot
            row_labels=gt_array[1],                               # add labels
            col_labels=gt_array[0],                               # add labels
            cmap="Greys",                                         # define the colormap, here we only want Black/White
            ax=figure_ax1)                                        # plot on ax1

    figure_ax2.plot(rmean, color="Black")                         # plot the rolling mean on ax2
    figure_ax2.get_xaxis().set_ticks([])                          # hide x-tick-labels
    figure_ax2.set_xlim(0,2690)                                   # set limits for the x-axis
    figure_ax2.set_ylim(0,0.5)                                    # set limits for the y-axis

    plt.show()                                                    # show the plot



### main function ############################################################################
def main():
    """ define the overall worflow"""
    args = cli_parser()                             # parse the commandline input
    aaf, gt_array = load_vcf(args.input_vcf)        # load the vcf, extract gt_matrix and alt. allele freq.
    rmean = make_sliding_aaf_mean(aaf)              # calculate alt. allele. freq sliding mean
    make_plots(rmean=rmean, gt_array = gt_array)    # plot aaf and heatmap


### execute code ############################################################################
main()
