import matplotlib                        # plotting library
import matplotlib.pyplot as plt          # plotting library
plt.switch_backend('TkAgg')              # switch plotting backend to enable X11 forwarding
import numpy as np                       # math
import argparse                          # library to make the command line interface pretty


## function to create the CLI ###############################################################
def cli_parser():                                                                           # "def" is the start of every function definition.
                                                                                            # here we define the function cli_parser()
    '''Parse command line input.'''                                                         # description of the function
    parser_main = argparse.ArgumentParser(prog='make_delta_af.py',                          # define name of the program
                                          description="plot position/delta-af as scatter")  # define description
    parser_main.add_argument("-a","--listA",                                                 # define the long and short flag for input A
                             help="path/to/list_with_allele_frequencies",                   # define the associated help msg.
                             required = True)                                               # make it a required argument
    parser_main.add_argument("-b","--listB",                                                     # define the long and short flag for input B
                             help="path/to/list_with_allele_frequencies",                   # define the associated help msg.
                             required = True)                                               # make it a required argument
    args = parser_main.parse_args()                                                         # execute parsing
    return args                                                                             # return resulting arguments


## load allele frequency lists ##############################################################
def load_aaf_list(listfile):                                                                # function to load the allele-frequency list
    """ load list into pos:aaf dictionary """
    with open(listfile, "r") as handle:                                                     #  open the file
        lines =handle.read().split("\n")                                                    # read the file, split the text into chunks at line-endings "\n"
        split_lines = []                                                                    # create empty list for split lines
        for line in lines:                                                                  # for each line do:
            split_lines.append(line.split(" "))                                             # split the line at spaces into the filed chromosome, position, aaf
        listdict = dict()                                                                   # create empty dictionary
        for sline in split_lines:                                                           # for each split line:
            if len(sline)==3:                                                               # if the line is 3 fields long:
                listdict[int(sline[1])] = float(sline[2])                                   # create dictionary with field 1 (pos) as key, and field 2 (aaf) as item
    return listdict                                                                         # return the dictionary

def load_aaf_list_V2(listfile):                                                             # same as the first one, just shorter and less readable. for fun.
    """ load list into pos:aaf dictionary """                                               # using nested dict and list comprehension
    return  {int(j[1]):float(j[2]) for j in [line.split(" ") for line in open(listfile, "r").read().split("\n")] if len(j)==3 }


## get allele-frequency difference #########################################################
def make_diff_dict(dictA, dictB):
    """ take the absolute allele-frequency between the two lists for each position """
    diff_d = {}                                                                            # init empty dict, same as "diff_d = dict()"
    for key in dictA.keys():                                                               # for each key (position) in dict A:
        diff = abs(dictA[key] - dictB[key])                                                # take the absolute of the difference between the delta af between the two list of samples.
        diff_d[key]=diff                                                                   # add it to the diff_dictionary as a pos:diff pair.
    return diff_d                                                                          # return dictionary


## subset allele-frequency difference ######################################################
def subset_d(d, thresh=0.9):                                                               # function to filter all positions with a high allele-frequency, default thresh is 0.9
    """ subsets the diff_dict """
    return {key:item for key, item in d.items() if item>thresh}                            # take all key:value pairs if the value is above the threshold


## define workflow and plot ################################################################
def main():                                                                                      # function that defines the workflow
    """ define workflow """
    args = cli_parser()                                                                          # parse commandline arguments
    A = load_aaf_list(args.listA)                                                                # load the first list using the long-form function
    B = load_aaf_list_V2(args.listB)                                                             # doing the same thing for the second list with the short-form function
    ddict = make_diff_dict(dictA=A, dictB=B)                                                     # get allele-frequency differences
    sddict = subset_d(d=ddict, thresh=0.9)                                                       # get all variants that have a allele-frequency difference above 0.9
    pos, delta = [key for key, item in ddict.items()], [item for key, item in ddict.items()]     # get data into list form
    pos2, delta2 = [key for key, item in sddict.items()], [item for key, item in sddict.items()] # get subset-data into list form
    fig, ax = plt.subplots(ncols=1, nrows=1, figsize=(30,5))                                     # create plot, one 30:5 subplot
    ax.scatter(pos,delta, color="Black", alpha=0.3, s=10)                                        # plot a scatter for all variants in black with an alpha ( see-throughness) of 0.3
    ax.scatter(pos2,delta2, color="Red", alpha=1, s=20,)                                         # plot a scatter for all variants above 0.9 in red with an alpha of 1  (solid)
    for i, txt in enumerate(pos2):                                                               # for each position in the subset and their index:
        ax.annotate(txt, (pos2[i], delta2[i]))                                                   # annotate the point with the position
    ax.get_xaxis().set_ticks([])                                                                 # dont plot the x-tick labels because it takes aaaaaaaaaaages.
    ax.set_xlim(0,4411212)                                                                       # set limits for the x-axis
    ax.set_ylabel("allele-frequency difference between group A and B")                           # set y-axis label
    plt.show()                                                                                   # display the plot

## execute workflow and plot ###############################################################
if __name__ == "__main__":                                                                 # ignore this
    main()                                                                                 # execute workflow
