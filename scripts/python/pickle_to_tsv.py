#!/usr/bin/env python3 
import os
import pickle
import argparse

import numpy as np

from check_consensus_units import get_consensus_unit_new

def cla_parser():
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "-p", "--pickle", type=str, required=True, help="Path to pickle file to turn into tsv"
    )
    parser.add_argument(
        "-o", "--outdir", type=str, required=False, help="Path to output directory where tsv file will be written"
    )
    parser.add_argument(
        "-i", "--identifier", type=str, required=False, help="If specified, this string will be added as the first field in each row in the output tsv.\
            Useful to e.g. associate sequence identifiers with each repeat in the output when comparing multiple sequences."
    )

    return parser.parse_args()

# def get_consensus_unit(repeat):
#     """ Determine the consensus unit/motif for a repeata from a given multiple sequence alignment

#     Parameters
#     repeat (Repeat):    A TRAL tandem repeat for which a consensus unit will be determined

#     Returns
#     consensus_unit (str):
#                         The consensus unit for the provided repeat
#     """
#     nt_tuple = ("T", "C", "A", "G") # This is the order the nts appear in in a TRAL Repeat object's msaTDN attribute
#     consensus_unit = ""
    
#     for position in repeat.msaTDN: # iterate over msa columns in msaTDN
#         # get the index of the most occuring nt in the msa column, and select consensus nt
#         # in case of tie: take the first one in the order given in nt_tuple
#         consensus_nt_idx = np.asarray(position == np.amax(position)).nonzero()[0][0]
#         consensus_unit += nt_tuple[consensus_nt_idx]
    
#     return consensus_unit

def write_file(tr_list, destination, score, identifier=None):
    """ Write all tandem repeats from a list of Repeats to a specified output file (.tsv format)

    Parameters
    tr_list (list):
                    A basic Python list (NOT RepeatList) of Repeat instances to be written to file
    destination (str):
                    Name of file to be made
    """
    # Sort TRs based on order of appearance in the input sequence
    sorted_trs = sorted(tr_list, key=lambda x: x.begin)

    header = "\t".join(["begin",                        
                        "l_effective",
                        "n_effective",
                        "repeat_region_length",
                        "score",
                        "pvalue",
                        "divergence",
                        "consensus_unit",
                        "msa\n"])
    if identifier:
        header = "seq\t" + header
    with open(destination, "w") as f:
        f.write(header)
        for tr in sorted_trs:
            line = [str(i) for i in [
                    tr.begin,                    
                    tr.l_effective,
                    tr.n_effective,
                    tr.repeat_region_length,
                    tr.d_score[score],
                    tr.d_pvalue[score],
                    tr.d_divergence[score],
                    get_consensus_unit_new(",".join(tr.msa)),
                    ",".join(tr.msa)]]
            if identifier:
                line = [identifier] + line
            f.write("\t".join(line) + "\n")

def main():
    args = cla_parser()
    if not os.path.isfile(args.pickle) or not args.pickle.endswith(".pickle"):
        raise ValueError("Path specified to --pickle should point to file with '.pickle' extension")    
    
    if not args.outdir:
        output_path = args.pickle.replace("pickle", "tsv")
    else:
        output_path = os.path.join(args.outdir, args.pickle.split(os.sep)[-1].replace("pickle", "tsv"))

    with open(args.pickle, 'rb') as f:
        repeat_list = pickle.load(f)
    # for repeat in repeat_list.repeats:
    #     print(repeat)
    #     print("consensus: ", get_consensus_unit(repeat))
    #     print()
    write_file(tr_list=repeat_list.repeats, destination=output_path, score="phylo_gap01", identifier=args.identifier)


if __name__ == "__main__":
    main()
