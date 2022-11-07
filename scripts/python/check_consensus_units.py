from collections import Counter

import numpy as np

def get_consensus_unit_new(msa: str) -> str:    
    """ Given a comma-separated set of units representing a multiple sequence
    alignment (from DB), determine the consensus unit at each column of the alignment and 
    return the consensus unit. Insertion columns where the number of 
    insertions >= number of nucleotides are skipped. For each column, most common nt
    is selected as consensus, in case of tie: pick one at random.

    Parameters
    msa (str):  String representation of a multiple sequence alignment, with units
                delimited by commas
    
    Returns  
    consensus_unit (str): 
                A consensus unit derived from the provided msa
    """    
    # Convert msa into list of lists, convert into np.Array() and transpose
    msa_matrix_t = np.array([list(unit) for unit in msa.split(",")]).transpose()
    consensus_unit = ""

    for msa_col in msa_matrix_t:
        # if half or more of the msa column are gap entries, skip column
        if np.count_nonzero(msa_col == '-') >= 0.5 * len(msa_col):
            continue

        # discard gap entries, get most common nt
        msa_col = msa_col[msa_col != '-']
        consensus_unit += Counter(msa_col).most_common(1)[0][0] # most_common() will return e.g. [('A': 5)]

    return consensus_unit


def main():
    input_file = "/cfs/earth/scratch/verb/projects/dopamine_transporter_repeats/results/human/tsv_refined_p0.05_d0.05/merged_p0.05_d0.05.tsv"
    output_file = "/cfs/earth/scratch/verb/projects/dopamine_transporter_repeats/results/human/tsv_refined_p0.05_d0.05/merged_corrected_consensus_p0.05_d0.05.tsv"
    
    with open(output_file, 'w') as o:
        with open(input_file, 'r') as f:
            o.write(next(f)) # write header line to output file
            for i, line in enumerate(f):
                # Column 8 is the consensus unit, 9 is the msa
                line_split = line.strip().split("\t")
                consensus_unit = get_consensus_unit_new(msa=line_split[9])
                line_split[8] = consensus_unit

                out_line = '\t'.join(line_split)
                o.write(f"{out_line}\n")

    # with open(input_file, 'r') as f:
    #     next(f) # skip header line
    #     for i, line in enumerate(f):
    #         if 1 >= 15:
    #             exit()
    #         # Column 8 is the consensus unit, 9 is the msa
    #         line_split = line.strip().split("\t")
    #         consensus_unit = get_consensus_unit_new(msa=line_split[9])

    #         if not consensus_unit == line_split[8]:
    #                 print(i + 2)
    #                 print(line_split[8])
    #                 print(consensus_unit)
    #                 print()
    #                 exit()


if __name__ == "__main__":
    main()
