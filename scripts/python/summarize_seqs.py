#!/usr/bin/env python3
import os
import argparse

from Bio import SeqIO

def cla_parser():
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--dir", "-d", type=str, required=True, help="Directory where sequences (fasta format) are stored"
    )

    return parser.parse_args()

def get_fasta_files(seq_dir, extensions={"fna", "fasta", "fa"}):
    if not isinstance(extensions, (list, set, tuple)):
        raise ValueError("Specify extensions in a list, set or tuple")
    dir_contents = [os.path.join(seq_dir, any_path) for any_path in os.listdir(seq_dir)]

    return list(filter(lambda x: any(x.endswith(i) for i in extensions), dir_contents))

def summarize_seqs(seq_dir):
    fasta_files = get_fasta_files(seq_dir)
    
    for fasta_file in fasta_files:
        for record in SeqIO.parse(fasta_file, "fasta"):
            print(record.id)


def main():
    args = cla_parser()

    summarize_seqs(args.dir)

if __name__ == "__main__":
    main()
