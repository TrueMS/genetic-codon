#!/usr/bin/env python3

import argparse
from collections import defaultdict, OrderedDict
import pandas as pd


def parse_arguments():
    parser = argparse.ArgumentParser(description='Merge identical protein sequences and sum their counts')
    parser.add_argument('-pro', '--protein', required=True, help='Input protein FASTA file')
    parser.add_argument('-t', '--table', required=True, help='Input ASV count table')
    parser.add_argument('-o', '--output', required=True, help='Output merged count table')
    return parser.parse_args()


def read_fasta(filepath):
    """Read FASTA file and return OrderedDict of {seq_id: sequence}"""
    sequences = OrderedDict()
    with open(filepath, 'r') as f:
        seq_id = None
        seq = []
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if seq_id:
                    sequences[seq_id] = ''.join(seq)
                seq_id = line[1:].split()[0]  # Get ID before first space
                seq = []
            else:
                seq.append(line.upper())
        if seq_id:
            sequences[seq_id] = ''.join(seq)
    return sequences


def find_sequence_groups(sequences):
    """Group sequence IDs by identical sequences"""
    seq_to_ids = defaultdict(list)
    for seq_id, seq in sequences.items():
        seq_to_ids[seq].append(seq_id)

    # Create mapping from all IDs to their representative (first occurrence)
    id_to_representative = {}
    representative_ids = []

    for seq, ids in seq_to_ids.items():
        representative = ids[0]  # First ID as representative
        representative_ids.append(representative)
        for seq_id in ids:
            id_to_representative[seq_id] = representative

    return id_to_representative, representative_ids


def merge_counts(table_file, id_to_representative, representative_ids):
    """Read count table and merge counts for identical sequences"""
    # Read the table
    df = pd.read_csv(table_file, sep='\t', index_col=0)

    # Create new dataframe for merged counts
    merged_df = pd.DataFrame(0, index=representative_ids, columns=df.columns)

    # Sum counts for each group
    for asv_id in df.index:
        if asv_id in id_to_representative:
            representative = id_to_representative[asv_id]
            merged_df.loc[representative] += df.loc[asv_id]

    # Remove rows with all zeros
    merged_df = merged_df.loc[(merged_df != 0).any(axis=1)]

    return merged_df


def main():
    args = parse_arguments()

    # Read protein sequences
    print(f"Reading protein sequences from {args.protein}")
    sequences = read_fasta(args.protein)
    print(f"Found {len(sequences)} sequences")

    # Find identical sequences and group them
    print("Identifying identical sequences...")
    id_to_representative, representative_ids = find_sequence_groups(sequences)

    # Count unique sequences
    unique_count = len(representative_ids)
    print(f"Found {unique_count} unique sequences")

    # Show merging information
    merge_count = 0
    for seq_id, rep in id_to_representative.items():
        if seq_id != rep:
            merge_count += 1
    if merge_count > 0:
        print(f"Merging {merge_count} duplicate sequences into their representatives")

    # Read and merge count table
    print(f"Reading count table from {args.table}")
    merged_df = merge_counts(args.table, id_to_representative, representative_ids)

    # Save results
    merged_df.to_csv(args.output, sep='\t')
    print(f"Merged count table saved to {args.output}")
    print(f"Final table contains {len(merged_df)} representative sequences")


if __name__ == "__main__":
    main()