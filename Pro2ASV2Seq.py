#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
from collections import defaultdict
import re


def translate_dna(dna_seq):
    """Translate DNA sequence to protein sequence in all 6 reading frames"""
    codon_table = {
        'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
        'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
        'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*',
        'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W',
        'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
        'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
        'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
        'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
        'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
        'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
        'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
        'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
        'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
        'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
        'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
        'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G'
    }

    translations = []
    dna_seq = dna_seq.upper()

    # Forward 3 reading frames
    for frame in range(3):
        protein = ''
        for i in range(frame, len(dna_seq) - 2, 3):
            codon = dna_seq[i:i + 3]
            protein += codon_table.get(codon, 'X')
        translations.append(protein)

    # Reverse complement
    complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
    rev_comp = ''.join([complement.get(base, 'N') for base in reversed(dna_seq)])

    # Reverse complement 3 reading frames
    for frame in range(3):
        protein = ''
        for i in range(frame, len(rev_comp) - 2, 3):
            codon = rev_comp[i:i + 3]
            protein += codon_table.get(codon, 'X')
        translations.append(protein)

    return translations


def read_fasta(fasta_file):
    """Read FASTA file and return sequence information dictionary"""
    sequences = {}
    current_id = None
    current_seq = []
    current_sample = None
    current_size = None

    with open(fasta_file, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if current_id:
                    # Save previous sequence
                    seq_str = ''.join(current_seq)
                    sequences[current_id] = {
                        'sequence': seq_str,
                        'sample': current_sample,
                        'size': current_size
                    }

                # Parse new sequence header
                header = line[1:]

                # Extract sequence ID (first field before semicolon)
                current_id = header.split(';')[0]

                # Extract sample and size information
                sample_match = re.search(r'sample=([^;]+)', header)
                size_match = re.search(r'size=(\d+)', header)

                current_sample = sample_match.group(1) if sample_match else 'Unknown'
                current_size = int(size_match.group(1)) if size_match else 1

                current_seq = []
            else:
                current_seq.append(line)

        # Save last sequence
        if current_id:
            seq_str = ''.join(current_seq)
            sequences[current_id] = {
                'sequence': seq_str,
                'sample': current_sample,
                'size': current_size
            }

    return sequences


def process_data(table_file, seq_file, output_file):
    """Main processing function"""
    # Read FASTA sequences
    print("Reading FASTA sequence file...")
    sequences = read_fasta(seq_file)
    print(f"Total sequences read: {len(sequences)}")

    # Build protein sequence index
    print("Translating DNA sequences and building protein index...")
    protein_to_dna = defaultdict(lambda: {
        'dna_seqs': set(),
        'samples': set(),
        'total_copies': 0
    })

    for seq_id, seq_info in sequences.items():
        dna_seq = seq_info['sequence']
        sample = seq_info['sample']
        size = seq_info['size']

        # Translate to 6 reading frames
        translations = translate_dna(dna_seq)

        # Build index for each translation
        for translation in translations:
            if translation:
                protein_to_dna[translation]['dna_seqs'].add(seq_id)
                protein_to_dna[translation]['samples'].add(sample)
                protein_to_dna[translation]['total_copies'] += size

    print(f"Unique protein sequences generated: {len(protein_to_dna)}")

    # Read and process table file
    print("Processing table file...")
    with open(table_file, 'r') as f_in, open(output_file, 'w') as f_out:
        # Read header
        header = f_in.readline().strip()

        # Add tab prefix to header and new columns
        new_header = '\t' + header + '\tNum_DNA_Variants\tDNA_Samples\tTotal_DNA_Copies\n'
        f_out.write(new_header)

        # Process each line
        processed = 0
        matched = 0

        for line in f_in:
            line = line.strip()
            if not line:
                continue

            fields = line.split('\t')
            protein_seq = fields[-3]  # Sequence column (third to last)

            # Find matching DNA sequence information
            if protein_seq in protein_to_dna:
                info = protein_to_dna[protein_seq]
                num_variants = len(info['dna_seqs'])
                samples = ','.join(sorted(info['samples']))
                total_copies = info['total_copies']
                matched += 1
            else:
                num_variants = 0
                samples = ''
                total_copies = 0

            # Write new line without tab prefix
            new_line = f"{line}\t{num_variants}\t{samples}\t{total_copies}\n"
            f_out.write(new_line)

            processed += 1
            if processed % 1000 == 0:
                print(f"Processed {processed} lines...")

        print(f"Processing complete. Total lines processed: {processed}")
        print(f"Matched protein sequences: {matched}")


def main():
    parser = argparse.ArgumentParser(
        description='Translate DNA sequences to protein and count DNA variants for each protein'
    )
    parser.add_argument('-table', required=True, help='Input table file (tab-separated)')
    parser.add_argument('-seq', required=True, help='Input FASTA sequence file')
    parser.add_argument('-o', required=True, help='Output file name')

    args = parser.parse_args()

    process_data(args.table, args.seq, args.o)
    print(f"Results saved to {args.o}")


if __name__ == '__main__':
    main()