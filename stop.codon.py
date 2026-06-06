#!/usr/bin/env python3
"""
Detect stop codons in protein sequences and locate them in the corresponding DNA sequences
Infer possible amino acids at stop codon positions through homologous protein alignment
"""

import argparse
import sys
from collections import Counter


def parse_fasta(filename):
    """
    Parse a FASTA format file, return a dictionary {sequence_name: sequence_content}
    """
    sequences = {}
    current_name = None
    current_seq = []

    try:
        with open(filename, 'r') as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue

                if line.startswith('>'):
                    # Save the previous sequence
                    if current_name:
                        sequences[current_name] = ''.join(current_seq)
                    # Start a new sequence
                    current_name = line[1:].split()[0]  # Only take the part before the first space as the sequence name
                    current_seq = []
                else:
                    current_seq.append(line.upper())

            # Save the last sequence
            if current_name:
                sequences[current_name] = ''.join(current_seq)

    except FileNotFoundError:
        print(f"Error: File not found {filename}")
        sys.exit(1)
    except Exception as e:
        print(f"Error reading file {filename}: {e}")
        sys.exit(1)

    return sequences


def get_dna_name_from_protein(protein_name):
    """
    Get the corresponding DNA sequence name from the protein sequence name
    The DNA name is identical to the protein name
    """
    return protein_name


def translate_dna(dna_seq, frame_start):
    """
    Translate a DNA sequence
    frame_start: 0, 1, or 2, indicating which base to start from (0-based)
    """
    codon_table = {
        'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
        'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
        'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
        'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
        'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
        'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
        'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
        'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
        'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*',
        'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
        'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
        'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
        'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W',
        'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
        'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
        'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G'
    }

    protein = []
    for i in range(frame_start, len(dna_seq) - 2, 3):
        codon = dna_seq[i:i + 3]
        if len(codon) == 3:
            aa = codon_table.get(codon, 'X')
            protein.append(aa)

    return ''.join(protein)


def find_best_frame(protein_seq, dna_seq):
    """
    Determine the correct reading frame by aligning the protein sequence
    Only return the reading frame if there is an exact match, otherwise return None
    """
    for frame_start in [0, 1, 2]:
        translated = translate_dna(dna_seq, frame_start)

        # Check for an exact match
        if translated == protein_seq:
            return frame_start

    # No exact-match reading frame found
    return None


def find_stop_codons_with_positions(protein_seq, dna_seq, frame_start):
    """
    Find stop codons in the protein sequence, return the corresponding DNA codons and positions
    Also return the position of the stop codon in the protein sequence (used for homologous alignment)
    """
    results = []

    # Find all positions of '*' in the protein sequence
    for i, aa in enumerate(protein_seq):
        if aa == '*':
            # Calculate the start position of this stop codon in the DNA (1-based)
            dna_position_0based = frame_start + i * 3
            dna_position_1based = dna_position_0based + 1  # Convert to 1-based

            # Extract the corresponding codon
            if dna_position_0based + 3 <= len(dna_seq):
                codon = dna_seq[dna_position_0based:dna_position_0based + 3]
                results.append((codon, dna_position_1based, i))  # Add the position in the protein sequence

    return results


def find_sequence_in_alignment(target_seq, alignment_sequences):
    """
    Find the sequence in the alignment file that matches the target sequence
    Return the sequence name and the aligned sequence
    """
    # Remove stop codons from the target sequence for comparison
    target_no_stop = target_seq.replace('*', '')

    for seq_name, aligned_seq in alignment_sequences.items():
        # Remove gaps and stop codons from the aligned sequence for comparison
        seq_no_gap = aligned_seq.replace('-', '').replace('*', '')

        # Check if it is the same sequence (allowing differences in sequence names)
        if seq_no_gap == target_no_stop:
            return seq_name, aligned_seq

        # Also try partial matching (if the sequence names are similar)
        if any(part in seq_name for part in seq_name.split('_')) or \
                any(part in seq_name for part in seq_name.split('_')):
            if seq_no_gap == target_no_stop:
                return seq_name, aligned_seq

    return None, None


def map_position_to_alignment(protein_seq, aligned_seq, position):
    """
    Map a position in the protein sequence to a position in the aligned sequence
    protein_seq: the original protein sequence (may contain stop codons)
    aligned_seq: the aligned sequence (may contain gaps)
    position: the position in the original protein sequence (0-based)
    Return: the position in the aligned sequence (0-based)
    """
    # Compute the mapping relationship
    protein_index = 0
    alignment_index = 0

    while alignment_index < len(aligned_seq) and protein_index <= position:
        if aligned_seq[alignment_index] != '-':
            if protein_index == position:
                return alignment_index
            protein_index += 1
        alignment_index += 1

    return None


def analyze_homologous_positions_with_alignment(protein_seq, protein_name,
                                                alignment_sequences,
                                                stop_positions,
                                                max_predictions=3):
    """
    Analyze the amino acid distribution of homologous sequences at stop codon positions through the alignment file
    protein_seq: the protein sequence containing stop codons
    protein_name: the name of the protein sequence
    alignment_sequences: the dictionary of homologous sequence alignments
    stop_positions: the list of stop codon positions in the protein sequence (0-based)
    max_predictions: the maximum number of predictions to return
    """
    position_predictions = {}

    # Find the corresponding sequence in the alignment file
    query_name, query_aligned = find_sequence_in_alignment(protein_seq, alignment_sequences)

    if query_aligned is None:
        print(f"    Warning: Could not find a match for sequence {protein_name} in the homologous alignment file")
        return {pos: "Could not find corresponding sequence in the alignment file" for pos in stop_positions}

    print(f"    Found matching sequence in the alignment file: {query_name}")

    for pos in stop_positions:
        # Map the position in the original sequence to the position in the aligned sequence
        aligned_pos = map_position_to_alignment(protein_seq, query_aligned, pos)

        if aligned_pos is None:
            position_predictions[pos] = "Position mapping failed"
            continue

        # Collect the amino acids of all sequences at this alignment position
        amino_acids = []
        for seq_name, seq in alignment_sequences.items():
            # Skip the query sequence itself
            if seq_name == query_name:
                continue

            if aligned_pos < len(seq):
                aa = seq[aligned_pos]
                # Ignore gaps, stop codons, and unknown amino acids
                if aa not in ['-', '*', 'X']:
                    amino_acids.append(aa)

        if amino_acids:
            # Count amino acid frequencies
            aa_counter = Counter(amino_acids)
            total = len(amino_acids)

            # Get the most common amino acids and their proportions
            predictions = []
            for aa, count in aa_counter.most_common(max_predictions):
                percentage = (count / total) * 100
                predictions.append(f"{aa}({percentage:.1f}%)")

            position_predictions[pos] = ','.join(predictions) if predictions else "No prediction"
        else:
            position_predictions[pos] = "No valid homologous sequences"

    return position_predictions


def main():
    parser = argparse.ArgumentParser(description='Detect stop codons in protein sequences and infer possible amino acids')
    parser.add_argument('-g', '--dna', required=True, help='DNA sequence file (FASTA format)')
    parser.add_argument('-p', '--protein', required=True, help='Protein sequence file (FASTA format)')
    parser.add_argument('-pro', '--homologous', help='Homologous protein alignment file (FASTA format)')
    parser.add_argument('-m', '--max-predictions', type=int, default=3,
                        help='Maximum number of inferred amino acids to provide (default 3)')
    parser.add_argument('-o', '--output', required=True, help='Output file')

    args = parser.parse_args()

    # Read the sequence files
    print("Reading DNA sequence file...")
    dna_sequences = parse_fasta(args.dna)
    print(f"Read {len(dna_sequences)} DNA sequences")

    print("Reading protein sequence file...")
    protein_sequences = parse_fasta(args.protein)
    print(f"Read {len(protein_sequences)} protein sequences")

    # Read the homologous protein alignment file (if provided)
    homologous_sequences = {}
    if args.homologous:
        print("Reading homologous protein alignment file...")
        homologous_sequences = parse_fasta(args.homologous)
        print(f"Read {len(homologous_sequences)} homologous sequences")

    # Check for stop codons
    results = []
    has_stop = False
    skipped = 0  # Count protein sequences skipped because no matching DNA was found

    print("\nDetecting stop codons...")
    for protein_name, protein_seq in protein_sequences.items():
        # Get the corresponding DNA sequence name
        dna_name = get_dna_name_from_protein(protein_name)

        # Only process sequences that have a corresponding DNA sequence; otherwise skip silently
        if dna_name not in dna_sequences:
            skipped += 1
            continue

        dna_seq = dna_sequences[dna_name]

        # Find the correct reading frame
        frame_start = find_best_frame(protein_seq, dna_seq)

        if frame_start is None:
            # No correct reading frame found; skip this sequence
            print(f"  Skip {protein_name}: no correct reading frame found")
            skipped += 1
            continue

        frame_number = frame_start + 1  # Convert to 1-based for display
        print(f"  Processing {protein_name}, using reading frame {frame_number}")

        # Find stop codons
        if '*' in protein_seq:
            stop_codon_info = find_stop_codons_with_positions(protein_seq, dna_seq, frame_start)

            if stop_codon_info:
                has_stop = True

                # If there are homologous sequences, perform prediction
                position_predictions = {}
                if homologous_sequences:
                    stop_positions = [info[2] for info in stop_codon_info]
                    position_predictions = analyze_homologous_positions_with_alignment(
                        protein_seq, protein_name, homologous_sequences,
                        stop_positions, args.max_predictions
                    )

                # Record each stop codon
                for codon, dna_position, protein_position in stop_codon_info:
                    prediction = position_predictions.get(protein_position, "No prediction") if homologous_sequences else "-"
                    results.append((protein_name, dna_name, codon, dna_position, protein_position + 1, prediction))
                    print(f"    Found stop codon: {codon} (DNA position: {dna_position}, protein position: {protein_position + 1})")
                    if prediction != "-" and prediction != "No prediction" and "No" not in prediction:
                        print(f"      Possible amino acids: {prediction}")
        else:
            print(f"    No stop codon found")

    # Output the results
    with open(args.output, 'w') as f:
        if not has_stop:
            f.write("No stop codon exists\n")
            print("\nResult: No stop codon exists")
        else:
            # Write the stop codon information
            f.write("## Stop Codon Detection Results ##\n")
            if homologous_sequences:
                f.write("Protein Sequence Name\tDNA Sequence Name\tStop Codon\tDNA Start Position\tProtein Position\tInferred Amino Acid\n")
            else:
                f.write("Protein Sequence Name\tDNA Sequence Name\tStop Codon\tDNA Start Position\tProtein Position\n")

            for result in results:
                if homologous_sequences:
                    protein_name, dna_name, codon, dna_pos, protein_pos, prediction = result
                    f.write(f"{protein_name}\t{dna_name}\t{codon}\t{dna_pos}\t{protein_pos}\t{prediction}\n")
                else:
                    protein_name, dna_name, codon, dna_pos, protein_pos, _ = result
                    f.write(f"{protein_name}\t{dna_name}\t{codon}\t{dna_pos}\t{protein_pos}\n")

            print(f"\nFound {len(results)} stop codons in total")

        if skipped:
            print(f"Skipped {skipped} protein sequences without a matching DNA sequence")

        print(f"Results saved to {args.output}")


if __name__ == "__main__":
    main()