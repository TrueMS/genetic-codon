#!/usr/bin/env python3

import argparse
from collections import Counter


def parse_arguments():
    parser = argparse.ArgumentParser(description='Calculate character frequency at each position in aligned sequences')
    parser.add_argument('-i', '--input', required=True, help='Input aligned sequence file (FASTA format)')
    parser.add_argument('-o', '--output', required=True, help='Output file path')
    return parser.parse_args()


def read_fasta(filepath):
    sequences = []
    with open(filepath, 'r') as f:
        seq = []
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if seq:
                    sequences.append(''.join(seq))
                    seq = []
            else:
                seq.append(line.upper())
        if seq:
            sequences.append(''.join(seq))
    return sequences


def calculate_frequencies(sequences):
    alignment_length = len(sequences[0])
    num_sequences = len(sequences)

    # Get all unique characters
    all_chars = set()
    for seq in sequences:
        all_chars.update(seq)
    all_chars = sorted(all_chars)

    # Calculate frequencies for each position
    position_stats = []
    dominant_chars = []
    dominant_freqs = []

    for pos in range(alignment_length):
        column = [seq[pos] for seq in sequences]
        counts = Counter(column)

        # Store frequencies for this position
        pos_freqs = {}

        # First, calculate * and - frequencies based on total sequences
        star_count = counts.get('*', 0)
        dash_count = counts.get('-', 0)
        star_freq = star_count / num_sequences
        dash_freq = dash_count / num_sequences

        # Calculate non-* and non- count for normalization
        excluded_count = star_count + dash_count
        valid_count = num_sequences - excluded_count

        for char in all_chars:
            if char == '*':
                pos_freqs[char] = star_freq
            elif char == '-':
                pos_freqs[char] = dash_freq
            else:
                # For other characters, calculate frequency relative to valid (non-* and non-) total
                if valid_count > 0:
                    pos_freqs[char] = counts.get(char, 0) / valid_count
                else:
                    pos_freqs[char] = 0.0

        position_stats.append(pos_freqs)

        # Find dominant character (among non-* and non- characters)
        valid_counts = {k: v for k, v in counts.items() if k not in ['*', '-']}
        if valid_counts:
            max_char = max(valid_counts, key=valid_counts.get)
            if valid_count > 0:
                dominant_freqs.append(valid_counts[max_char] / valid_count)
            else:
                dominant_freqs.append(0.0)
        else:
            # If only * and - exist, choose the one with higher count
            if star_count >= dash_count:
                max_char = '*'
                dominant_freqs.append(star_freq)
            else:
                max_char = '-'
                dominant_freqs.append(dash_freq)

        dominant_chars.append(max_char)

    return all_chars, position_stats, dominant_chars, dominant_freqs


def save_results(output_path, all_chars, position_stats, dominant_chars, dominant_freqs):
    with open(output_path, 'w') as f:
        # Write header
        header = '\t' + '\t'.join([f'Position_{i + 1}' for i in range(len(position_stats))])
        f.write(header + '\n')

        # Write character frequencies
        for char in all_chars:
            row = char
            for pos_freq in position_stats:
                row += f'\t{pos_freq[char]:.4f}'
            f.write(row + '\n')

        # Write dominant character row
        row = 'Dominant_Character\t' + '\t'.join(dominant_chars)
        f.write(row + '\n')

        # Write dominant frequency row
        row = 'Dominant_Frequency\t' + '\t'.join([f'{freq:.4f}' for freq in dominant_freqs])
        f.write(row + '\n')


def main():
    args = parse_arguments()

    # Read sequences
    sequences = read_fasta(args.input)
    print(f"Read {len(sequences)} sequences, length: {len(sequences[0])}")

    # Calculate frequencies
    all_chars, position_stats, dominant_chars, dominant_freqs = calculate_frequencies(sequences)

    # Save results
    save_results(args.output, all_chars, position_stats, dominant_chars, dominant_freqs)
    print(f"Results saved to: {args.output}")


if __name__ == "__main__":
    main()