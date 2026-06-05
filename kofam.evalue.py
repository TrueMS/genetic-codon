import argparse
import sys
import os
from collections import defaultdict

def process_line(line, remove_leading_tab=False):
    line = line.strip().lstrip('*')
    parts = line.split(None, 5)
    result = '\t'.join(parts)
    if remove_leading_tab and result.startswith('\t'):
        result = result[1:]
    return result + '\n'

def filter_file(input_file, e_value, n_hits):
    try:
        # Generate output file name
        base_name = os.path.splitext(os.path.basename(input_file))[0]
        output_file = f"{base_name}_{e_value:.0e}_top{n_hits}.out"

        gene_hits = defaultdict(list)

        with open(input_file, 'r') as infile:
            for line in infile:
                processed_line = process_line(line)
                fields = processed_line.strip().split('\t')
                
                if len(fields) >= 5:
                    try:
                        gene_id = fields[0]
                        fifth_column = float(fields[4])
                        if fifth_column <= e_value:
                            gene_hits[gene_id].append((fifth_column, processed_line))
                    except ValueError:
                        continue

        with open(output_file, 'w') as outfile:
            for gene_id in gene_hits:
                sorted_hits = sorted(gene_hits[gene_id], key=lambda x: x[0])
                top_hits = sorted_hits[:n_hits]
                for _, line in top_hits:
                    outfile.write(line)

        print(f"Filtered results saved to: {output_file}")
    except IOError as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)

def main():
    parser = argparse.ArgumentParser(description="Filter tab-separated file based on the fifth column value and keep top N hits per gene.")
    parser.add_argument('-i', '--input', required=True, help="Path to the input tab-separated file")
    parser.add_argument('-e', '--evalue', type=float, required=True, help="E-value threshold (e.g., 1e-5)")
    parser.add_argument('-n', '--nhits', type=int, default=1, help="Number of top hits to keep per gene (default: 1)")

    args = parser.parse_args()

    filter_file(args.input, args.evalue, args.nhits)

if __name__ == "__main__":
    main()
