#!/usr/bin/env python3
"""
Count stop codon usage frequencies in Prodigal-predicted gene files
"""

import argparse
from pathlib import Path
from Bio import SeqIO


def count_stop_codons(gene_file):
    """
    Count stop codon usage in a gene file
    Gene sequences output by Prodigal include stop codons
    """
    stop_codons = {'TAA': 0, 'TAG': 0, 'TGA': 0}
    
    try:
        for record in SeqIO.parse(str(gene_file), 'fasta'):
            seq_str = str(record.seq).upper()
            if len(seq_str) >= 3:
                last_codon = seq_str[-3:]
                if last_codon in stop_codons:
                    stop_codons[last_codon] += 1
    except Exception as e:
        print(f"Warning: error reading file {gene_file}: {e}")
        return None
    
    # Calculate frequencies
    total = sum(stop_codons.values())
    if total == 0:
        frequencies = {'TAA': 0.0, 'TAG': 0.0, 'TGA': 0.0}
    else:
        frequencies = {codon: count / total for codon, count in stop_codons.items()}
    
    return {
        'counts': stop_codons,
        'frequencies': frequencies,
        'total': total
    }


def main():
    parser = argparse.ArgumentParser(
        description='Count stop codon usage frequencies in Prodigal-predicted gene files'
    )
    parser.add_argument(
        '-i', '--input',
        required=True,
        help='Input folder path containing Prodigal-predicted gene files'
    )
    parser.add_argument(
        '-o', '--output',
        required=True,
        help='Output file path'
    )
    
    args = parser.parse_args()
    
    input_dir = Path(args.input)
    output_file = Path(args.output)
    
    if not input_dir.exists():
        print(f"Error: input folder {input_dir} does not exist")
        return
    
    if not input_dir.is_dir():
        print(f"Error: {input_dir} is not a directory")
        return
    
    # Supported file extensions
    supported_extensions = ['.fa', '.fasta', '.fna', '.ffn']
    
    # Get all gene files
    gene_files = []
    for ext in supported_extensions:
        gene_files.extend(input_dir.glob(f'*{ext}'))
        gene_files.extend(input_dir.glob(f'*{ext.upper()}'))
    
    gene_files = sorted(set(gene_files))
    
    if not gene_files:
        print(f"Warning: no gene files found in {input_dir}")
        return
    
    print(f"Found {len(gene_files)} gene files")
    
    # Analyze each file
    results = []
    for filepath in gene_files:
        print(f"Processing: {filepath.name}")
        result = count_stop_codons(filepath)
        if result is not None:
            results.append({
                'filename': filepath.name,
                **result
            })
    
    if not results:
        print("Error: no gene files were processed successfully")
        return
    
    # Write the results file
    with open(output_file, 'w', encoding='utf-8') as f:
        header = ['Genome', 'TAA_count', 'TAG_count', 'TGA_count',
                  'TAA_frequency', 'TAG_frequency', 'TGA_frequency', 'Total_genes']
        f.write('\t'.join(header) + '\n')
        
        for result in results:
            row = [
                result['filename'],
                str(result['counts']['TAA']),
                str(result['counts']['TAG']),
                str(result['counts']['TGA']),
                f"{result['frequencies']['TAA']:.4f}",
                f"{result['frequencies']['TAG']:.4f}",
                f"{result['frequencies']['TGA']:.4f}",
                str(result['total'])
            ]
            f.write('\t'.join(row) + '\n')
    
    print(f"Results saved to: {output_file}")
    print(f"Processed {len(results)} gene files in total")


if __name__ == '__main__':
    main()
