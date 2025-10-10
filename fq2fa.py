#!/usr/bin/env python3
import os
import argparse
from Bio import SeqIO


def convert_fastq_to_fasta(input_folder, rename_titles=False, suffix=""):
    """
    Convert fastq files in a folder to fasta format

    Parameters:
    input_folder -- input folder path
    rename_titles -- whether to rename sequence names
    suffix -- file suffix, used to generate new sequence names
    """
    # Ensure input folder exists
    if not os.path.isdir(input_folder):
        print(f"Error: Input folder {input_folder} does not exist")
        return

    # Get all files in the folder
    files = os.listdir(input_folder)

    for file in files:
        input_path = os.path.join(input_folder, file)

        # Skip subfolders
        if os.path.isdir(input_path):
            continue

        # Generate output filename (change extension to .fasta)
        base_name = os.path.splitext(file)[0]
        output_path = os.path.join(input_folder, f"{base_name}.fasta")

        # If sequence names need to be renamed
        if rename_titles:
            # Get filename prefix (remove suffix)
            if suffix and file.endswith(suffix):
                name_prefix = file[:-(len(suffix))]
            else:
                name_prefix = base_name

            # Counter for sequence numbering
            counter = 0

            # Open input and output files
            with open(output_path, 'w') as output_handle:
                for record in SeqIO.parse(input_path, "fastq"):
                    counter += 1
                    # Modify sequence name
                    record.id = f"{name_prefix}_{counter}"
                    record.description = ""
                    # Write sequence in fasta format
                    output_handle.write(f">{record.id}\n{record.seq}\n")

            print(f"Converted {file} to {os.path.basename(output_path)} and renamed sequences")

        # Don't rename sequence names, just convert format
        else:
            SeqIO.convert(input_path, "fastq", output_path, "fasta")
            print(f"Converted {file} to {os.path.basename(output_path)}")


def main():
    # Create command-line argument parser
    parser = argparse.ArgumentParser(description="Convert fastq files to fasta format")
    parser.add_argument("-i", "--input", required=True, help="Input folder path")
    parser.add_argument("-title", "--rename_titles", action="store_true", help="Whether to rename sequence names")
    parser.add_argument("-pfix", "--suffix", default="", help="File suffix, used to generate new sequence names")

    # Parse command-line arguments
    args = parser.parse_args()

    # Execute conversion
    convert_fastq_to_fasta(args.input, args.rename_titles, args.suffix)


if __name__ == "__main__":
    main()