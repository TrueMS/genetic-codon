#!/usr/bin/env python3
"""
Batch processing of fasta files for sequence dereplication and denoising, generating ASV table
"""

import os
import sys
import argparse
import subprocess
import shutil
from pathlib import Path
from collections import defaultdict
import hashlib
import traceback


def run_command(cmd, description=""):
    """Run shell command"""
    print(f"Running: {description if description else cmd}")
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"Command failed: {cmd}")
        print(f"Error output: {result.stderr}")
        print(f"Standard output: {result.stdout}")
        return False
    print(f"Command completed successfully")
    return True


def process_single_file(input_file, output_dir, threads=20, keep_temp=False):
    """Process a single fasta file"""

    # Get base name (remove path and extension)
    base_name = Path(input_file).stem

    # Create independent temporary directory for each file
    temp_dir = Path(output_dir) / "temp" / base_name
    temp_dir.mkdir(parents=True, exist_ok=True)

    derep_file = temp_dir / f"dereplicated.{base_name}.fa"
    derep_sample_file = temp_dir / f"dereplicated_with_sample.{base_name}.fasta"
    final_file = Path(output_dir) / f"unoise.alpha2.{base_name}.fa"

    print(f"\n{'=' * 60}")
    print(f"Processing file: {input_file}")
    print(f"Base name: {base_name}")
    print(f"Temporary directory: {temp_dir}")
    print(f"{'=' * 60}")

    try:
        # Step 1: vsearch dereplication
        print(f"\nStep 1: Dereplicating sequences")
        print(f"Input: {input_file}")
        print(f"Output: {derep_file}")

        cmd1 = f"vsearch --derep_fulllength {input_file} --output {derep_file} --sizeout --threads {threads}"
        if not run_command(cmd1, "Dereplicating sequences"):
            print(f"❌ Error: Dereplication failed - {input_file}")
            return None

        # Check if dereplicated file was generated
        if not derep_file.exists():
            print(f"❌ Error: Dereplicated file not generated - {derep_file}")
            return None
        else:
            file_size = derep_file.stat().st_size
            print(f"✓ Dereplicated file generated: {derep_file} (size: {file_size} bytes)")

        # Step 2: Add sample information
        print(f"\nStep 2: Adding sample information")
        print(f"Input: {derep_file}")
        print(f"Output: {derep_sample_file}")

        with open(derep_file, 'r') as infile, open(derep_sample_file, 'w') as outfile:
            seq_count = 0
            for line in infile:
                if line.startswith('>'):
                    seq_count += 1
                    header = line.strip()
                    sample_name = base_name
                    new_header = f"{header};sample={sample_name}\n"
                    outfile.write(new_header)
                else:
                    outfile.write(line)

        # Check if sample file was generated
        if not derep_sample_file.exists():
            print(f"❌ Error: Sample file not generated - {derep_sample_file}")
            return None
        else:
            file_size = derep_sample_file.stat().st_size
            print(f"✓ Sample file generated: {derep_sample_file} (size: {file_size} bytes, sequences: {seq_count})")

        # Step 3: vsearch clustering and denoising
        print(f"\nStep 3: Clustering and denoising")
        print(f"Input: {derep_sample_file}")
        print(f"Output: {final_file}")

        cmd3 = f"vsearch --cluster_unoise {derep_sample_file} --sizein --sizeout --centroids {final_file} --threads {threads} --unoise_alpha 2"
        if not run_command(cmd3, "Clustering and denoising"):
            print(f"❌ Error: Clustering and denoising failed - {input_file}")
            return None

        # Check if final file was generated
        if not final_file.exists():
            print(f"❌ Error: Final file not generated - {final_file}")
            return None
        else:
            file_size = final_file.stat().st_size
            print(f"✓ Final file generated: {final_file} (size: {file_size} bytes)")

        print(f"\n✅ Successfully completed processing: {final_file}")

        # List temporary directory contents
        if keep_temp:
            print(f"\nTemporary files kept in: {temp_dir}")
            temp_files = list(temp_dir.iterdir())
            for tf in temp_files:
                print(f"  - {tf.name} ({tf.stat().st_size} bytes)")

        return final_file

    except Exception as e:
        print(f"\n❌ Exception occurred while processing file {input_file}: {str(e)}")
        traceback.print_exc()
        return None

    finally:
        # Clean up only when not keeping temporary files
        if not keep_temp:
            print(f"\nCleaning up temporary files...")
            try:
                if derep_file.exists():
                    os.remove(derep_file)
                    print(f"  Deleted: {derep_file.name}")
                if derep_sample_file.exists():
                    os.remove(derep_sample_file)
                    print(f"  Deleted: {derep_sample_file.name}")
                # Delete temporary directory if empty
                if temp_dir.exists() and not any(temp_dir.iterdir()):
                    temp_dir.rmdir()
                    print(f"  Deleted empty directory: {temp_dir.name}")
            except Exception as e:
                print(f"Error cleaning up temporary files: {e}")


def parse_fasta_with_metadata(file_path):
    """Parse fasta file with sample and size information"""
    sequences = {}
    current_id = None
    current_seq = []
    sample_info = {}

    with open(file_path, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if current_id:
                    sequences[current_id] = ''.join(current_seq)

                # Parse header: >AS_450;sample=AS;size=217
                parts = line[1:].split(';')
                current_id = parts[0]

                # Extract sample and size information
                metadata = {}
                for part in parts[1:]:
                    if '=' in part:
                        key, value = part.split('=', 1)
                        metadata[key] = value

                sample_info[current_id] = {
                    'sample': metadata.get('sample', 'unknown'),
                    'size': int(metadata.get('size', 0))
                }
                current_seq = []
            else:
                current_seq.append(line)

        if current_id:
            sequences[current_id] = ''.join(current_seq)

    return sequences, sample_info


def generate_asv_table(output_dir):
    """Generate ASV table and sequence file"""

    print(f"\n{'=' * 60}")
    print("Generating ASV table...")
    print(f"{'=' * 60}")

    # Collect all processed files
    processed_files = list(Path(output_dir).glob("unoise.alpha2.*.fa"))

    if not processed_files:
        print("❌ Error: No processed files found")
        return

    print(f"Found {len(processed_files)} processed files for ASV table generation")
    for pf in processed_files:
        print(f"  - {pf.name}")

    # Store all sequences and their abundance in each sample
    asv_sequences = {}  # {seq_hash: sequence}
    asv_abundance = defaultdict(lambda: defaultdict(int))  # {seq_hash: {sample: count}}

    # Read all files
    all_samples = set()

    for file_path in processed_files:
        print(f"\nReading file: {file_path.name}")
        try:
            sequences, sample_info = parse_fasta_with_metadata(file_path)
            print(f"  Found {len(sequences)} sequences")

            for seq_id, sequence in sequences.items():
                info = sample_info[seq_id]
                sample = info['sample']
                size = info['size']
                all_samples.add(sample)

                # Use sequence hash as unique identifier
                seq_hash = hashlib.md5(sequence.encode()).hexdigest()

                # Store sequence
                if seq_hash not in asv_sequences:
                    asv_sequences[seq_hash] = sequence

                # Accumulate abundance
                asv_abundance[seq_hash][sample] += size

        except Exception as e:
            print(f"❌ Error reading file {file_path}: {e}")
            traceback.print_exc()
            continue

    if not asv_sequences:
        print("❌ Error: No sequences read")
        return

    # Generate ASV IDs
    asv_mapping = {}
    asv_counter = 1
    for seq_hash in sorted(asv_sequences.keys()):
        asv_id = f"ASV_{asv_counter:04d}"
        asv_mapping[seq_hash] = asv_id
        asv_counter += 1

    # Write ASV sequence file
    asv_seq_file = Path(output_dir) / "ASV.seq.fa"
    with open(asv_seq_file, 'w') as f:
        for seq_hash, asv_id in sorted(asv_mapping.items(), key=lambda x: x[1]):
            f.write(f">{asv_id}\n")
            sequence = asv_sequences[seq_hash]
            # Wrap at 80 characters per line
            for i in range(0, len(sequence), 80):
                f.write(sequence[i:i + 80] + '\n')

    print(f"\n✓ ASV sequence file saved: {asv_seq_file}")

    # Write ASV abundance table
    asv_table_file = Path(output_dir) / "ASV.table.txt"
    sorted_samples = sorted(all_samples)

    with open(asv_table_file, 'w') as f:
        # Write header
        f.write("ASV_ID\t" + "\t".join(sorted_samples) + "\n")

        # Write data
        for seq_hash, asv_id in sorted(asv_mapping.items(), key=lambda x: x[1]):
            row = [asv_id]
            for sample in sorted_samples:
                count = asv_abundance[seq_hash].get(sample, 0)
                row.append(str(count))
            f.write("\t".join(row) + "\n")

    print(f"✓ ASV abundance table saved: {asv_table_file}")

    # Print statistics
    print(f"\n📊 Statistics:")
    print(f"  Total ASVs: {len(asv_mapping)}")
    print(f"  Number of samples: {len(all_samples)}")
    print(f"  Sample list: {', '.join(sorted_samples)}")


def main():
    parser = argparse.ArgumentParser(description='Batch processing of fasta files for sequence dereplication and denoising')
    parser.add_argument('-i', '--input', required=True, help='Input folder path')
    parser.add_argument('-o', '--output', required=True, help='Output folder path')
    parser.add_argument('-t', '--threads', type=int, default=20, help='Number of threads (default: 20)')
    parser.add_argument('--keep-temp', action='store_true', help='Keep intermediate files')

    args = parser.parse_args()

    # Check input directory
    input_dir = Path(args.input)
    if not input_dir.exists():
        print(f"❌ Error: Input directory does not exist: {input_dir}")
        sys.exit(1)

    # Create output directory
    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Find all fasta/fa/fq files
    input_files = []
    for ext in ['*.fasta', '*.fa', '*.fq']:
        input_files.extend(input_dir.glob(ext))

    if not input_files:
        print(f"❌ Error: No fasta files found in {input_dir}")
        sys.exit(1)

    print(f"\n📁 Found {len(input_files)} files to process:")
    for f in input_files:
        print(f"  - {f.name}")

    # Process each file
    processed_files = []
    failed_files = []

    for i, input_file in enumerate(input_files, 1):
        print(f"\n[{i}/{len(input_files)}] Starting processing...")
        result = process_single_file(input_file, output_dir, args.threads, args.keep_temp)
        if result:
            processed_files.append(result)
            print(f"✅ [{i}/{len(input_files)}] Success")
        else:
            failed_files.append(input_file)
            print(f"❌ [{i}/{len(input_files)}] Failed")

    # Print processing summary
    print(f"\n{'=' * 60}")
    print(f"Processing Summary:")
    print(f"{'=' * 60}")
    print(f"  ✅ Success: {len(processed_files)} files")
    print(f"  ❌ Failed: {len(failed_files)} files")

    if failed_files:
        print(f"\nFailed files:")
        for f in failed_files:
            print(f"  - {f.name}")

    if processed_files:
        print(f"\nSuccessful files:")
        for f in processed_files:
            print(f"  - {f.name}")

    # Generate ASV table
    if processed_files:
        generate_asv_table(output_dir)
    else:
        print("\n⚠️ Warning: No successfully processed files, cannot generate ASV table")

    # Check temporary directory status
    if args.keep_temp:
        temp_dir = output_dir / "temp"
        if temp_dir.exists():
            print(f"\n📁 Temporary files kept in: {temp_dir}")
            subdirs = list(temp_dir.iterdir())
            print(f"  Contains {len(subdirs)} subdirectories")
            for subdir in subdirs:
                if subdir.is_dir():
                    files = list(subdir.iterdir())
                    print(f"    - {subdir.name}: {len(files)} files")

    print("\n🎉 Processing complete!")


if __name__ == "__main__":
    main()