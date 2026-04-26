#!/usr/bin/env python3
"""
Batch processing of fasta files for sequence dereplication and denoising, generating ASV table and unchim
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


def process_single_file(input_file, output_dir, threads=20, minsize=8, keep_temp=False):
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
    print(f"Minsize threshold: {minsize}")
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
        print(f"\nStep 3: Clustering and denoising (minsize={minsize})")
        print(f"Input: {derep_sample_file}")
        print(f"Output: {final_file}")

        cmd3 = f"vsearch --cluster_unoise {derep_sample_file} --sizein --sizeout --centroids {final_file} --threads {threads} --unoise_alpha 2 --minsize {minsize}"
        if not run_command(cmd3, f"Clustering and denoising with minsize={minsize}"):
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
    parser.add_argument('-m', '--minsize', type=int, default=8, help='Minimum abundance for UNOISE clustering (default: 8)')
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
    
    print(f"\n⚙️ Processing parameters:")
    print(f"  - Threads: {args.threads}")
    print(f"  - Minsize: {args.minsize}")
    print(f"  - Keep temporary files: {args.keep_temp}")

    # Process each file
    processed_files = []
    failed_files = []

    for i, input_file in enumerate(input_files, 1):
        print(f"\n[{i}/{len(input_files)}] Starting processing...")
        result = process_single_file(input_file, output_dir, args.threads, args.minsize, args.keep_temp)
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
    main()#!/usr/bin/env python3
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


def process_single_file(input_file, output_dir, threads=20, minsize=8, keep_temp=False):
    """Process a single fasta file"""

    # Get base name (remove path and extension)
    base_name = Path(input_file).stem

    # Create independent temporary directory for each file
    temp_dir = Path(output_dir) / "temp" / base_name
    temp_dir.mkdir(parents=True, exist_ok=True)

    derep_file = temp_dir / f"dereplicated.{base_name}.fa"
    derep_sample_file = temp_dir / f"dereplicated_with_sample.{base_name}.fasta"
    unoise_file = temp_dir / f"unoise.alpha2.{base_name}.fa"
    
    # Chimera detection output files
    chimera_aln_file = temp_dir / f"chimaln.{base_name}.txt"
    chimera_file = temp_dir / f"chimeras.{base_name}.fa"
    chimera_info_file = temp_dir / f"chiminfo.{base_name}.txt"
    nonchimera_file = Path(output_dir) / f"nonchimeras.{base_name}.fa"

    print(f"\n{'=' * 60}")
    print(f"Processing file: {input_file}")
    print(f"Base name: {base_name}")
    print(f"Temporary directory: {temp_dir}")
    print(f"Minsize threshold: {minsize}")
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
        print(f"\nStep 3: Clustering and denoising (minsize={minsize})")
        print(f"Input: {derep_sample_file}")
        print(f"Output: {unoise_file}")

        cmd3 = f"vsearch --cluster_unoise {derep_sample_file} --sizein --sizeout --centroids {unoise_file} --threads {threads} --unoise_alpha 2 --minsize {minsize}"
        if not run_command(cmd3, f"Clustering and denoising with minsize={minsize}"):
            print(f"❌ Error: Clustering and denoising failed - {input_file}")
            return None

        # Check if UNOISE file was generated
        if not unoise_file.exists():
            print(f"❌ Error: UNOISE file not generated - {unoise_file}")
            return None
        else:
            file_size = unoise_file.stat().st_size
            print(f"✓ UNOISE file generated: {unoise_file} (size: {file_size} bytes)")

        # Step 4: Chimera detection using uchime3_denovo
        print(f"\nStep 4: Chimera detection")
        print(f"Input: {unoise_file}")
        print(f"Output (non-chimeras): {nonchimera_file}")
        print(f"Output (chimeras): {chimera_file}")
        print(f"Output (alignments): {chimera_aln_file}")
        print(f"Output (info): {chimera_info_file}")

        cmd4 = (f"vsearch --uchime3_denovo {unoise_file} "
                f"--uchimealns {chimera_aln_file} "
                f"--chimeras {chimera_file} "
                f"--uchimeout {chimera_info_file} "
                f"--fasta_score "
                f"--nonchimeras {nonchimera_file} "
                f"--threads {threads} "
                f"--sizein --sizeout")
        
        if not run_command(cmd4, "Detecting chimeras with uchime3_denovo"):
            print(f"❌ Error: Chimera detection failed - {input_file}")
            return None

        # Check if non-chimera file was generated
        if not nonchimera_file.exists():
            print(f"❌ Error: Non-chimera file not generated - {nonchimera_file}")
            return None
        else:
            file_size = nonchimera_file.stat().st_size
            print(f"✓ Non-chimera file generated: {nonchimera_file} (size: {file_size} bytes)")

        # Count sequences in each file for statistics
        def count_sequences(file_path):
            """Count sequences in fasta file"""
            if not file_path.exists():
                return 0
            count = 0
            with open(file_path, 'r') as f:
                for line in f:
                    if line.startswith('>'):
                        count += 1
            return count

        unoise_count = count_sequences(unoise_file)
        chimera_count = count_sequences(chimera_file) if chimera_file.exists() else 0
        nonchimera_count = count_sequences(nonchimera_file)

        print(f"\n📊 Chimera detection statistics:")
        print(f"  UNOISE sequences: {unoise_count}")
        print(f"  Chimeras detected: {chimera_count}")
        print(f"  Non-chimeras: {nonchimera_count}")
        print(f"  Chimera percentage: {chimera_count / unoise_count * 100:.2f}%" if unoise_count > 0 else "  N/A")

        print(f"\n✅ Successfully completed processing: {nonchimera_file}")

        # List temporary directory contents
        if keep_temp:
            print(f"\nTemporary files kept in: {temp_dir}")
            temp_files = list(temp_dir.iterdir())
            for tf in temp_files:
                print(f"  - {tf.name} ({tf.stat().st_size} bytes)")

        return nonchimera_file

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
                if unoise_file.exists():
                    os.remove(unoise_file)
                    print(f"  Deleted: {unoise_file.name}")
                # Keep chimera detection files in temp if they exist
                # They will be useful for quality control
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

    # Collect all non-chimera files (updated pattern)
    processed_files = list(Path(output_dir).glob("nonchimeras.*.fa"))

    if not processed_files:
        print("❌ Error: No non-chimera files found")
        return

    print(f"Found {len(processed_files)} non-chimera files for ASV table generation")
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


def generate_chimera_summary(output_dir):
    """Generate chimera detection summary report"""
    
    print(f"\n{'=' * 60}")
    print("Generating chimera detection summary...")
    print(f"{'=' * 60}")
    
    temp_dir = Path(output_dir) / "temp"
    if not temp_dir.exists():
        print("⚠️ Warning: Temporary directory not found, cannot generate chimera summary")
        return
    
    summary_file = Path(output_dir) / "chimera_detection_summary.txt"
    
    with open(summary_file, 'w') as f:
        f.write("Chimera Detection Summary\n")
        f.write("=" * 80 + "\n\n")
        
        total_unoise = 0
        total_chimeras = 0
        total_nonchimeras = 0
        
        # Find all sample subdirectories
        sample_dirs = [d for d in temp_dir.iterdir() if d.is_dir()]
        
        for sample_dir in sorted(sample_dirs):
            sample_name = sample_dir.name
            
            unoise_file = sample_dir / f"unoise.alpha2.{sample_name}.fa"
            chimera_file = sample_dir / f"chimeras.{sample_name}.fa"
            
            def count_seqs(file_path):
                if not file_path.exists():
                    return 0
                count = 0
                with open(file_path, 'r') as fh:
                    for line in fh:
                        if line.startswith('>'):
                            count += 1
                return count
            
            unoise_count = count_seqs(unoise_file)
            chimera_count = count_seqs(chimera_file)
            nonchimera_count = unoise_count - chimera_count
            
            total_unoise += unoise_count
            total_chimeras += chimera_count
            total_nonchimeras += nonchimera_count
            
            f.write(f"Sample: {sample_name}\n")
            f.write(f"  UNOISE sequences: {unoise_count}\n")
            f.write(f"  Chimeras detected: {chimera_count}\n")
            f.write(f"  Non-chimeras: {nonchimera_count}\n")
            if unoise_count > 0:
                f.write(f"  Chimera percentage: {chimera_count / unoise_count * 100:.2f}%\n")
            f.write("\n")
        
        f.write("=" * 80 + "\n")
        f.write("Overall Summary\n")
        f.write("=" * 80 + "\n")
        f.write(f"Total samples: {len(sample_dirs)}\n")
        f.write(f"Total UNOISE sequences: {total_unoise}\n")
        f.write(f"Total chimeras detected: {total_chimeras}\n")
        f.write(f"Total non-chimeras: {total_nonchimeras}\n")
        if total_unoise > 0:
            f.write(f"Overall chimera percentage: {total_chimeras / total_unoise * 100:.2f}%\n")
    
    print(f"✓ Chimera detection summary saved: {summary_file}")


def main():
    parser = argparse.ArgumentParser(description='Batch processing of fasta files for sequence dereplication and denoising')
    parser.add_argument('-i', '--input', required=True, help='Input folder path')
    parser.add_argument('-o', '--output', required=True, help='Output folder path')
    parser.add_argument('-t', '--threads', type=int, default=20, help='Number of threads (default: 20)')
    parser.add_argument('-m', '--minsize', type=int, default=8, help='Minimum abundance for UNOISE clustering (default: 8)')
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
    
    print(f"\n⚙️ Processing parameters:")
    print(f"  - Threads: {args.threads}")
    print(f"  - Minsize: {args.minsize}")
    print(f"  - Keep temporary files: {args.keep_temp}")

    # Process each file
    processed_files = []
    failed_files = []

    for i, input_file in enumerate(input_files, 1):
        print(f"\n[{i}/{len(input_files)}] Starting processing...")
        result = process_single_file(input_file, output_dir, args.threads, args.minsize, args.keep_temp)
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
        
        # Generate chimera detection summary
        if args.keep_temp:
            generate_chimera_summary(output_dir)
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