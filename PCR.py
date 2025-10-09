# Require both upstream and downstream primers to match

#!/usr/bin/env python3

import argparse
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import csv
from tqdm import tqdm
import multiprocessing
import os
import shutil

IUPAC_CODES = {
    'A': ['A'],
    'C': ['C'],
    'G': ['G'],
    'T': ['T'],
    'R': ['A', 'G'],
    'Y': ['C', 'T'],
    'S': ['G', 'C'],
    'W': ['A', 'T'],
    'K': ['G', 'T'],
    'M': ['A', 'C'],
    'B': ['C', 'G', 'T'],
    'D': ['A', 'G', 'T'],
    'H': ['A', 'C', 'T'],
    'V': ['A', 'C', 'G'],
    'N': ['A', 'C', 'G', 'T']
}


def parse_args():
    parser = argparse.ArgumentParser(description="Extract sequences between primers from a FASTA file.")
    parser.add_argument("input_file", help="Input file in FASTA format")
    parser.add_argument("forward_primers", help="Forward primer sequences")
    parser.add_argument("reverse_primers", help="Reverse primer sequences")
    parser.add_argument("output_file", help="Output file in FASTA format")
    parser.add_argument("log_file", help="Tab-separated log file for processing details")
    parser.add_argument("alignment_file", help="File to store alignment results")
    parser.add_argument("--forward_mismatches", type=int, default=3,
                        help="Maximum number of mismatches allowed in forward primer alignment")
    parser.add_argument("--reverse_mismatches", type=int, default=3,
                        help="Maximum number of mismatches allowed in reverse primer alignment")
    parser.add_argument("-n", "--num_processes", type=int, default=4,
                        help="Number of processes to use for parallel processing")
    return parser.parse_args()


def is_match(primer_base, seq_base):
    return seq_base.upper() in IUPAC_CODES.get(primer_base.upper(), [primer_base.upper()])


def find_all_best_matches(sequence, primer, max_mismatches):
    best_score = float('inf')
    best_positions = []
    for i in range(len(sequence) - len(primer) + 1):
        score = sum(0 if is_match(p, s) else 1 for p, s in zip(primer, sequence[i:i + len(primer)]))
        if score <= max_mismatches:
            if score < best_score:
                best_score = score
                best_positions = [i]
            elif score == best_score:
                best_positions.append(i)
    return best_positions, best_score


def visualize_alignment(sequence, primer, position, is_start=True):
    if position == -1:
        return f"Primer: {primer}\nSeq:    {'?' * len(primer)}\n        (Not found or too many mismatches)"

    if is_start and position == 1:
        return f"Primer: {primer}\nSeq:    (Assumed to be before the sequence start)"

    if not is_start and position == len(sequence) + 1:
        return f"Primer: {primer}\nSeq:    (Assumed to be after the sequence end)"

    # Calculate the start and end positions of the display region
    start = max(0, position - 1 - 5)  # Additional 5 bases upstream
    end = min(len(sequence), position - 1 + len(primer) + 5)  # Additional 5 bases downstream

    # Extract the complete display sequence
    full_seq_segment = sequence[start:end]

    # Calculate the offset for primer alignment position
    offset = position - 1 - start

    # Build display lines
    seq_display = full_seq_segment
    match_line = ' ' * offset  # Spaces before alignment

    # Calculate match status
    mismatches = 0
    for i, (p, s) in enumerate(zip(primer, sequence[position - 1:position - 1 + len(primer)])):
        if is_match(p, s):
            match_line += '|'
        else:
            match_line += ' '
            mismatches += 1

    # Add spaces after alignment
    match_line += ' ' * (len(full_seq_segment) - len(match_line))

    # Build primer display line
    primer_line = ' ' * offset + primer + ' ' * (len(full_seq_segment) - len(primer) - offset)

    # Add position marker
    position_info = f"Position: {position}"

    # Build return result
    result = (
        f"{position_info}\n"
        f"Sequence:  {seq_display}\n"
        f"           {match_line}\n"
        f"Primer:    {primer_line}\n"
        f"Mismatches: {mismatches}\n"
    )

    return result


def extract_sequence(record, forward_primers, reverse_primers, forward_mismatches, reverse_mismatches):
    sequence = str(record.seq)
    rev_comp_sequence = str(record.seq.reverse_complement())

    def process_single_sequence(seq, is_original=True):
        # Find all matching positions for forward and reverse primers
        forward_matches = []
        for i, primer in enumerate(forward_primers):
            positions, score = find_all_best_matches(seq, primer, forward_mismatches)
            for pos in positions:
                forward_matches.append((pos + 1, i, score))
        forward_matches.sort()

        reverse_matches = []
        for i, primer in enumerate(reverse_primers):
            positions, score = find_all_best_matches(seq, primer, reverse_mismatches)
            for pos in positions:
                reverse_matches.append((pos + 1, i, score))
        reverse_matches.sort()

        # Generate all possible amplicons
        amplicons = set()

        # Modified: Only when both forward and reverse primers are found, an amplicon is considered to be generated
        if forward_matches and reverse_matches:
            for f_pos, f_index, f_score in forward_matches:
                for r_pos, r_index, r_score in reverse_matches:
                    if r_pos > f_pos:
                        end_pos = r_pos + len(reverse_primers[r_index])
                        amplicons.add((f_pos, end_pos, f_index, r_index, f_score, r_score))

        alignment_results = []
        logs = []
        extracted_records = []

        for i, (start, end, f_index, r_index, f_score, r_score) in enumerate(sorted(amplicons)):
            f_primer = forward_primers[f_index]
            r_primer = reverse_primers[r_index]

            # Add sequence identifier and boundary information
            alignment_results.append(f"\n{'=' * 60}\n")
            alignment_results.append(
                f"Sequence: {record.id}_{'rev_comp_' if not is_original else ''}amplicon_{i + 1}\n")
            alignment_results.append(f"Sequence length: {len(seq)}\n")

            # Add forward primer alignment results
            alignment_results.append("\nForward primer alignment:\n")
            f_alignment = visualize_alignment(seq, f_primer, start)
            alignment_results.append(f"{f_alignment}\n")

            # Add reverse primer alignment results
            alignment_results.append("\nReverse primer alignment:\n")
            r_alignment = visualize_alignment(seq, r_primer, end - len(r_primer))
            alignment_results.append(f"{r_alignment}\n")

            alignment_results.append(f"{'=' * 60}\n")

            extracted_seq = seq[max(0, start - 1):min(len(seq), end - 1)]

            if is_original:
                log = {
                    "Sequence_ID": f"{record.id}_amplicon_{i + 1}",
                    "Original_Length": len(seq),
                    "Forward_Primer_Index": f_index + 1,
                    "Forward_Primer_Position": start,
                    "Forward_Primer_Mismatches": f_score,
                    "Reverse_Primer_Index": r_index + 1,
                    "Reverse_Primer_Position": end - len(r_primer),
                    "Reverse_Primer_Mismatches": r_score,
                    "Extracted_Length": len(extracted_seq),
                    "Is_Reverse_Complement": "No"
                }
            else:
                log = {
                    "Sequence_ID": f"{record.id}_rev_comp_amplicon_{i + 1}",
                    "Original_Length": len(seq),
                    "Forward_Primer_Index": f_index + 1,
                    "Forward_Primer_Position": len(seq) - start + 1,
                    "Forward_Primer_Mismatches": f_score,
                    "Reverse_Primer_Index": r_index + 1,
                    "Reverse_Primer_Position": len(seq) - (end - len(r_primer)) + 1,
                    "Reverse_Primer_Mismatches": r_score,
                    "Extracted_Length": len(extracted_seq),
                    "Is_Reverse_Complement": "Yes"
                }
            logs.append(log)

            if extracted_seq:
                extracted_record = SeqRecord(
                    Seq(extracted_seq),
                    id=f"{record.id}_{'rev_comp_' if not is_original else ''}amplicon_{i + 1}",
                    description=f"F{f_index + 1}_R{r_index + 1}"
                )
                extracted_records.append(extracted_record)

        return extracted_records, logs, alignment_results

    original_records, original_logs, original_alignments = process_single_sequence(sequence, is_original=True)
    rev_comp_records, rev_comp_logs, rev_comp_alignments = process_single_sequence(rev_comp_sequence, is_original=False)

    return (original_records + rev_comp_records,
            original_logs + rev_comp_logs,
            original_alignments + rev_comp_alignments)


def process_file(input_file, forward_primers, reverse_primers, forward_mismatches, reverse_mismatches, progress_queue):
    extracted_records = []
    logs = []
    alignments = []
    for record in SeqIO.parse(input_file, "fasta"):
        record_extracted, record_logs, record_alignments = extract_sequence(record, forward_primers, reverse_primers,
                                                                            forward_mismatches, reverse_mismatches)
        extracted_records.extend(record_extracted)
        logs.extend(record_logs)
        alignments.extend(record_alignments)
        progress_queue.put(1)
    return extracted_records, logs, alignments


def split_fasta(input_file, num_files, tmp_dir):
    print(f"Reading input file: {input_file}")
    records = list(SeqIO.parse(input_file, "fasta"))
    total_sequences = len(records)
    sequences_per_file = total_sequences // num_files

    print(f"Total sequences: {total_sequences}")
    print(f"Splitting into {num_files} files")

    split_files = []
    for i in range(num_files):
        start = i * sequences_per_file
        end = start + sequences_per_file if i < num_files - 1 else total_sequences

        temp_file = os.path.join(tmp_dir, f"temp_{i}.fasta")
        SeqIO.write(records[start:end], temp_file, "fasta")
        split_files.append(temp_file)
        print(f"Created temporary file: {temp_file}")

    return split_files, total_sequences


def main():
    args = parse_args()

    forward_primers = args.forward_primers.split(',')
    reverse_primers = args.reverse_primers.split(',')

    tmp_dir = os.path.join(os.getcwd(), "tmp")
    os.makedirs(tmp_dir, exist_ok=True)
    print(f"Created temporary directory: {tmp_dir}")

    print("Splitting input file...")
    split_files, total_sequences = split_fasta(args.input_file, args.num_processes, tmp_dir)
    print("File splitting complete.")

    print("Starting parallel processing...")
    pool = multiprocessing.Pool(processes=args.num_processes)
    manager = multiprocessing.Manager()
    progress_queue = manager.Queue()

    results = []
    for file in split_files:
        result = pool.apply_async(process_file,
                                  (file, forward_primers, reverse_primers,
                                   args.forward_mismatches, args.reverse_mismatches,
                                   progress_queue))
        results.append(result)

    with tqdm(total=total_sequences, desc="Processing sequences") as pbar:
        processed_sequences = 0
        while processed_sequences < total_sequences:
            processed = progress_queue.get()
            processed_sequences += processed
            pbar.update(processed)

    extracted_records = []
    all_logs = []
    all_alignments = []
    for result in results:
        chunk_records, chunk_logs, chunk_alignments = result.get()
        extracted_records.extend(chunk_records)
        all_logs.extend(chunk_logs)
        all_alignments.extend(chunk_alignments)

    pool.close()
    pool.join()

    print("Parallel processing complete.")
    print("Writing results...")

    with open(args.log_file, 'w', newline='') as log_file:
        fieldnames = ["Original_Sequence_Name", "Sequence_ID", "Original_Length", "Forward_Primer_Index",
                      "Forward_Primer_Position", "Forward_Primer_Mismatches",
                      "Reverse_Primer_Index", "Reverse_Primer_Position",
                      "Reverse_Primer_Mismatches", "Extracted_Length", "Is_Reverse_Complement", "Extracted>0"]
        writer = csv.DictWriter(log_file, fieldnames=fieldnames, delimiter='\t')

        writer.writeheader()
        for log in all_logs:
            original_sequence_name = log["Sequence_ID"].split("_amplicon_")[0].split("_rev_comp")[0]
            log["Original_Sequence_Name"] = original_sequence_name
            log["Extracted>0"] = 'T' if log["Extracted_Length"] > 0 else 'F'
            writer.writerow(log)

    print(f"Processing logs have been written to {args.log_file}")

    filtered_records = [record for record in extracted_records
                        if "N" not in record.description]

    SeqIO.write(filtered_records, args.output_file, "fasta")
    print(f"Extracted {len(filtered_records)} sequences have been written to {args.output_file}")

    with open(args.alignment_file, 'w') as alignment_file:
        for alignment in all_alignments:
            alignment_file.write(alignment)
    print(f"Alignment results have been written to {args.alignment_file}")

    print("Cleaning up temporary files...")
    shutil.rmtree(tmp_dir)
    print(f"Deleted temporary directory: {tmp_dir}")

    print("Process completed successfully.")


if __name__ == "__main__":
    main()
