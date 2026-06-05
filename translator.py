import argparse
from Bio import SeqIO
import subprocess

def modify_fasta_headers(input_file, output_file):
    intermediate_file = "temp_translated.fasta"
    
    cmd = f"seqkit translate -T 11 -f 6 -j 1 {input_file} > {intermediate_file}"
    subprocess.run(cmd, shell=True)
    
    with open(output_file, 'w') as out:
        for idx, record in enumerate(SeqIO.parse(intermediate_file, 'fasta')):
            header = record.description

            parts = header.split(maxsplit=1)
            
            frame_num = (idx % 6) + 1
            
            if len(parts) > 1:
                new_header = f"{parts[0]}_{frame_num} {parts[1]}"
            else:
                new_header = f"{parts[0]}_{frame_num}"
                
            out.write(f">{new_header}\n{str(record.seq)}\n")
            
    subprocess.run(f"rm {intermediate_file}", shell=True)

def main():
    parser = argparse.ArgumentParser(description='Modify headers of translated sequences')
    parser.add_argument('-i', '--input', required=True, help='Input gene sequences file')
    parser.add_argument('-o', '--output', required=True, help='Output protein sequences file')
    
    args = parser.parse_args()
    
    modify_fasta_headers(args.input, args.output)

if __name__ == "__main__":
    main()