# genetic-codon
This project stores the scripts used in the article "A hidden layer of stop-codon recoding in the anaerobic methanotrophic archaeon"

# Usage examples

## Perform quality control on the raw amplicons (single-end sequencing) and convert the FASTQ format to FASTA format.

python /data/cailab/script.sh/batch_sickle.py -i 0.data
  · Quality control
mkdir 1.sickle && mv 0.data/*sickle 1.sickle

python /path/to/fq2fa.py -i 1.sickle -title -pfix .fq.sickle
  · Coverting FastQ to Fasta
cp 1.sickle/*fasta 2.fasta/

## mcrASVtable.py 
python /path/to/mcrASVtable.py -i 2.fasta/ -o 3.outputchim -t 20 --keep-temp -m 3
  · Batch processing of fasta files for sequence dereplication and denoising, generating ASV table and unchim
  
## PCR.py 
python /path/to/PCR.py mcr.nucl.ref.fa GGAACAGATATCGTRTGYGA AACTAYGCHATGAACGTAGG mcr.nucl.amplicon.fa pcr.process alignment_results --forward_mismatches 3 --reverse_mismatches 3 -n 15
  · Extract the primer-flanked fragments from the mcr.nucl.ref.fa sequence to mcr.nucl.amplicon.fa.
  
## Coding.py
python /path/to/Coding.py -i dereplicated.mcr.pro.DRNA.mafft -o coding.stats
  · Calculate character frequency at each position in aligned sequences

## stop_codon_usage.py
python /data/cailab/script.sh/stop_codon_usage.py -i mcrA -o stopcodon.ANME2dmcrA
  · Calculate the usage frequency of each stop codon in the sequences.

## AlnView.py
python .\AlnView.py -i .\Ostar.Pro.mafft -diff -n 90
  · Visualize the aligned sequences.
  
## Pro2ASV2Seq
python .\Pro2ASV2Seq.py -table .\ASV.mcr.pro.table.DNA.with_sequences.txt -seq .\DNAsamples.nonchimeras.fasta -o DNA.Pro.table.Seq.txt
  · Translate DNA sequences to protein and count DNA variants for each protein
