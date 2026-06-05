# genetic-codon
This project stores the scripts used in the article "A hidden layer of stop-codon recoding in the anaerobic methanotrophic archaeon"

# Usage examples

## Perform quality control on the raw amplicons (single-end sequencing) and convert the FASTQ format to FASTA format.

  * Quality control
    
  `python /data/cailab/script.sh/batch_sickle.py -i 0.data`

 * Coverting FastQ to Fasta
   
  `mkdir 1.sickle && mv 0.data/*sickle 1.sickle`

`python /path/to/fq2fa.py -i 1.sickle -title -pfix .fq.sickle`
 
`cp 1.sickle/*fasta 2.fasta/`

## mcrASVtable.py 

* Batch processing of fasta files for sequence dereplication and denoising, generating ASV table and unchim

`python /path/to/mcrASVtable.py -i 2.fasta/ -o 3.outputchim -t 20 --keep-temp -m 3`


## Protable.py

* Merge identical protein sequences and sum their counts

`python /path/to/Protable.py -pro ASV.mcr.pro.forabun -t ASV.mcr.table.txt -o ASV.mcr.pro.table.txt`

## PCR.py 

 * Extract the primer-flanked fragments from the mcr.nucl.ref.fa sequence to mcr.nucl.amplicon.fa

`python /path/to/PCR.py mcr.nucl.ref.fa GGAACAGATATCGTRTGYGA AACTAYGCHATGAACGTAGG mcr.nucl.amplicon.fa pcr.process alignment_results --forward_mismatches 3 --reverse_mismatches 3 -n 15`
 
## Coding.py

  * Calculate character frequency at each position in aligned sequences
    
`python /path/to/Coding.py -i dereplicated.mcr.pro.DRNA.mafft -o coding.stats`

## stop_codon_usage.py

* Calculate the usage frequency of each stop codon in the sequences
  
`python /path/to/stop_codon_usage.py -i mcrA -o stopcodon.ANME2dmcrA`

## AlnView.py

python /path/to/AlnView.py -i Ostar.Pro.mafft -diff -n 90

  · Visualize the aligned sequences.
  
## Pro2ASV2Seq.py

python /path/to/Pro2ASV2Seq.py -table ASV.mcr.pro.table.DNA.with_sequences.txt -seq DNAsamples.nonchimeras.fasta -o DNA.Pro.table.Seq.txt

  · Translate DNA sequences to protein and count DNA variants for each protein
