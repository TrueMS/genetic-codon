# genetic-codon
This project stores the scripts used in the article "A hidden layer of stop-codon recoding in the anaerobic methanotrophic archaeon"

# Usage examples

## Perform quality control on the raw amplicons (single-end sequencing) and convert the FASTQ format to FASTA format.

  * Quality control
    
```
python /data/cailab/script.sh/batch_sickle.py -i 0.data
```

 * Coverting FastQ to Fasta
   
```
mkdir 1.sickle && mv 0.data/*sickle 1.sickle
python /path/to/fq2fa.py -i 1.sickle -title -pfix .fq.sickle
cp 1.sickle/*fasta 2.fasta/
```

## mcrASVtable.py 

* Batch processing of fasta files for sequence dereplication and denoising, generating ASV table and unchim

```
python /path/to/mcrASVtable.py -i 2.fasta/ -o 3.outputchim -t 20 --keep-temp -m 3
```

* Directory structure of the output-folder `3.outputchim/`
```
3.outputchim/
├── nonchimeras.SampleA.fa      # Sample A: denoised & chimera-filtered per-sample ASV sequences (main output)
├── nonchimeras.SampleB.fa      # Sample B: same as above
├── ASV.seq.fa                  # [Core Result 1] Consolidated representative FASTA of all ASVs across samples
├── ASV.table.txt               # [Core Result 2] ASV-by-sample abundance matrix (tab-delimited)
├── chimera_detection_summary.txt # [Only generated with --keep-temp] Global chimera statistics summary across all samples
└── temp/                       # [Only generated with --keep-temp] Intermediate working directory, separated by individual sample subfolders
    ├── SampleA/
    │   ├── dereplicated.SampleA.fa         # Raw dereplicated sequences
    │   ├── dereplicated_with_sample.SampleA.fasta # Dereplicated sequences tagged with sample ID in FASTA headers
    │   ├── unoise.alpha2.SampleA.fa        # UNOISE-denoised sequences prior to chimera removal
    │   ├── chimaln.SampleA.txt             # UCHIME chimera alignment details
    │   ├── chimeras.SampleA.fa             # Detected putative chimeric sequences
    │   └── chiminfo.SampleA.txt            # Chimera scoring and parent-source trace information
    └── SampleB/...
 ```

### translator.py

* Translate the nucleotide sequences in ASV.seq.fa into protein sequences across all six reading frames

```
python /path/to/translator.py -i ASV.seq.fa -o ASV.pro
```

### Identity McrA amplicons from the six sets of translation products

```
# The overall goal of the subsequent steps is to screen mcrA sequences from the amplified ASV.seq.fa and compile them into ASV.mcr.seq.fa.
# Specifically, all translated sequences are scanned against the KEGG K00399.hmm profile; hits with an E-value ≤ 100 are defined as McrA.

/data/cailab/database.HQ/kofam/kofamscan/bin/kofam_scan-1.3.0/exec_annotation -o kofam.ASV.mcr.out -c /data/cailab/database.HQ/kofam/kofamscan/bin/config.mcrA.yml --cpu 15 -E 1e-5 ASV.pro
python /data/cailab/script.sh/kofam.evalue.py -i kofam.ASV.mcr.out -e 1e-100 -n 1
cut -f1 kofam.ASV.mcr_1e-100_top1.out | sed 's/_[^_]*$//' | sort -u > ASV.mcr.seqname
seqkit grep -f ASV.mcr.seqname ASV.seq.fa > ASV.mcr.seq.fa
```

## Protable.py

* Merge identical protein sequences and sum their counts
  
```
python /path/to/Protable.py -pro ASV.mcr.pro.forabun -t ASV.mcr.table.txt -o ASV.mcr.pro.table.txt
```

## PCR.py 

 * Extract the primer-flanked fragments from the mcr.nucl.ref.fa sequence to mcr.nucl.amplicon.fa
   
```
python /path/to/PCR.py mcr.nucl.ref.fa GGAACAGATATCGTRTGYGA AACTAYGCHATGAACGTAGG mcr.nucl.amplicon.fa pcr.process alignment_results --forward_mismatches 3 --reverse_mismatches 3 -n 15
```
 
## Coding.py

  * Calculate character frequency at each position in aligned sequences
    
```
python /path/to/Coding.py -i dereplicated.mcr.pro.DRNA.mafft -o coding.stats
```

## stop_codon_usage.py

* Calculate the usage frequency of each stop codon in the sequences
  
```
python /path/to/stop_codon_usage.py -i mcrA -o stopcodon.ANME2dmcrA
```

## AlnView.py

* Visualize the aligned sequences.

```
python /path/to/AlnView.py -i Ostar.Pro.mafft -diff -n 90
```

## Pro2ASV2Seq.py

 * Translate DNA sequences to protein and count DNA variants for each protein

```
python /path/to/Pro2ASV2Seq.py -table ASV.mcr.pro.table.DNA.with_sequences.txt -seq DNAsamples.nonchimeras.fasta -o DNA.Pro.table.Seq.txt
```

