# genetic-codon
This project stores the scripts used in the article "A hidden layer of stop-codon recoding in the anaerobic methanotrophic archaeon"

Please email heq22@mail.ustc.edu.cn for any questions about this project.

# Usage examples

## 1. Perform quality control on the raw amplicons (single-end sequencing) and convert the FASTQ format to FASTA format.

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

## 2. Generate mcrA ASV table

### 2.1 mcrASVtable.py 

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

### 2.2 translator.py

* Translate the nucleotide sequences in ASV.seq.fa into protein sequences across all six reading frames

```
python /path/to/translator.py -i ASV.seq.fa -o ASV.pro
```

### 2.3 Identity McrA amplicons from the six sets of translation products

```
# The overall goal of the subsequent steps is to screen mcrA sequences from the amplified ASV.seq.fa and compile them into ASV.mcr.seq.fa.
# Specifically, all translated sequences are scanned against the KEGG K00399.hmm profile; hits with an E-value ≤ 100 are defined as McrA.

/data/cailab/database.HQ/kofam/kofamscan/bin/kofam_scan-1.3.0/exec_annotation -o kofam.ASV.mcr.out -c /data/cailab/database.HQ/kofam/kofamscan/bin/config.mcrA.yml --cpu 15 -E 1e-5 ASV.pro
python /path/to/kofam.evalue.py -i kofam.ASV.mcr.out -e 1e-100 -n 1
cut -f1 kofam.ASV.mcr_1e-100_top1.out | sed 's/_[^_]*$//' | sort -u > ASV.mcr.seqname
seqkit grep -f ASV.mcr.seqname ASV.seq.fa > ASV.mcr.seq.fa
```

### 2.4 Generate an ASV table retaining only mcrA sequences

```
awk 'NR==FNR{if(/^>/){gsub(/^>/,"");ids[$0]=1}}NR!=FNR{if(FNR==1||$1 in ids)print}' ASV.mcr.seq.fa ASV.table.txt > ASV.mcr.table.txt
```

## 3. Extract McrA sequences containing stop codons (*, star)

```
cut -f1 kofam.ASV.mcr_1e-100_top1.out > ASV.pro.seqname
seqkit grep -f ASV.pro.seqname ASV.pro > ASV.mcr.pro
```

```
awk 'BEGIN{seq="";header=""}$0~/^>/{if(length(seq)>0){if(seq~/\*/){print header>"mcr.pro.star.fa";print seq>"mcr.pro.star.fa"}else{print header>"mcr.pro.nostar.fa";print seq>"mcr.pro.nostar.fa"}}header=$0;seq="";next}{seq=seq$0}END{if(length(seq)>0){if(seq~/\*/){print header>"mcr.pro.star.fa";print seq>"mcr.pro.star.fa"}else{print header>"mcr.pro.nostar.fa";print seq>"mcr.pro.nostar.fa"}}}' ASV.mcr.pro

# We observed multiple * residues within several McrA protein sequences.
# Further inspection revealed their amplicons carry a single-base deletion relative to canonical mcrA amplicons, which causes a frameshift and erroneous downstream translation.
# This phenomenon of single-base deletion is not addressed in the present study.
# We extract sequences with a single * from mcr.pro.star.fa to ASV.mcr.pro.ostar

awk 'BEGIN{s="";h=""}$0~/^>/{if(length(s)>0){n=gsub(/\*/,"",s);if(n==1){print h>"ASV.mcr.pro.ostar";print s>"ASV.mcr.pro.ostar"}}h=$0;s="";next}{s=s$0}END{if(length(s)>0){n=gsub(/\*/,"",s);if(n==1){print h>"ASV.mcr.pro.ostar";print s>"ASV.mcr.pro.ostar"}}}' mcr.pro.star.fa

seqkit rmdup -s ASV.mcr.pro.ostar > ASV.mcr.pro.ostar.dmp # Keep one representative for identical protein sequences after degeneracy collapsing
```


## 4. Protable.py

* Merge identical protein sequences and sum their counts
  
```
cp ASV.mcr.pro ASV.mcr.pro.forabun
sed -i '/^>/s/_[^_]*$//' ASV.mcr.pro.forabun
seqkit rmdup -s ASV.mcr.pro.forabun > dereplicated.ASV.mcr.pro.forabun # Keep one representative for identical protein sequences after degeneracy collapsing
python /path/to/Protable.py -pro ASV.mcr.pro.forabun -t ASV.mcr.table.txt -o ASV.mcr.pro.table.txt # The Protable.py can also rmdup
```

## 5. PCR.py 

 * Extract the primer-flanked fragments from the mcr.nucl.ref.fa sequence to mcr.nucl.amplicon.fa
   
```
# mcr.nucl.ref.fa contains mcrA genes retrived from public available genomes
python /path/to/PCR.py mcr.nucl.ref.fa GGAACAGATATCGTRTGYGA AACTAYGCHATGAACGTAGG mcr.nucl.amplicon.fa pcr.process alignment_results --forward_mismatches 3 --reverse_mismatches 3 -n 15

# Sequences from mcr.nucl.amplicon.fa were translated and combined with stop-codon-containing McrA sequences identified in this study for phylogenetic tree construction.
# KnownAndOstar.trimal is the trimmed alignment file and KnownAndOstar.trimal.treefile is the resulting tree file, both stored under the Fig.1 directory of this project.
```
 
## 6. Coding predictions

  * Calculate character frequency at each position in aligned sequences
    
```
# Position_number in the output file is the index within the aligned protein sequence
# clade1.pro.mafft and clade2.pro.mafft are stored under the Coding directory of this project
python /path/to/Coding.py -i clade1.pro.mafft -o clade1.coding.stats
```

* Predict the reassignment for stop codon

```
# Protein Position in the output file is the index within the unaligned protein sequence. So it may be small than Coding.py output
# Ostar.nucl.fa and clade1.pro are stored under the Coding directory of this project
python .\stop.codon.py -g .\Ostar.nucl.fa -p .\clade1.pro -pro .\clade1.pro.mafft -o reassignment.clade1.out
```

## 7. stop_codon_usage.py

* Calculate the usage frequency of each stop codon in the sequences
  
```
python /path/to/stop_codon_usage.py -i mcrA -o stopcodon.ANME2dmcrA
```

## 8. AlnView.py

* Visualize the aligned sequences.

```
# Ostar.Pro.mafft is stored under the AlnView directory of this project
python /path/to/AlnView.py -i Ostar.Pro.mafft -diff -n 90
```

## 9. Pro2ASV2Seq.py

 * Translate DNA sequences to protein and count DNA variants for each protein

```
# ASV.mcr.pro.table.DNA.with_sequences.txt and DNAsamples.nonchimeras.fasta are stored under the CountVarients directory of this project
python /path/to/Pro2ASV2Seq.py -table ASV.mcr.pro.table.DNA.with_sequences.txt -seq DNAsamples.nonchimeras.fasta -o DNA.Pro.table.Seq.txt
```

