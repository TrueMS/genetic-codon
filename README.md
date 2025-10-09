# genetic-codon
# This project stores the scripts used in the article "Unexpectedly Diverse Deviations of Stop Codons in the Methyl-Coenzyme M Reductase"

# Usage examples
# ASVtable.py
python /path/to/ASVtable.py -i 1.fasta -o 2.output -t 20 --keep-temp

# PCR.py
python /path/to/PCR.py mcr.nucl.ref.fa GGAACAGATATCGTRTGYGA AACTAYGCHATGAACGTAGG mcr.nucl.amplicon.fa pcr.process alignment_results --forward_mismatches 3 --reverse_mismatches 3 -n 15

# Coding.py
python /path/to/Coding.py -i dereplicated.mcr.pro.DRNA.mafft -o coding.stats
