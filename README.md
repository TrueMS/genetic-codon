# genetic-codon
This project stores the scripts used in the article "A hidden layer of stop-codon recoding in the anaerobic methanotrophic archaeon"

# Usage examples

python /data/cailab/script.sh/batch_sickle.py -i 0.data

mkdir 1.sickle && mv 0.data/*sickle 1.sickle

python /data/cailab/script.sh/fq2fa.py -i 1.sickle -title -pfix .fq.sickle

cp 1.sickle/*fasta 2.fasta/

## mcrASVtable.py
python /path/to/mcrASVtable.py -i 2.fasta/ -o 3.outputchim -t 20 --keep-temp -m 3

## PCR.py
python /path/to/PCR.py mcr.nucl.ref.fa GGAACAGATATCGTRTGYGA AACTAYGCHATGAACGTAGG mcr.nucl.amplicon.fa pcr.process alignment_results --forward_mismatches 3 --reverse_mismatches 3 -n 15

## Coding.py
python /path/to/Coding.py -i dereplicated.mcr.pro.DRNA.mafft -o coding.stats
