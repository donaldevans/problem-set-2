#! /usrs/ecdnv/bin bash

# Use BEDtools intersect to identify the size of the largest overlap
# between CTCF and H3K4me3 locations.

datasets=$HOME/MOLB7621/data-sets
tfbs_bed=$datasets/bed/encode.tfbs.chr22.bed.gz
histone_bed=$datasets/bed/encode.h3k4me3.hela.chr22.bed.gz
hg19_chr22=$datasets/fasta/hg19.chr22.fa
ctcf=$datasets/bedtools/ctcf.hela.chr22.bg
tss=$datasets/bed/tss.hg19.chr22.bed.gz
gzcat $tfbs_bed | awk '$4 == "CTCF"' > ctcf-peaks.bed
hg19=$datasets/bedtools/hg19.genome

answer_1=$(bedtools intersect -a ctcf-peaks.bed -b $histone_bed -wo \
    | awk '{print $NF}' \
    | sort -nr | head -n1)

echo "answer-1: $answer_1" > answers.yml

# Use BEDtools to calculate the GC content of nucleotides 19,000,000 to
# 19,000,500 on chr22 of hg19 genome build. Report the GC content as a
# fraction (e.g., 0.50).

echo -e "chr22\t19000000\t19000500" > region.bed

answer_2=$(bedtools nuc -fi $hg19_chr22 -bed region.bed \
    | grep -v "^#" \
    | cut -f5) 
echo "answer-2: $answer_2" >> answers.yml

# Use BEDtools to identify the length of the CTCF ChIP-seq peak (i.e.,
# interval) that has the largest mean signal in ctcf.hela.chr22.bg.gz.

answer_3=$(bedtools map -a ctcf-peaks.bed -b $ctcf -c 4 -o mean \
    | sort -k5 \
    | tail -n1 \
    | awk '{print $3 -$2}')
echo "answer-3: $answer_3" >> answers.yml



# Use BEDtools to identify the gene promoter (defined as 1000 bp upstream
# of a TSS) with the highest median signal in ctcf.hela.chr22.bg.gz.
# Report the gene name (e.g., 'ABC123')

answer_4=$(bedtools flank -l 1000 -r 0 -s -i $tss -g $hg19 \
    | bedtools sort -i - \
    | bedtools map -a - -b $ctcf -c 4 -o median \
    | sort -k7n \
    | tail -n1 \
    | cut -f4)
echo "answer-4: $answer_4" >> answers.yml

#bedtools flank -i $tss -g $hg19 -l 1000 -r 0 -s | head

# Use BEDtools to identify the longest interval on chr22 that is not
# covered by genes.hg19.bed.gz. Report the interval like chr1:100-500.

answer_5=$(bedtools intersect -v -a)
echo "answer-5: $answer_5" >> answers.yml


# Use one or more BEDtools that we haven't covered in class. Be creative.
