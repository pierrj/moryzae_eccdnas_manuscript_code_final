#!/bin/bash
#SBATCH --job-name=ecc_count_comparisons_final
#SBATCH --partition=savio2
#SBATCH --qos=savio_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=24:00:00
#SBATCH --mail-user=pierrj@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --output=/global/home/users/pierrj/slurm_stdout/slurm-%j.out
#SBATCH --error=/global/home/users/pierrj/slurm_stderr/slurm-%j.out
#MIT License
#
#Copyright (c) 2021 Pierre Michel Joubert
#
#Permission is hereby granted, free of charge, to any person obtaining a copy
#of this software and associated documentation files (the "Software"), to deal
#in the Software without restriction, including without limitation the rights
#to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
#copies of the Software, and to permit persons to whom the Software is
#furnished to do so, subject to the following conditions:
#
#The above copyright notice and this permission notice shall be included in all
#copies or substantial portions of the Software.
#
#THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
#IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
#AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
#OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
#SOFTWARE.

cd /global/scratch/users/pierrj/eccDNA/rerunning_stuff_final/comparative/count_comparison_re

if [ -f "ecc_count_comparison.txt" ]; then
    rm ecc_count_comparison.txt
fi

while read -r ORGANISM_NAME GENE_GFF SAMPLE_MAPFILE ECC_DIR GENOME_FILE COPIA_FILE GYPSY_FILE CONTIGNAMES_FILE; do
    # filter genome to contigs of interest, same as those used for ecc_caller
    samtools faidx ${GENOME_FILE} $(cat ${CONTIGNAMES_FILE} | tr "\n" " ") > ${GENOME_FILE}.filtered
    GENOME_FILE=$(realpath ${GENOME_FILE}.filtered)
    samtools faidx ${GENOME_FILE}
    cut -f1,2 ${GENOME_FILE}.fai > ${GENOME_FILE}.chromsizes
    GENOME_CHROMSIZES=$(realpath ${GENOME_FILE}.chromsizes)
    awk -v OFS='\t' '{print $1, 1, $2}' ${GENOME_CHROMSIZES} > ${GENOME_CHROMSIZES}.forintersect
    GENOME_CHROMSIZES_FOR_INTERSECT=$(realpath ${GENOME_CHROMSIZES}.forintersect)
    # filter gene locations, copia locations and gypsy locations for only contigs of interest
    bedtools intersect -header -u -a ${GENE_GFF} -b ${GENOME_CHROMSIZES_FOR_INTERSECT} > ${GENE_GFF}.filtered
    GENE_GFF=$(realpath ${GENE_GFF}.filtered)
    bedtools intersect -header -u -a ${COPIA_FILE} -b ${GENOME_CHROMSIZES_FOR_INTERSECT} > ${COPIA_FILE}.filtered
    COPIA_FILE=$(realpath ${COPIA_FILE}.filtered)
    bedtools intersect -header -u -a ${GYPSY_FILE} -b ${GENOME_CHROMSIZES_FOR_INTERSECT} > ${GYPSY_FILE}.filtered
    GYPSY_FILE=$(realpath ${GYPSY_FILE}.filtered)
    cd $ECC_DIR
    while read sample;
    do
        cd ${sample}
        echo $sample 1>&2
        # deal with overlaps in annotation
        awk '{if ($3 == "exon") print $0}' ${GENE_GFF} | sort -k1,1 -k4,4n | bedtools merge -i - > tmp.exons.bed
        sort -k1,1 -k2,2n ${COPIA_FILE} | bedtools merge -i - > tmp.copia.bed
        sort -k1,1 -k2,2n ${GYPSY_FILE} | bedtools merge -i - > tmp.gypsy.bed
        # no lowq eccdnas
        awk ' $5 != "lowq" ' ${sample}.ecc_caller_out.details.txt > ${sample}.ecc_caller_out.details.nolowq.txt
        # if 50 percent of sequence is made up of exons it is a genic eccdna otherwise it is noncoding
        # middle awk line deals with multiple hits from bedtools intersect per eccdna
        bedtools intersect -wao -a ${sample}.ecc_caller_out.details.nolowq.txt -b tmp.exons.bed | \
                    awk -v OFS='\t' '{a[$1"\t"$2"\t"$3] += $(NF)} END{for (i in a) print i, a[i]}' | \
                    awk -v OFS='\t' '{ if ($4/($3-$2) >= 0.5) {print $0}}' | sort -k1,1 -k2,2n > tmp.genic_eccs
        bedtools intersect -wao -a ${sample}.ecc_caller_out.details.nolowq.txt -b tmp.exons.bed | \
                    awk -v OFS='\t' '{a[$1"\t"$2"\t"$3] += $(NF)} END{for (i in a) print i, a[i]}' | \
                    awk -v OFS='\t' '{ if ($4/($3-$2) < 0.5) {print $0}}' | sort -k1,1 -k2,2n > tmp.noncoding_eccs
        # if 90 percent of sequence is made up of ltr sequence that is an ltr eccdna
        # middle awk line deals with multiple hits from bedtools intersect per eccdna
        bedtools intersect -wao -a ${sample}.ecc_caller_out.details.nolowq.txt -b tmp.copia.bed | \
            awk -v OFS='\t' '{a[$1"\t"$2"\t"$3] += $(NF)} END{for (i in a) print i, a[i]}' | \
            awk -v OFS='\t' '{ if ($4/($3-$2) >= 0.9) {print $0}}' | sort -k1,1 -k2,2n > tmp.copia_eccs
        bedtools intersect -wao -a ${sample}.ecc_caller_out.details.nolowq.txt -b tmp.gypsy.bed | \
            awk -v OFS='\t' '{a[$1"\t"$2"\t"$3] += $(NF)} END{for (i in a) print i, a[i]}' | \
            awk -v OFS='\t' '{ if ($4/($3-$2) >= 0.9) {print $0}}' | sort -k1,1 -k2,2n > tmp.gypsy_eccs
        # if not ltr-eccdnas, sort by size
        bedtools intersect -wao -a ${sample}.ecc_caller_out.details.nolowq.txt -b tmp.copia.bed | \
            awk -v OFS='\t' '{a[$1"\t"$2"\t"$3] += $(NF)} END{for (i in a) print i, a[i]}' | \
            awk -v OFS='\t' '{ if ($4/($3-$2) < 0.9) {print $0}}' | \
            bedtools intersect -wao -a - -b tmp.gypsy.bed | \
            awk -v OFS='\t' '{a[$1"\t"$2"\t"$3] += $(NF)} END{for (i in a) print i, a[i]}' | \
            awk -v OFS='\t' '{ if ($4/($3-$2) < 0.9) {print $0}}'| sort -k1,1 -k2,2n | awk '$3-$2 <= 400' > tmp.micro_dnas
        bedtools intersect -wao -a ${sample}.ecc_caller_out.details.nolowq.txt -b tmp.copia.bed | \
            awk -v OFS='\t' '{a[$1"\t"$2"\t"$3] += $(NF)} END{for (i in a) print i, a[i]}' | \
            awk -v OFS='\t' '{ if ($4/($3-$2) < 0.9) {print $0}}' | \
            bedtools intersect -wao -a - -b tmp.gypsy.bed | \
            awk -v OFS='\t' '{a[$1"\t"$2"\t"$3] += $(NF)} END{for (i in a) print i, a[i]}' | \
            awk -v OFS='\t' '{ if ($4/($3-$2) < 0.9) {print $0}}'| sort -k1,1 -k2,2n | awk '$3-$2 > 400' > tmp.large_eccdnas
        # count each category
        ecc_count_genic=$(wc -l tmp.genic_eccs | awk '{print $1}')
        ecc_count_noncoding=$(wc -l tmp.noncoding_eccs | awk '{print $1}')
        ecc_count_copia=$(wc -l tmp.copia_eccs | awk '{print $1}')
        ecc_count_gypsy=$(wc -l tmp.gypsy_eccs | awk '{print $1}')
        ecc_count_micro=$(wc -l tmp.micro_dnas | awk '{print $1}')
        ecc_count_large=$(wc -l tmp.large_eccdnas | awk '{print $1}')
        ecc_count_all=$(wc -l ${sample}.ecc_caller_out.details.nolowq.txt | awk '{print $1}')
        ## this is in number of mapped reads
        read_count=$(samtools view -c -F 4 -F 2048 no_secondary.filtered.sorted.${sample}.bam | awk '{print $1}') ## mapped reads
        genome_size=$(awk 'BEGIN{sum=0;}{if($0 !~ /^>/){sum+=length($0);}}END{print sum;}' ${GENOME_FILE} ) ## basepairs
        coding_bp=$(awk '{SUM+=$3-$2}END{print SUM}' tmp.exons.bed)
        copia_bp=$(awk '{SUM+=$3-$2}END{print SUM}' tmp.copia.bed)
        gypsy_bp=$(awk '{SUM+=$3-$2}END{print SUM}' tmp.gypsy.bed)
        # put it all in a table
        echo -e $sample'\t'$ecc_count_genic'\t'$ecc_count_noncoding'\t'$ecc_count_copia\
            '\t'$ecc_count_gypsy'\t'$ecc_count_micro'\t'$ecc_count_large'\t'$ecc_count_all\
            '\t'$read_count'\t'$genome_size\
            '\t'$coding_bp'\t'$copia_bp'\t'$gypsy_bp >> /global/scratch/users/pierrj/eccDNA/rerunning_stuff_final/comparative/count_comparison_re/ecc_count_comparison.txt
        rm tmp.exons.bed
        rm tmp.genic_eccs
        rm tmp.noncoding_eccs
        rm tmp.copia_eccs
        rm tmp.gypsy_eccs
        rm tmp.micro_dnas
        rm tmp.large_eccdnas
        cd ..
    done < $SAMPLE_MAPFILE
done < mapfile_ecc_count

cd /global/scratch/users/pierrj/eccDNA/rerunning_stuff_final/comparative/count_comparison_re

paste ecc_count_comparison.txt ecc_count_metadata.txt > ecc_count_comparison_w_metadata.txt