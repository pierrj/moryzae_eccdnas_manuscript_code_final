#!/bin/bash
#SBATCH --job-name=effector_on_eccdna_figure
#SBATCH --account=fc_kvkallow
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


cd /global/scratch/users/pierrj/eccDNA/magnaporthe_pureculture/illumina

# looking at avrpita3
effector_name=AvrPita3
effector_chrom=MQOP01000016.1
effector_start=368546
effector_end=370280
effector_gene_name=gene_11229 ## no t1 please

## get gff

GFF_FILE=/global/scratch/users/pierrj/references/guy11_fungap_out_12_28_20.gff3

# get location of effector from gff
awk -v g=$effector_gene_name '$9 ~ g' ${GFF_FILE} > /global/scratch/users/pierrj/eccDNA/pipeline_tests/gviz_effector_ecc/${effector_name}.gff3

## get sample of interest
while read sample; do
cd ${sample}
    echo ${sample}
    awk -v c=$effector_chrom -v s=$effector_start -v e=$effector_end '$1 == c && $2 < s && $3 > e' ${sample}.ecc_caller_out.details.nolowq.txt
cd ..
done < mapfile

## pick sample_target manually
sample_target=G3_1A

# pick a region bordering effectors
region_start=$(($effector_start-4000))
region_end=$(($effector_end+4000))

## get region mapped reads

cd ${sample_target}

samtools index no_secondary.filtered.sorted.${sample_target}.bam

## get total reads in region, sort and index
samtools view -bh no_secondary.filtered.sorted.${sample_target}.bam "${effector_chrom}:${region_start}-${region_end}" > \
    /global/scratch/users/pierrj/eccDNA/pipeline_tests/gviz_effector_ecc/${effector_name}.${sample_target}.totalreads.bam

samtools sort /global/scratch/users/pierrj/eccDNA/pipeline_tests/gviz_effector_ecc/${effector_name}.${sample_target}.totalreads.bam > \
    /global/scratch/users/pierrj/eccDNA/pipeline_tests/gviz_effector_ecc/${effector_name}.${sample_target}.totalreads.sorted.bam

samtools index /global/scratch/users/pierrj/eccDNA/pipeline_tests/gviz_effector_ecc/${effector_name}.${sample_target}.totalreads.sorted.bam

## get eccs in region

awk -v c=$effector_chrom -v s=$region_start -v e=$region_end '$1 == c && $2 > s && $3 < e' ${sample_target}.ecc_caller_out.details.nolowq.txt > \
    /global/scratch/users/pierrj/eccDNA/pipeline_tests/gviz_effector_ecc/${effector_name}.${sample_target}.eccs.bed


cd /global/scratch/users/pierrj/eccDNA/pipeline_tests/gviz_effector_ecc/

SAMPLE=${effector_name}.${sample_target}
FILTERED_BAMFILE=${effector_name}.${sample_target}.totalreads.bam
GENOME_CHROMSIZES=/global/scratch/users/pierrj/references/guy11_genome_baoetal2017.chromsizes



## get junction split reads, taken directly from ecc_caller_mapq0.sh

samtools view -f 81 -F 4 ${FILTERED_BAMFILE} > tmp.reverseread1.${SAMPLE}.sam
splitread_file="reverseread1.${SAMPLE}.sam"
awk -v OFS='\t' '{a=gensub(/^([0-9]+)M.*[HS]$/, "\\1", "", $6); b=gensub(/.*[HS]([0-9]+)M$/, "\\1", "", $6); if (a !~ /[DMIHS]/ && int(a) > 19 ) print $0, 1; else if (b !~ /[DMIHS]/ && int(b) > 19) print $0, 2}' tmp.${splitread_file} > tmp.qualityfiltered.${splitread_file}
awk 'NR==FNR{a[$1, $3]++; next} a[$1, $3]==2' tmp.qualityfiltered.${splitread_file} tmp.qualityfiltered.${splitread_file} | sort -k1,1 -k18,18n > tmp.samechromosome.exactlytwice.qualityfiltered.${splitread_file}
awk -v OFS='\t' '{
    prev=$0; f4=$4; f1=$1
    getline 
    if ($1 == f1 && f4 > $4) {
        print prev
        print $0
    }
}' tmp.samechromosome.exactlytwice.qualityfiltered.${splitread_file} | awk -v OFS='\t' '{$NF=""; print $0}' | sed 's/[ \t]\+$//' > tmp.oriented.samechromosome.exactlytwice.qualityfiltered.${splitread_file}

samtools view -f 145 -F 4 ${FILTERED_BAMFILE} > tmp.reverseread2.${SAMPLE}.sam
splitread_file="reverseread2.${SAMPLE}.sam"
awk -v OFS='\t' '{a=gensub(/^([0-9]+)M.*[HS]$/, "\\1", "", $6); b=gensub(/.*[HS]([0-9]+)M$/, "\\1", "", $6); if (a !~ /[DMIHS]/ && int(a) > 19 ) print $0, 1; else if (b !~ /[DMIHS]/ && int(b) > 19) print $0, 2}' tmp.${splitread_file} > tmp.qualityfiltered.${splitread_file}
awk 'NR==FNR{a[$1, $3]++; next} a[$1, $3]==2' tmp.qualityfiltered.${splitread_file} tmp.qualityfiltered.${splitread_file} | sort -k1,1 -k18,18n > tmp.samechromosome.exactlytwice.qualityfiltered.${splitread_file}
awk -v OFS='\t' '{
    prev=$0; f4=$4; f1=$1
    getline 
    if ($1 == f1 && f4 > $4) {
        print prev
        print $0
    }
}' tmp.samechromosome.exactlytwice.qualityfiltered.${splitread_file} | awk -v OFS='\t' '{$NF=""; print $0}' | sed 's/[ \t]\+$//' > tmp.oriented.samechromosome.exactlytwice.qualityfiltered.${splitread_file}

samtools view -f 65 -F 20 ${FILTERED_BAMFILE} > tmp.forwardread1.${SAMPLE}.sam
splitread_file="forwardread1.${SAMPLE}.sam"
awk -v OFS='\t' '{a=gensub(/^([0-9]+)M.*[HS]$/, "\\1", "", $6); b=gensub(/.*[HS]([0-9]+)M$/, "\\1", "", $6); if (a !~ /[DMIHS]/ && int(a) > 19 ) print $0, 1; else if (b !~ /[DMIHS]/ && int(b) > 19) print $0, 2}' tmp.${splitread_file} > tmp.qualityfiltered.${splitread_file}
awk 'NR==FNR{a[$1, $3]++; next} a[$1, $3]==2' tmp.qualityfiltered.${splitread_file} tmp.qualityfiltered.${splitread_file} | sort -k1,1 -k18,18n > tmp.samechromosome.exactlytwice.qualityfiltered.${splitread_file}
awk -v OFS='\t' '{
    prev=$0; f4=$4; f1=$1
    getline 
    if ($1 == f1 && f4 > $4) {
        print prev
        print $0
    }
}' tmp.samechromosome.exactlytwice.qualityfiltered.${splitread_file} | awk -v OFS='\t' '{$NF=""; print $0}' | sed 's/[ \t]\+$//' > tmp.oriented.samechromosome.exactlytwice.qualityfiltered.${splitread_file}

samtools view -f 129 -F 20 ${FILTERED_BAMFILE} > tmp.forwardread2.${SAMPLE}.sam
splitread_file="forwardread2.${SAMPLE}.sam"
awk -v OFS='\t' '{a=gensub(/^([0-9]+)M.*[HS]$/, "\\1", "", $6); b=gensub(/.*[HS]([0-9]+)M$/, "\\1", "", $6); if (a !~ /[DMIHS]/ && int(a) > 19 ) print $0, 1; else if (b !~ /[DMIHS]/ && int(b) > 19) print $0, 2}' tmp.${splitread_file} > tmp.qualityfiltered.${splitread_file}
awk 'NR==FNR{a[$1, $3]++; next} a[$1, $3]==2' tmp.qualityfiltered.${splitread_file} tmp.qualityfiltered.${splitread_file} | sort -k1,1 -k18,18n > tmp.samechromosome.exactlytwice.qualityfiltered.${splitread_file}
awk -v OFS='\t' '{
    prev=$0; f4=$4; f1=$1
    getline 
    if ($1 == f1 && f4 > $4) {
        print prev
        print $0
    }
}' tmp.samechromosome.exactlytwice.qualityfiltered.${splitread_file} | awk -v OFS='\t' '{$NF=""; print $0}' | sed 's/[ \t]\+$//' > tmp.oriented.samechromosome.exactlytwice.qualityfiltered.${splitread_file}

# I did some merging of reads in the past but this ended up being detrimental
# these next two chunks should be removed eventually
# currently both of these should contain 0 reads
samtools view -f 16 -F 5 ${FILTERED_BAMFILE} > tmp.reversemerged.${SAMPLE}.sam
splitread_file="reversemerged.${SAMPLE}.sam"
awk -v OFS='\t' '{a=gensub(/^([0-9]+)M.*[HS]$/, "\\1", "", $6); b=gensub(/.*[HS]([0-9]+)M$/, "\\1", "", $6); if (a !~ /[DMIHS]/ && int(a) > 19 ) print $0, 1; else if (b !~ /[DMIHS]/ && int(b) > 19) print $0, 2}' tmp.${splitread_file} > tmp.qualityfiltered.${splitread_file}
awk 'NR==FNR{a[$1, $3]++; next} a[$1, $3]==2' tmp.qualityfiltered.${splitread_file} tmp.qualityfiltered.${splitread_file} | sort -k1,1 -k18,18n > tmp.samechromosome.exactlytwice.qualityfiltered.${splitread_file}
awk -v OFS='\t' '{
    prev=$0; f4=$4; f1=$1
    getline 
    if ($1 == f1 && f4 > $4) {
        print prev
        print $0
    }
}' tmp.samechromosome.exactlytwice.qualityfiltered.${splitread_file} | awk -v OFS='\t' '{$NF=""; print $0}' | sed 's/[ \t]\+$//' > tmp.oriented.samechromosome.exactlytwice.qualityfiltered.${splitread_file}

samtools view -F 21 ${FILTERED_BAMFILE} > tmp.forwardmerged.${SAMPLE}.sam
splitread_file="forwardmerged.${SAMPLE}.sam"
awk -v OFS='\t' '{a=gensub(/^([0-9]+)M.*[HS]$/, "\\1", "", $6); b=gensub(/.*[HS]([0-9]+)M$/, "\\1", "", $6); if (a !~ /[DMIHS]/ && int(a) > 19 ) print $0, 1; else if (b !~ /[DMIHS]/ && int(b) > 19) print $0, 2}' tmp.${splitread_file} > tmp.qualityfiltered.${splitread_file}
awk 'NR==FNR{a[$1, $3]++; next} a[$1, $3]==2' tmp.qualityfiltered.${splitread_file} tmp.qualityfiltered.${splitread_file} | sort -k1,1 -k18,18n > tmp.samechromosome.exactlytwice.qualityfiltered.${splitread_file}
awk -v OFS='\t' '{
    prev=$0; f4=$4; f1=$1
    getline 
    if ($1 == f1 && f4 > $4) {
        print prev
        print $0
    }
}' tmp.samechromosome.exactlytwice.qualityfiltered.${splitread_file} | awk -v OFS='\t' '{$NF=""; print $0}' | sed 's/[ \t]\+$//' > tmp.oriented.samechromosome.exactlytwice.qualityfiltered.${splitread_file}

# putting them all back together
cat tmp.oriented.samechromosome.exactlytwice.qualityfiltered.reverseread1.${SAMPLE}.sam \
    tmp.oriented.samechromosome.exactlytwice.qualityfiltered.reverseread2.${SAMPLE}.sam \
    tmp.oriented.samechromosome.exactlytwice.qualityfiltered.forwardread1.${SAMPLE}.sam \
    tmp.oriented.samechromosome.exactlytwice.qualityfiltered.forwardread2.${SAMPLE}.sam \
    tmp.oriented.samechromosome.exactlytwice.qualityfiltered.reversemerged.${SAMPLE}.sam \
    tmp.oriented.samechromosome.exactlytwice.qualityfiltered.forwardmerged.${SAMPLE}.sam > tmp.oriented.samechromosome.exactlytwice.qualityfiltered.all.${SAMPLE}.sam

# verify that the match lengths on either side match up so that only proper split reads are looked on
python /global/home/users/pierrj/git/python/filter_for_match_lengths.py tmp.oriented.samechromosome.exactlytwice.qualityfiltered.all.${SAMPLE}.sam tmp.match_length_filtered.oriented.samechromosome.exactlytwice.qualityfiltered.all.${SAMPLE}.sam

### COMBINE OUTPUT WITH OUTPUT FROM MULTI MAPPING SPLIT READS BELOW ###


### get opposite facing read pairs, taken from ecc_caller_mapq0.sh

# get outward facing read pairs using sam flags
# convert to bed file
# fix names for filtering
# filter to appearing only exactly twice, meaning that only complete read pairs are present
samtools view ${FILTERED_BAMFILE} | awk '{ if (($2 == 81 || $2 == 83 || $2 == 145 || $2 == 147 ) && $9 > 0) print $0 ; else if (($2 == 97 || $2 == 99 || $2 == 161 || $2 == 163) && $9 <0) print $0}' | cat <(samtools view -H ${FILTERED_BAMFILE}) - | samtools view -b -h - > tmp.outwardfacing.${SAMPLE}.bam
bedtools bamtobed -i tmp.outwardfacing.${SAMPLE}.bam | sort -k 4,4 > tmp.outwardfacing.${SAMPLE}.bed
mv tmp.outwardfacing.${SAMPLE}.bed tmp.outwardfacing.${SAMPLE}.bed.old
awk 'BEGIN {OFS="\t"}; {print $1,$2,$3,substr($4, 1, length($4)-2),$5,$6}' tmp.outwardfacing.${SAMPLE}.bed.old > tmp.outwardfacing.${SAMPLE}.bed.old.trimmed
awk 'NR==FNR{a[$4]++; next} a[$4]==2' tmp.outwardfacing.${SAMPLE}.bed.old.trimmed tmp.outwardfacing.${SAMPLE}.bed.old.trimmed > outwardfacing.${SAMPLE}.bed

## get outward facing pe reads
## convert to bedpe format and then to bam file
awk -v OFS='\t' '{
    f1=$1; f2=$2; f3=$3; f4=$4; f5=$5; f6=$6
    getline 
    {
        print f1, f2, f3, $1, $2, $3, $4, $5, f6, $6
    }
}' outwardfacing.${SAMPLE}.bed > outwardfacing.${SAMPLE}.bedpe

## deal with bedtools older version bug
source activate /global/scratch/users/pierrj/conda_envs/bedtools

/global/scratch/users/pierrj/conda_envs/bedtools/bin/bedtools bedpetobam -i outwardfacing.${SAMPLE}.bedpe -g ${GENOME_CHROMSIZES} > outwardfacing.${SAMPLE}.bam

samtools sort outwardfacing.${SAMPLE}.bam > outwardfacing.${SAMPLE}.sorted.bam

samtools index outwardfacing.${SAMPLE}.sorted.bam

## get pacbio bam reads, index and sort

samtools sort /global/scratch/users/pierrj/eccDNA/rerunning_stuff_final/pipeline_qc/pacbio/${sample_target}_aligned.bam > \
    /global/scratch/users/pierrj/eccDNA/rerunning_stuff_final/pipeline_qc/pacbio/${sample_target}_aligned.sorted.bam

samtools index /global/scratch/users/pierrj/eccDNA/rerunning_stuff_final/pipeline_qc/pacbio/${sample_target}_aligned.sorted.bam

samtools view -h /global/scratch/users/pierrj/eccDNA/rerunning_stuff_final/pipeline_qc/pacbio/${sample_target}_aligned.sorted.bam \
        "${effector_chrom}:${region_start}-${region_end}" > \
        ${effector_name}.${sample_target}.pacbioreads.bam

samtools sort ${effector_name}.${sample_target}.pacbioreads.bam > ${effector_name}.${sample_target}.pacbioreads.sorted.bam

samtools index ${effector_name}.${sample_target}.pacbioreads.sorted.bam

## pacbio splitreads specifically, index and sort

samtools view -F 4 ${effector_name}.${sample_target}.pacbioreads.bam | awk '{print $2}' > ${effector_name}.${sample_target}_aligned_samflags

bedtools bamtobed -cigar -i /global/scratch/users/pierrj/eccDNA/pipeline_tests/gviz_effector_ecc/${effector_name}.${sample_target}.pacbioreads.bam \
    | paste - ${effector_name}.${sample_target}_aligned_samflags > ${effector_name}.${sample_target}_aligned.bed

python /global/home/users/pierrj/git/python/ecc_caller_pacbio_splitreadout.py ${effector_name}.${sample_target}_aligned.bed \
    ${effector_name}.${sample_target}_pacbio_eccs.bed 20 50000 ${effector_name}.${sample_target}_pacbio_splitreads.bed

/global/scratch/users/pierrj/conda_envs/bedtools/bin/bedtools bedtobam -i ${effector_name}.${sample_target}_pacbio_splitreads.bed -g ${GENOME_CHROMSIZES} > ${effector_name}.${sample_target}_pacbio_splitreads.bam

samtools sort ${effector_name}.${sample_target}_pacbio_splitreads.bam > ${effector_name}.${sample_target}_pacbio_splitreads.sorted.bam

samtools index ${effector_name}.${sample_target}_pacbio_splitreads.sorted.bam


### to deal with getting multi mapping split reads not being in junction split reads ##

## get genome file with everything but region of interest masked
## get read names
## extract read names from fastq file
## remap with BWA MEM 
## call junction split reads as before

GENOME_FASTA=/global/scratch/users/pierrj/references/guy11_genome_baoetal2017.fasta

module load seqtk

echo ${effector_chrom} | seqtk subseq ${GENOME_FASTA} - > effector_chrom.fasta

end_of_chrom=$(grep ${effector_chrom} $GENOME_CHROMSIZES | cut -f 2)

echo -e ${effector_chrom}"\t"1"\t"${region_start} > to_mask
echo -e ${effector_chrom}"\t"${region_end}"\t"${end_of_chrom} >> to_mask

/global/scratch/users/pierrj/conda_envs/bedtools/bin/bedtools maskfasta -fi effector_chrom.fasta -bed to_mask -fo effector_chrom_masked.fasta

samtools view -q 1 ${FILTERED_BAMFILE} | cut -f1 | sort | uniq > readnames_no_mapq0
samtools view ${FILTERED_BAMFILE} | cut -f1 | sort | uniq > readnames_w_mapq0

cat readnames_no_mapq0 readnames_w_mapq0 | sort | uniq -u > readnames_only_mapq0


seqtk subseq /global/scratch/users/pierrj/eccDNA/magnaporthe_pureculture/illumina/${sample_target}/${sample_target}_R1.fastq readnames_only_mapq0 > mapq0_reads_R1.fastq
seqtk subseq /global/scratch/users/pierrj/eccDNA/magnaporthe_pureculture/illumina/${sample_target}/${sample_target}_R2.fastq readnames_only_mapq0 > mapq0_reads_R2.fastq

module load cutadapt
module load bwa

THREADS=${SLURM_NTASKS}
READONE=mapq0_reads_R1.fastq
READTWO=mapq0_reads_R2.fastq
SAMPLE=${effector_name}.${sample_target}.remapped

cutadapt -j ${THREADS} -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT --nextseq-trim=20 -o tmp.trimmed.seqprep.${SAMPLE}_R1.fastq -p tmp.trimmed.seqprep.${SAMPLE}_R2.fastq ${READONE} ${READTWO}

bwa index effector_chrom_masked.fasta

GENOME_DB=effector_chrom_masked.fasta

bwa mem -q -a -t ${THREADS} ${GENOME_DB} tmp.trimmed.seqprep.${SAMPLE}_R1.fastq tmp.trimmed.seqprep.${SAMPLE}_R2.fastq -o tmp.seqprep.trimmed.${SAMPLE}_bwamem.sam

samtools view -S -b tmp.seqprep.trimmed.${SAMPLE}_bwamem.sam > ${SAMPLE}.mergedandpe.bwamem.bam

samtools sort -@ ${THREADS} ${SAMPLE}.mergedandpe.bwamem.bam > ${SAMPLE}.sorted.mergedandpe.bwamem.bam
samtools index ${SAMPLE}.sorted.mergedandpe.bwamem.bam

samtools view -b -F 256 ${SAMPLE}.sorted.mergedandpe.bwamem.bam > no_secondary.filtered.sorted.${SAMPLE}.bam

FILTERED_BAMFILE=no_secondary.filtered.sorted.${effector_name}.${sample_target}.remapped.bam

## get junction split reads, taken directly from ecc_caller_mapq0.sh

samtools view -f 81 -F 4 ${FILTERED_BAMFILE} > tmp.reverseread1.${SAMPLE}.sam
splitread_file="reverseread1.${SAMPLE}.sam"
awk -v OFS='\t' '{a=gensub(/^([0-9]+)M.*[HS]$/, "\\1", "", $6); b=gensub(/.*[HS]([0-9]+)M$/, "\\1", "", $6); if (a !~ /[DMIHS]/ && int(a) > 19 ) print $0, 1; else if (b !~ /[DMIHS]/ && int(b) > 19) print $0, 2}' tmp.${splitread_file} > tmp.qualityfiltered.${splitread_file}
awk 'NR==FNR{a[$1, $3]++; next} a[$1, $3]==2' tmp.qualityfiltered.${splitread_file} tmp.qualityfiltered.${splitread_file} | sort -k1,1 -k18,18n > tmp.samechromosome.exactlytwice.qualityfiltered.${splitread_file}
awk -v OFS='\t' '{
    prev=$0; f4=$4; f1=$1
    getline 
    if ($1 == f1 && f4 > $4) {
        print prev
        print $0
    }
}' tmp.samechromosome.exactlytwice.qualityfiltered.${splitread_file} | awk -v OFS='\t' '{$NF=""; print $0}' | sed 's/[ \t]\+$//' > tmp.oriented.samechromosome.exactlytwice.qualityfiltered.${splitread_file}

samtools view -f 145 -F 4 ${FILTERED_BAMFILE} > tmp.reverseread2.${SAMPLE}.sam
splitread_file="reverseread2.${SAMPLE}.sam"
awk -v OFS='\t' '{a=gensub(/^([0-9]+)M.*[HS]$/, "\\1", "", $6); b=gensub(/.*[HS]([0-9]+)M$/, "\\1", "", $6); if (a !~ /[DMIHS]/ && int(a) > 19 ) print $0, 1; else if (b !~ /[DMIHS]/ && int(b) > 19) print $0, 2}' tmp.${splitread_file} > tmp.qualityfiltered.${splitread_file}
awk 'NR==FNR{a[$1, $3]++; next} a[$1, $3]==2' tmp.qualityfiltered.${splitread_file} tmp.qualityfiltered.${splitread_file} | sort -k1,1 -k18,18n > tmp.samechromosome.exactlytwice.qualityfiltered.${splitread_file}
awk -v OFS='\t' '{
    prev=$0; f4=$4; f1=$1
    getline 
    if ($1 == f1 && f4 > $4) {
        print prev
        print $0
    }
}' tmp.samechromosome.exactlytwice.qualityfiltered.${splitread_file} | awk -v OFS='\t' '{$NF=""; print $0}' | sed 's/[ \t]\+$//' > tmp.oriented.samechromosome.exactlytwice.qualityfiltered.${splitread_file}

samtools view -f 65 -F 20 ${FILTERED_BAMFILE} > tmp.forwardread1.${SAMPLE}.sam
splitread_file="forwardread1.${SAMPLE}.sam"
awk -v OFS='\t' '{a=gensub(/^([0-9]+)M.*[HS]$/, "\\1", "", $6); b=gensub(/.*[HS]([0-9]+)M$/, "\\1", "", $6); if (a !~ /[DMIHS]/ && int(a) > 19 ) print $0, 1; else if (b !~ /[DMIHS]/ && int(b) > 19) print $0, 2}' tmp.${splitread_file} > tmp.qualityfiltered.${splitread_file}
awk 'NR==FNR{a[$1, $3]++; next} a[$1, $3]==2' tmp.qualityfiltered.${splitread_file} tmp.qualityfiltered.${splitread_file} | sort -k1,1 -k18,18n > tmp.samechromosome.exactlytwice.qualityfiltered.${splitread_file}
awk -v OFS='\t' '{
    prev=$0; f4=$4; f1=$1
    getline 
    if ($1 == f1 && f4 > $4) {
        print prev
        print $0
    }
}' tmp.samechromosome.exactlytwice.qualityfiltered.${splitread_file} | awk -v OFS='\t' '{$NF=""; print $0}' | sed 's/[ \t]\+$//' > tmp.oriented.samechromosome.exactlytwice.qualityfiltered.${splitread_file}

samtools view -f 129 -F 20 ${FILTERED_BAMFILE} > tmp.forwardread2.${SAMPLE}.sam
splitread_file="forwardread2.${SAMPLE}.sam"
awk -v OFS='\t' '{a=gensub(/^([0-9]+)M.*[HS]$/, "\\1", "", $6); b=gensub(/.*[HS]([0-9]+)M$/, "\\1", "", $6); if (a !~ /[DMIHS]/ && int(a) > 19 ) print $0, 1; else if (b !~ /[DMIHS]/ && int(b) > 19) print $0, 2}' tmp.${splitread_file} > tmp.qualityfiltered.${splitread_file}
awk 'NR==FNR{a[$1, $3]++; next} a[$1, $3]==2' tmp.qualityfiltered.${splitread_file} tmp.qualityfiltered.${splitread_file} | sort -k1,1 -k18,18n > tmp.samechromosome.exactlytwice.qualityfiltered.${splitread_file}
awk -v OFS='\t' '{
    prev=$0; f4=$4; f1=$1
    getline 
    if ($1 == f1 && f4 > $4) {
        print prev
        print $0
    }
}' tmp.samechromosome.exactlytwice.qualityfiltered.${splitread_file} | awk -v OFS='\t' '{$NF=""; print $0}' | sed 's/[ \t]\+$//' > tmp.oriented.samechromosome.exactlytwice.qualityfiltered.${splitread_file}

# I did some merging of reads in the past but this ended up being detrimental
# these next two chunks should be removed eventually
# currently both of these should contain 0 reads
samtools view -f 16 -F 5 ${FILTERED_BAMFILE} > tmp.reversemerged.${SAMPLE}.sam
splitread_file="reversemerged.${SAMPLE}.sam"
awk -v OFS='\t' '{a=gensub(/^([0-9]+)M.*[HS]$/, "\\1", "", $6); b=gensub(/.*[HS]([0-9]+)M$/, "\\1", "", $6); if (a !~ /[DMIHS]/ && int(a) > 19 ) print $0, 1; else if (b !~ /[DMIHS]/ && int(b) > 19) print $0, 2}' tmp.${splitread_file} > tmp.qualityfiltered.${splitread_file}
awk 'NR==FNR{a[$1, $3]++; next} a[$1, $3]==2' tmp.qualityfiltered.${splitread_file} tmp.qualityfiltered.${splitread_file} | sort -k1,1 -k18,18n > tmp.samechromosome.exactlytwice.qualityfiltered.${splitread_file}
awk -v OFS='\t' '{
    prev=$0; f4=$4; f1=$1
    getline 
    if ($1 == f1 && f4 > $4) {
        print prev
        print $0
    }
}' tmp.samechromosome.exactlytwice.qualityfiltered.${splitread_file} | awk -v OFS='\t' '{$NF=""; print $0}' | sed 's/[ \t]\+$//' > tmp.oriented.samechromosome.exactlytwice.qualityfiltered.${splitread_file}

samtools view -F 21 ${FILTERED_BAMFILE} > tmp.forwardmerged.${SAMPLE}.sam
splitread_file="forwardmerged.${SAMPLE}.sam"
awk -v OFS='\t' '{a=gensub(/^([0-9]+)M.*[HS]$/, "\\1", "", $6); b=gensub(/.*[HS]([0-9]+)M$/, "\\1", "", $6); if (a !~ /[DMIHS]/ && int(a) > 19 ) print $0, 1; else if (b !~ /[DMIHS]/ && int(b) > 19) print $0, 2}' tmp.${splitread_file} > tmp.qualityfiltered.${splitread_file}
awk 'NR==FNR{a[$1, $3]++; next} a[$1, $3]==2' tmp.qualityfiltered.${splitread_file} tmp.qualityfiltered.${splitread_file} | sort -k1,1 -k18,18n > tmp.samechromosome.exactlytwice.qualityfiltered.${splitread_file}
awk -v OFS='\t' '{
    prev=$0; f4=$4; f1=$1
    getline 
    if ($1 == f1 && f4 > $4) {
        print prev
        print $0
    }
}' tmp.samechromosome.exactlytwice.qualityfiltered.${splitread_file} | awk -v OFS='\t' '{$NF=""; print $0}' | sed 's/[ \t]\+$//' > tmp.oriented.samechromosome.exactlytwice.qualityfiltered.${splitread_file}

# putting them all back together
cat tmp.oriented.samechromosome.exactlytwice.qualityfiltered.reverseread1.${SAMPLE}.sam \
    tmp.oriented.samechromosome.exactlytwice.qualityfiltered.reverseread2.${SAMPLE}.sam \
    tmp.oriented.samechromosome.exactlytwice.qualityfiltered.forwardread1.${SAMPLE}.sam \
    tmp.oriented.samechromosome.exactlytwice.qualityfiltered.forwardread2.${SAMPLE}.sam \
    tmp.oriented.samechromosome.exactlytwice.qualityfiltered.reversemerged.${SAMPLE}.sam \
    tmp.oriented.samechromosome.exactlytwice.qualityfiltered.forwardmerged.${SAMPLE}.sam > tmp.oriented.samechromosome.exactlytwice.qualityfiltered.all.${SAMPLE}.sam

# verify that the match lengths on either side match up so that only proper split reads are looked on
python /global/home/users/pierrj/git/python/filter_for_match_lengths.py tmp.oriented.samechromosome.exactlytwice.qualityfiltered.all.${SAMPLE}.sam tmp.match_length_filtered.oriented.samechromosome.exactlytwice.qualityfiltered.all.${SAMPLE}.sam



## combine non-multimapped with multi-mapped split reads

SAMPLE=${effector_name}.${sample_target}

cat tmp.match_length_filtered.oriented.samechromosome.exactlytwice.qualityfiltered.all.${SAMPLE}.sam \
    tmp.match_length_filtered.oriented.samechromosome.exactlytwice.qualityfiltered.all.${SAMPLE}.remapped.sam > \
    tmp.match_length_filtered.oriented.samechromosome.exactlytwice.qualityfiltered.all.${SAMPLE}.both.sam


# converting to bed file
samtools view -b -h <(cat <(samtools view -H ${FILTERED_BAMFILE}) tmp.match_length_filtered.oriented.samechromosome.exactlytwice.qualityfiltered.all.${SAMPLE}.both.sam) > tmp.oriented.samechromosome.exactlytwice.qualityfiltered.all.${SAMPLE}.bam

cp tmp.oriented.samechromosome.exactlytwice.qualityfiltered.all.${SAMPLE}.bam junction_splitreads.${SAMPLE}.bam

samtools sort junction_splitreads.${SAMPLE}.bam > junction_splitreads.${SAMPLE}.sorted.bam

samtools index junction_splitreads.${SAMPLE}.sorted.bam