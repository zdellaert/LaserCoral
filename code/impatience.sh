#test trimmed V5 impatient

salloc -c 6 --mem 64G

cd /home/zdellaert_uri_edu/zdellaert_uri_edu-shared/output_WGBS/test

module load uri/main
module load Bismark/0.23.1-foss-2021b
module load bowtie2/2.5.2

cp /home/zdellaert_uri_edu/LaserCoral/references/Pocillopora_acuta_HIv2.assembly.fasta .
bismark_genome_preparation --verbose --parallel 4 ./

# Run Bismark align
bismark \
    -genome ./ \
    -p 4 \
    -u 10000 \
    -score_min L,0,-1.0 \
    -1 ../../data_WGBS/trimmed_V5/trimmed_V5_trimmed_V3_LCM_24_S5_R1_001.fastq.gz \
    -2 ../../data_WGBS/trimmed_V5/trimmed_V5_trimmed_V3_LCM_24_S5_R2_001.fastq.gz \
    -o ./ \
    2> "trimmed_V5_LCM_24-bismark_summary.txt"


# Results: Mapping efficincey increased from 53% to 54.7%

# Run Bismark align
bismark \
    -genome ./ \
    -p 4 \
    -u 10000 \
    -score_min L,0,-1.0 \
    -1 ../../data_WGBS/trimmed_V5/trimmed_V5_trimmed_V3_LCM_11_S3_R1_001.fastq.gz \
    -2 ../../data_WGBS/trimmed_V5/trimmed_V5_trimmed_V3_LCM_11_S3_R2_001.fastq.gz \
    -o ./ \
    2> "trimmed_V5_LCM_11-bismark_summary.txt"


# Results: Mapping efficincey decreased from 28% to 26.1%