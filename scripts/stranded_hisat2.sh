#!/bin/bash
#SBATCH -t 120:00:00
#SBATCH --nodes=1 --ntasks-per-node=20
#SBATCH --mem=200GB
#SBATCH --export=NONE
#SBATCH --error=../scripts/outs_errs/"%x_error.%j" #write out slurm error reports
#SBATCH --output=../scripts/outs_errs/"%x_output.%j" #write out any program outpus
#SBATCH --mail-type=BEGIN,END,FAIL #email you when job starts, stops and/or fails
#SBATCH --mail-user=zdellaert@uri.edu #your email to send notifications
#SBATCH -D /data/putnamlab/zdellaert/LaserCoral/output_RNA #set working directory

# load modules needed
module load HISAT2/2.2.1-gompi-2022a #Alignment to reference genome: HISAT2
module load SAMtools/1.16.1-GCC-11.3.0 #Preparation of alignment for assembly: SAMtools

# index the reference genome, will write to a directory called Pacuta_ref
hisat2-build -f ../references/Pocillopora_acuta_HIv2.assembly.fasta ./Pacuta_ref

echo "Reference genome indexed. Starting alingment" $(date)

# make the output directory if it does not exist (-p checks for this)
mkdir -p stranded_hisat2

# call the oligo-trimmed sequences into an array
array=(../data_RNA/trimmed_oligo*_R1_001.fastq.gz)

# align the files to the indexed genome using hisat2
for read1 in ${array[@]}; do
    
    # extract the sample name of R1 files as LCM_##
    sample_name=$(basename $read1 | sed 's/.*trimmed_\([A-Za-z0-9]*_[0-9]*\).*/\1/')
    
    # Define the corresponding reverse read file (R2)
    read2=${read1/_R1_/_R2_}
    
    # perform alignment
    hisat2 -p 16 --time --dta -q --rna-strandness RF -x Pacuta_ref -1 ${read1} -2 ${read2} -S stranded_hisat2/${sample_name}.sam
    echo "${sample_name} aligned!"

    # sort the sam file into a bam file
    samtools sort -@ 8 -o stranded_hisat2/${sample_name}.bam stranded_hisat2/${sample_name}.sam
    echo "${sample_name} bam-ified!"
    
    # index bam file , creating a .bai file which is nice for viewing in IGB
    samtools index stranded_hisat2/${sample_name}.bam stranded_hisat2/${sample_name}.bai
    
    # remove sam file to save disk space
    rm stranded_hisat2/${sample_name}.sam
done

# move the reference index files into the stranded_hisat2 directory
mv Pacuta_ref.* stranded_hisat2/

#  Calculate mapping percentages
for i in stranded_hisat2/*.bam; do
    echo "${i}" >> stranded_hisat2/mapped_reads_counts_Pacuta
    samtools flagstat ${i} | grep "mapped (" >> stranded_hisat2/mapped_reads_counts_Pacuta
done
