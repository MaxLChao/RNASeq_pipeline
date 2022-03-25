#!/bin/bash -l
#SBATCH --partition=angsd_class
#SBATCH --ntasks=1
#SBATCH --job-name=alignment
#SBATCH --time=24:00:00
#SBATCH --mem=30G
#SBATCH --ntasks=1 --cpus-per-task=6 #threads
#SBATCH --mail-user=mlc333@cornell.edu
#SBATCH --mail-type=END,FAIL


# script to fastp, run STAR, and run samtools to sort these files. Then summarize via
# multiqc and use feature counts

# set up locs
THREADS=6
samps=(/athena/angsd/scratch/mac4017/Final_Project/Fastqs/zmp_ph12[0-2]_*/)
gDir="/athena/angsd/scratch/mac4017/programs/STAR/annotations/GRCz11"
fastqc_dir="/athena/angsd/scratch/mac4017/Final_Project/Processed/fastqc"
bamqc_dir="/athena/angsd/scratch/mac4017/Final_Project/Processed/bamqc"
Processed="/athena/angsd/scratch/mac4017/Final_Project/Processed"
# set up loads
spack load star@2.7.0
spack load fastqc
spack load subread
spack load -r py-multiqc
spack load -r trimgalore


cd $Processed
# here we generate` fastqc output + generate alignment, run bamqc
for i in ${samps[@]}; do
   filedir=$(echo ${i} | cut -d / -f8)
   mkdir $filedir
   cd $filedir
   zcat ${i}/*_1.fastq.gz | gzip > ${filedir}_merged_1.fastq.gz
   zcat ${i}/*_2.fastq.gz | gzip > ${filedir}_merged_2.fastq.gz
   files_to_align=$(ls *_merged_*)
   trim_galore --illumina --paired --fastqc_args "--extract --outdir ../fastqc/"  ${files_to_align}
   rm $files_to_align
   files_to_align=$(ls *fq.gz)
   STAR --runMode alignReads --runThreadN 6 --genomeDir $gDir/STARindex --readFilesIn $files_to_align --readFilesCommand zcat --outFileNamePrefix ./${files_to_align%%.fq.gz}_STAR_ --outSAMtype BAM SortedByCoordinate --outFilterMultimapNmax 4
   bam=$(ls *.bam)
   /softlib/apps/EL7/BamQC/bin/bamqc --dir /scratch -o $bamqc_dir $bam
   cd ..
   echo $filedir sample aligned, $i of ${#samps[@]}
done

# here we summarize and generate the multiqc output
multiqc $fastqc_dir
bams=$(ls zmp_ph12*/*.bam)
featureCounts -a $gDir/Danio_rerio.GRCz11.105.gtf -T 6 -o zebra_counts.txt ${bams[@]}

#
