#!/bin/bash
#SBATCH -p long
#SBATCH --job-name=HEPG2_rna_seq
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=alhi9177@colorado.edu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=6gb
#SBATCH --time=20:00:00
#SBATCH --output=nextflow.out
#SBATCH --error=nextflow.err

pwd; hostname; date
echo "Here we go You've requested $SLURM_CPUS_ON_NODE core."

module load singularity/3.1.1

nextflow run nf-core/rnaseq -r 1.4.2 \
-resume \
-profile singularity \
--reads '/scratch/Shares/rinnclass/CLASS_2022/alhi9177/bchm5631_alex/CLASS_2022/class_exeRcises/analysis/17_API_to_RNA/fastq/*{_read1,_read2}.fastq.gz' \
--fasta /scratch/Shares/rinn/genomes/Homo_sapiens/Gencode/v32/GRCh38.p13.genome.fa \
--gtf /scratch/Shares/rinn/genomes/Homo_sapiens/Gencode/v32/gencode.v32.annotation.gtf \
--pseudo_aligner salmon \
--gencode \
--email alhi9177@colorado.edu \
-c nextflow.config

date
/scratch/Shares/rinnclass/CLASS_2022/alhi9177/bchm5631_alex/CLASS_2022/class_exeRcises/analysis/18/rnaseq
