#!/bin/bash
#SBATCH -p nbi-short
#SBATCH --mem=1G
#SBATCH -o /jic/scratch/groups/Philippa-Borrill/Heidi/fastqc_results/slurm_output/%x.%N.%j.out
#SBATCH -e /jic/scratch/groups/Philippa-Borrill/Heidi/fastqc_results/slurm_output/%x.%N.%j.err
#SBATCH --mail-type=end,fail
#SBATCH --mail-user=town@nbi.ac.uk

source fastqc-0.11.8

cd /jic/research-groups/Philippa-Borrill/raw_data/example_RNA_seq/fastq/

fastqc -o /jic/scratch/groups/Philippa-Borrill/Heidi/fastqc_results *.fastq.gz

