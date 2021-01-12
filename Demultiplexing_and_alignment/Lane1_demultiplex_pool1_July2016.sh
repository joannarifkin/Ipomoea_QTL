#!/bin/bash
#
#SBATCH --mem=3000
#SBATCH --job-name=L_1_P_1_demultiplex
#SBATCH -p rausherlab 
#SBATCH --mail-user=jlr42
#SBATCH --mail-type=ALL
python John_Kelly_script_8_base_bc.py /work/rausher/jlr42/raw_reads/Data_order_2701_October_2015/unzipped/Pool_1_ATCACG/all_R1.fastq /work/rausher/jlr42/raw_reads/Data_order_2701_October_2015/unzipped/Pool_1_ATCACG/all_R2.fastq ATCACG
