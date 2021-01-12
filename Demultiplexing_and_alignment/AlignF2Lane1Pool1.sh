while read i
do

echo "Processing sample $i ..." 

ngm -t 10 -r /ohta/joanna.rifkin/Genomes/Ipomoea/Lacunosa/BioNano/EXP_REFINEFINAL1_bppAdjust_cmap_lacunosa_genome_070416_bionano_fasta_NGScontigs_HYBRID_SCAFFOLD.fasta -b -1 /ohta/joanna.rifkin/RawReads/IpomoeaF2/RADSeq_data_JR_15-16/Demultiplexed_reads/Lane_1/Pool1/r1.F2_sample_$i\.fq.gz -2 /ohta/joanna.rifkin/RawReads/IpomoeaF2/RADSeq_data_JR_15-16/Demultiplexed_reads/Lane_1/Pool1/r2.F2_sample_$i\.fq.gz -o /ohta/joanna.rifkin/Alignments/NGM/Ipomoea/2019/Lane1/$i.bam 

done < /ohta/joanna.rifkin/Alignments/NGM/Ipomoea/Lane1Pool1List.txt