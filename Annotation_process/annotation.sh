##build database for repeat identification
BuildDatabase -name lacu_final -engine ncbi ../../Ipomoea_Lep-Anchor_final.fasta

##identify de novo repeatness
RepeatModeler -database lacu_final -engine ncbi -pa 10

##identify universal repeatness
RepeatMasker -lib lacu_final-families.fa -e ncbi -pa 10 -dir . ../../Ipomoea_Lep-Anchor_final.fasta

##make index for RNAseq alignment
STAR  --runMode genomeGenerate \
    --genomeDir index \
    --runThreadN 11 \
    --genomeFastaFiles ../../Ipomoea_Lep-Anchor_final.fasta \
    --genomeSAindexNbases 13

##align RNAseq against genome
##only shows file names for one group, exists three others
STAR \
    --genomeDir /hpchome/rausherlab/gc158/braker2/star/index \
    --runThreadN 7 \
    --readFilesIn KS_anther_S44_L001_R1_001.fastq.gz.P.qtrim.gz,KS_anther_S44_L002_R1_001.fastq.gz.P.qtrim.gz,KS_anther_S81_L006_R1_001.fastq.gz.P.qtrim.gz,KS_leaf_S28_L001_R1_001.fastq.gz.P.qtrim.gz,KS_leaf_S28_L002_R1_001.fastq.gz.P.qtrim.gz,KS_leaf_S65_L006_R1_001.fastq.gz.P.qtrim.gz,KS_pistil_S60_L001_R1_001.fastq.gz.P.qtrim.gz,KS_pistil_S60_L002_R1_001.fastq.gz.P.qtrim.gz,KS_pistil_S97_L006_R1_001.fastq.gz.P.qtrim.gz,KS_pollen_S105_L006_R1_001.fastq.gz.P.qtrim.gz,KS_pollen_S68_L001_R1_001.fastq.gz.P.qtrim.gz,KS_pollen_S68_L002_R1_001.fastq.gz.P.qtrim.gz,KS_root_S40_L001_R1_001.fastq.gz.P.qtrim.gz,KS_root_S40_L002_R1_001.fastq.gz.P.qtrim.gz,KS_root_S77_L006_R1_001.fastq.gz.P.qtrim.gz,KS_self_S101_L006_R1_001.fastq.gz.P.qtrim.gz,KS_self_S64_L001_R1_001.fastq.gz.P.qtrim.gz,KS_self_S64_L002_R1_001.fastq.gz.P.qtrim.gz KS_anther_S44_L001_R2_001.fastq.gz.P.qtrim.gz,KS_anther_S44_L002_R2_001.fastq.gz.P.qtrim.gz,KS_anther_S81_L006_R2_001.fastq.gz.P.qtrim.gz,KS_leaf_S28_L001_R2_001.fastq.gz.P.qtrim.gz,KS_leaf_S28_L002_R2_001.fastq.gz.P.qtrim.gz,KS_leaf_S65_L006_R2_001.fastq.gz.P.qtrim.gz,KS_pistil_S60_L001_R2_001.fastq.gz.P.qtrim.gz,KS_pistil_S60_L002_R2_001.fastq.gz.P.qtrim.gz,KS_pistil_S97_L006_R2_001.fastq.gz.P.qtrim.gz,KS_pollen_S105_L006_R2_001.fastq.gz.P.qtrim.gz,KS_pollen_S68_L001_R2_001.fastq.gz.P.qtrim.gz,KS_pollen_S68_L002_R2_001.fastq.gz.P.qtrim.gz,KS_root_S40_L001_R2_001.fastq.gz.P.qtrim.gz,KS_root_S40_L002_R2_001.fastq.gz.P.qtrim.gz,KS_root_S77_L006_R2_001.fastq.gz.P.qtrim.gz,KS_self_S101_L006_R2_001.fastq.gz.P.qtrim.gz,KS_self_S64_L001_R2_001.fastq.gz.P.qtrim.gz,KS_self_S64_L002_R2_001.fastq.gz.P.qtrim.gz \
    --readFilesCommand zcat \
    --outFileNamePrefix ks_P_ \
    --outSAMtype BAM SortedByCoordinate \
    --outBAMsortingThreadN 7 \
    --outSAMstrandField intronMotif

##merge alignments
stringtie /hpchome/rausherlab/gc158/braker2/star/ks_P_Aligned.sortedByCoord.out.bam \
         -o ks_p.gtf \
         -p 22 
stringtie --merge -c 0 -F 0.1 -T 0.1 -g 300 -f 0 -m 30 -o lacu_tie.gtf -p 10 mergelist

##extract fasta
gffread -w lacu_final.trans.fasta -g /hpchome/rausherlab/gc158/genome/Ipomoea_Lep-Anchor_final.fasta lacu_tie.gtf

##braker annotation
module load Perl/5.30.0
export PERL5LIB=/opt/apps/rhel7/perl-5.30.0/lib/site_perl/5.30.0
braker.pl --cores 11 \
          --species=lacunosa_final \
          --genome=repeat/Ipomoea_Lep-Anchor_final.fasta.masked \
          --softmasking \
          --bam=star/kent_P_Aligned.sortedByCoord.out.bam,star/ks_P_Aligned.sortedByCoord.out.bam,star/pitx_P_Aligned.sortedByCoord.out.bam,star/KOm_P.bam \
          --gff3 \
          --UTR=on

##measure chracteristics
agat_sp_statistics.pl --gff ../braker2/braker/braker_utr.gff3 \
        -g /hpchome/rausherlab/gc158/genome/Ipomoea_Lep-Anchor_final.fasta \
        -o braker_utr.gff3.stat

##run with maker for the 1st time
##maker_opt file
##the re-annotation section and gene-prediction section are subject to change during iteration
	#-----Genome (these are always required)
	genome=/hpchome/rausherlab/gc158/genome/Ipomoea_Lep-Anchor_final.fasta.masked #genome sequence (fasta file or fasta embeded in GFF3 file)
	organism_type=eukaryotic #eukaryotic or prokaryotic. Default is eukaryotic

	#-----Re-annotation Using MAKER Derived GFF3
	maker_gff= #MAKER derived GFF3 file
	est_pass=0 #use ESTs in maker_gff: 1 = yes, 0 = no
	altest_pass=0 #use alternate organism ESTs in maker_gff: 1 = yes, 0 = no
	protein_pass=0 #use protein alignments in maker_gff: 1 = yes, 0 = no
	rm_pass=0 #use repeats in maker_gff: 1 = yes, 0 = no
	model_pass=0 #use gene models in maker_gff: 1 = yes, 0 = no
	pred_pass=0 #use ab-initio predictions in maker_gff: 1 = yes, 0 = no
	other_pass=0 #passthrough anyything else in maker_gff: 1 = yes, 0 = no
	#-----EST Evidence (for best results provide a file for at least one)
	est=lacunosa.trans.fasta #set of ESTs or assembled mRNA-seq in fasta format
	altest= #EST/cDNA sequence file in fasta format from an alternate organism
	est_gff= #aligned ESTs or mRNA-seq from an external GFF3 file
	altest_gff= #aligned ESTs from a closly relate species in GFF3 format

	#-----Protein Homology Evidence (for best results provide a file for at least one)
	protein=  #protein sequence file in fasta format (i.e. from mutiple oransisms)
	protein_gff=  #aligned protein homology evidence from an external GFF3 file

	#-----Repeat Masking (leave values blank to skip repeat masking)
	model_org= #select a model organism for RepBase masking in RepeatMasker
	rmlib=lacu_final-families.fa #provide an organism specific repeat library in fasta format for RepeatMasker
	repeat_protein=/opt/maker/data/te_proteins.fasta #provide a fasta file of transposable element proteins for RepeatRunner
	rm_gff= #pre-identified repeat elements from an external GFF3 file
	prok_rm=0 #forces MAKER to repeatmask prokaryotes (no reason to change this), 1 = yes, 0 = no
	softmask=1 #use soft-masking rather than hard-masking in BLAST (i.e. seg and dust filtering)

	#-----Gene Prediction
	snaphmm= #SNAP HMM file
	gmhmm= #GeneMark HMM file
	augustus_species= #Augustus gene prediction species model
	fgenesh_par_file= #FGENESH parameter file
	pred_gff=../braker2/braker/braker_utr.gff3 #ab-initio predictions from an external GFF3 file
	model_gff= #annotated gene models from an external GFF3 file (annotation pass-through)
	est2genome=0 #infer gene predictions directly from ESTs, 1 = yes, 0 = no
	in2genome=0 #infer predictions from protein homology, 1 = yes, 0 = no
	trna=0 #find tRNAs with tRNAscan, 1 = yes, 0 = no
	snoscan_rrna= #rRNA file to have Snoscan find snoRNAs
	unmask=0 #also run ab-initio prediction programs on unmasked sequence, 1 = yes, 0 = no

	#-----Other Annotation Feature Types (features MAKER doesn't recognize)
	other_gff= #extra features to pass-through to final MAKER generated GFF3 file

	#-----External Application Behavior Options
	alt_peptide=C #amino acid used to replace non-standard amino acids in BLAST databases
	cpus=1 #max number of cpus to use in BLAST and RepeatMasker (not for MPI, leave 1 when using MPI)

	#-----MAKER Behavior Options
	max_dna_len=100000 #length for dividing up contigs into chunks (increases/decreases memory usage)
	min_contig=1 #skip genome contigs below this length (under 10kb are often useless)

	pred_flank=200 #flank for extending evidence clusters sent to gene predictors
	pred_stats=1 #report AED and QI statistics for all predictions as well as models
	AED_threshold=1 #Maximum Annotation Edit Distance allowed (bound by 0 and 1)
	min_protein=0 #require at least this many amino acids in predicted proteins
	alt_splice=0 #Take extra steps to try and find alternative splicing, 1 = yes, 0 = no
	always_complete=0 #extra steps to force start and stop codons, 1 = yes, 0 = no
	map_forward=0 #map names and attributes forward from old GFF3 genes, 1 = yes, 0 = no
	keep_preds=0 #Concordance threshold to add unsupported gene prediction (bound by 0 and 1)

	split_hit=20000 #length for the splitting of hits (expected max intron size for evidence alignments)
	single_exon=0 #consider single exon EST evidence when generating annotations, 1 = yes, 0 = no
	single_length=250 #min length required for single exon ESTs if 'single_exon is enabled'
	correct_est_fusion=0 #limits use of ESTs in annotation to avoid fusion genes

	tries=2 #number of times to try a contig if there is a failure for some reason
	clean_try=0 #remove all data from previous run before retrying, 1 = yes, 0 = no
	clean_up=1 #removes theVoid directory with individual analysis files, 1 = yes, 0 = no
	TMP= #specify a directory other than the system default temporary directory for temporary files

##train SNAP
maker2zff ../iterate/maker1.gff
fathom genome.ann genome.dna -gene-stats > 1.stats
fathom genome.ann genome.dna -validate > 1.validate
fathom genome.ann genome.dna -categorize 300
fathom uni.ann uni.dna -export 300 -plus
forge export.ann export.dna
hmm-assembler.pl maker1 params > maker1.hmm

##evaluate AED score
AED_cdf_generator.pl -b 0.025 maker1.gff

##remove duplication and miss-matched sequences due to the duplication of RNAseq
agat_convert_sp_gxf2gxf.pl \
        -g maker1.gff \
        -o maker1.sort.gff

#
#calculate statistics, re-train SNAP, re-annotate with MAKER using the codes above
#compare and select one final version for downstream analysis
#
