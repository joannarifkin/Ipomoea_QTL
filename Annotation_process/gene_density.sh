#get gene list
grep '[[:blank:]]gene[[:blank:]]' data/lacunosa.gff | cut -f 1,4,5 | awk '{print "chr"$1"\t"$2"\t"$3}' > genes.bed

#make 50kb windows
cut -d ' ' -f 3,6 chrlength.txt | tr ' ' '\t' > genome.bed
bedtools makewindows -g genome.bed -w 50000 > lacu.window

#calculate density
bedtools coverage -a lacu.window -b genes.bed| cut -f 1-4 > number.txt

