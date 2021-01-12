while read i
do

echo "Processing sample $i ..." 

java -jar /ohta/joanna.rifkin/picard/build/libs/picard.jar FixMateInformation \
       I=/ohta/joanna.rifkin/Alignments/NGM/Ipomoea/2019/Lane1/$i.bam 
       
java -jar /ohta/joanna.rifkin/picard/build/libs/picard.jar CleanSam \
INPUT=/ohta/joanna.rifkin/Alignments/NGM/Ipomoea/2019/Lane1/$i.bam \
OUTPUT=/ohta/joanna.rifkin/Alignments/NGM/Ipomoea/2019/Lane1/tmp_$i\_cleaned.bam

java -jar /ohta/joanna.rifkin/picard/build/libs/picard.jar SortSam \
 INPUT=/ohta/joanna.rifkin/Alignments/NGM/Ipomoea/2019/Lane1/tmp_$i\_cleaned.bam \
 OUTPUT=/ohta/joanna.rifkin/Alignments/NGM/Ipomoea/2019/Lane1/tmp_$i\_cleaned_sorted.bam \
 SORT_ORDER=coordinate CREATE_INDEX=true \
 VALIDATION_STRINGENCY=LENIENT \
 TMP_DIR=`pwd`/tmp 

 
java -jar /ohta/joanna.rifkin/picard/build/libs/picard.jar AddOrReplaceReadGroups \
      VALIDATION_STRINGENCY=LENIENT \
      I=/ohta/joanna.rifkin/Alignments/NGM/Ipomoea/2019/Lane1/tmp_$i\_cleaned_sorted.bam \
      O=/ohta/joanna.rifkin/Alignments/NGM/Ipomoea/2019/Lane1/mergeready_$i\.bam \
      RGID=1$i \
      RGLB=lib1$i \
      RGPL=illumina \
      RGPU=unit1$i \
      RGSM=sample$i
	  
	  
done < /ohta/joanna.rifkin/Alignments/NGM/Ipomoea/Lane1List.txt
