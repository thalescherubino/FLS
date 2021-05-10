screen -r carabica

cd /media/newhd

aws s3 sync s3://f20ftsusat13211 ./


#index the genome
/home/carabica/bin/STAR/source/STAR runThreadN 24 --runMode genomeGenerate --genomeDir ./ --sjdbGTFfile ./ GCF_003713225.1_Cara_1.0_genomic.gff --sjdbGTFtagExonParentTranscript Parent --sjdbGTFfeatureExon exon --genomeFastaFiles ./GCF_003713225.1_Cara_1.0_genomic.fna --genomeSAindexNbases 13

cd ~/RNAseq/analysis

allPath=/media/newhd/F20FTSUSAT1321_COFayozT/

for library in `ls /media/newhd/F20FTSUSAT1321_COFayozT/*1.fq.gz`
do
  echo $library
pigz -p 8 -d -f -k -c $library > `basename -s .gz $library` &

#Add a 10 seconds pause just to avoid problems while decompressing the second file
sleep 10s
echo Resumming

pigz -p 8 -d -f -k -c $allPath`basename -s 1.fq.gz $library`2.fq.gz > `basename -s 1.fq.gz $library`2.fq

#aling paired end reads

~/bin/STAR/source/STAR --runThreadN 8 --genomeDir ~/genome/ --readFilesIn `basename -s .gz $library` `basename -s 1.fq.gz $library`2.fq --outFileNamePrefix  `basename -s 1.fq.gz $library` --outSAMtype BAM Unsorted

rm `basename -s .gz $library` `basename -s 1.fq.gz $library`2.fq

#keep only the mapped reads in bam file
#samtools view -F 4 -@ 8 `basename -s 1.fq.gz $library`Aligned.out.bam > `basename -s 1.fq.gz $library`onlyMapped.bam
#not a fuking good idea the last command, for some reason the resulting file was file 3 times bigger than the original

#Mark dup with picard

#The program can take either coordinate-sorted or query-sorted inputs, however the behavior is slightly different. When the input is coordinate-sorted, unmapped mates of mapped records and supplementary/secondary alignments are not marked as duplicates. However, when the input is query-sorted (actually query-grouped), then unmapped mates and secondary/supplementary reads are not excluded from the duplication test and can be marked as duplicate reads.


java -XX:ParallelGCThreads=8 -jar ~/bin/picard.jar SortSam I=`basename -s 1.fq.gz $library`Aligned.out.bam O=`basename -s 1.fq.gz $library`sorted.bam SORT_ORDER=queryname

rm `basename -s 1.fq.gz $library`Aligned.out.bam

java -XX:ParallelGCThreads=8 -jar ~/bin/picard.jar MarkDuplicates I=`basename -s 1.fq.gz $library`sorted.bam O=`basename -s 1.fq.gz $library`sorted.rmdup.bam REMOVE_DUPLICATES=true M=`basename -s 1.fq.gz $library`rmdupReport.txt

done
