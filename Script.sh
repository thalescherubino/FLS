screen -r carabica

cd /media/newhd

aws s3 sync s3://f20ftsusat1321 ./

cd

/home/carabica/bin/STAR/source/STAR runThreadN 24 --runMode genomeGenerate --genomeDir ./ --sjdbGTFfile ./ GCF_003713225.1_Cara_1.0_genomic.gff --sjdbGTFtagExonParentTranscript Parent --sjdbGTFfeatureExon exon --genomeFastaFiles ./GCF_003713225.1_Cara_1.0_genomic.fna --genomeSAindexNbases 13

allPath=/media/newhd/teste/

for library in `ls /media/newhd/teste/*1.fq.gz`
do
  echo $library
pigz -p 8 -d -f -k -c $library > `basename -s .gz $library` &

#Add a 10 seconds pause just to avoid problems while decompressing the second file
sleep 10s
echo Resumming

pigz -p 8 -d -f -k -c $allPath`basename -s 1.fq.gz $library`2.fq.gz > `basename -s 1.fq.gz $library`2.fq

done
