screen -r carabica

cd /media/newhd

aws s3 sync s3://f20ftsusat1321 ./

cd

/home/carabica/bin/STAR/source/STAR runThreadN 24 --runMode genomeGenerate --genomeDir ./ --sjdbGTFfile ./ GCF_003713225.1_Cara_1.0_genomic.gff --sjdbGTFtagExonParentTranscript Parent --sjdbGTFfeatureExon exon --genomeFastaFiles ./GCF_003713225.1_Cara_1.0_genomic.fna --genomeSAindexNbases 13

for library in 'ls /media/newhd/teste/*1_fq.gz'
do
  echo $library
pigz -p 8 -d -f -k $library > `basename -s .gz $library`

pigz -p 8 -d -f -k `basename -s 1_fq.gz $library`2_fq.gz > `basename -s 1_fq.gz $library`2_fq `

done
