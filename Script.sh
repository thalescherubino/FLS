screen -r carabica

cd /media/newhd

aws s3 sync s3://f20ftsusat13211 ./


#index the genome
/home/carabica/bin/STAR/source/STAR runThreadN 24 --runMode genomeGenerate --genomeDir ./ --sjdbGTFfile ./ GCF_003713225.1_Cara_1.0_genomic.gff --sjdbGTFtagExonParentTranscript Parent --sjdbGTFfeatureExon exon --genomeFastaFiles ./GCF_003713225.1_Cara_1.0_genomic.fna --genomeSAindexNbases 13

cd ~/RNAseq/analysis

allPath=/media/newhd/F20FTSUSAT1321_COFayozT/

for library in `ls /media/newhd/F20FTSUSAT1321_COFayozT/*1.fq.gz`
do

if [ -f "`basename -s 1.fq.gz $library`counts.txt" ]; then
echo "`basename -s 1.fq.gz $library`counts.txt exists."
else
echo "`basename -s 1.fq.gz $library`counts.txt does not exist."

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

rm `basename -s 1.fq.gz $library`sorted.bam

/usr/bin/htseq-count  -a 10 -t exon -i Parent -f bam --stranded=no `basename -s 1.fq.gz $library`sorted.rmdup.bam ../../genome/compatible.gff > `basename -s 1.fq.gz $library`counts.txt &
fi
done

#The Differential expression analysis stats here

library("edgeR")
targets <- readTargets()
names <- targets$description
matrix_input <- readDGE(targets, comment.char = "!")

#remove meta Tags
MetaTags <- grep("^__", rownames(matrix_input))
matrix_input <- matrix_input[-MetaTags, ]


reads_before <- sum(matrix_input$counts)
#remove low expressed genes
rnaseqmatrix <- matrix_input$counts
rnaseqmatrix <- rnaseqmatrix[rowMeans(rnaseqmatrix) >=35,]

sum(rnaseqmatrix)/reads_before

conditions = matrix_input$samples[,2]
analysis_matrix <- DGEList(counts = rnaseqmatrix,group = conditions)
colnames(analysis_matrix$counts) <- names
design <- model.matrix(~0+group, data=analysis_matrix$samples)
colnames(design) <- levels(analysis_matrix$samples$group)


#NORMALIZATIONS
analysis_matrix <- calcNormFactors(analysis_matrix)

#To estimate common dispersion:
analysis_matrix <- estimateGLMCommonDisp(analysis_matrix, design)
#To estimate trended dispersions:
analysis_matrix <- estimateGLMTrendedDisp(analysis_matrix, design)
#To estimate tagwise dispersions:
analysis_matrix <- estimateGLMTagwiseDisp(analysis_matrix, design)

#Fit a negative binomial generalized log-linear model to the read counts for each gene.
fit <- glmFit(analysis_matrix,design)

pdf(file = "edgeR_BCV.pdf",h=8,w=8)
plotBCV(analysis_matrix)
dev.off()

system("open edgeR_BCV.pdf")

#An MDS plots shows distances, in terms of biological coeficient of variation (BCV) - An MDS plot shows the relative similarities of the samples.

cols1 <- colorRampPalette(c("#c24a3a","#522c28","#7bc437"))(48)

pdf(file = "edgeR_MDS.pdf",h=15,w=15)
plotMDS(analysis_matrix,col = cols1,cex=1.5)
dev.off()
system("open edgeR_MDS.pdf")


#samples_trees
cpm.matrix <- cpm(analysis_matrix,normalized.lib.sizes=F)
colnames(cpm.matrix) <- names
t.cpm.matrix <- t(cpm.matrix)
sampleTree <- hclust(dist(t.cpm.matrix), method = "average");

# Plot the sample tree: Open a graphic output window of size 12 by 9 inches
# The user should change the dimensions if the window is too large or too small.
#sizeGrWindow(12,9)
pdf(file = "sampleClustering.pdf", width = 12, height = 9);
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
cex.axis = 1.5, cex.main = 2)
# Plot a line to show the cut
#abline(h = 15, col = "red")
dev.off()

system("open sampleClustering.pdf")


#analyse densisties before normalization - input
pdf("Densities_input.pdf")
plotDensities( log(matrix_input$counts), legend = "topright")
dev.off()

#analyse densisties low expression filter
pdf("Densities_low_expression_fiter.pdf")
plotDensities( log(rnaseqmatrix), legend = "topright")
dev.off()

#analyse densisties after normalization
pdf("Densities_normalization.pdf")
plotDensities( log(cpm.matrix), legend = "topright")
dev.off()

pdf("Densities_log_cpm_fitted_norm.pdf")
plotDensities(log(cpm(fit$fitted.values, normalized.lib.sizes=TRUE)), legend = "topright")
dev.off()

#histogram of densities log10
pdf("Log10_histogram_normilized.pdf", h=10,w=10)
hist(log(cpm.matrix+1,10), col=gray.colors(19, start = 0.9, end = 0.3))
dev.off()

#histogram of densities Log2
pdf("Log2_histogram_normilized.pdf", h=10,w=10)
hist(log(cpm.matrix+1,2), col=gray.colors(19, start = 0.9, end = 0.3))
dev.off()

#histogram of densities no log
pdf("histogram_normilized.pdf", h=10,w=10)
hist(cpm.matrix, col=gray.colors(19, start = 0.9, end = 0.3))
dev.off()

# Load library for pheatmap
cpm.matrix.corr <- cor(cpm.matrix, method="spearman",use="pairwise.complete.obs")

library("pheatmap")
library(gplots)
pdf(file = "sampleClusteringHeatmap.pdf", width = 50, height = 50);
pheatmap(cpm.matrix.corr, fontsize=50,cellwidth=50, cellheight=50, treeheight_col= 450, treeheight_row=450,angle_col=90, legend=F, cex.legend=1.2)
#heatmap.2(cpm.matrix.corr*10000,trace = "none",margins = c(5, 11),keysize =1 , key.title="",col ="bluered",density.info="none")
dev.off()

system("open sampleClusteringHeatmap.pdf")

system("open *pdf")
