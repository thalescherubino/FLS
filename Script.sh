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

#from the compatble gff extract usefull information about the quantified transcripts

 awk '{if ($3 != "exon") print $0}' compatible.gff > description.gff

awk -F"\t" '{print $3"\t"$9}' description.gff > temp

cat temp | tr ";" "\t" > temp2

awk -F"\t" '{for(i=5;i<=NF;i++){if($i~/^product=/){a=$i}} print $2"\t"$1"\t"a}' temp2 > temp3

sed -i"" 's/ID=//g' temp3

sed -i"" 's/product=//g' temp3

mv  temp3 C.a.transcriptAnnot.tab

rm temp temp2 description.gff

##################################################
#The Differential expression analysis stats here#
################################################


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


#get gene data
genData <- read.delim("C.a.transcriptAnnot.tab",sep="\t",header=F,fill=T)

colnames(genData) <- c("transcriptName", "type","description")

#######################
#######CONTRASTS#######
#######################


#AcaFBudMar_20_x_OeiFBudMar_20

lrt_AcaFBudMar_20_x_OeiFBudMar_20<- glmLRT(fit, contrast = c(-1,0,0,0,1,0,0,0,0,0,0,0))
tTags_lrt_AcaFBudMar_20_x_OeiFBudMar_20 <- topTags(lrt_AcaFBudMar_20_x_OeiFBudMar_20, n= NULL)
keep_sig <- matrix(tTags_lrt_AcaFBudMar_20_x_OeiFBudMar_20$table$FDR <= 0.05 & abs(tTags_lrt_AcaFBudMar_20_x_OeiFBudMar_20$table$logFC) >= 1)

sig_kept_AcaFBudMar_20_x_OeiFBudMar_20<- tTags_lrt_AcaFBudMar_20_x_OeiFBudMar_20$table[keep_sig[,1],]

sig_kept_AcaFBudMar_20_x_OeiFBudMar_20 <- sig_kept_AcaFBudMar_20_x_OeiFBudMar_20[order(sig_kept_AcaFBudMar_20_x_OeiFBudMar_20$logFC),]

sig_kept_AcaFBudMar_20_x_OeiFBudMar_20 <-cbind( genData[match(rownames(sig_kept_AcaFBudMar_20_x_OeiFBudMar_20),genData$transcriptName),],sig_kept_AcaFBudMar_20_x_OeiFBudMar_20)

write.table(sig_kept_AcaFBudMar_20_x_OeiFBudMar_20, file = "AcaFBudMar_20_x_OeiFBudMar_20.tab", row.names = T, quote = F, sep = "\t")

sig_kept_AcaFBudMar_20_x_OeiFBudMar_20<- tTags_lrt_AcaFBudMar_20_x_OeiFBudMar_20$table[keep_sig[,1],]

sig_kept_AcaFBudMar_20_x_OeiFBudMar_20 <- sig_kept_AcaFBudMar_20_x_OeiFBudMar_20[order(sig_kept_AcaFBudMar_20_x_OeiFBudMar_20$logFC),]

###############################
#AcaFBudMar_20_x_SemFBudMar_20

lrt_AcaFBudMar_20_x_SemFBudMar_20<- glmLRT(fit, contrast = c(-1,0,0,0,0,0,0,0,1,0,0,0))
tTags_lrt_AcaFBudMar_20_x_SemFBudMar_20 <- topTags(lrt_AcaFBudMar_20_x_SemFBudMar_20, n= NULL)
keep_sig <- matrix(tTags_lrt_AcaFBudMar_20_x_SemFBudMar_20$table$FDR <= 0.05 & abs(tTags_lrt_AcaFBudMar_20_x_SemFBudMar_20$table$logFC) >= 1)

sig_kept_AcaFBudMar_20_x_SemFBudMar_20<- tTags_lrt_AcaFBudMar_20_x_SemFBudMar_20$table[keep_sig[,1],]

sig_kept_AcaFBudMar_20_x_SemFBudMar_20 <- sig_kept_AcaFBudMar_20_x_SemFBudMar_20[order(sig_kept_AcaFBudMar_20_x_SemFBudMar_20$logFC),]

sig_kept_AcaFBudMar_20_x_SemFBudMar_20 <-cbind( genData[match(rownames(sig_kept_AcaFBudMar_20_x_SemFBudMar_20),genData$transcriptName),],sig_kept_AcaFBudMar_20_x_SemFBudMar_20)

write.table(sig_kept_AcaFBudMar_20_x_SemFBudMar_20, file = "AcaFBudMar_20_x_SemFBudMar_20.tab", row.names = T, quote = F, sep = "\t")

sig_kept_AcaFBudMar_20_x_SemFBudMar_20<- tTags_lrt_AcaFBudMar_20_x_SemFBudMar_20$table[keep_sig[,1],]

sig_kept_AcaFBudMar_20_x_SemFBudMar_20 <- sig_kept_AcaFBudMar_20_x_SemFBudMar_20[order(sig_kept_AcaFBudMar_20_x_SemFBudMar_20$logFC),]

###############################
#OeiFBudMar_20_x_SemFBudMar_20

lrt_OeiFBudMar_20_x_SemFBudMar_20<- glmLRT(fit, contrast = c(0,0,0,0,-1,0,0,0,1,0,0,0))
tTags_lrt_OeiFBudMar_20_x_SemFBudMar_20 <- topTags(lrt_OeiFBudMar_20_x_SemFBudMar_20, n= NULL)
keep_sig <- matrix(tTags_lrt_OeiFBudMar_20_x_SemFBudMar_20$table$FDR <= 0.05 & abs(tTags_lrt_OeiFBudMar_20_x_SemFBudMar_20$table$logFC) >= 1)

sig_kept_OeiFBudMar_20_x_SemFBudMar_20<- tTags_lrt_OeiFBudMar_20_x_SemFBudMar_20$table[keep_sig[,1],]

sig_kept_OeiFBudMar_20_x_SemFBudMar_20 <- sig_kept_OeiFBudMar_20_x_SemFBudMar_20[order(sig_kept_OeiFBudMar_20_x_SemFBudMar_20$logFC),]

sig_kept_OeiFBudMar_20_x_SemFBudMar_20 <-cbind( genData[match(rownames(sig_kept_OeiFBudMar_20_x_SemFBudMar_20),genData$transcriptName),],sig_kept_OeiFBudMar_20_x_SemFBudMar_20)

write.table(sig_kept_OeiFBudMar_20_x_SemFBudMar_20, file = "OeiFBudMar_20_x_SemFBudMar_20.tab", row.names = T, quote = F, sep = "\t")

sig_kept_OeiFBudMar_20_x_SemFBudMar_20<- tTags_lrt_OeiFBudMar_20_x_SemFBudMar_20$table[keep_sig[,1],]

sig_kept_OeiFBudMar_20_x_SemFBudMar_20 <- sig_kept_OeiFBudMar_20_x_SemFBudMar_20[order(sig_kept_OeiFBudMar_20_x_SemFBudMar_20$logFC),]

###############################
#OeiFBudMar_20_x_OeiSamJan_21

lrt_OeiFBudMar_20_x_OeiSamJan_21<- glmLRT(fit, contrast = c(0,0,0,0,-1,0,0,1,0,0,0,0))
tTags_lrt_OeiFBudMar_20_x_OeiSamJan_21 <- topTags(lrt_OeiFBudMar_20_x_OeiSamJan_21, n= NULL)
keep_sig <- matrix(tTags_lrt_OeiFBudMar_20_x_OeiSamJan_21$table$FDR <= 0.05 & abs(tTags_lrt_OeiFBudMar_20_x_OeiSamJan_21$table$logFC) >= 1)

sig_kept_OeiFBudMar_20_x_OeiSamJan_21<- tTags_lrt_OeiFBudMar_20_x_OeiSamJan_21$table[keep_sig[,1],]

sig_kept_OeiFBudMar_20_x_OeiSamJan_21 <- sig_kept_OeiFBudMar_20_x_OeiSamJan_21[order(sig_kept_OeiFBudMar_20_x_OeiSamJan_21$logFC),]

sig_kept_OeiFBudMar_20_x_OeiSamJan_21 <-cbind( genData[match(rownames(sig_kept_OeiFBudMar_20_x_OeiSamJan_21),genData$transcriptName),],sig_kept_OeiFBudMar_20_x_OeiSamJan_21)

write.table(sig_kept_OeiFBudMar_20_x_OeiSamJan_21, file = "OeiFBudMar_20_x_OeiSamJan_21.tab", row.names = T, quote = F, sep = "\t")

sig_kept_OeiFBudMar_20_x_OeiSamJan_21<- tTags_lrt_OeiFBudMar_20_x_OeiSamJan_21$table[keep_sig[,1],]

sig_kept_OeiFBudMar_20_x_OeiSamJan_21 <- sig_kept_OeiFBudMar_20_x_OeiSamJan_21[order(sig_kept_OeiFBudMar_20_x_OeiSamJan_21$logFC),]

###############################
#AcaFBudMar_20_x_AcaSamJan_21

lrt_AcaFBudMar_20_x_AcaSamJan_21<- glmLRT(fit, contrast = c(-1,0,0,1,0,0,0,0,0,0,0,0))
tTags_lrt_AcaFBudMar_20_x_AcaSamJan_21 <- topTags(lrt_AcaFBudMar_20_x_AcaSamJan_21, n= NULL)
keep_sig <- matrix(tTags_lrt_AcaFBudMar_20_x_AcaSamJan_21$table$FDR <= 0.05 & abs(tTags_lrt_AcaFBudMar_20_x_AcaSamJan_21$table$logFC) >= 1)

sig_kept_AcaFBudMar_20_x_AcaSamJan_21<- tTags_lrt_AcaFBudMar_20_x_AcaSamJan_21$table[keep_sig[,1],]

sig_kept_AcaFBudMar_20_x_AcaSamJan_21 <- sig_kept_AcaFBudMar_20_x_AcaSamJan_21[order(sig_kept_AcaFBudMar_20_x_AcaSamJan_21$logFC),]

sig_kept_AcaFBudMar_20_x_AcaSamJan_21 <-cbind( genData[match(rownames(sig_kept_AcaFBudMar_20_x_AcaSamJan_21),genData$transcriptName),],sig_kept_AcaFBudMar_20_x_AcaSamJan_21)

write.table(sig_kept_AcaFBudMar_20_x_AcaSamJan_21, file = "AcaFBudMar_20_x_AcaSamJan_21.tab", row.names = T, quote = F, sep = "\t")

sig_kept_AcaFBudMar_20_x_AcaSamJan_21<- tTags_lrt_AcaFBudMar_20_x_AcaSamJan_21$table[keep_sig[,1],]

sig_kept_AcaFBudMar_20_x_AcaSamJan_21 <- sig_kept_AcaFBudMar_20_x_AcaSamJan_21[order(sig_kept_AcaFBudMar_20_x_AcaSamJan_21$logFC),]


lrt_SemFBudMar_20_x_SemSamJan_21<- glmLRT(fit, contrast = c(0,0,0,0,0,0,0,0,-1,0,0,1))
tTags_lrt_SemFBudMar_20_x_SemSamJan_21 <- topTags(lrt_SemFBudMar_20_x_SemSamJan_21, n= NULL)
keep_sig <- matrix(tTags_lrt_SemFBudMar_20_x_SemSamJan_21$table$FDR <= 0.05 & abs(tTags_lrt_SemFBudMar_20_x_SemSamJan_21$table$logFC) >= 1)

sig_kept_SemFBudMar_20_x_SemSamJan_21<- tTags_lrt_SemFBudMar_20_x_SemSamJan_21$table[keep_sig[,1],]

sig_kept_SemFBudMar_20_x_SemSamJan_21 <- sig_kept_SemFBudMar_20_x_SemSamJan_21[order(sig_kept_SemFBudMar_20_x_SemSamJan_21$logFC),]

sig_kept_SemFBudMar_20_x_SemSamJan_21 <-cbind( genData[match(rownames(sig_kept_SemFBudMar_20_x_SemSamJan_21),genData$transcriptName),],sig_kept_SemFBudMar_20_x_SemSamJan_21)

write.table(sig_kept_SemFBudMar_20_x_SemSamJan_21, file = "SemFBudMar_20_x_SemSamJan_21.tab", row.names = T, quote = F, sep = "\t")

sig_kept_SemFBudMar_20_x_SemSamJan_21<- tTags_lrt_SemFBudMar_20_x_SemSamJan_21$table[keep_sig[,1],]

sig_kept_SemFBudMar_20_x_SemSamJan_21 <- sig_kept_SemFBudMar_20_x_SemSamJan_21[order(sig_kept_SemFBudMar_20_x_SemSamJan_21$logFC),]

###############################
#AcaSamJan_21_x_OeiSamJan_21

lrt_AcaSamJan_21_x_OeiSamJan_21<- glmLRT(fit, contrast = c(0,0,0,-1,0,0,0,1,0,0,0,0))
tTags_lrt_AcaSamJan_21_x_OeiSamJan_21 <- topTags(lrt_AcaSamJan_21_x_OeiSamJan_21, n= NULL)
keep_sig <- matrix(tTags_lrt_AcaSamJan_21_x_OeiSamJan_21$table$FDR <= 0.05 & abs(tTags_lrt_AcaSamJan_21_x_OeiSamJan_21$table$logFC) >= 1)

sig_kept_AcaSamJan_21_x_OeiSamJan_21<- tTags_lrt_AcaSamJan_21_x_OeiSamJan_21$table[keep_sig[,1],]

sig_kept_AcaSamJan_21_x_OeiSamJan_21 <- sig_kept_AcaSamJan_21_x_OeiSamJan_21[order(sig_kept_AcaSamJan_21_x_OeiSamJan_21$logFC),]

sig_kept_AcaSamJan_21_x_OeiSamJan_21 <-cbind( genData[match(rownames(sig_kept_AcaSamJan_21_x_OeiSamJan_21),genData$transcriptName),],sig_kept_AcaSamJan_21_x_OeiSamJan_21)

write.table(sig_kept_AcaSamJan_21_x_OeiSamJan_21, file = "AcaSamJan_21_x_OeiSamJan_21.tab", row.names = T, quote = F, sep = "\t")

sig_kept_AcaSamJan_21_x_OeiSamJan_21<- tTags_lrt_AcaSamJan_21_x_OeiSamJan_21$table[keep_sig[,1],]

sig_kept_AcaSamJan_21_x_OeiSamJan_21 <- sig_kept_AcaSamJan_21_x_OeiSamJan_21[order(sig_kept_AcaSamJan_21_x_OeiSamJan_21$logFC),]

###############################
#AcaSamJan_21_x_SemSamJan_21

lrt_AcaSamJan_21_x_SemSamJan_21<- glmLRT(fit, contrast = c(0,0,0,-1,0,0,0,0,0,0,0,1))
tTags_lrt_AcaSamJan_21_x_SemSamJan_21 <- topTags(lrt_AcaSamJan_21_x_SemSamJan_21, n= NULL)
keep_sig <- matrix(tTags_lrt_AcaSamJan_21_x_SemSamJan_21$table$FDR <= 0.05 & abs(tTags_lrt_AcaSamJan_21_x_SemSamJan_21$table$logFC) >= 1)

sig_kept_AcaSamJan_21_x_SemSamJan_21<- tTags_lrt_AcaSamJan_21_x_SemSamJan_21$table[keep_sig[,1],]

sig_kept_AcaSamJan_21_x_SemSamJan_21 <- sig_kept_AcaSamJan_21_x_SemSamJan_21[order(sig_kept_AcaSamJan_21_x_SemSamJan_21$logFC),]

sig_kept_AcaSamJan_21_x_SemSamJan_21 <-cbind( genData[match(rownames(sig_kept_AcaSamJan_21_x_SemSamJan_21),genData$transcriptName),],sig_kept_AcaSamJan_21_x_SemSamJan_21)

write.table(sig_kept_AcaSamJan_21_x_SemSamJan_21, file = "AcaSamJan_21_x_SemSamJan_21.tab", row.names = T, quote = F, sep = "\t")

sig_kept_AcaSamJan_21_x_SemSamJan_21<- tTags_lrt_AcaSamJan_21_x_SemSamJan_21$table[keep_sig[,1],]

sig_kept_AcaSamJan_21_x_SemSamJan_21 <- sig_kept_AcaSamJan_21_x_SemSamJan_21[order(sig_kept_AcaSamJan_21_x_SemSamJan_21$logFC),]

###############################
#OeiSamJan_21_x_SemSamJan_21

lrt_OeiSamJan_21_x_SemSamJan_21<- glmLRT(fit, contrast = c(0,0,0,0,0,0,0,-1,0,0,0,1))
tTags_lrt_OeiSamJan_21_x_SemSamJan_21 <- topTags(lrt_OeiSamJan_21_x_SemSamJan_21, n= NULL)
keep_sig <- matrix(tTags_lrt_OeiSamJan_21_x_SemSamJan_21$table$FDR <= 0.05 & abs(tTags_lrt_OeiSamJan_21_x_SemSamJan_21$table$logFC) >= 1)

sig_kept_OeiSamJan_21_x_SemSamJan_21<- tTags_lrt_OeiSamJan_21_x_SemSamJan_21$table[keep_sig[,1],]

sig_kept_OeiSamJan_21_x_SemSamJan_21 <- sig_kept_OeiSamJan_21_x_SemSamJan_21[order(sig_kept_OeiSamJan_21_x_SemSamJan_21$logFC),]

sig_kept_OeiSamJan_21_x_SemSamJan_21 <-cbind( genData[match(rownames(sig_kept_OeiSamJan_21_x_SemSamJan_21),genData$transcriptName),],sig_kept_OeiSamJan_21_x_SemSamJan_21)

write.table(sig_kept_OeiSamJan_21_x_SemSamJan_21, file = "OeiSamJan_21_x_SemSamJan_21.tab", row.names = T, quote = F, sep = "\t")

sig_kept_OeiSamJan_21_x_SemSamJan_21<- tTags_lrt_OeiSamJan_21_x_SemSamJan_21$table[keep_sig[,1],]

sig_kept_OeiSamJan_21_x_SemSamJan_21 <- sig_kept_OeiSamJan_21_x_SemSamJan_21[order(sig_kept_OeiSamJan_21_x_SemSamJan_21$logFC),]

###############################
#AcaLeafAug_20_x_OeiLeafAug_20

lrt_AcaLeafAug_20_x_OeiLeafAug_20<- glmLRT(fit, contrast = c(0,-1,0,0,0,1,0,0,0,0,0,0))
tTags_lrt_AcaLeafAug_20_x_OeiLeafAug_20 <- topTags(lrt_AcaLeafAug_20_x_OeiLeafAug_20, n= NULL)
keep_sig <- matrix(tTags_lrt_AcaLeafAug_20_x_OeiLeafAug_20$table$FDR <= 0.05 & abs(tTags_lrt_AcaLeafAug_20_x_OeiLeafAug_20$table$logFC) >= 1)

sig_kept_AcaLeafAug_20_x_OeiLeafAug_20<- tTags_lrt_AcaLeafAug_20_x_OeiLeafAug_20$table[keep_sig[,1],]

sig_kept_AcaLeafAug_20_x_OeiLeafAug_20 <- sig_kept_AcaLeafAug_20_x_OeiLeafAug_20[order(sig_kept_AcaLeafAug_20_x_OeiLeafAug_20$logFC),]

sig_kept_AcaLeafAug_20_x_OeiLeafAug_20 <-cbind( genData[match(rownames(sig_kept_AcaLeafAug_20_x_OeiLeafAug_20),genData$transcriptName),],sig_kept_AcaLeafAug_20_x_OeiLeafAug_20)

write.table(sig_kept_AcaLeafAug_20_x_OeiLeafAug_20, file = "AcaLeafAug_20_x_OeiLeafAug_20.tab", row.names = T, quote = F, sep = "\t")

sig_kept_AcaLeafAug_20_x_OeiLeafAug_20<- tTags_lrt_AcaLeafAug_20_x_OeiLeafAug_20$table[keep_sig[,1],]

sig_kept_AcaLeafAug_20_x_OeiLeafAug_20 <- sig_kept_AcaLeafAug_20_x_OeiLeafAug_20[order(sig_kept_AcaLeafAug_20_x_OeiLeafAug_20$logFC),]

###############################
#AcaLeafAug_20_x_SemLeafAug_20

lrt_AcaLeafAug_20_x_SemLeafAug_20<- glmLRT(fit, contrast = c(0,-1,0,0,0,0,0,0,0,1,0,0))
tTags_lrt_AcaLeafAug_20_x_SemLeafAug_20 <- topTags(lrt_AcaLeafAug_20_x_SemLeafAug_20, n= NULL)
keep_sig <- matrix(tTags_lrt_AcaLeafAug_20_x_SemLeafAug_20$table$FDR <= 0.05 & abs(tTags_lrt_AcaLeafAug_20_x_SemLeafAug_20$table$logFC) >= 1)

sig_kept_AcaLeafAug_20_x_SemLeafAug_20<- tTags_lrt_AcaLeafAug_20_x_SemLeafAug_20$table[keep_sig[,1],]

sig_kept_AcaLeafAug_20_x_SemLeafAug_20 <- sig_kept_AcaLeafAug_20_x_SemLeafAug_20[order(sig_kept_AcaLeafAug_20_x_SemLeafAug_20$logFC),]

sig_kept_AcaLeafAug_20_x_SemLeafAug_20 <-cbind( genData[match(rownames(sig_kept_AcaLeafAug_20_x_SemLeafAug_20),genData$transcriptName),],sig_kept_AcaLeafAug_20_x_SemLeafAug_20)

write.table(sig_kept_AcaLeafAug_20_x_SemLeafAug_20, file = "AcaLeafAug_20_x_SemLeafAug_20.tab", row.names = T, quote = F, sep = "\t")

sig_kept_AcaLeafAug_20_x_SemLeafAug_20<- tTags_lrt_AcaLeafAug_20_x_SemLeafAug_20$table[keep_sig[,1],]

sig_kept_AcaLeafAug_20_x_SemLeafAug_20 <- sig_kept_AcaLeafAug_20_x_SemLeafAug_20[order(sig_kept_AcaLeafAug_20_x_SemLeafAug_20$logFC),]

###############################
#OeiLeafAug_20_x_SemLeafAug_20

lrt_OeiLeafAug_20_x_SemLeafAug_20<- glmLRT(fit, contrast = c(0,0,0,0,0,-1,0,0,0,1,0,0))
tTags_lrt_OeiLeafAug_20_x_SemLeafAug_20 <- topTags(lrt_OeiLeafAug_20_x_SemLeafAug_20, n= NULL)
keep_sig <- matrix(tTags_lrt_OeiLeafAug_20_x_SemLeafAug_20$table$FDR <= 0.05 & abs(tTags_lrt_OeiLeafAug_20_x_SemLeafAug_20$table$logFC) >= 1)

sig_kept_OeiLeafAug_20_x_SemLeafAug_20<- tTags_lrt_OeiLeafAug_20_x_SemLeafAug_20$table[keep_sig[,1],]

sig_kept_OeiLeafAug_20_x_SemLeafAug_20 <- sig_kept_OeiLeafAug_20_x_SemLeafAug_20[order(sig_kept_OeiLeafAug_20_x_SemLeafAug_20$logFC),]

sig_kept_OeiLeafAug_20_x_SemLeafAug_20 <-cbind( genData[match(rownames(sig_kept_OeiLeafAug_20_x_SemLeafAug_20),genData$transcriptName),],sig_kept_OeiLeafAug_20_x_SemLeafAug_20)

write.table(sig_kept_OeiLeafAug_20_x_SemLeafAug_20, file = "OeiLeafAug_20_x_SemLeafAug_20.tab", row.names = T, quote = F, sep = "\t")

sig_kept_OeiLeafAug_20_x_SemLeafAug_20<- tTags_lrt_OeiLeafAug_20_x_SemLeafAug_20$table[keep_sig[,1],]

sig_kept_OeiLeafAug_20_x_SemLeafAug_20 <- sig_kept_OeiLeafAug_20_x_SemLeafAug_20[order(sig_kept_OeiLeafAug_20_x_SemLeafAug_20$logFC),]

###############################
#AcaLeafMar_20_x_OeiLeafMar_20

lrt_AcaLeafMar_20_x_OeiLeafMar_20<- glmLRT(fit, contrast = c(0,0,-1,0,0,0,1,0,0,0,0,0))
tTags_lrt_AcaLeafMar_20_x_OeiLeafMar_20 <- topTags(lrt_AcaLeafMar_20_x_OeiLeafMar_20, n= NULL)
keep_sig <- matrix(tTags_lrt_AcaLeafMar_20_x_OeiLeafMar_20$table$FDR <= 0.05 & abs(tTags_lrt_AcaLeafMar_20_x_OeiLeafMar_20$table$logFC) >= 1)

sig_kept_AcaLeafMar_20_x_OeiLeafMar_20<- tTags_lrt_AcaLeafMar_20_x_OeiLeafMar_20$table[keep_sig[,1],]

sig_kept_AcaLeafMar_20_x_OeiLeafMar_20 <- sig_kept_AcaLeafMar_20_x_OeiLeafMar_20[order(sig_kept_AcaLeafMar_20_x_OeiLeafMar_20$logFC),]

sig_kept_AcaLeafMar_20_x_OeiLeafMar_20 <-cbind( genData[match(rownames(sig_kept_AcaLeafMar_20_x_OeiLeafMar_20),genData$transcriptName),],sig_kept_AcaLeafMar_20_x_OeiLeafMar_20)

write.table(sig_kept_AcaLeafMar_20_x_OeiLeafMar_20, file = "AcaLeafMar_20_x_OeiLeafMar_20.tab", row.names = T, quote = F, sep = "\t")

sig_kept_AcaLeafMar_20_x_OeiLeafMar_20<- tTags_lrt_AcaLeafMar_20_x_OeiLeafMar_20$table[keep_sig[,1],]

sig_kept_AcaLeafMar_20_x_OeiLeafMar_20 <- sig_kept_AcaLeafMar_20_x_OeiLeafMar_20[order(sig_kept_AcaLeafMar_20_x_OeiLeafMar_20$logFC),]

###############################
#AcaLeafMar_20_x_SemLeafMar_20

lrt_AcaLeafMar_20_x_SemLeafMar_20<- glmLRT(fit, contrast = c(0,0,-1,0,0,0,0,0,0,0,1,0))
tTags_lrt_AcaLeafMar_20_x_SemLeafMar_20 <- topTags(lrt_AcaLeafMar_20_x_SemLeafMar_20, n= NULL)
keep_sig <- matrix(tTags_lrt_AcaLeafMar_20_x_SemLeafMar_20$table$FDR <= 0.05 & abs(tTags_lrt_AcaLeafMar_20_x_SemLeafMar_20$table$logFC) >= 1)

sig_kept_AcaLeafMar_20_x_SemLeafMar_20<- tTags_lrt_AcaLeafMar_20_x_SemLeafMar_20$table[keep_sig[,1],]

sig_kept_AcaLeafMar_20_x_SemLeafMar_20 <- sig_kept_AcaLeafMar_20_x_SemLeafMar_20[order(sig_kept_AcaLeafMar_20_x_SemLeafMar_20$logFC),]

sig_kept_AcaLeafMar_20_x_SemLeafMar_20 <-cbind( genData[match(rownames(sig_kept_AcaLeafMar_20_x_SemLeafMar_20),genData$transcriptName),],sig_kept_AcaLeafMar_20_x_SemLeafMar_20)

write.table(sig_kept_AcaLeafMar_20_x_SemLeafMar_20, file = "AcaLeafMar_20_x_SemLeafMar_20.tab", row.names = T, quote = F, sep = "\t")

sig_kept_AcaLeafMar_20_x_SemLeafMar_20<- tTags_lrt_AcaLeafMar_20_x_SemLeafMar_20$table[keep_sig[,1],]

sig_kept_AcaLeafMar_20_x_SemLeafMar_20 <- sig_kept_AcaLeafMar_20_x_SemLeafMar_20[order(sig_kept_AcaLeafMar_20_x_SemLeafMar_20$logFC),]

###############################
#OeiLeafMar_20_x_SemLeafMar_20

lrt_OeiLeafMar_20_x_SemLeafMar_20<- glmLRT(fit, contrast = c(0,0,0,0,0,0,-1,0,0,0,1,0))
tTags_lrt_OeiLeafMar_20_x_SemLeafMar_20 <- topTags(lrt_OeiLeafMar_20_x_SemLeafMar_20, n= NULL)
keep_sig <- matrix(tTags_lrt_OeiLeafMar_20_x_SemLeafMar_20$table$FDR <= 0.05 & abs(tTags_lrt_OeiLeafMar_20_x_SemLeafMar_20$table$logFC) >= 1)

sig_kept_OeiLeafMar_20_x_SemLeafMar_20<- tTags_lrt_OeiLeafMar_20_x_SemLeafMar_20$table[keep_sig[,1],]

sig_kept_OeiLeafMar_20_x_SemLeafMar_20 <- sig_kept_OeiLeafMar_20_x_SemLeafMar_20[order(sig_kept_OeiLeafMar_20_x_SemLeafMar_20$logFC),]

sig_kept_OeiLeafMar_20_x_SemLeafMar_20 <-cbind( genData[match(rownames(sig_kept_OeiLeafMar_20_x_SemLeafMar_20),genData$transcriptName),],sig_kept_OeiLeafMar_20_x_SemLeafMar_20)

write.table(sig_kept_OeiLeafMar_20_x_SemLeafMar_20, file = "OeiLeafMar_20_x_SemLeafMar_20.tab", row.names = T, quote = F, sep = "\t")

sig_kept_OeiLeafMar_20_x_SemLeafMar_20<- tTags_lrt_OeiLeafMar_20_x_SemLeafMar_20$table[keep_sig[,1],]

sig_kept_OeiLeafMar_20_x_SemLeafMar_20 <- sig_kept_OeiLeafMar_20_x_SemLeafMar_20[order(sig_kept_OeiLeafMar_20_x_SemLeafMar_20$logFC),]

###############################
#AcaLeafMar_20_x_AcaLeafAug_20

lrt_AcaLeafMar_20_x_AcaLeafAug_20<- glmLRT(fit, contrast = c(0,1,-1,0,0,0,0,0,0,0,0,0))
tTags_lrt_AcaLeafMar_20_x_AcaLeafAug_20 <- topTags(lrt_AcaLeafMar_20_x_AcaLeafAug_20, n= NULL)
keep_sig <- matrix(tTags_lrt_AcaLeafMar_20_x_AcaLeafAug_20$table$FDR <= 0.05 & abs(tTags_lrt_AcaLeafMar_20_x_AcaLeafAug_20$table$logFC) >= 1)

sig_kept_AcaLeafMar_20_x_AcaLeafAug_20<- tTags_lrt_AcaLeafMar_20_x_AcaLeafAug_20$table[keep_sig[,1],]

sig_kept_AcaLeafMar_20_x_AcaLeafAug_20 <- sig_kept_AcaLeafMar_20_x_AcaLeafAug_20[order(sig_kept_AcaLeafMar_20_x_AcaLeafAug_20$logFC),]

sig_kept_AcaLeafMar_20_x_AcaLeafAug_20 <-cbind( genData[match(rownames(sig_kept_AcaLeafMar_20_x_AcaLeafAug_20),genData$transcriptName),],sig_kept_AcaLeafMar_20_x_AcaLeafAug_20)

write.table(sig_kept_AcaLeafMar_20_x_AcaLeafAug_20, file = "AcaLeafMar_20_x_AcaLeafAug_20.tab", row.names = T, quote = F, sep = "\t")

sig_kept_AcaLeafMar_20_x_AcaLeafAug_20<- tTags_lrt_AcaLeafMar_20_x_AcaLeafAug_20$table[keep_sig[,1],]

sig_kept_AcaLeafMar_20_x_AcaLeafAug_20 <- sig_kept_AcaLeafMar_20_x_AcaLeafAug_20[order(sig_kept_AcaLeafMar_20_x_AcaLeafAug_20$logFC),]

###############################
#OeiLeafMar_20_x_OeiLeafAug_20

lrt_OeiLeafMar_20_x_OeiLeafAug_20<- glmLRT(fit, contrast = c(0,0,0,0,0,1,-1,0,0,0,0,0))
tTags_lrt_OeiLeafMar_20_x_OeiLeafAug_20 <- topTags(lrt_OeiLeafMar_20_x_OeiLeafAug_20, n= NULL)
keep_sig <- matrix(tTags_lrt_OeiLeafMar_20_x_OeiLeafAug_20$table$FDR <= 0.05 & abs(tTags_lrt_OeiLeafMar_20_x_OeiLeafAug_20$table$logFC) >= 1)

sig_kept_OeiLeafMar_20_x_OeiLeafAug_20<- tTags_lrt_OeiLeafMar_20_x_OeiLeafAug_20$table[keep_sig[,1],]

sig_kept_OeiLeafMar_20_x_OeiLeafAug_20 <- sig_kept_OeiLeafMar_20_x_OeiLeafAug_20[order(sig_kept_OeiLeafMar_20_x_OeiLeafAug_20$logFC),]

sig_kept_OeiLeafMar_20_x_OeiLeafAug_20 <-cbind( genData[match(rownames(sig_kept_OeiLeafMar_20_x_OeiLeafAug_20),genData$transcriptName),],sig_kept_OeiLeafMar_20_x_OeiLeafAug_20)

write.table(sig_kept_OeiLeafMar_20_x_OeiLeafAug_20, file = "OeiLeafMar_20_x_OeiLeafAug_20.tab", row.names = T, quote = F, sep = "\t")

sig_kept_OeiLeafMar_20_x_OeiLeafAug_20<- tTags_lrt_OeiLeafMar_20_x_OeiLeafAug_20$table[keep_sig[,1],]

sig_kept_OeiLeafMar_20_x_OeiLeafAug_20 <- sig_kept_OeiLeafMar_20_x_OeiLeafAug_20[order(sig_kept_OeiLeafMar_20_x_OeiLeafAug_20$logFC),]

###############################
#SemLeafMar_20_x_SemLeafAug_20

lrt_SemLeafMar_20_x_SemLeafAug_20<- glmLRT(fit, contrast = c(0,0,0,0,0,0,0,0,0,1,-1,0))
tTags_lrt_SemLeafMar_20_x_SemLeafAug_20 <- topTags(lrt_SemLeafMar_20_x_SemLeafAug_20, n= NULL)
keep_sig <- matrix(tTags_lrt_SemLeafMar_20_x_SemLeafAug_20$table$FDR <= 0.05 & abs(tTags_lrt_SemLeafMar_20_x_SemLeafAug_20$table$logFC) >= 1)

sig_kept_SemLeafMar_20_x_SemLeafAug_20<- tTags_lrt_SemLeafMar_20_x_SemLeafAug_20$table[keep_sig[,1],]

sig_kept_SemLeafMar_20_x_SemLeafAug_20 <- sig_kept_SemLeafMar_20_x_SemLeafAug_20[order(sig_kept_SemLeafMar_20_x_SemLeafAug_20$logFC),]

sig_kept_SemLeafMar_20_x_SemLeafAug_20 <-cbind( genData[match(rownames(sig_kept_SemLeafMar_20_x_SemLeafAug_20),genData$transcriptName),],sig_kept_SemLeafMar_20_x_SemLeafAug_20)

write.table(sig_kept_SemLeafMar_20_x_SemLeafAug_20, file = "SemLeafMar_20_x_SemLeafAug_20.tab", row.names = T, quote = F, sep = "\t")

sig_kept_SemLeafMar_20_x_SemLeafAug_20<- tTags_lrt_SemLeafMar_20_x_SemLeafAug_20$table[keep_sig[,1],]

sig_kept_SemLeafMar_20_x_SemLeafAug_20 <- sig_kept_SemLeafMar_20_x_SemLeafAug_20[order(sig_kept_SemLeafMar_20_x_SemLeafAug_20$logFC),]

#########################################
#####GO enrichment############################
##########################################
awk '{if ($2 != "") print $0 }' GO_terms_NCBI_AnnotationOnePerRow.txt > GO.C.a.OnePerRow.txt

for folder in `ls -d */`
do
  echo $folder
  cd $folder/DOWN
  sed -i "" 's/rna-//g' rnaIDs.tab
  for rna in `awk '{print $1}' rnaIDs.tab`
  do
    grep $rna ../../namesEquivalences.tab >> temp
  done
  for prot in `awk '{print $2}' temp`
  do
    grep $prot ../../GO.C.a.OnePerRow.txt | tee -a GOterms.tab
  done

  cd ../UP

  sed -i "" 's/rna-//g' rnaIDs.tab

  for rna in `awk '{print $1}' rnaIDs.tab`
  do
    grep $rna ../../namesEquivalences.tab >> temp
  done
  for prot in `awk '{print $2}' temp`
  do
    grep $prot ../../GO.C.a.OnePerRow.txt | tee -a GOterms.tab
  done
  cd ../../
done
