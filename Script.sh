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
################# Alignment Stats ###############
################################################


alignmentStats <- read.delim("alignmentStats.tab",row.names=1)

alignmentStats$effectiveCounts <- fit$samples$lib.size

cols <- c("#4859c2", "#6b38b8","#39abad","#40b84a","#b06042")

meanMatrix <- apply(alignmentStats,2,mean)
sdMatrix <- apply(alignmentStats,2,sd)


pdf("AlignmentData.pdf",h=10,w=17)
barplot <- barplot(meanMatrix,beside=T,col=cols, yaxt='n',xaxt='n',ylim=c(0,max(meanMatrix +sdMatrix)*1.164),main ="Mean Number of Fragments Per Sample",cex.main=2.3,width=1.2,space=c(.2,1))

title(ylab="Millions of RNA fragments",cex.lab=1.5,line=2.7)

axis(2,seq(0,60000000,5000000),labels=F)
mtext(side=2,text=seq(0,60,5),outer=F,las=2,line=.8,at=seq(0,60000000,5000000),cex=1.5)
mtext(side=1,text=c("Sequenced And QC","Uniquelly Mapped","Multi Mappers","Unmapped","Effective Counts"),outer=F,line=2,at=barplot,cex=2)

arrows(x0 = barplot, y0 = meanMatrix - sdMatrix, x1 = barplot, y1=meanMatrix + sdMatrix ,code=3,angle=90,length=0.05,col="black",lwd=3.5)
dev.off()

system("open AlignmentData.pdf")

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


################
#General Plots##
################

sig <- list()

sig <- ls(pattern="sig_kept")

dedupKeepedNames <- c()

for (DE in sig){
dedupKeepedNames <- c(dedupKeepedNames,rownames(get(DE)))

}

dedupKeepedNames <- unique(dedupKeepedNames)

#awk -F"\t" '$3=="exon"{ID=substr($9, length($9)-14, 37); L[ID]+=$5-$4+1}END{for(i in L){print i"\t"L[i]}}' compatible.gff > processed.geneLenght.txt

#cat processed.geneLenght.txt | sed 's/=/rna-/g' > geneLenths.txt

genes.len <- read.table("geneLenths.txt",header = F)
analysis_matrix$genes <- genes.len
colnames(analysis_matrix$genes) <- c("names", "lengths")

DE.RPKM.matrix <- rpkm(analysis_matrix,normalized.lib.sizes=T)

DE.RPKM.matrix <- rpkm(analysis_matrix[intersect(dedupKeepedNames, rownames(DE.RPKM.matrix)),],normalized.lib.sizes=T)



OeiFBudMar_20.means <- as.matrix(rowMeans(DE.RPKM.matrix[,c(1:4)]))
OeiLeafMar_20.means <- as.matrix(rowMeans(DE.RPKM.matrix[,c(5:8)]))
AcaFBudMar_20.means <- as.matrix(rowMeans(DE.RPKM.matrix[,c(9:12)]))
AcaLeafMar_20.means <- as.matrix(rowMeans(DE.RPKM.matrix[,c(13:16)]))
SemFBudMar_20.means <- as.matrix(rowMeans(DE.RPKM.matrix[,c(17:20)]))
SemLeafMar_20.means <- as.matrix(rowMeans(DE.RPKM.matrix[,c(21:24)]))
OeiLeafAug_20.means <- as.matrix(rowMeans(DE.RPKM.matrix[,c(25:28)]))
AcaLeafAug_20.means <- as.matrix(rowMeans(DE.RPKM.matrix[,c(29:32)]))
SemLeafAug_20.means <- as.matrix(rowMeans(DE.RPKM.matrix[,c(33:36)]))
OeiSamJan_21.means <- as.matrix(rowMeans(DE.RPKM.matrix[,c(37:40)]))
AcaSamJan_21.means <- as.matrix(rowMeans(DE.RPKM.matrix[,c(41:44)]))
SemSamJan_21.means <- as.matrix(rowMeans(DE.RPKM.matrix[,c(45:48)]))

per_condition_means <- cbind(AcaFBudMar_20.means,
SemFBudMar_20.means,
OeiFBudMar_20.means,
AcaSamJan_21.means,
SemSamJan_21.means,
OeiSamJan_21.means,
AcaLeafMar_20.means,
SemLeafMar_20.means,
OeiLeafMar_20.means,
AcaLeafAug_20.means,
SemLeafAug_20.means,
OeiLeafAug_20.means)

colnames(per_condition_means) <- c("AcaFBudMar_20",
"SemFBudMar_20",
"OeiFBudMar_20",
"AcaSamJan_21",
"SemSamJan_21",
"OeiSamJan_21",
"AcaLeafMar_20",
"SemLeafMar_20",
"OeiLeafMar_20",
"AcaLeafAug_20",
"SemLeafAug_20",
"OeiLeafAug_20")


blue.red <- colorRampPalette(c("white","#d93e23"))(8)


library(gplots)
pdf("heatmap_log_absolute_Expression_DE.pdf", he=10,wi=15)
heatmap.2(log(per_condition_means+1,10),dendrogram = "both",trace = "none",key = T, Colv=T,Rowv=T,margins = c(5, 2),offsetRow =100, key.xlab="Log 10 RPKM + 1", key.ylab="", keysize =.9 , key.title="", cexRow=1, cexCol=1,srtCol=20,adjCol=c(0.5,1.5), offsetCol = .3, col=blue.red, denscol="green", densadj = 5,colsep=c(0:ncol(per_condition_means)),sepcolor="black",sepwidth=c(0.01, 0.01))
dev.off()

system("open heatmap_log_absolute_Expression_DE.pdf")


OeiFBudMar_20.ste <- as.matrix(apply(DE.RPKM.matrix[,c(1:4)],1,sd)/sqrt(4))
OeiLeafMar_20.ste <- as.matrix(apply(DE.RPKM.matrix[,c(5:8)],1,sd)/sqrt(4))
AcaFBudMar_20.ste <- as.matrix(apply(DE.RPKM.matrix[,c(9:12)],1,sd)/sqrt(4))
AcaLeafMar_20.ste <- as.matrix(apply(DE.RPKM.matrix[,c(13:16)],1,sd)/sqrt(4))
SemFBudMar_20.ste <- as.matrix(apply(DE.RPKM.matrix[,c(17:20)],1,sd)/sqrt(4))
SemLeafMar_20.ste <- as.matrix(apply(DE.RPKM.matrix[,c(21:24)],1,sd)/sqrt(4))
OeiLeafAug_20.ste <- as.matrix(apply(DE.RPKM.matrix[,c(25:28)],1,sd)/sqrt(4))
AcaLeafAug_20.ste <- as.matrix(apply(DE.RPKM.matrix[,c(29:32)],1,sd)/sqrt(4))
SemLeafAug_20.ste <- as.matrix(apply(DE.RPKM.matrix[,c(33:36)],1,sd)/sqrt(4))
OeiSamJan_21.ste <- as.matrix(apply(DE.RPKM.matrix[,c(37:40)],1,sd)/sqrt(4))
AcaSamJan_21.ste <- as.matrix(apply(DE.RPKM.matrix[,c(41:44)],1,sd)/sqrt(4))
SemSamJan_21.ste <- as.matrix(apply(DE.RPKM.matrix[,c(45:48)],1,sd)/sqrt(4))

per_condition_sde <-

cbind(AcaFBudMar_20.ste,
SemFBudMar_20.ste,
OeiFBudMar_20.ste,
AcaSamJan_21.ste,
SemSamJan_21.ste,
OeiSamJan_21.ste,
AcaLeafMar_20.ste,
SemLeafMar_20.ste,
OeiLeafMar_20.ste,
AcaLeafAug_20.ste,
SemLeafAug_20.ste,
OeiLeafAug_20.ste)

colnames(per_condition_sde) <- c("AcaFBudMar_20",
"SemFBudMar_20",
"OeiFBudMar_20",
"AcaSamJan_21",
"SemSamJan_21",
"OeiSamJan_21",
"AcaLeafMar_20",
"SemLeafMar_20",
"OeiLeafMar_20",
"AcaLeafAug_20",
"SemLeafAug_20",
"OeiLeafAug_20")


cols <- c(rep("#1c5414"),rep("#731c23"))



cols <- colorRampPalette(c("#c24a3a","#522c28","#7bc437"))(12)


i=10
system("mkdir barplot")
for (gene in rownames(per_condition_means)){
print(gene)
limitOfY <- max(per_condition_means[gene,]+per_condition_sde[gene,]/sqrt(1))*1.20
pdf(paste0("./barplot/barplot",gene,".pdf"),h=5,w=20)
barplot <- barplot(per_condition_means[gene,],beside=T,col=
, yaxt='n',xaxt='n',ylim=c(0,limitOfY),main =paste("Expression of gene",gene),cex.main=1,width=1.2)

title(ylab="FPKM",cex.lab=1.5,line=2.3)

#legend("topleft",legend=colnames(per_condition_sde),col = cols,pch=19,cex=1.5)


if(limitOfY <= 10 ){
i=1
}
if(limitOfY > 10 & limitOfY <= 50 ){
i=3
}
if( limitOfY > 50 & limitOfY <= 100 ){
i=10
}
if( limitOfY > 100 & limitOfY <= 1000 ){
i=50
}
if( limitOfY > 1000 & limitOfY <= 10000 ){
i=100
}
if( limitOfY > 10000 & limitOfY <= 100000 ){
i=1000
}

axis(2,seq(0,limitOfY,i),labels=F)
mtext(side=2,text=seq(0,limitOfY,i),outer=F,las=2,line=.8,at=seq(0,limitOfY,i),cex=.8)


mtext(side=1,text=colnames(per_condition_sde),outer=F,line=1.2,at=barplot,cex=1,srt=30)

arrows(x0 = barplot, y0 = per_condition_means[gene,] - per_condition_sde[gene,], x1 = barplot, y1=per_condition_means[gene,] + per_condition_sde[gene,] ,code=3,angle=90,length=0.05,col="black",lwd=1.3)

dev.off()

}

system("mkdir volcanoes")
#l <- list()

l <- ls(pattern="tTags_lrt_")

#Volcano plots
for (comp in l){
print(comp)
FC <-  get(comp)$table$logFC
names(FC) <- rownames(get(comp)$table)
print(max(FC))

FDR <- get(comp)$table$FDR
DATA <- cbind(FC,FDR)
print(-log(min(FDR),10))

pdf(paste0("./volcanoes/",comp,".pdf"),h=5,w=5)
par(mar = c(4.1, 4.1, 4.1, 2.1))

plot(-log(DATA[,2],10) ~ DATA[,1],xlim=c(-16,16),ylim=c(0,150),cex=.1,col="grey",yaxt='n',xaxt='n',ylab="",xlab="")

title(main=gsub("tTags_lrt_","" ,comp),ylab="-log10(FDR)",cex.lab=.9,line=2.3,xlab="log2FC")

legend("topleft",legend=c("Up-regulated","Not DE","Down-regulated"),col = c("#d93e23","grey","#236cd9" ),pch=19,cex=.6)

axis(2,seq(0,150,10),labels=F)
mtext(side=2,text=seq(0,150,10),outer=F,las=2,line=.8,at=seq(0,150,10),cex=.8)

text(x=-12.5,y=-log(0.0000001,10)-.6,labels="FDR = 0.05",las=2, col = "#0e4f34",cex=.9)

axis(1,seq(-16,16,4),labels=F)
mtext(side=1,text=seq(-16,16,4),outer=F,las=1,line=.8,at=seq(-16,16,4),cex=.8)

axis(1,0,labels=F)
mtext(side=1,text=0,outer=F,las=1,line=.8,at=0,cex=.8)

points(-log(DATA[,2][DATA[,2] < 0.05 & DATA[,1] >= 1],10) ~DATA[,1][DATA[,2] < 0.05 & DATA[,1] >= 1], col = "#d93e23",cex=.1)

points(-log(DATA[,2][DATA[,2] < 0.05 & DATA[,1] <= -1],10) ~ DATA[,1][DATA[,2] < 0.05 & DATA[,1] <= -1], col = "#236cd9",cex=.1)

abline(v= 1, col = "#0e4f34",lty = 2)
abline(v= -1, col = "#0e4f34",lty = 2)
abline(h= -log(0.05,10), col = "#0e4f34",lty = 2)

#lightBlue "#8ce4ff"


#text(-log(tail(DATA[order(DATA[,1]),2]),10) ~ tail(DATA[order(DATA[,1]),1]),labels= names(tail(DATA[order(DATA[,1]),1])),cex=.25,pos=1)

dev.off()
}




system("mkdir MA")

for (comp in l){
print(comp)
FC <-  get(comp)$table$logFC
names(FC) <- rownames(get(comp)$table)

FDR <- get(comp)$table$FDR
lCPM <- get(comp)$table$logCPM
DATA <- cbind(FC,lCPM,FDR)

print(max(lCPM))


pdf(paste0("./MA/",comp,".pdf"),h=5,w=5)
par(mar = c(4.1, 4.1, 4.1, 2.1))

plot(DATA[,1] ~ DATA[,2] ,xlim=c(0,12),ylim=c(-16,16),cex=.1,col="grey",yaxt='n',xaxt='n',ylab="",xlab="")

title(main=gsub("tTags_lrt_","" ,comp),ylab="log2FC",cex.lab=.9,line=2.3,xlab="logCounts (CPM)")

#legend("topleft",legend=c("Up-regulated","Not DE","Down-regulated"),col = c("#d93e23","grey","#236cd9" ),pch=19,cex=.6)

axis(1,seq(0,12,2),labels=F)
mtext(side=1,text=seq(0,12,2),outer=F,las=1,line=.8,at=seq(0,12,2),cex=.8)

#text(x=-8.5,y=-log(0.05,10) +1.5,labels="FDR = 0.05",las=2, col = "#0e4f34",cex=.9)

axis(2,seq(-16,16,2),labels=F)
mtext(side=2,text=seq(-16,16,2),outer=F,las=1,line=.8,at=seq(-16,16,2),cex=.8)

#axis(1,0,labels=F)
#mtext(side=1,text=0,outer=F,las=1,line=.8,at=0,cex=.8)

points(DATA[,1][DATA[,3] < 0.05 & DATA[,1] >= 1] ~DATA[,2][DATA[,3] < 0.05 & DATA[,1] >= 1], col = "#d93e23",cex=.1)

points(DATA[,1][DATA[,3] < 0.05 & DATA[,1] <= -1] ~DATA[,2][DATA[,3] < 0.05 & DATA[,1] <= -1], col = "#236cd9",cex=.1)

#abline(v= 1, col = "#0e4f34",lty = 2)
#abline(v= -1, col = "#0e4f34",lty = 2)
#abline(h= -log(0.05,10), col = "#0e4f34",lty = 2)

#lightBlue "#8ce4ff"

#DOWN NAMES
#text(DATA[DownInside,1] ~ DATA[DownInside,2],labels= DownInside,cex=.25,pos=1)

dev.off()
}


#######################
#####Heatmap selected##
#######################

#DE.RPKM.matrix

#FloralIdentityLeaf.tab

FloralIdentityLeaves <- read.table("FloralIdentityLeaf.tab", sep = "\t")

FloralIdentityLeavesDescription <- genData[ genData$transcriptName %in% FloralIdentityLeaves$V1,]

leavesRPKM <- DE.RPKM.matrix[,c(grep("Leaf",colnames(DE.RPKM.matrix)))]

selectedGenes <- leavesRPKM[rownames(leavesRPKM) %in% FloralIdentityLeaves$V1,]

rownames(selectedGenes) <- paste(FloralIdentityLeavesDescription[match(rownames(selectedGenes),FloralIdentityLeavesDescription$transcriptName),1],FloralIdentityLeavesDescription[match(rownames(selectedGenes),FloralIdentityLeavesDescription$transcriptName),3], sep = "-")


rownames(selectedGenes) <- gsub( "rna-","",rownames(selectedGenes))

LogSelectedGenes <- log(selectedGenes+1,10)

LogSelectedGenes <- LogSelectedGenes[rowMeans(LogSelectedGenes) >= 2,]


pdf("FloralIdentityLeaf.pdf",h=8.08,w=8.08)
par(cex.main=1)
heatmap.2(LogSelectedGenes,
trace = "none",
key = T,
margins = c(10, 15),
offsetRow = .001,
 key.xlab="Log(10) FPKM + 1",
  key.ylab="",
keysize =1,
key.title="",
  cexRow=.457,
  srtRow=330,
    col=myCol,
    main="Floral Identity Genes in Leaves",
    density.info="density",
    cexCol=.7)
dev.off()


#FloralIdentitySAMandFBud.tab

myCol <- colorRampPalette(c("white","blue", "red"))(32)

FloralIdentitySAMandFBud <- read.table("FloralIdentitySAMandFBud.tab", sep = "\t")

FloralIdentitySAMandFBudDescription <- genData[ genData$transcriptName %in% FloralIdentitySAMandFBud$V1,]

FloralIdentitySAMandFBudRPKM <- DE.RPKM.matrix[,c(grep("FBud|Sam",colnames(DE.RPKM.matrix)))]

selectedGenes <- FloralIdentitySAMandFBudRPKM[rownames(FloralIdentitySAMandFBudRPKM) %in% FloralIdentitySAMandFBudDescription$transcriptName,]

rownames(selectedGenes) <- paste(FloralIdentitySAMandFBudDescription[match(rownames(selectedGenes),FloralIdentitySAMandFBudDescription$transcriptName),1],FloralIdentitySAMandFBudDescription[match(rownames(selectedGenes),FloralIdentitySAMandFBudDescription$transcriptName),3], sep = "-")

rownames(selectedGenes) <- gsub( "rna-","",rownames(selectedGenes))

LogSelectedGenes <- log(selectedGenes+1,10)

LogSelectedGenes <- LogSelectedGenes[rowMeans(LogSelectedGenes) >= 2,]

##########################
#FloralIdentitySAMandFBud#
##########################

pdf("FloralIdentitySAMandFBud.pdf",h=8.08,w=8.08)
par(cex.main=1)
heatmap.2(LogSelectedGenes,
trace = "none",
key = T,
margins = c(10, 15),
offsetRow = .001,
 key.xlab="Log(10) FPKM + 1",
  key.ylab="",
keysize =1,
key.title="",
  cexRow=.457,
  srtRow=330,
    col=myCol,
    main="Floral Identity Genes in SAM and FBuds",
    density.info="density",
    cexCol=.7)
dev.off()

###################
#EthyleneLeaf.tab#
#################

EthyleneLeaves <- read.table("EthyleneLeaf.tab", sep = "\t")

EthyleneLeavesDescription <- genData[ genData$transcriptName %in% EthyleneLeaves$V1,]

EthyleneleavesRPKM <- DE.RPKM.matrix[,c(grep("Leaf",colnames(DE.RPKM.matrix)))]

selectedGenes <- leavesRPKM[rownames(EthyleneleavesRPKM) %in% EthyleneLeaves$V1,]

rownames(selectedGenes) <- paste(EthyleneLeavesDescription[match(rownames(selectedGenes),EthyleneLeavesDescription$transcriptName),1],EthyleneLeavesDescription[match(rownames(selectedGenes),EthyleneLeavesDescription$transcriptName),3], sep = "-")

rownames(selectedGenes) <- gsub( "rna-","",rownames(selectedGenes))

LogSelectedGenes <- log(selectedGenes+1,10)

LogSelectedGenes <- LogSelectedGenes[rowMeans(LogSelectedGenes) >= 2,]


pdf("EthyleneLeaf.pdf",h=8.08,w=8.08)
par(cex.main=1)
heatmap.2(LogSelectedGenes,
trace = "none",
key = T,
margins = c(10, 15),
offsetRow = .001,
 key.xlab="Log(10) FPKM + 1",
  key.ylab="",
keysize =1,
key.title="",
  cexRow=.457,
  srtRow=330,
    col=myCol,
    main="Ethylene Genes in Leaves",
    density.info="density",
    cexCol=.7)
dev.off()


#########################
#EthyleneSAMandFBud.tab#
#######################

EthyleneSAMandFBud <- read.table("EthyleneSAMandFBud.tab", sep = "\t")

EthyleneSAMandFBudDescription <- genData[ genData$transcriptName %in% EthyleneSAMandFBud$V1,]

EthyleneSAMandFBudRPKM <- DE.RPKM.matrix[,c(grep("FBud|Sam",colnames(DE.RPKM.matrix)))]

selectedGenes <- SAMandFBudRPKM[rownames(EthyleneSAMandFBudRPKM) %in% EthyleneSAMandFBud$V1,]

rownames(selectedGenes) <- paste(EthyleneSAMandFBudDescription[match(rownames(selectedGenes),EthyleneSAMandFBudDescription$transcriptName),1],EthyleneSAMandFBudDescription[match(rownames(selectedGenes),EthyleneSAMandFBudDescription$transcriptName),3], sep = "-")

rownames(selectedGenes) <- gsub( "rna-","",rownames(selectedGenes))

LogSelectedGenes <- log(selectedGenes+1,10)

LogSelectedGenes <- LogSelectedGenes[rowMeans(LogSelectedGenes) >= 2,]


pdf("EthyleneSAMandFBud.pdf",h=8.08,w=8.08)
par(cex.main=1)
heatmap.2(LogSelectedGenes,
trace = "none",
key = T,
margins = c(10, 15),
offsetRow = .001,
 key.xlab="Log(10) FPKM + 1",
  key.ylab="",
keysize =1,
key.title="",
  cexRow=.65,
  srtRow=330,
    col=myCol,
    main="Ethylene Genes in SAM and FBuds",
    density.info="density",
    cexCol=.6)
dev.off()




########################
#CarbohydratesLeaf.tab#
######################

CarbohydratesLeaves <- read.table("CarbohydratesLeaf.tab", sep = "\t")

CarbohydratesLeavesDescription <- genData[ genData$transcriptName %in% CarbohydratesLeaves$V1,]

CarbohydratesleavesRPKM <- DE.RPKM.matrix[,c(grep("Leaf",colnames(DE.RPKM.matrix)))]

selectedGenes <- leavesRPKM[rownames(CarbohydratesleavesRPKM) %in% CarbohydratesLeaves$V1,]

rownames(selectedGenes) <- paste(CarbohydratesLeavesDescription[match(rownames(selectedGenes),CarbohydratesLeavesDescription$transcriptName),1],CarbohydratesLeavesDescription[match(rownames(selectedGenes),CarbohydratesLeavesDescription$transcriptName),3], sep = "-")

rownames(selectedGenes) <- gsub( "rna-","",rownames(selectedGenes))

LogSelectedGenes <- log(selectedGenes+1,10)

LogSelectedGenes <- LogSelectedGenes[rowMeans(LogSelectedGenes) >= 2,]


pdf("CarbohydratesLeaf.pdf",h=8.08,w=8.08)
par(cex.main=1)
heatmap.2(LogSelectedGenes,
trace = "none",
key = T,
margins = c(10, 15),
offsetRow = .001,
 key.xlab="Log(10) FPKM + 1",
  key.ylab="",
keysize =1,
key.title="",
  cexRow=.457,
  srtRow=330,
    col=myCol,
    main="Carbohydrate Genes in Leaves",
    density.info="density",
    cexCol=.7)
dev.off()





########################
#CarbohydratesSAMandFBud.tab#
######################

CarbohydratesSAMandFBud <- read.table("CarbohydratesSAMandFBud.tab", sep = "\t")

CarbohydratesSAMandFBudDescription <- genData[ genData$transcriptName %in% CarbohydratesSAMandFBud$V1,]

CarbohydratesSAMandFBudRPKM <- DE.RPKM.matrix[,c(grep("FBud|Sam",colnames(DE.RPKM.matrix)))]

selectedGenes <- SAMandFBudRPKM[rownames(CarbohydratesSAMandFBudRPKM) %in% CarbohydratesSAMandFBud$V1,]

rownames(selectedGenes) <- paste(CarbohydratesSAMandFBudDescription[match(rownames(selectedGenes),CarbohydratesSAMandFBudDescription$transcriptName),1],CarbohydratesSAMandFBudDescription[match(rownames(selectedGenes),CarbohydratesSAMandFBudDescription$transcriptName),3], sep = "-")

rownames(selectedGenes) <- gsub( "rna-","",rownames(selectedGenes))

LogSelectedGenes <- log(selectedGenes+1,10)

LogSelectedGenes <- LogSelectedGenes[rowMeans(LogSelectedGenes) >= 2,]


pdf("CarbohydratesSAMandFBud.pdf",h=8.08,w=8.08)
par(cex.main=1)
heatmap.2(LogSelectedGenes,
trace = "none",
key = T,
margins = c(10, 15),
offsetRow = .001,
 key.xlab="Log(10) FPKM + 1",
  key.ylab="",
keysize =1,
key.title="",
  cexRow=.457,
  srtRow=330,
    col=myCol,
    main="Carbohydrate Genes in SAM and FBuds",
    density.info="density",
    cexCol=.6)
dev.off()


########################
#RealFloralIdentitySAMandFBud.tab#
######################

RealFloralIdentitySAMandFBud <- read.table("RealFloralIdentitySAMandFBud.tab", sep = "\t")

RealFloralIdentitySAMandFBudDescription <- genData[ genData$transcriptName %in% RealFloralIdentitySAMandFBud$V1,]

RealFloralIdentitySAMandFBudRPKM <- DE.RPKM.matrix[,c(grep("FBud|Sam",colnames(DE.RPKM.matrix)))]

selectedGenes <- SAMandFBudRPKM[rownames(RealFloralIdentitySAMandFBudRPKM) %in% RealFloralIdentitySAMandFBud$V1,]

rownames(selectedGenes) <- paste(RealFloralIdentitySAMandFBudDescription[match(rownames(selectedGenes),RealFloralIdentitySAMandFBudDescription$transcriptName),1],RealFloralIdentitySAMandFBudDescription[match(rownames(selectedGenes),RealFloralIdentitySAMandFBudDescription$transcriptName),3], sep = "-")

rownames(selectedGenes) <- gsub( "rna-","",rownames(selectedGenes))

LogSelectedGenes <- log(selectedGenes+1,10)

LogSelectedGenes <- LogSelectedGenes[rowMeans(LogSelectedGenes) >= 1,]


pdf("RealFloralIdentitySAMandFBud.pdf",h=8.08,w=8.08)
par(cex.main=1)
heatmap.2(LogSelectedGenes,
trace = "none",
key = T,
margins = c(10, 15),
offsetRow = .001,
 key.xlab="Log(10) FPKM + 1",
  key.ylab="",
keysize =1,
key.title="",
  cexRow=.457,
  srtRow=330,
    col=myCol,
    main="Floral Identity Genes in SAM and FBuds",
    density.info="density",
    cexCol=.6)
dev.off()

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
