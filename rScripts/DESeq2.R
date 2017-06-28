###################################################
## DESeq2 differential bound expression analysis ##
###################################################

library('DESeq2')

path = '/ps/imt/e/20141009_AG_Bauer_peeyush_re_analysis/further_analysis/results/RNAseq/NT2D1+-RA+KO/count'
count_data = read.csv(file = paste(path,'/NT2D1_PRMT6_KO+-ATRA_for_analysis.txt', 
                                      sep = ''), header = TRUE, row.names=1, sep = '\t')
pheno_data = read.csv(file = paste(path,'/pheno_-RA.txt',
                                   sep = ''), header = TRUE, row.names=1, sep = '\t')
print(colnames(count_data))
count_data = count_data[,-c(1)]

## Remove additional samples from calculation
print(colnames(count_data))
count_data = count_data[,c(1,3,4,5)]
print(colnames(count_data))

count_data <- count_data[ rowSums(count_data) >= 10, ]
count_data[count_data == 0] <- 1

##Matrix input
dds <- DESeqDataSetFromMatrix(countData = count_data,
                              colData = pheno_data,
                              design = ~ condition)
## ----prefilter-----------------------------------------------------------
#dds <- dds[ rowSums(counts(dds)) > 1, ]
#dds <- dds[ rowSums(counts(dds)) > 9, ]

## ----relevel-------------------------------------------------------------
dds$condition <- relevel(dds$condition, ref="CT")

## ----deseq---------------------------------------------------------------
dds <- DESeq(dds)
res <- results(dds)
## ----MA, fig.width=4.5, fig.height=4.5-----------------------------------
plotMA(res, main="MA Plot", ylim=c(-4,4))#, col = ifelse(res$log2FoldChange>=2 | res$log2FoldChange<=-2, "red3", "gray32"))

res <- as.data.frame(res)

## ----normalize counts----------------------------------------------------
normalizedcounts = t(t(counts(dds))/sizeFactors(dds))
normalizedcounts = as.data.frame(normalizedcounts)
write.table(normalizedcounts,
            file=paste(outpath,"/Norm_tag_counts.txt",sep = ''), sep = '\t')


colname = colnames(normalizedcounts)
length(colname)
for(i in colname){
  print(i)
  res[i] = normalizedcounts[i]
}
#res = res[res$log2FoldChange >= 1 | res$log2FoldChange <= -1,]
#res = res[res$pvalue <= 0.05,]

## ----resOrder------------------------------------------------------------
resOrdered <- res[order(res$padj),]

## ----sumRes--------------------------------------------------------------
summary(res)

## ----sumRes01------------------------------------------------------------
sum(res$padj < 0.1, na.rm=TRUE)

## ----export, eval=FALSE--------------------------------------------------
outpath = '/ps/imt/e/20141009_AG_Bauer_peeyush_re_analysis/further_analysis/results/RNAseq/NT2D1_KO_-ATRA/DESeq2'
write.table(res,
           file=paste(outpath,"/Deseq2_PRMT6_KO_-RA_removed_CT2_KO3_full.tsv",sep = ''), sep = '\t')

## Including factors

ddsMF<-dds

## ----replaceDesign-------------------------------------------------------
design(ddsMF) <- formula(~ day + condition)
ddsMF <- DESeq(ddsMF)

## ----multiResults--------------------------------------------------------
resMF <- results(ddsMF)
head(resMF)
 
## ----rlogAndVST----------------------------------------------------------
rld <- rlog(dds)
## ----figPCA, dev="pdf", fig.width=5, fig.height=3------------------------
plotPCA(rld, intgroup=c("condition", "sample"), ntop = 1000)
data = plotPCA(rld, intgroup=c("condition", "sample"), ntop = 2000, returnData = TRUE)
percentVar = round(100*attr(pcData, "percentVar"))



################################################
## DESeq analysis for differential expression ##
################################################

library('DESeq')

path = '/ps/imt/e/20141009_AG_Bauer_peeyush_re_analysis/further_analysis/results/RNAseq/NT2D1_KO_-ATRA/DESeq'
outpath='/ps/imt/e/20141009_AG_Bauer_peeyush_re_analysis/further_analysis/results/RNAseq/NT2D1_KO_-ATRA/DESeq'
countData = read.csv(file = paste(path,'/TopHat2_NT2D1_PRMT6_KO_E9vsB6.txt', 
                                      sep = ''), header = TRUE, row.names=1, sep = '\t')
phenodata = read.csv(file = paste(path,'/pheno.txt', 
                                   sep = ''), header = TRUE, row.names=1, sep = '\t')
countData = countData[,-c(1,2,3,4,5)]
countData <- countData[ rowSums(countData) >= 10, ]
countData[countData == 0] <- 1
head(countData)
colnames(countData)
cds = newCountDataSet(countData, phenodata$condition)
cds = estimateSizeFactors(cds)
sizeFactors(cds)
normalizedCounts <- t( t(counts(cds)) / sizeFactors(cds) )
head(normalizedCounts)
head( counts( cds, normalized=TRUE ) )
normalizedCounts = counts( cds, normalized=TRUE )

cds = estimateDispersions(cds)
plotDispEsts(cds)
head( fData(cds) )

vsdFull = varianceStabilizingTransformation( cds )
print(plotPCA(vsdFull, intgroup=c("condition", "sample")))

res = nbinomTest(cds, 'GFP', 'KO')
res1 = res
res = res1
plotMA(res, main="MA Plot", ylim=c(-8,8), col = ifelse(res$log2FoldChange>=2 | res$log2FoldChange<=-2, "red3", "gray32"))
colname = colnames(normalizedCounts)
length(colname)
for(i in 1:length(colname)){
  print(i)
  res[colname[i]] = normalizedCounts[,i]
}
res = res[res$log2FoldChange >= 1 | res$log2FoldChange <= -1,]
write.table(as.data.frame(res),
            file=paste(outpath,"/TopHat2_DESeq_FC_NT2D1_PRMT6_KO_results.txt",sep = ''), sep = '\t')

####### HeatMap ##############
res1 = res
res = res1
res = res[res$log2FoldChange >= 1.5 | res$log2FoldChange <= -1.5,]
res = res[res$padj <= 0.05,]

colnames(res)
res = res[order(res$log2FoldChange, decreasing = TRUE),]
res = res[order(res$HL60_2_10_1, decreasing = TRUE),]
df = res[, -c(1,2,3,4,5,6,7,8)]
df = log(df,2)
max(df)
mat = as.matrix(df)
my_palette <- colorRampPalette(c("white","wheat","ightcoral"))(n = 299)
# (optional) defines the color breaks manually for a "skewed" color transition

col_breaks = c(seq(0,5,length=100),
               seq(5.1,9.9,length=100),
               seq(10,15,length=100))

########## creating images ############
map_path = paste(outpath,'/heatmap.png', sep = "")
print(map_path)
png(map_path, 
    width = 12*300,        # 5 x 300 pixels
    height = 8*300,
    res = 300,            # 300 pixels per inch
    pointsize = 8)
#svg(path, width = 10, height = 8)
#print(path)
################################

heatmap.2(mat, 
          #cellnote = mat,  # same data set for cell labels
          main = "DEGenes log2(count)", # heat map title
          #xlab = "Differential tag densities",
          #ylab = "Peaks",
          labRow = "",
          labCol = colnames(df),
          #cexRow = 0.7,           # Changes the size of col and row font size
          cexCol = 1,
          notecol="black",      # change font color of cell labels to black
          density.info="none",  # turns off density plot inside color legend
          trace="none",         # turns off trace lines inside the heat map
          margins =c(6,3),     # widens margins around plot
          col=my_palette,       # use on color palette defined earlier
          breaks=col_breaks,    # enable color transition at specified limits
          dendrogram="none",     # only draw a row dendrogram
          #Colv= TRUE,           #cluster column
          #Rowv = TRUE,
          keysize = 0.5,
          na.color = 'black'
)
dev.off()

         