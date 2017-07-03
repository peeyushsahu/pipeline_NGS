library('gplots')

setwd('/ps/imt/e/HL60_Christene/further_analysis/RNA_seq/cuffDiff')
cuff<-readCufflinks()
cuff
gene.diff<-diffData(genes(cuff))
head(gene.diff)

gene_id sample_1 sample_2 status  value_1  value_2 log2_fold_change test_stat p_value  q_value significant
1 XLOC_000001 SKI_EGFP   SKI_KO NOTEST 0.000000 0.000000         0.000000  0.000000 1.00000 1.000000          no
2 XLOC_000002 SKI_EGFP   SKI_KO NOTEST 0.000000 0.000000         0.000000  0.000000 1.00000 1.000000          no
3 XLOC_000003 SKI_EGFP   SKI_KO NOTEST 0.000000 0.000000         0.000000  0.000000 1.00000 1.000000          no
4 XLOC_000004 SKI_EGFP   SKI_KO NOTEST 0.000000 0.000000         0.000000  0.000000 1.00000 1.000000          no
5 XLOC_000005 SKI_EGFP   SKI_KO NOTEST 0.000000 0.000000         0.000000  0.000000 1.00000 1.000000          no
6 XLOC_000006 SKI_EGFP   SKI_KO     OK 0.830004 0.646382        -0.360731 -0.796726 0.17045 0.577002          no

gene.matrix<-fpkmMatrix(genes(cuff))
head(gene.matrix)

gene.rep.matrix<-repFpkmMatrix(genes(cuff))
head(gene.rep.matrix)

SKI_EGFP_0 SKI_EGFP_1 SKI_EGFP_2 SKI_KO_0 SKI_KO_1 SKI_KO_2
XLOC_000001   0.000000   0.000000   0.000000 0.000000 0.000000 0.000000
XLOC_000002   0.000000   0.000000   0.000000 0.000000 0.000000 0.000000
XLOC_000003   0.000000   0.000000   0.000000 0.000000 0.000000 0.000000
XLOC_000004   0.000000   0.000000   0.000000 0.000000 0.000000 0.000000
XLOC_000005   0.000000   0.000000   0.000000 0.000000 0.000000 0.000000
XLOC_000006   0.680907   0.812414   0.996692 0.884393 0.634819 0.419935


## Selecting genes with diff. expression
gene.count.matrix<-repCountMatrix(genes(cuff))
head(gene.count.matrix)

gene.diff1 = gene.diff[gene.diff$log2_fold_change >= 1.5 | gene.diff$log2_fold_change <= -1.5,]
gene.diff1 = gene.diff1[gene.diff1$significant == 'yes',]

myGeneIds = gene.diff1$gene_id
myGenes<-getGenes(cuff,myGeneIds)
h<-csHeatmap(myGenes,cluster='both', labRow = FALSE) #, replicates = TRUE
h

## Filtering count fpkm data using DIFFGene ids

df4heatmap = gene.rep.matrix[rownames(gene.rep.matrix) %in% myGeneIds,] + 1
max(df4heatmap); min(df4heatmap)
df4heatmap1 = log10(df4heatmap)
max(df4heatmap1); min(df4heatmap1)

my_palette <- colorRampPalette(c("green", 'black', "red"))(n = 299)
# (optional) defines the color breaks manually for a "skewed" color transition

col_breaks = c(seq(0, 0.3,length=100),
               seq(0.31, 0.4,length=100),
               seq(0.41, 3, length=100))

########## creating images ############
map_path = paste(getwd(), '/FPKM_heatmap.png', sep = "")
print(map_path)
png(map_path, 
    width = 12*300,        # 5 x 300 pixels
    height = 10*300,
    res = 300,            # 300 pixels per inch
    pointsize = 8)
#svg(path, width = 10, height = 8)
#print(path)
################################

heatmap.2(as.matrix(df4heatmap1), 
          #cellnote = mat,  # same data set for cell labels
          main = "FPKM values for Diff. Genes", # heat map title
          xlab = "Samples",
          ylab = "Genes",
          labRow = "",
          labCol = colnames(df4heatmap1),
          #cexRow = 0.7,           # Changes the size of col and row font size
          cexCol = 1,
          notecol="black",      # change font color of cell labels to black
          density.info="none",  # turns off density plot inside color legend
          trace="none",         # turns off trace lines inside the heat map
          margins =c(6,4),     # widens margins around plot
          col=my_palette,       # use on color palette defined earlier
          breaks=col_breaks,    # enable color transition at specified limits
          dendrogram="col",     # only draw a row dendrogram
          Colv= TRUE,           #cluster column
          Rowv = FALSE,
          keysize = 0.5
)
dev.off()







