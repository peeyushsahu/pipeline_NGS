
############### plot single S-shaped tssHeatmap #################
library('RColorBrewer')
path = '/ps/imt/e/HL60_Christene/further_analysis/heatmap_tss_plots'
name = '/region_40000_interval_200_S_HmapWRTTss_SKI_GFP3_norm'
pkmat = read.csv(file = paste(path, name, '.tsv', sep = ""), header = TRUE, sep = '\t')

colpersample = 100

sort_pkmat = pkmat[order(pkmat$Next.Transcript.tss.distance, decreasing = TRUE),]
mat = as.matrix(sort_pkmat[,-c(1,2,3,4,5,6,7,8,9,10)])

my_palette <- colorRampPalette(c("white","red"))(n = 199)
# (optional) defines the color breaks manually for a "skewed" color transition

col_breaks = c(seq(0,2,length=100),
               seq(2.1,3,length=100))

########## creating images ############
map_path = paste(path, name, '.png', sep = "")
print(map_path)
png(map_path, width = 6*300, height = 8*300, res = 300, pointsize = 8)

################################
## colsep col number
nSam = ceiling(dim(mat)[2]/colpersample) + 1

#### heatmap3 #################
heatmap.3(mat,
          main = "Tag_density_HeatMap", # heat map title
          xlab = "Distribution of peaks over 40kb+-",
          ylab = "Peaks",
          col=my_palette,
          breaks=col_breaks,
          dendrogram="none",
          Colv= NA,           #cluster column
          Rowv = NA,
          labRow = "",
          labCol = c('-20kb', rep('',98),'TSS',rep('',99),'+20kb'),
          cexCol = 1.5,
          #RowSideColors = as.matrix(t(rowsidecolor)),
          #colsep = c(1,100),
          rowsep = c(0, dim(mat)[1]),
          #sepwidth = c(0.1,0.5),
          sepcolor = "#4d4d4d",
          keysize = 0.5
)
dev.off()


######### plot two combined tssHeatmaps #########
library('RColorBrewer')
path = '/ps/imt/e/HL60_Christene/further_analysis/heatmap_tss_plots/'

name1 = 'region_40000_interval_200_S_HmapWRTTss_SKI_GFP3_norm'
df1 = read.csv(file = paste(path, name1, '.tsv', sep = ""), header = TRUE, sep = '\t')
df1 = df1[order(df1$Next.Transcript.tss.distance, decreasing = TRUE),]

name2 = 'region_40000_interval_200_S_HmapWRTTss_H3K4me1_EN_norm'
df2 = read.csv(file = paste(path, name2, '.tsv', sep = ""), header = TRUE, sep = '\t')
df2 = df2[order(df2$Next.Transcript.tss.distance, decreasing = TRUE),]

mat1 = as.matrix(df1[,-c(1,2,3,4,5,6,7,8,9,10)])
mat2 = as.matrix(df2[,-c(1,2,3,4,5,6,7,8,9,10)])

quan1 = quantile(mat1)
quan2 = quantile(mat2)
print(quan1)
print(quan2)

mat1 <- ifelse(mat1<1,0,1)
mat2 <- ifelse(mat2<2,0,3)
mat <- mat1 + mat2
colpersample = 100

my_palette <- colorRampPalette(c("white", "red3", "blue", "yellow"))(n = 399)
# (optional) defines the color breaks manually for a "skewed" color transition

col_breaks = c(seq(0,0,length=100),
               seq(0.1,2,length=100),
               seq(2.1,3.5,length=100),
               seq(3.6,4,length=100))

########## creating images ############
map_path = paste(path, name1,'_',name2, '.png', sep = "")
print(map_path)
png(map_path, width = 6*300, height = 8*300, res = 300, pointsize = 8)

################################
## colsep col number
nSam = ceiling(dim(mat)[2]/colpersample) + 1

#### heatmap3 #################
heatmap.3(mat,
          main = "Tag_density_HeatMap", # heat map title
          xlab = "Distribution of peaks over 10000bp+-",
          ylab = "Peaks",
          col=my_palette,
          breaks=col_breaks,
          dendrogram="none",
          Colv= NA,           #cluster column
          Rowv = NA,
          labRow = "",
          labCol = c('-10kb', rep('',98),'TSS',rep('',99),'+10kb'),
          cexCol = 1.5,
          #RowSideColors = as.matrix(t(rowsidecolor)),
          #colsep = c(1,100),
          rowsep = c(0, dim(mat)[1]),
          #sepwidth = c(0.1,0.5),
          sepcolor = "#4d4d4d",
          keysize = 0.5
)
dev.off()


