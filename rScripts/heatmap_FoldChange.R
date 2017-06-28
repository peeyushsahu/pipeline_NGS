library('gplots')
library('RColorBrewer')

## Joining multiple dataframe using key column
path = '/ps/imt/e/20141009_AG_Bauer_peeyush_re_analysis/further_analysis/PRMT6_KO_analysis'
df1 = read.csv(file = paste(path,'/H3R2ame2_E9_H3R2ame2_B5.1_H3R2me2a_B6.2.txt', 
                                   sep = ''), header = TRUE, row.names=1, sep = '\t')

df2 = read.csv(file = paste(path,'/PRMT6_DB_E9_B6.1.txt', 
                            sep = ''), header = TRUE, sep = '\t')

colnames(df1); colnames(df2); #colnames(df3); colnames(df4)
name = "merge_dataframe"

df2 = df2[df2$log2FC_E9_vs_B5.1_H3K4me3 >= 1 | df2$log2FC_E9_vs_B5.1_H3K4me3 <= -1 | df2$log2FC_WT_vs_B5.1_H3K4me3 >= 1 | df2$log2FC_WT_vs_B5.1_H3K4me3 <= -1, ]
df1 = df1[df1$log2FC_E9_vs_B5.1 > 2 | df1$log2FC_E9_vs_B5.1 < -2 | df1$log2FC_WT_B5.1 > 2 | df1$log2FC_WT_B5.1 < -2, ]

## Merge dataframes
df11 = merge( df1, df2, by ='Next.transcript.gene.name', all.x = T, sort= F)
#df12 = merge( df11, df1, by ='Next.transcript.gene.name', all.x = T, sort= F)
#df13 = merge( df12, df1, by ='Next.transcript.gene.name', all.x = T, sort= F)

## Select only columns with norm
#df = df11[,grep("log", colnames(df11))]
df = df11
## remove rows with NA from dataframe 
df = df[complete.cases(df),]
df['Overlap'] = 0
for(row in 1:dim(df)[1]){
  if(as.character(df[row,'chr.x']) == as.character(df[row,'chr.y'])){
    if(max(as.numeric(df[row,'start.x']),as.numeric(df[row,'start.y'])) < min(as.numeric(df[row,'stop.x']), as.numeric(df[row,'stop.y']))){
      df[row, 'Overlap'] = 1
      print(max(as.numeric(df[row,'start.x']),as.numeric(df[row,'start.y'])))
      print(min(as.numeric(df[row,'stop.x']), as.numeric(df[row,'stop.y'])))
    }
    else{
      df[row, 'Overlap'] = 0
    }
  }
}
df = df[df$Overlap == 1,]
write.table(df, file = '/ps/imt/e/20141009_AG_Bauer_peeyush_re_analysis/further_analysis/PRMT6_KO_analysis/overlapping/common_DE_H3R2me2a_B6_PRMT6_B6.txt', sep='\t', row.names = FALSE)

RowNames = df11$Next.transcript.gene.name

mat = as.matrix(df)
my_palette <- colorRampPalette(c("red3","white","blue"))(n = 299)
# (optional) defines the color breaks manually for a "skewed" color transition

col_breaks = c(seq(-5,-2.1,length=100),
               seq(-2,2,length=100),
               seq(2.1,5,length=100))

########## creating images ############
map_path = paste(path, name, '_P6_K4_heatmap.png', sep = "")
print(map_path)
png(map_path, 
    width = 10*300,        # 5 x 300 pixels
    height = 8*300,
    res = 300,            # 300 pixels per inch
    pointsize = 8)
#svg(path, width = 10, height = 8)
#print(path)
################################

heatmap.2(mat, 
          #cellnote = mat,  # same data set for cell labels
          main = "Tag_density_HeatMap", # heat map title
          xlab = "Differential tag densities",
          ylab = "Peaks",
          labRow = "",
          labCol = colnames(df),
          #cexRow = 0.7,           # Changes the size of col and row font size
          #cexCol = 0.8,
          notecol="black",      # change font color of cell labels to black
          density.info="none",  # turns off density plot inside color legend
          trace="none",         # turns off trace lines inside the heat map
          margins =c(2,4),     # widens margins around plot
          col=my_palette,       # use on color palette defined earlier
          breaks=col_breaks,    # enable color transition at specified limits
          dendrogram="row",     # only draw a row dendrogram
          Colv= FALSE,           #cluster column
          Rowv = TRUE,
          keysize = 0.5,
          na.color = 'black'
)
dev.off()



#########################################################################################
##Kmeans of data
cluster = kmeans(df, 10, iter.max = 1000, nstart = 1, algorithm = "Forgy")
df['cluster'] = cluster$cluster
df = df[order(df$cluster),]

## rowsidecol for showing clustering
selcol <- colorRampPalette(brewer.pal(9,"Set1"))
hight = as.list(table(df$cluster))
color = selcol(length(hight))
rowsidecolor = c()
for(i in 1:length(hight)){
  rowsidecolor = append(rowsidecolor, rep(color[i], hight[i]))
}

mat = as.matrix(df[,-c(6)])
my_palette <- colorRampPalette(c("red3","white","blue"))(n = 299)
# (optional) defines the color breaks manually for a "skewed" color transition

col_breaks = c(seq(-5,-2.1,length=100),
               seq(-2,2,length=100),
               seq(2.1,5,length=100))

########## creating images ############
map_path = paste(path, name, '_P6_K4_Kmeans_heatmap.png', sep = "")
print(map_path)
png(map_path, 
    width = 10*300,        # 5 x 300 pixels
    height = 8*300,
    res = 300,            # 300 pixels per inch
    pointsize = 8)
#svg(path, width = 10, height = 8)
#print(path)
################################

heatmap.2(mat, 
          #cellnote = mat,  # same data set for cell labels
          main = "Tag_density_HeatMap", # heat map title
          xlab = "Differential tag densities",
          ylab = "Peaks",
          labRow = "",
          labCol = colnames(df),
          #cexRow = 0.7,           # Changes the size of col and row font size
          cexCol = 0.8,
          notecol="black",      # change font color of cell labels to black
          density.info="none",  # turns off density plot inside color legend
          trace="none",         # turns off trace lines inside the heat map
          margins =c(3,3),     # widens margins around plot
          col=my_palette,       # use on color palette defined earlier
          breaks=col_breaks,    # enable color transition at specified limits
          dendrogram="none",     # only draw a row dendrogram
          Colv= FALSE,           #cluster column
          Rowv = FALSE,
          keysize = 0.5,
          na.color = 'black',
          RowSideColors = rowsidecolor
)
dev.off()



#########################################################################################
## Load dataframe with foldchanges

path = '/ps/imt/e/20141009_AG_Bauer_peeyush_re_analysis/further_analysis/differential'
name = 'PRMT6_seq6_PRMT6_KO_E.9_PRMT6_KO_B6.2_PRMT6_KO_10.8_PRMT6_KO_B5.1'
diff_count = read.csv(file = paste(path,'/PRMT6_seq6_PRMT6_KO_E.9_PRMT6_KO_B6.2_PRMT6_KO_10.8_PRMT6_KO_B5.1.csv', 
                                      sep = ''), header = TRUE, row.names=1, sep = ',')
colNames = colnames(diff_count)
diff_count1 = diff_count[diff_count$log2FC_E9_vs_B6.2 > 1 | diff_count$log2FC_E9_vs_B6.2 < -1 | 
                           diff_count$log2FC_E9_10.8 < -1 | diff_count$log2FC_E9_10.8 > 1 |
                         diff_count$log2FC_E9_vs_B5.1 < -1 | diff_count$log2FC_E9_vs_B5.1 > 1,]
mat = as.matrix(diff_count1[,c(21,22,23)])

my_palette <- colorRampPalette(c("blue","white","red3"))(n = 299)
# (optional) defines the color breaks manually for a "skewed" color transition

col_breaks = c(seq(-5,-2.1,length=100),
               seq(-2,2,length=100),
               seq(2.1,5,length=100))

########## creating images ############
map_path = paste(path, name, '_heatmap.png', sep = "")
print(map_path)
png(map_path, 
    width = 10*300,        # 5 x 300 pixels
    height = 8*300,
    res = 300,            # 300 pixels per inch
    pointsize = 8)
#svg(path, width = 10, height = 8)
#print(path)
################################

heatmap.2(mat, 
          #cellnote = mat,  # same data set for cell labels
          main = "Tag_density_HeatMap", # heat map title
          xlab = "Differential tag densities",
          ylab = "Peaks",
          labRow = "",
          labCol = colNames[21:23],
          #cexRow = 0.7,           # Changes the size of col and row font size
          cexCol = 1,
          notecol="black",      # change font color of cell labels to black
          density.info="none",  # turns off density plot inside color legend
          trace="none",         # turns off trace lines inside the heat map
          margins =c(5,5),     # widens margins around plot
          col=my_palette,       # use on color palette defined earlier
          breaks=col_breaks,    # enable color transition at specified limits
          dendrogram="both",     # only draw a row dendrogram
          Colv= TRUE,           #cluster column
          Rowv = TRUE,
          keysize = 0.5
)
dev.off()
