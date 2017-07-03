#library('gplots')
library('RColorBrewer')

path = '/ps/imt/e/HL60_Christene/further_analysis/overlapping_plots/HL60_SKI_GFP_P3,B4-HL60-anti-Runx,HL60_CEBPB_WT/all14019_HL60_SKI_GFP_P3 vs HL60_IgG_GFP_P7 filtered/norm'
name = '/tagcountDF_all_norm'
pkmat = read.csv(file = paste(path, name, '.txt', sep = ""), header = TRUE, sep = '\t')
reclustering = FALSE
manualOrdering = FALSE
selectSamples = FALSE
Length = TRUE
colpersample = 60
##############################################################################

if ('cluster' %in% names(pkmat)){
  #sort_pkmat = pkmat
  sort_pkmat = pkmat[order(pkmat$cluster, decreasing = FALSE), ]
  #sort_pkmat = pkmat[order(pkmat$tagcount, decreasing = FALSE), ]
  #sort_pkmat = pkmat[order(pkmat$group, -pkmat$tagcount), ]
  #sort_pkmat = sort_pkmat[sort_pkmat$cluster %in% c(11,12,13,14,161,1),]
  mat = as.matrix(sort_pkmat[,-c(1,2,3,4,5,6,7,8,9)])
  #mat = as.matrix(sort_pkmat[,-c(1,2,3,4,5,6,7,8,9,10,11,12)])
  condition = "_sort_cluster"
  
  ## manually order clusters
  if (manualOrdering){
    desired_order = c(2,3,4,5,6,7,8,21,0,23,24,25,22)
    pkmat$cluster <- factor( as.character(pkmat$cluster), levels=desired_order )
    # Re-order the data.frame
    sort_pkmat <- pkmat[order(pkmat$cluster),]
    #mat = as.matrix(sort_pkmat[,-c(1,2,3,4,5,6,7,8,9,10,11,12)])
    mat = as.matrix(sort_pkmat[,-c(1,2,3,4,5,6,7,8,9)])
    if (selectSamples){
      mat1 = mat[,1:60]
      mat2 = mat[,121:240]
      mat = cbind2(mat1,mat2)
      print('here')
    }
    condition = "_manual_sort_cluster"
  }
  ## rowsidecol for showing clustering
  selcol <- colorRampPalette(brewer.pal(9,"Set1"))
  hight = as.list(table(sort_pkmat$cluster))
  color = selcol(length(hight))
  rowsidecolor = c()
  for(i in 1:length(hight)){
    rowsidecolor = append(rowsidecolor, rep(color[i], hight[i]))
  }
} 
###########################################################################

if (!'cluster' %in% names(pkmat)){
    sort_pkmat = pkmat
    #sort_pkmat = pkmat[pkmat$cluster %in% c(8),]
    #sort_pkmat = sort_pkmat[,-c(12)]
    clus_df = sort_pkmat[,-c(1,2,3,4,5,6,7,8,9,10)]
    cluster = kmeans(clus_df, 9, iter.max = 10000, nstart = 1, algorithm = "Forgy")
    sort_pkmat['cluster'] = cluster$cluster
    sort_pkmat = sort_pkmat[order(sort_pkmat$cluster),]
    mat = as.matrix(sort_pkmat[,-c(1,2,3,4,5,6,7,8,9,10)]) #,6,7,8,9,10,11
    condition = "_R-Cluster"
    ## Write clustring dataframe
    write.table(sort_pkmat, paste(path, name, condition, '_.txt', sep = ""), sep = '\t')
    ## rowsidecol for showing clustering
    selcol <- colorRampPalette(brewer.pal(9,"Set1"))
    hight = as.list(table(sort_pkmat$cluster))
    color = selcol(length(hight))
    rowsidecolor = c()
    for(i in 1:length(hight)){
      rowsidecolor = append(rowsidecolor, rep(color[i], hight[i]))
    }
}

#######################################################################

if (reclustering){
  whichcluster = 0
  #sort_pkmat = pkmat[order(pkmat$cluster, decreasing = FALSE),]
  clus_pkmat = sort_pkmat[sort_pkmat$cluster %in% c(whichcluster),]
  colInd = grep("^cluster$", colnames(clus_pkmat))
  clus_pkmat = clus_pkmat[,-c(colInd)]
  clus_df = clus_pkmat[,-c(1,2,3,4,5,6,7,8,9,10,11,12)]
  
  ## Clustering & joining cluster values in df
  cluster = kmeans(clus_df, 5, iter.max = 10000, nstart = 1, algorithm = "Forgy")
  clus_pkmat['cluster'] = cluster$cluster + (10*(whichcluster+1))
  clus_pkmat = clus_pkmat[order(clus_pkmat$cluster),]
  clus_pkmat$cluster = as.factor(clus_pkmat$cluster)
  
  sort_pkmat = sort_pkmat[ which(sort_pkmat$cluster != whichcluster ), ]
  sort_pkmat <- rbind(sort_pkmat, clus_pkmat)
  
  ## Write clustring dataframe
  condition = paste("_ReClusterwid",whichcluster)
  write.table(sort_pkmat, paste(path, name, condition, '_.txt', sep = ""), sep = '\t')
  sort_pkmat = sort_pkmat[order(sort_pkmat$cluster, decreasing = FALSE),]
  ## preparing mat for heatmap
  mat = as.matrix(sort_pkmat[,-c(1,2,3,4,5,6,7,8,9,10,11,12)])

  ## rowsidecol for showing clustering
  selcol <- colorRampPalette(brewer.pal(9,"Set1"))
  hight = as.list(table(sort_pkmat$cluster))
  color = selcol(length(hight))
  rowsidecolor = c()
  for(i in 1:length(hight)){
    rowsidecolor = append(rowsidecolor, rep(color[i], hight[i]))
  }
}

#######################################################################
if(Length){
  #sort_pkmat = pkmat
  sort_pkmat = pkmat[order(pkmat$tagcount, decreasing = TRUE),]
  #sort_pkmat = pkmat[order(pkmat$Next.Transcript.tss.distance, decreasing = TRUE),]
  mat = as.matrix(sort_pkmat[,-c(1,2,3,4,5,6,7,8,9,10,11)])
  if (selectSamples){
    mat1 = mat[,1:180]
    mat2 = mat[,241:300]
    mat = cbind2(mat1,mat2)
    print('here')
    condition = "removed_B5"
  }
  condition = "_sort_count"
  rowsidecolor = c()
}

my_palette <- colorRampPalette(c("white","red","red3"))(n = 299)
# (optional) defines the color breaks manually for a "skewed" color transition
if(length(grep('norm', name)) != 0){ #"norm" %in% unlist(strsplit(path, '/'))
  print('norm')
  col_breaks = c(seq(0,2,length=100),
                 seq(2.1,4,length=100),
                 seq(4.1,8,length=100))
} else
  {
  col_breaks = c(seq(0,10,length=100),
                 seq(10.1,15,length=100),
                 seq(15.1,30,length=100))
}

########## creating images ############
map_path = paste(path, name, condition, '.png', sep = "")
print(map_path)
png(map_path, width = 6*300, height = 8*300, res = 300, pointsize = 8)
#svg(map_path, width = 10, height = 8)
#print(path)
################################
## colsep col number
nSam = ceiling(dim(mat)[2]/colpersample) + 1
ColSep = c(1)
for(i in 2:nSam){
  if(i == 2){ColSep[i] = ColSep[i-1] + colpersample-1}
  else{ColSep[i] = ColSep[i-1] + colpersample}
}

#### heatmap3 #################
heatmap.3(mat,
         main = "Tag_density_HeatMap", # heat map title
         xlab = "Distribution of peaks over 3000bp+-",
         ylab = "Peaks",
         col=my_palette,
         breaks=col_breaks,
         dendrogram="none",
         Colv= NA,           #cluster column
         Rowv = NA,
         labRow = "",
         labCol = "",
         #RowSideColors = as.matrix(t(rowsidecolor)),
         colsep = ColSep,
         rowsep = c(0, dim(mat)[1]),
         sepwidth = c(0.1,0.1),
         sepcolor = "#4d4d4d",
         keysize = 0.5
         )
dev.off()


