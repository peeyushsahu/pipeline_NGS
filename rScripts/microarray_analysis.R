#install the core bioconductor packages, if not already installed
source("http://bioconductor.org/biocLite.R")
biocLite()

# install additional bioconductor libraries, if not already installed
biocLite("GEOquery")
biocLite("affy")
biocLite("gcrma")
biocLite("hugene10stv1cdf")
biocLite("hugene10stv1probe")
biocLite("hugene10stprobeset.db")
biocLite("hugene10sttranscriptcluster.db")

# for annotation of Affy-U133 chip
source("https://bioconductor.org/biocLite.R")
biocLite("hgu133plus2.db")

#Load the necessary libraries
library(GEOquery)
library(affy)
library(gcrma)
library(hugene10stv1cdf)
library(hugene10stv1probe)
library(hugene10stprobeset.db)
library(hugene10sttranscriptcluster.db)
library(hgu133plus2.db)
library(pd.hg.u133.plus.2)

# Now the data loading will take place

## You can download data directly from GEO
#Set working directory for download
setwd("/path/to/dir/to/save_data")

#Download the CEL file package for this dataset (by GSE - Geo series id)
getGEOSuppFiles("GSE27447")

#Unpack the CEL files
setwd("/path/to/dir/to/save_data/GSE27447")
untar("GSE27447_RAW.tar", exdir="data")

# Now if you already have .CEL files
setwd("/ps/imt/e/HL60_Christene/Haferlach_data_GSE13159/GSE13159_RAW")
cels = list.files(".", pattern = "CEL")
sapply(paste(getwd(), cels, sep="/"), gunzip)
cels = list.files(".", pattern = "CEL")

## reading metadata file
meta.data = read.csv('/ps/imt/e/HL60_Christene/Haferlach_data_GSE13159/GSE13159_metadata.txt', sep='\t', stringsAsFactors = FALSE)
# choose which sample to take
cels = c()
for(i in 1:dim(meta.data)[1]){
  non_leu = grepl('Non-leukemia', meta.data[i,2])
  aml = grepl('AML', meta.data[i,2])
  all = grepl('ALL', meta.data[i,2])
  if (non_leu == TRUE || aml == TRUE) { 
    cels[length(cels)+1] = paste(meta.data[i,1],'.CEL', sep = "")
  }
}

## reading .cel files
raw.data=ReadAffy(verbose=TRUE, filenames=cels, cdfname="HG-U133_Plus_2") #From bioconductor

#perform RMA normalization (I would normally use GCRMA but it did not work with this chip)
#data.rma.norm=rma(raw.data)
data.gcrma.norm=gcrma(raw.data)

#Get the important stuff out of the data - the expression estimates for each array
gcrma=exprs(data.gcrma.norm)

#Format values to 5 decimal places
gcrma=format(gcrma, digits=5)


## Quality control checks

# load colour libraries
library(RColorBrewer)
# set colour palette
cols <- brewer.pal(8, "Set1")
# plot a boxplot of unnormalised intensity values
boxplot(raw.data, col=cols)
# plot a boxplot of normalised intensity values, affyPLM is required to interrogate celfiles.gcrma
library(affyPLM)
boxplot(data.gcrma.norm, col=cols)
# the boxplots are somewhat skewed by the normalisation algorithm
# and it is often more informative to look at density plots
# Plot a density vs log intensity histogram for the unnormalised data
hist(raw.data, col=cols)
# Plot a density vs log intensity histogram for the normalised data
hist(data.gcrma.norm, col=cols)


#Map probe sets to gene symbols or other annotations
#To see all available mappings for this platform
#ls("package:hugene10stprobeset.db") #Annotations at the exon probeset level
#ls("package:hugene10sttranscriptcluster.db") #Annotations at the transcript-cluster level (more gene-centric view)
ls("package:hgu133plus2.db") #Annotations at the transcript-cluster level (more gene-centric view)

#Extract probe ids, entrez symbols, and entrez ids
probes=row.names(gcrma)
Symbols = unlist(mget(probes, hgu133plus2SYMBOL, ifnotfound=NA))
Entrez_IDs = unlist(mget(probes, hgu133plus2ENTREZID, ifnotfound=NA))

#Combine gene annotations with raw data
gcrma=cbind(probes,Symbols,Entrez_IDs,gcrma)

#Write RMA-normalized, mapped data to file
write.table(gcrma, file = "/ps/imt/e/HL60_Christene/haferlach_gcrma_healthy_aml_hgu133.txt", quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)













