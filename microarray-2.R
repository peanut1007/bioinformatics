# Web Document 11.2. Pevsner from GEO2R
# GEO2R commands

# Version info: R 2.14.1, Biobase 2.15.3, GEOquery 2.23.2, limma 3.10.1
# R scripts generated  Thu Apr 3 13:47:04 EDT 2014

################################################################
#   Differential expression analysis with limma
################################################################

### To load a library in R you will need to first install it, e.g.:
# source("http://bioconductor.org/biocLite.R")
# biocLite("Biobase")
# biocLite("GEOquery")
# biocLite("limma")
library(Biobase)
library(GEOquery)
library(limma)
### You can then get information about various functions in these  
### packages, e.g. > limmaUsersGuide()
# load series and platform data from GEO
### Note that the getGEO command of the GEOquery library is useful to 
### extract GEO datasets

gset <- getGEO("GSE1397", GSEMatrix =TRUE)
class(gset) # list
if (length(gset) > 1) idx <- grep("GPL96", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]  # idxth list element

####
class(gset)# gset has 22283 features, 28 samples
# [1] "ExpressionSet"
# attr(,"package")
# [1] "Biobase

gset[idx] # has 1 feature, 28 samples
class(gset[idx]) # gset[1] seems to be the first feature of the array
# [1] "ExpressionSet"
# attr(,"package")
# [1] "Biobase
####

fvarLabels(gset)
featureNames(gset)
fvarMetadata(gset)
sampleNames(gset)

# make proper column names to match toptable 
fvarLabels(gset) <- make.names(fvarLabels(gset))

# group names for all samples
### Here the object sml will concatenate (abbreviated c) the samples.
### The trisomy 13 samples and controls (n=6) are marked "X".
sml <- c("G0","G0","G0","G0","G1","G1","G1","G1","G1","G1","G1","G0","G0","G0","G1","G1","G0","G0","G1","G1","G0","G0","X","X","X","X","X","X");

# eliminate samples marked as "X"
sel <- which(sml != "X")
sml <- sml[sel]
gset <- gset[ ,sel] # slice of gset all rows, columns that are case

# log2 transform
### We will discuss the rationale for log2 transformation below
ex <- exprs(gset) # ?expression levels ?
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
#      0%       25%       50%       75%       99%      100% 
#   0.100    22.900    73.800   182.500  3927.062 35504.000 

LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0) ||
  (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)
if (LogC) { ex[which(ex <= 0)] <- NaN   # won't take the log of zero or neg (downregulated)
exprs(gset) <- log2(ex) }

# set up the data and proceed with analysis
fl <- as.factor(sml)
gset$description <- fl
design <- model.matrix(~ description + 0, gset)
### model.matrix (from the stats package) creates a design matrix as ### specified.
colnames(design) <- levels(fl)
fit <- lmFit(gset, design) 
### we will discuss lmFit when we perform analyses with limma (below). 
### lmFit (from the limma package) fits a linear model to the log-
### transformed expression values for each probe in a series of 
### arrays.
cont.matrix <- makeContrasts(G1-G0, levels=design)
### makeContrasts determines fold change between groups of samples
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2, 0.01)
### eBayes (from the limma package) uses empirical Bayes statistics to ### determine differential expression. For usage, details, references, ### and examples use 
### > ?eBayes
tT <- topTable(fit2, adjust="fdr", sort.by="B", number=250)
### topTable (from the limma package) extracts a table of the top-
### ranked genes from a linear model fit that has been processed by ### eBayes.

# load NCBI platform annotation
gpl <- annotation(gset)
platf <- getGEO(gpl, AnnotGPL=TRUE)
ncbifd <- data.frame(attr(dataTable(platf), "table"))

# replace original platform annotation
tT <- tT[setdiff(colnames(tT), setdiff(fvarLabels(gset), "ID"))]
tT <- merge(tT, ncbifd, by="ID")
tT <- tT[order(tT$P.Value), ]  # restore correct order

tT <- subset(tT, select=c("ID","adj.P.Val","P.Value","t","B","logFC","Gene.symbol","Gene.title","Chromosome.location"))
write.table(tT, file=stdout(), row.names=F, sep="\t")


# Boxplot for selected GEO samples
library(Biobase)
library(GEOquery)
# load series and platform data from GEO
gset <- getGEO("GSE1397", GSEMatrix =TRUE)
if (length(gset) > 1) idx <- grep("GPL96", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]
# group names for all samples in a series
sml <- c("G0","G0","G0","G0","G1","G1","G1","G1","G1","G1","G1","G0","G0","G0","G1","G1","G0","G0","G1","G1","G0","G0","X","X","X","X","X","X")
# eliminate samples marked as "X"
sel <- which(sml != "X")
sml <- sml[sel]
gset <- gset[ ,sel]
# order samples by group
ex <- exprs(gset)[ , order(sml)]
sml <- sml[order(sml)]
fl <- as.factor(sml)
labels <- c("TS21","euploid")
# set parameters and draw the plot
palette(c("#dfeaf4","#f4dfdf", "#AABBCC"))
dev.new(width=4+dim(gset)[[2]]/5, height=6)
par(mar=c(2+round(max(nchar(sampleNames(gset)))/2),4,2,1))
title <- paste ("GSE1397", '/', annotation(gset), " selected samples", sep ='')
boxplot(ex, boxwex=0.6, notch=T, main=title, outline=FALSE, las=2, col=fl)
legend("topleft", labels, fill=palette(), bty="n")

