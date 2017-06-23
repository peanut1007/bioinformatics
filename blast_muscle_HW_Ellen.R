# from page bottom of:
# http://www.molgen.ua.ac.be/bioinfo/acourse/Bioconductor-by-examples.txt

# Some web services based packages
# ================================

#install.packages("RCurl")
#install.packages("rentrez")
#source("https://bioconductor.org/biocLite.R")
#biocLite() 
#biocLite("Biostrings")
library(Biostrings)
library(RCurl)
library(rentrez)
seq1 <- entrez_fetch(db = "sequences", rettype = 'fasta', id = "1139595") #human
seq2 <- entrez_fetch(db = "sequences", rettype = 'fasta', id = "332368596") #domestic goose
seq3 <- entrez_fetch(db = "sequences", rettype = 'fasta', id = "7544130") # junglefowl 
seq4 <- entrez_fetch(db = "sequences", rettype = 'fasta', id = "311909432") # westernClawedFrog 
seq5 <- entrez_fetch(db = "sequences", rettype = 'fasta', id = "62825355") # zebraFish
seq6 <- entrez_fetch(db = "sequences", rettype = 'fasta', id = "7021524") # rhesusMonkey
seq7 <- entrez_fetch(db = "sequences", rettype = 'fasta', id = "2760950") # houseMouse
seq8 <- entrez_fetch(db = "sequences", rettype = 'fasta', id = "1510129") # NorwayRat
seq9 <- entrez_fetch(db = "sequences", rettype = 'fasta', id = "66730862") # dog
seq10 <- entrez_fetch(db = "sequences", rettype = 'fasta', id = "57864232") # cow
seq11 <- entrez_fetch(db = "sequences", rettype = 'fasta', id = "677308066") # bar-tailed trogon


sequences <- c(seq1, seq2, seq3, seq4, seq5, seq6, seq7, seq8, seq9, seq10, seq11)

# # get a sequence by id (using NCBI eutil efetch)
getprotein = function(id) {
  seq=getForm("http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi", db="protein", id=id, retmode="text", rettype="fasta", binary=F)
  s=strsplit(seq[1],"\\n")[[1]]
  seq=AAString(paste(s[2:length(s)],collapse=""))
}
seq=getprotein("15718680")

#Here is the FASTA sequence for one of the isoforms of human P53 protein:
# seq <- paste("1 meepqsdpsv epplsqetfs dlwkllpenn vlsplpsqam ddlmlspddi eqwftedpgp
#  61 deaprmpeaa ppvapapaap tpaapapaps wplsssvpsq ktyqgsygfr lgflhsgtak
# 121 svtctyspal nkmfcqlakt cpvqlwvdst pppgtrvram aiykqsqhmt evvrrcphhe
# 181 rcsdsdglap pqhlirvegn lrveylddrn tfrhsvvvpy eppevgsdct tihynymcns
# 241 scmggmnrrp iltiitleds sgnllgrnsf evrvcacpgr drrteeenlr kkgephhelp
# 301 pgstkralpn ntssspqpkk kpldgeyftl qdqtsfqken c")

# blast a sequence using the NCBI REST interface
ncbi_blast = function(seq,database="nr",expect=10,program="blastp",hitlistsize=100) {
  # based on the NCBI REST interface (http://www.ncbi.nlm.nih.gov/blast/Doc/urlapi.html)
  job = as.character(getForm("http://blast.ncbi.nlm.nih.gov/Blast.cgi", QUERY=seq, DATABASE=database,
                             HITLIST_SIZE=hitlistsize, FILTER="L", EXPECT=expect, FORMAT_TYPE="Text", PROGRAM=program, CLIENT="web",
                             SERVICE="plain", NCBI_GI="on", PAGE="Nucleotides", CMD="Put"))
  m = regexpr("RID = ([^\n]+)",job)
  id = substring(job,m[1]+6,m[1]+attributes(m)$match.length-1)
  m = regexpr("RTOE = ([^\n]+)",job)
  rtoe = substring(job,m[1]+7,m[1]+attributes(m)$match.length-1)
  result = "Status=WAITING"
  while (length(grep("Status=WAITING",result))) {
    Sys.sleep(1)
    result = as.character(getForm("http://blast.ncbi.nlm.nih.gov/Blast.cgi", RID=id, FORMAT_TYPE="Text", CMD="Get"))
  }
  getForm("http://blast.ncbi.nlm.nih.gov/Blast.cgi", RID=id, CMD="Delete")
  return(result)
}
b = ncbi_blast(seq,database="nr",program="blastp") 
#b = ncbi_blast(seq,database="nr",program="PSI-BLAST") 
head(b)
cat(b)

# multiple sequence alignment
# EBI REST service for multiple sequence alignment using muscle (http://www.ebi.ac.uk/Tools/webservices/services/msa/muscle_rest#resulttypes)
muscle = function(sequences) {
  # convert sequences to fasta
  fasta=paste(sequences,collapse="\n")
  #fasta=paste(paste(">",names(sequences),sep=""),as.character(sequences),sep="\n",collapse="\n")
  id=postForm("http://www.ebi.ac.uk/Tools/services/rest/muscle/run/",style="post",email="pfrancislyon@usfca.edu",format="fasta",sequence=fasta)
  status = ""
  while (status != "FINISHED") {
    Sys.sleep(1)
    status=getURL(sprintf("http://www.ebi.ac.uk/Tools/services/rest/muscle/status/%s",id))
  }
  resulttype="aln-fasta"
  ali=getURL(sprintf("http://www.ebi.ac.uk/Tools/services/rest/muscle/result/%s/%s", id, resulttype))
}


#You may want to fetch sequences for MSA with MUSCLE using this sort of command for each of your 11 sequences 
seq1 <- entrez_fetch(db = "sequences", rettype = 'fasta', id = "1139595") #human
seq2 <- entrez_fetch(db = "sequences", rettype = 'fasta', id = "332368596") #domestic goose
seq3 <- entrez_fetch(db = "sequences", rettype = 'fasta', id = "7544130") # junglefowl 
seq4 <- entrez_fetch(db = "sequences", rettype = 'fasta', id = "311909432") # westernClawedFrog 
seq5 <- entrez_fetch(db = "sequences", rettype = 'fasta', id = "62825355") # zebraFish
seq6 <- entrez_fetch(db = "sequences", rettype = 'fasta', id = "7021524") # rhesusMonkey
seq7 <- entrez_fetch(db = "sequences", rettype = 'fasta', id = "2760950") # houseMouse
seq8 <- entrez_fetch(db = "sequences", rettype = 'fasta', id = "1510129") # NorwayRat
seq9 <- entrez_fetch(db = "sequences", rettype = 'fasta', id = "66730862") # dog
seq10 <- entrez_fetch(db = "sequences", rettype = 'fasta', id = "57864232") # cow
seq11 <- entrez_fetch(db = "sequences", rettype = 'fasta', id = "677308066") # bar-tailed trogon


sequences <- c(seq1, seq2, seq3, seq4, seq5, seq6, seq7, seq8, seq9, seq10, seq11)

#The line beginning with "fasta=" has been commented out and replaced with a line
# to use when seq is created as a vector of sequences retrieved from Entrez DB
muscle = function(sequences) {
  # convert sequences to fasta
  fasta=paste(sequences,collapse="\n")
  #fasta=paste(paste(">",names(sequences),sep=""),as.character(sequences),sep="\n",collapse="\n")
  id=postForm("http://www.ebi.ac.uk/Tools/services/rest/muscle/run/",style="post",email="pfrancislyon@usfca.edu",format="fasta",sequence=fasta)
  status = ""
  while (status != "FINISHED") {
    Sys.sleep(1)
    status=getURL(sprintf("http://www.ebi.ac.uk/Tools/services/rest/muscle/status/%s",id))
  }
  resulttype="aln-fasta"
  ali=getURL(sprintf("http://www.ebi.ac.uk/Tools/services/rest/muscle/result/%s/%s", id, resulttype))
}

# write to file
result = muscle(sequences)
sink("muscle_align.txt")
cat(result)
sink()

## create an object from file of aligned sequences
msa <- readAAMultipleAlignment(filepath ="muscle_align.txt", format="fasta") 

## See a detailed pager view
if (interactive()) {
  detail(msa)
}

## calculate frequencies
alphabetFrequency(msa)
consensusMatrix(msa, baseOnly=TRUE)[, 84:90] 

## cluster the alignments
sdist <- stringDist(as(msa,"AAStringSet"), method="levenshtein")
clust <- hclust(sdist, method = "single")
plot(clust, hang = -1)
fourgroups <- cutree(clust, 4)
fourgroups

# ## if you have a phylips alignment write out the alignement object to Phylip format
# write.phylip(x = msa, filepath = tempfile("foo.txt",tempdir()))


autoMasked <- maskGaps(msa, min.fraction=0.5, min.block.width=4)
autoMasked


# code chunks from Biostrings documentation that comes with the Biostrings package:
### http://127.0.0.1:18596/library/Biostrings/doc/MultipleAlignments.R
###################################################
### code chunk number 14: cluster
###################################################
sdist <- stringDist(as(msa,"AAStringSet"), method="hamming")
clust <- hclust(sdist, method = "single")
pdf(file="badTree.pdf")
plot(clust, hang = -1)
dev.off()


###################################################
### code chunk number 15: cluster2
###################################################
sdist <- stringDist(as(autoMasked,"AAStringSet"), method="hamming")  ## DNAStringSet vs AAStringSet
clust <- hclust(sdist, method = "single")
pdf(file="goodTree.pdf")
plot(clust, hang = -1)
dev.off()
fourgroups <- cutree(clust, 4)
fourgroups


###################################################
### code chunk number 16: fastaExample (eval = FALSE)
###################################################
AAStr = as(msa, "AAStringSet")
writeXStringSet(AAStr, file="myFile.fa")

