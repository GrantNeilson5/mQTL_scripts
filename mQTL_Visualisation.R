#### READ IN ALL CHR FILES ####
library(ggplot2)
library(dplyr)
setwd('/gpfs/ts0/scratch/gn261/Nimhams/mQTL/Output/')

Frontal_mqtl <- read.table("/gpfs/ts0/scratch/gn261/Nimhams/mQTL/Output/Frontal_mQTL.txt")

mydir = "/gpfs/ts0/scratch/gn261/Nimhams/mQTL/Output/"
myfiles = list.files(pattern="*.txt")
myfiles = myfiles[-1]

#Read multiple files from each chr
for(i in 1:22) {
  chr <- paste("Ind_chr",i,".txt", sep = "")
  chr_ob <- strsplit(chr, ".txt") #removes the .txt taken from file name
  tmp <- read.table(paste("Ind_chr",i,".txt", sep=""), stringsAsFactors = F, header = T) #read in the file as a dataframe
  chr_ob <- chr_ob[[1]] #get the dataframe from the list
  assign(chr_ob, tmp) #assigns the file name tmp to the object chr_ob which is the dataframe
  chr_ob <- NULL
  tmp <- NULL}
rm(chr, chr_ob, i, mydir, myfiles, tmp)

#now we rbind all chr as a dataframe in a loop
all_chr <- Ind_chr1
for (i in 2:22){
  tmp <- paste("Ind_chr",i, sep = "")
  tmp <- get(tmp)
  all_chr  <- rbind(all_chr, tmp)
}



all_chr$bonf_corrected_p.value <- p.adjust(all_chr$p.value, method = "bonferroni")
#### SNP and DNAm Site location visualisation #####
##Add CPG annotation data
epicManifest <- read.csv("/gpfs/ts0/scratch/and202/NIMHAMS/epicManifest_hg38.csv", 
                         stringsAsFactors = F, header = T, sep =" ")



meth.all <- cbind(all_chr, epicManifest[match(all_chr$gene, epicManifest$IlmnID),
                                       c("CHR", "MAPINFO")])


#Recode X and Y chromosomes and remove snps that do not have chromosome locations
meth.all <- meth.all[which(meth.all$CHR != "X"),]
meth.all <- meth.all[which(meth.all$CHR != "Y"),]
meth.all<-meth.all[which(!is.na(meth.all$CHR)),] #5416 to 5416
meth.all$CHR<-as.numeric(as.character(meth.all$CHR))

all_chr <- meth.all
meth.all <- select(meth.all, SNP, CHR, MAPINFO)

##Add SNP annotation
genotypeManifest <- read.csv("/gpfs/ts0/scratch/gn261/Nimhams/GSA-24v2-0_A2.csv")

all_chr$SNP <- gsub("_A", "", all_chr$SNP)
all_chrSNP <- gsub("_T", "", all_chr$SNP)
all_chr$SNP <- gsub("_G", "", all_chr$SNP)
all_chr$SNP <- gsub("_C", "", all_chr$SNP)

snps.all <- cbind(all_chr, genotypeManifest[match(all_chr$SNP, genotypeManifest$Name), 
                                           c("Chr", "MapInfo")])

all_chr <- cbind(all_chr, genotypeManifest[match(all_chr$SNP, genotypeManifest$Name), 
                                            c("Chr", "MapInfo")])

snps.all <-snps.all[which(!is.na(snps.all$Chr)),] #8071 to 5632
snps.all$Chr<-as.numeric(as.character(snps.all$Chr))

snps.all <- select(snps.all, SNP, Chr, MapInfo )

chrEnd.meth<-vector(length = 22)
chrEnd.snps<-vector(length = 22)
for(i in 1:22){
  chrEnd.meth[i]<-max(meth.all[which(meth.all$CHR == i),3])
  chrEnd.snps[i]<-max(snps.all[which(snps.all[,2] == i),3])
}

chrStart.meth<-vector(length = 22)
chrStart.snps<-vector(length = 22)
for(i in 1:22){
  chrStart.meth[i]<-min(meth.all[which(meth.all$CHR == i),3])
  chrStart.snps[i]<-min(snps.all[which(snps.all[,2] == i),3])
}

chrSize.snps<-chrEnd.snps-chrStart.snps
chrSize.meth<-chrEnd.meth-chrStart.meth

snp.coord<-vector(length = nrow(snps.all))
meth.coord<-vector(length = nrow(meth.all))
for(i in 1:22){
  snp.coord[which(snps.all[,2] == i)]<-(snps.all[which(snps.all[,2] == i),3]-chrStart.snps[i])/chrSize.snps[i]+i-1
  meth.coord[which(meth.all[,2] == i)]<-(meth.all[which(meth.all[,2] == i),3]-chrStart.meth[i])/chrSize.meth[i]+i-1
}

snps.all<-cbind(snps.all, snp.coord)
meth.all<-cbind(meth.all, meth.coord)

points.y<-snps.all[match(gsub("_.", "", all_chr$SNP), snps.all[,1]),4]
points.x<-meth.all[match(all_chr$SNP, meth.all$SNP),4]

### take absolute values of effect size
logP<-abs(all_chr$beta)

### see source data for resulting data points.


logP_col_palette<-cbind(seq(from = 0, to = 0.61, by = 0.01) , colorRampPalette(c("tan1", "orange",  "red", "red3", "red4",  "darkred"))(62))

pdf("MQTL.pdf")
layout(matrix(c(1,2), ncol = 1), height = c(1,0.1))
op <- par(oma=c(5,7,1,1))
par(op)

points.col<-logP_col_palette[match(round(logP,2), logP_col_palette[,1]),2]
plot(points.x, points.y, pch = 15, cex =0.75, col = points.col, type = "n", ylab = "Position of SNP (Chr)", xlab = "Position of DNA methylation site (Chr)", main = "",  axes = FALSE, xlim = c(0,22), ylim = c(0,22), xaxs = "i", yaxs = "i",cex.lab = 1, cex.axis = 1)
axis(1, at = c(-1,seq(0.5,21.5, by = 1), 23), c("", 1:22, ""), cex.axis = 1)
axis(2, at = c(-1,seq(0.5,21.5, by = 1), 23), c("", 1:22, ""), cex.axis = 1, las = 2)
points(points.x[length(points.x):1], points.y[length(points.x):1], pch = 15, cex =0.75, col = points.col[length(points.x):1])



### plot legend
par(op)
plot(as.numeric(logP_col_palette[,1]), rep(1, nrow(logP_col_palette)), col = logP_col_palette[,2], pch = 15, cex = 1, axes = FALSE, main  = "", ylab = "", ylim = c(0.9999,1.001),xlab = "", cex.lab = 1, cex.axis = 1, cex.main = 1, adj = 0)
axis(1, at = c(0,0.2,0.4,0.6),  c(0,0.2,0.4,0.6)*100, cex.axis = 1)
dev.off()


## calculate distance between SNP and methylation probe
dist<-abs(all_chr$MAPINFO)
dist[which(as.character(all_chr$Chr) != as.character(all_chr$CHR))]<--9


par(op)
pdf("Distnace_of_betweenmQTLS.pdf")
signedDist<-(all_chr$MAPINFO)/1000
plot(signedDist[which(dist != -9)], abs(all_chr$beta[which(dist != -9)])*100, main = "", xlab = "Distance (Mb)", ylab = "Effect size (% difference in DNA methylation per minor allele)", pch = 16, cex = 0.8, ylim = c(-0,40), xlim = c(-500,500), cex.axis = 1, cex.lab = 0.8, axes = FALSE, xaxs = "i", yaxs = "i")
axis(1, seq(-1, 1, 0.5), at = seq(-1000, 1000, 500), cex.axis = 1, cex.lab = 1)
axis(2, las = 2, cex.axis = 1, cex.lab = 1)
 dev.off()
