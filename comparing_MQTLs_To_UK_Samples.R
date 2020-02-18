## This script creates an allele frquency chart for different population based on selected number of SNPs

####      Read in files     ####
library(ggplot2)
library(dplyr)
library(tidyr)
library(stringi)
library(stringr)

setwd('/gpfs/ts0/scratch/gn261/Nimhams//mQTL/Output/')
mydir = "/gpfs/ts0/scratch/and202/NIMHAMS//mQTL/Output/"
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



##Add annotation files
epicManifest <- read.csv("/gpfs/ts0/scratch/and202/NIMHAMS/epicManifest_hg38.csv", 
                         stringsAsFactors = F, header = T, sep =" ")

genotypeManifest <- read.csv("/gpfs/ts0/scratch/gn261/Nimhams/Infinium_Global_Screening-24_v2.0_manifest.csv")

##Read in the european samples that had mQTLs run on
load("/gpfs/ts0/scratch/and202/UKDataset/AllAutosomeResults_mQTLs_PFC_Anno_with2PCs.Rdata")


####      Find SNPs of interest     ####
#check which cpgs are in both datasets
length(intersect(res.pfc.anno$gene, all_chr$gene)) #1731 cpgs

#Take the snps that do not intersect in both datasets
#SNPs only in India. Which snps not in european mqtl
unique_ind_snps = all_chr[-which(all_chr$SNP %in% res.pfc.anno$SNP),] #6804
#remove the a/c/G/T on the SNP name then read out as text.file
unique_ind_snps$SNP <- gsub("_A", "", unique_ind_snps$SNP)
unique_ind_snps$SNP <- gsub("_T", "", unique_ind_snps$SNP)
unique_ind_snps$SNP <- gsub("_G", "", unique_ind_snps$SNP)
unique_ind_snps$SNP <- gsub("_C", "", unique_ind_snps$SNP)
unique_ind_snps <- unique_ind_snps[order(unique_ind_snps$FDR),]
unique_ind_snps$SNP <- gsub("GSA.", "", unique_ind_snps$SNP)
ind_snps = unique_ind_snps$SNP

write.table(ind_snps, 'SNPs_only_in_india_mqtl.txt', row.names = F, quote = F, col.names =F)

#save the snps that are present in ind_snps. These snps will be used to load the allele frequency in different
#populations. The code to do this is: cat SNPs_only_in_india_mqtl.txt | while read R ; do wget -q -O - "https://www.ncbi.nlm.nih.gov/snp/${R}?download=frequency" | grep -E '^(#Study|1000Genomes)' | sed "s/^/${R}\t/" >> snp_freq_1000G.txt ; done

####      Visualise Allele Frequency In SA and Euro populations       #####

##read in allele frequency
ind_snp_freq = read.table("snp_freq_1000G.txt",
                          sep = "\t", fill = TRUE, stringsAsFactors = F) #36756
ind_snp_freq <- ind_snp_freq[!duplicated(ind_snp_freq), ]

# make an index eg. every 7th to be removed as there are blank rows and remove from dataframe
ind <- seq(1, nrow(ind_snp_freq), by=7)
ind_snp_freq = ind_snp_freq[-ind, ] #25932

rownames(ind_snp_freq) <- c(1:nrow(ind_snp_freq))
colnames(ind_snp_freq) <- c("SNP", "Study", "Population", "Group", "SampleSize", "RefAllele",
                            "AltAllele", "BioProjectID", "BiosSampleID")
ind_snp_freq<- ind_snp_freq[which(ind_snp_freq$Population == c("Europe" , "South Asian")),]
#create two tables for each population then combine them later
ind_snp_freq_euro <- ind_snp_freq[which(ind_snp_freq$Population == "Europe"),]
ind_snp_freq_sa <- ind_snp_freq[which(ind_snp_freq$Population == "South Asian"),]


##make combinations of allele frequency in respective tables. This will be joined to the table later on.
#Euro
ind_snp_freq_euro$Freq <- rep(NA, nrow(ind_snp_freq_euro))
ind_snp_freq_euro$Freq = paste( gsub("[^0-9.]",'',ind_snp_freq_euro$RefAllele),
                                gsub("[^0-9.]",'',ind_snp_freq_euro$AltAllele),
                                sep = "/")
#SA
ind_snp_freq_sa$Freq <- rep(NA, nrow(ind_snp_freq_sa))
ind_snp_freq_sa$Freq = paste( gsub("[^0-9.]",'',ind_snp_freq_sa$RefAllele),
                              gsub("[^0-9.]",'',ind_snp_freq_sa$AltAllele),
                              sep = "/")


##Now merge two population to one
ind_snp_freq_tbl <- merge(ind_snp_freq_euro,ind_snp_freq_sa, by = "SNP", all= TRUE)

ind_snp_freq_tbl$AlleleComb = rep(NA, nrow(ind_snp_freq_tbl))
ind_snp_freq_tbl$AlleleComb = paste(stri_extract_first_regex(ind_snp_freq_tbl$RefAllele.x, ".{1}"),
                                    stri_extract_first_regex(ind_snp_freq_tbl$AltAllele.x, ".{1}"),
                                    sep = "/")
colnames(ind_snp_freq_tbl)[colnames(ind_snp_freq_tbl) == 'Freq.x'] <- 'Europe'
colnames(ind_snp_freq_tbl)[colnames(ind_snp_freq_tbl) == 'Freq.y'] <- 'South.Asia'
rem = c("Population.x", "Population.y", "SampleSize.x", "SampleSize.y",
        "BioProjectID.x", "BiosSampleID.x", "BioProjectID.y", "BiosSampleID.y",
        "Group.x", "Group.y", "Study.y", "Study.x", "RefAllele.x", "RefAllele.y", 
        "AltAllele.x", "AltAllele.y")
ind_snp_freq_tbl = ind_snp_freq_tbl[, -which(names(ind_snp_freq_tbl) %in% rem)]


ind_snp_freq_tbl2<- cbind(ind_snp_freq_tbl, genotypeManifest[match(ind_snp_freq_tbl$SNP, genotypeManifest$Name), 
                                                             c("Chr", "MapInfo")])
ind_snp_freq_tbl2$Chr <- as.numeric(as.character(ind_snp_freq_tbl2$Chr))
ind_snp_freq_tbl2$Location <- rep(NA, nrow(ind_snp_freq_tbl2))
ind_snp_freq_tbl2$Location <- paste(ind_snp_freq_tbl2$Chr, ind_snp_freq_tbl2$MapInfo, sep = ":")
ind_snp_freq_tbl2 = ind_snp_freq_tbl2[, -which(names(ind_snp_freq_tbl2) %in% c("Chr", "MapInfo"))]

ind_snp_freq_tbl2 = ind_snp_freq_tbl2[c("SNP", "Location", "AlleleComb", "Europe", "South.Asia")]
colnames(ind_snp_freq_tbl2)[colnames(ind_snp_freq_tbl2) == 'AlleleComb'] <- 'Allele'
colnames(ind_snp_freq_tbl2)[colnames(ind_snp_freq_tbl2) == 'South.Asia'] <- 'South Asia'
rownames(ind_snp_freq_tbl2) <- NULL
pdf('Ind SNP frequency in 1KG.pdf', w = 8, h = 8)
grid.table(ind_snp_freq_tbl2[1:20,], rows = NULL)
dev.off()


###snps that have a higher frequncy in european comapred to indian

euro_freq <- str_split_fixed(ind_snp_freq_tbl2$Europe, "/", 2)
colnames(euro_freq) <- c("euro_major", "euro_minor")

SA_freq <-str_split_fixed(ind_snp_freq_tbl2$`South Asia`, "/", 2)
colnames(SA_freq) <- c("SA_major", "SA_minor")

ind_snp_freq_tbl3 <- cbind(ind_snp_freq_tbl2, euro_freq, SA_freq)

ind_snp_freq_tbl3$euro_minor <- as.numeric(as.character(ind_snp_freq_tbl3$euro_minor))
ind_snp_freq_tbl3$SA_minor <- as.numeric(as.character(ind_snp_freq_tbl3$SA_minor))

ind_snp_freq_tbl3$difference <- abs(ind_snp_freq_tbl3$euro_minor - ind_snp_freq_tbl3$SA_minor)

ind_snp_freq_tbl3 <- ind_snp_freq_tbl3[order(ind_snp_freq_tbl3$difference, decreasing = TRUE),]

ind_snp_freq_top_20 <- ind_snp_freq_tbl3[1:20,]

ind_snp_freq_top_20$SNP <- as.factor(ind_snp_freq_top_20$SNP)

##Making plotting data frame
myvars <- c("SNP", "euro_minor")
Euro_minor <- ind_snp_freq_top_20[myvars]
colnames(Euro_minor)[2] <- "Minor"

myvars1 <- c("SNP", "SA_minor")
SA_minors  <- ind_snp_freq_top_20[myvars1]
colnames(SA_minors)[2] <- "Minor"

plotting <- bind_rows(Euro_minor, SA_minors)

pdf("SNPs_with_largest_Absolute_Diff_in_MAF_Fronal.pdf")
ggplot(plotting, aes(SNP, Minor))+
  geom_point()+
  geom_line()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  ylab("Minor Allele Frequency")
dev.off()
  
pdf("Top 20 SNPs with largest difference in MAF Frontal.pdf")
par(las=2,mar=c(7, 4.1, 4.1, 4), xpd=TRUE)
plot(seq(1,20,1), ind_snp_freq_top_20$euro_minor, xaxt='n',ann=FALSE)
title(main= "Top 20 SNPs with largest difference in MAF (frontal)",
      xlab="", ylab="MAF")
for( i in 1:20){
  segments(x0 = i, y0 = ind_snp_freq_top_20$SA_minor[i],
           x1 = i ,y1 = ind_snp_freq_top_20$euro_minor[i],
           col= 'grey', lwd = 2)
}
points(ind_snp_freq_top_20$SA_minor,pch=16,col="red")
points(ind_snp_freq_top_20$euro_minor,pch=16,col="blue")
legend("topright",inset=c(-0.1,0),c("SA", "EU"),cex=.8,col=c("red","blue"),pch=c(16,16))
axis(1, at=1:20, labels=ind_snp_freq_top_20$SNP)
dev.off()
  
               

#### Make a table using the freq imputed from hapmap? ####
#read in frequency file and clean
frq = read.table("/gpfs/ts0/scratch/and202/NIMHAMS/data_filtered_frq.frq", header = T,
                 stringsAsFactors = F)
#add annotation and combine chr and bp in all_chr file
all_chr2 <- all_chr

all_chr2$SNP <- gsub("_A", "", all_chr2$SNP)
all_chr2SNP <- gsub("_T", "", all_chr2$SNP)
all_chr2$SNP <- gsub("_G", "", all_chr2$SNP)
all_chr2$SNP <- gsub("_C", "", all_chr2$SNP)
#Add in the chr and bp information to all_chr2
all_chr2 <- cbind(all_chr2, genotypeManifest[match(all_chr2$SNP, genotypeManifest$Name), 
                                             c("Chr", "MapInfo")])

all_chr2 <- all_chr2[which(!is.na(all_chr2$Chr)),] #remove rows which are NA
all_chr2$Chr = paste("chr", all_chr2$Chr, sep = "")
all_chr2$SNP1 = paste(all_chr2$Chr, all_chr2$MapInfo, sep = ":")
all_chr2$SNP1 = paste(all_chr2$SNP1, "SNP", sep = ":")
all_chr2 <- all_chr2[order(all_chr2$SNP1),]

length(intersect(all_chr2$SNP1, frq$SNP)) #30
int_snps <- all_chr2[which(all_chr2$SNP1 %in% frq$SNP),9] #35



int_snp_df = frq[frq$SNP %in% int_snps,]
int_snp_df2 = separate(data = int_snp_df, col =SNP, 
                       into = c("chr", "bp", "snp"), sep = ":")
int_snp_df2 = int_snp_df2[, -which(names(int_snp_df2) %in% c("chr", "snp", "NCHROBS"))]
pdf('Euro imputed MAF, unique SNPs to India.pdf', w = 8, h = 10, onefile = TRUE)
grid.table(int_snp_df2, rows = NULL)
dev.off()
