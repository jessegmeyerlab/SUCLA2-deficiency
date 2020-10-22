if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("UniProt.ws")
BiocManager::install("vctrs")
library(ggplot2)
library(tidyr)
library(ggpubr)
library(UniProt.ws)
### read in sites from Rardin et al
rsk<-read.delim("D:/NatureComms_succK/R/Park2013sirt5KO.txt", stringsAsFactors = F)
head(rsk)
nrow(rsk)

## map park proteins to genes as done for fibro data
up <- UniProt.ws(taxId=10090)  ## mouse taxID
taxId(up) <- 10090
keys <- unique(rsk$Protein)
nchar(keys)
#keys<-keys[nchar(keys)>1]

columns <- c("SEQUENCE","GENES")
kt <- "UNIPROTKB"

rsk_map_result <- select(up, keys, columns, kt)
head(rsk_map_result)

first_gene = function(x)unlist(strsplit(x, " "))[1]
rsk_map_result$firstgene = unlist(lapply(FUN=first_gene, rsk_map_result$GENES))
head(rsk_map_result)
rsk_map_result[1,]

## have first gene, now map those back to the original list
rsk$gene <- rsk_map_result$firstgene[match(rsk$Protein, rsk_map_result$UNIPROTKB)]

rsk$gene[229]
rsk$Protein[229]
rsk$Fasta.headers
#select(up, rsk$Protein[229], columns, kt)

### some genes are missing because they report trembl IDs

head(rsk)

rsk$gene_site <-paste(toupper(rsk$gene), rsk$Position, sep="_")

rsk$Sirt5log2FC<-log2(as.numeric(rsk$Ratio.Sirt5.KO.WT))

### make small subset dataframe
rsks<-rsk[,c("gene_site","Protein", "Sirt5log2FC")]
colnames(rsks)[2]<-"uniprot"

### read myotube sites
msk<-read.delim("D:/NatureComms_succK/R/myo_sitelvl_final_v2.txt", stringsAsFactors = F)
msk$firstgene<-  unlist(lapply(FUN=first_gene, msk$genenames))
msk$gene_site<- paste(msk$firstgene, msk$site, sep="_")

head(msk)
## remove those without significant FDR
nrow(msk)
msk<-msk[msk$FDR<0.05,]
nrow(msk)
### make small subset frame
msks<-msk[,c("gene_site","uniprot", "log2FC")]

## read fibroblast sites
fsk<-read.delim("D:/NatureComms_succK/R/fibro_sitelvl_final_v2.txt", stringsAsFactors = F)
fsk$firstgene<-  unlist(lapply(FUN=first_gene, fsk$genenames))
fsk$gene_site<- paste(fsk$firstgene, fsk$sites, sep="")

head(fsk)

nrow(fsk)
fsk<-fsk[fsk$FDR<0.05,]
nrow(fsk)
## how many unique genes?
length(unique(c(fsk$firstgene,msk$firstgene)))  ## 359

## how many unique uniprot?
length(unique(c(fsk$uniprot,msk$uniprot)))  ## 359

length(unique(c(fsk$firstgene)))

### filter each by FDR < 0.01

fsks<-fsk[,c("gene_site","uniprot", "log2FC")]

# merge the myo and fibro dataframes
m1 <- merge(msks, fsks, by=c("gene_site", "uniprot"), all = T)
m1

#### map the uniprot ids to get genes and sequences
uph <- UniProt.ws(taxId=9606)  ## mouse taxID
taxId(uph) <- 9606
keys <- unique(m1$uniprot)
columns <- c("SEQUENCE", "GENES")
kt <- "UNIPROTKB"

m1_map_result <- select(uph, keys, columns, kt)

m1_map_result$firstgene  <-  unlist(lapply(FUN=first_gene, m1_map_result$GENES))

head(fsk)

table(fsk$firstgene)
rsk$gene<-toupper(rsk$gene)
table(rsk$gene)
data.frame(table(rsk$gene))
merged_genecounts<- merge(data.frame(table(rsk$gene)), data.frame(table(fsk$firstgene)), by="Var1")
colnames(merged_genecounts)[2] <- "Rardin"
colnames(merged_genecounts)[3] <- "Fibro"




# GET EQUATION AND R-SQUARED AS STRING
# SOURCE: https://groups.google.com/forum/#!topic/ggplot2/1TgH-kG5XMA

lm_eqn <- function(df){
  m <- lm(Rardin ~ Fibro, df);
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2, 
                   list(a = format(unname(coef(m)[1]), digits = 2),
                        b = format(unname(coef(m)[2]), digits = 2),
                        r2 = format(summary(m)$r.squared, digits = 3)))
  as.character(as.expression(eq));
}


merged_genecounts
library(ggplot2)
ggplot(merged_genecounts, aes(x=Fibro, y=Rardin))+ 
  geom_point() +
  theme_Publication()+ geom_smooth(method = lm)+ 
  geom_text(x = 5, y = 20, label = lm_eqn(merged_genecounts), parse = TRUE)
  

p<-ggplot(merged_genecounts, aes(x=Fibro, y=Rardin))+ 
  geom_point() +
  theme_Publication()+ geom_smooth(method = lm)


p1 <- p 


sum= lm(Rardin ~ Fibro, merged_genecounts)
summary(sum)
ggplot2.scatterplot(data=merged_genecounts, xName=Rardin, yName=Fibro)

dev.off()


write.table(merged_genecounts, "D:/NatureComms_succK/R/outputs/merged_genecounts.txt", sep="\t", row.names = F, )
