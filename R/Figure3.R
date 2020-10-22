## Figure 3 Nature Communications paper - overlap between human and mouse succinylation sites


########## ggpublication theme #########3
#### https://rpubs.com/Koundy/71792
theme_Publication <- function(base_size=14, base_family="helvetica") {
  library(grid)
  library(ggthemes)
  (theme_foundation(base_size=base_size, base_family=base_family)
    + theme(plot.title = element_text(face = "bold",
                                      size = rel(1.2), hjust = 0.5),
            text = element_text(),
            panel.background = element_rect(colour = NA),
            plot.background = element_rect(colour = NA),
            panel.border = element_rect(colour = NA),
            axis.title = element_text(face = "bold",size = rel(1)),
            axis.title.y = element_text(angle=90,vjust =2),
            axis.title.x = element_text(angle=90, vjust = -0.2),
            axis.text = element_text(), 
            axis.line = element_line(colour="black"),
            axis.ticks = element_line(),
            panel.grid.major = element_line(colour="#f0f0f0"),
            panel.grid.minor = element_blank(),
            legend.key = element_rect(colour = NA),
            legend.position = "bottom",
            legend.direction = "horizontal",
            legend.key.size= unit(0.2, "cm"),
            legend.margin = unit(0.2, "cm"),
            legend.title = element_text(face="italic"),
            plot.margin=unit(c(15,7,7,7),"mm"),
            strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
            strip.text = element_text(face="bold")
    ))
}

scale_fill_Publication <- function(...){
  library(scales)
  discrete_scale("fill","Publication",manual_pal(values = c("#386cb0","#fdb462","#7fc97f","#ef3b2c","#662506","#a6cee3","#fb9a99","#984ea3","#ffff33")), ...)
  
}

scale_colour_Publication <- function(...){
  library(scales)
  discrete_scale("colour","Publication",manual_pal(values = c("#386cb0","#fdb462","#7fc97f","#ef3b2c","#662506","#a6cee3","#fb9a99","#984ea3","#ffff33")), ...)
  
}

###################################################33

#BiocManager::install("msa")
#BiocManager::install("UniProt.ws")
library(msa)
library(data.table)


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


### test sequence alignment by getting human and mouse version of one protein
human_testseq1<-m1_map_result[which(m1_map_result$firstgene=="ACAT1"),"SEQUENCE"]

mouse_testseq1<-rsk_map_result[which(toupper(rsk_map_result$firstgene)=="ACAT1"),"SEQUENCE"]
human_testseq1
mouse_testseq1

test_align<-pairwiseAlignment(human_testseq1, mouse_testseq1)

test_align
insertion(test_align)
deletion(test_align)
indel(test_align)

rsks[1,]


### psuedo code
#1 get shared genes
mouse_genes<-toupper(rsk$gene)
mouse_genes
human_genes<-fsk$firstgene
human_genes


### loop through shared genes
gshared  <-  intersect(human_genes, mouse_genes)
#x<-gshared[2]
rsk$alignedPosition<-rsk$Position
for(x in gshared){
  print(x)
  htmp<-m1_map_result[which(m1_map_result$firstgene==x),"SEQUENCE"]
  mtmp<-rsk_map_result[which(toupper(rsk_map_result$firstgene)==x),"SEQUENCE"]
  ta<-pairwiseAlignment(htmp, mtmp)
  tmpindel<-indel(ta)
  ## add the insertions to mouse #s
  tmp_gene_index<-which(toupper(rsk$gene)==x)
  # loop through the succinyl sites for that gene
  for(pos in tmp_gene_index){
    #print(pos)
    tmp_position<-rsk$Position[pos]
    #print(tmp_position)
    delta=0 ### keep track of change in position
    ## loop through insertions
    if(length(tmpindel@insertion[[1]])>0){
      print("contains insertions")
      for(i in 1:length(tmpindel@insertion[[1]])){
        y= tmpindel@insertion[[1]][i]
        if(y@start<tmp_position){
          delta=delta+y@width
        }
      }
    }
    #loop through deletions, take away width
    if(length(tmpindel@deletion[[1]])>0){
      print("contains deletions")
      for(i in 1:length(tmpindel@deletion[[1]])){
        y=tmpindel@deletion[[1]][i]
        if(y@start<tmp_position){
          delta=delta-y@width
        }
      }
    }
    rsk$alignedPosition[pos]<- tmp_position+delta
  }
}

head(rsk)

## how many sites have shifted?  ==118 
length(which(rsk$Position!=rsk$alignedPosition))

# what percent of all sites is that? ==11%
length(which(rsk$Position!=rsk$alignedPosition))/nrow(rsk)

### which sites have changed
which(rsk$Position!=rsk$alignedPosition)

### check some of those manually

rsk[which(rsk$Position!=rsk$alignedPosition)[1],]  ### aldolase, checks out

rsk[which(rsk$Position!=rsk$alignedPosition)[10],]  ### G3P, checks out



######## merged by unaligned site ######

rsk$gene_site <-paste(toupper(rsk$gene), rsk$Position, sep="_")
### make small subset dataframe
rsks<-rsk[,c("gene_site","Protein", "Sirt5log2FC")]
colnames(rsks)[2]<-"uniprot"

m_unaligned <- merge(fsks, rsks, by=c("gene_site"))

nrow(m_unaligned)

getgene=function(x){
  unlist(strsplit(x, "_"))[1]
}


genes_unal<-unlist(lapply(FUN=getgene, X=m_unaligned$gene_site))
length(unique(genes_unal))

############ merged by aligned site  ################
## update the mini dataframe
rsk$gene_site <-paste(toupper(rsk$gene), rsk$alignedPosition, sep="_")
### make small subset dataframe
rsks<-rsk[,c("gene_site","Protein", "Sirt5log2FC")]
colnames(rsks)[2]<-"uniprot"

malign <- merge(fsks, rsks, by=c("gene_site"))

nrow(malign)  ## 238 sites in both comparisons

genes_al<-unlist(lapply(FUN=getgene, X=malign$gene_site))
length(unique(genes_al))  ## 102 proteins shared

### number of unique genes in the Park et al dataset
length(unique(rsk$gene)) ## 457
##
### num uniprot ids in park et al
length(unique(rsk$Protein))
### number of unique genes in the human fibro
length(unique(fsk$firstgene)) # 353
##
install.packages("VennDiagram")
library(VennDiagram)

??VennDiagram
dev.off()
#svg("D:/NatureComms_succK/R/outputs/3A_protien_overlap.svg",width=5, height=5)
draw.pairwise.venn(353, 457, 102, c("SCL","SIRT5"), lwd=4, cex=c(3,3,3))
dev.off()
nrow(fsk) # 933 in human fibro
length(unique(rsk$gene_site)) ## 1080 sites in Park et al
dev.off()
#svg("D:/NatureComms_succK/R/outputs/3A_site_overlap.svg",width=5, height=5)
draw.pairwise.venn( 933,1080, 238, c("SCL","SIRT5"), lwd=4, cex=c(3,3,3))
dev.off()
#2 for x in gshared:
#  align sequences
#  compute new numbering for mouse sites
# 
tca_genes<-c("Cs", "Aco2", "Idh2", "Idh3a", "Idh3g", "Dld", "Dlst", "Suclg1", "Sucla2", 
             "Suclg2", "Sdha", "Fh", "Mdh")
malign$genes<-genes_al
head(malign)
tca_overlap <- malign[malign$genes %in% toupper(tca_genes),]

head(tca_overlap)
## convert to long format
tca1<-cbind(tca_overlap[,c("gene_site", "log2FC", "genes")], rep("SCL", nrow(tca_overlap)))
colnames(tca1)[4]<-"group"

tca2<-cbind(tca_overlap[,c("gene_site", "Sirt5log2FC", "genes")], rep("SIRT5", nrow(tca_overlap)))
colnames(tca2)[2]<-"log2FC"
colnames(tca2)[4]<-"group"

tca<-rbind(tca1, tca2)
typeof(tca$log2FC)

library(ggplot2)


p<-ggplot(tca, aes(x=genes, y=log2FC, fill=group))+
  geom_boxplot()+
  geom_point(pch = 21, position = position_jitterdodge(), size=2, color="black")
p





### set jitter val manually
set.seed(2016)
jitterVal <- runif(nrow(tca), max = 0.5)
jitterVal <- jitterVal * ifelse(tca$group == "SCL", -1, +1)



ggplot(tca, aes(x=as.numeric(as.factor(tca$genes))+jitterVal, y=log2FC, fill=group, color=gene_site)) +
  geom_point() +
  geom_line(aes(group = gene_site))+ 

  scale_x_continuous(breaks = c(1,2,3,4,5,6,7,8,9), labels = levels(factor(tca$genes))) +
  xlab(NULL) +
  expand_limits(x = c(0.5, 3.5))


ggplot(tca, aes(x=as.numeric(as.factor(tca$genes))+jitterVal, y=log2FC, fill=group, color=gene_site)) +
  geom_boxplot(aes(x))+
  geom_point() +
  geom_line(aes(group = gene_site))+ 
  scale_x_continuous(breaks = c(1,2,3,4,5,6,7,8,9), labels = levels(factor(tca$genes))) +
  xlab(NULL) +
  expand_limits(x = c(0.5, 3.5))


########## figure 3E #############
top_genes<-c("ATP5A1", "ACADVL", "DLD", "ENO1", "FLNA", "HADHA", "HSPD1",  "IDH2",
             "LMNA", "MDH2", "MYH9", "VIM")
top_prot <- malign[malign$genes %in% toupper(top_genes),]
top1<-cbind(top_prot[,c("gene_site", "log2FC", "genes")], rep("SCL", nrow(top_prot)))
colnames(top1)[4]<-"group"

top2<-cbind(top_prot[,c("gene_site", "Sirt5log2FC", "genes")], rep("SIRT5", nrow(top_prot)))
colnames(top2)[2]<-"log2FC"
colnames(top2)[4]<-"group"

top<-rbind(top1, top2)
unique(top$genes)
head(top)
### set jitter val manually
set.seed(2016)
jitterVal <- runif(nrow(top), max = 0.1)
jitterVal <- jitterVal * ifelse(top$group == "SCL", -5, +5)

#install.packages("ggthemes")


### supplemental with all the sites connected
ggplot(top, aes(x=as.numeric(as.factor(top$genes))+jitterVal, 
                y=log2FC, fill=group, color=gene_site,shape=group)) +
  geom_point(color="black") +
  geom_line(aes(group = gene_site))+ 
  
  scale_x_continuous(breaks = c(1,2,3,4,5,6,7,8,9,10,11), labels = levels(factor(top$genes))) +
  xlab(NULL) +
  expand_limits(x = c(0.5, 3.5))+theme_Publication()



### just boxplots FIGURE 3E
p<-ggplot(top, aes(x=genes, y=log2FC, fill=group))+
  geom_boxplot()+
  geom_point(pch = 21, position = position_jitterdodge(), size=2, color="black")
p +theme_Publication()


## write table for source data
write.table(top, "D:/NatureComms_succK/R/outputs/source_data_3e.txt", row.names = F, sep="\t", quote=F)



### what about all the sites?


## convert to long format
mall1<-cbind(malign[,c("gene_site", "log2FC", "genes")], rep("SCL", nrow(malign)))
colnames(mall1)[4]<-"group"

mall2<-cbind(malign[,c("gene_site", "Sirt5log2FC", "genes")], rep("SIRT5", nrow(malign)))
colnames(mall2)[2]<-"log2FC"
colnames(mall2)[4]<-"group"

mall<-rbind(mall1, mall2)

p<-ggplot(mall, aes(x=genes, y=log2FC, fill=group))+
  geom_boxplot()+
  geom_point(pch = 21, position = position_jitterdodge(), size=2, color="black")+
  theme(axis.text.x = element_text(angle = 90,size =6), 
        panel.background = element_rect(fill="white", 
                                        color="grey", 
                                        size=0.5, linetype="solid")) 
p

p<-ggplot(mall, aes(x=genes, y=log2FC, fill=group))+
  geom_boxplot()+
  geom_point(pch = 21, position = position_jitterdodge(), size=2, color="black")
p +theme_Publication()

######## check panel B, log2fold changes overall
y1<-na.omit(rsk$Sirt5log2FC)
rsk$gene_site
b<-data.frame(rsk$Sirt5log2FC, rep("SIRT5", length(rsk$Sirt5log2FC)), rsk$gene_site)
head(b)
b<-b[is.na(b[,1])==FALSE,]
colnames(b)[1]<-"log2FC"
colnames(b)[2]<-"group"
colnames(b)[3]<- "gene_site"
head(b)
b1<-data.frame(fsk$log2FC, rep("SCL", nrow(fsk)), fsk$gene_site)
colnames(b1)[1]<-"log2FC"             
colnames(b1)[2]<-"group"             
colnames(b1)[3]<-"gene_site"             

b<-rbind(b, b1 )
head(b)
nrow(b) ### should be 1908 rows
dev.off()
ggplot(b, aes(y=log2FC, x=factor(group, levels=c("SCL", "SIRT5")), fill=factor(group, levels=c("SCL", "SIRT5"))))+
  geom_boxplot() 


+theme_Publication()


boxplot(log2FC~group, b)

write.table(b, "D:/NatureComms_succK/R/outputs/source_data_3c.txt", sep="\t")
########### write table #########33

write.table(malign, "D:/NatureComms_succK/R/outputs/sites_overlap.txt", sep="\t")

?setDT
malign
setDT(malign)[,freq := .N, by=c("genes")]
head(malign)
malign_sorted<-malign[order(freq, decreasing=T),]
write.table(malign_sorted, "D:/NatureComms_succK/R/outputs/sites_overlap_sorted.txt", sep="\t")
