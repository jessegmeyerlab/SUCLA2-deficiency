
library(UniProt.ws)
library(stringr)
library(ggplot2)
sk<-read.delim("D:/NatureComms_succK/rardin_sK_site_quant.txt", stringsAsFactors = F)

head(sk)
nrow(sk)

sk$uniprot_site<-paste(sk$Accessions, sk$Site, sep = "_")

# is there only one peptide per site?
length(unique(sk$uniprot_site)) == length(sk$uniprot_site)


tca_genes<-c("Cs", "Aco2", "Idh2", "Idh3a", "Idh3g", "Dld", "Dlst", "Suclg1", "Sucla2", 
             "Suclg2", "Sdha", "Fh", "Mdh")


only_tca<- sk[which(sk$Gene.Name %in% tca_genes),]

nrow(only_tca)
only_tca$Normalized.Ratio.KO.WT
only_tca$Gene.Name
only_tca$uniprot_site
p<-ggplot(only_tca, 
          aes(x=Gene.Name, y=log2(Normalized.Ratio.KO.WT)))+
  geom_boxplot(size=2)
p

tlm<-lm(only_tca$Normalized.Ratio.KO.WT ~ only_tca$Gene.Name)
tmpanova<-Anova(tlm)

summary(fm1<-aov(log2(Normalized.Ratio.KO.WT) ~ Gene.Name, data=only_tca))
write.table(file="D:/NatureComms_succK/R/outputs/comment2_tukeyHSD.txt", 
            sep="\t", col.names = T, row.names = T,
            data.frame(TukeyHSD(fm1, "Gene.Name", ordered=TRUE)$Gene.Name))
### print table 
data.frame(TukeyHSD(fm1, "Gene.Name", ordered=TRUE)$Gene.Name)
plot(TukeyHSD(fm1, "Gene.Name", ordered=TRUE))

## save the table of comparisons for philipp 
