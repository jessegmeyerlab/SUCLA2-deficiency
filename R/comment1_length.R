
library(UniProt.ws)
library(stringr)
sk<-read.delim("D:/NatureComms_succK/Rardin_succinyl_sites.txt", stringsAsFactors = F)

head(sk)
nrow(sk)
sku<-sk[match(unique(sk$uniprot_site), sk$uniprot_site),]  ## keep only first occurance of a site
nrow(sku)
length(unique(sku$uniprot_site))

up <- UniProt.ws(taxId=10090)
taxId(up) <- 10090
columns(up)
keytypes(up)

keys <- unique(sku$Accession..)
nchar(keys)
keys<-keys[nchar(keys)>1]
onekey<-sku$Accession..[2]
onekey
columns <- c("LENGTH","SEQUENCE","GENES")
kt <- "UNIPROTKB"

res <- select(up, keys, columns, kt)


res


## get number of sites per protein
site_list<-list()
nsites<-c()
for (x in unique(sku$Accession..)){
  if(nchar(x)>1){
    nsites<-c(nsites,length(which(sku$Accession..==x)))
  }
}
  

unique(sku$Accession..)[1:5]
res$UNIPROTKB[1:5]
res$nsites<-nsites

res$nsites_per_len<-nsites/as.numeric(res$LENGTH)
num_k<-c()

for(x in res$SEQUENCE){
  num_k <- c(num_k, str_count(x, pattern="K") )
}

num_k
length(num_k)
length(nsites)

res$num_k<-nsites/num_k
dev.off()
plot(num_k, as.numeric(res$LENGTH))
nsites
length(res$LENGTH)
length(unique(sku$Accession..))

first_gene = function(x)unlist(strsplit(x, " "))[1]
res$firstgene = unlist(lapply(FUN=first_gene, res$GENES))
head(res)






### add protein quantities to the df
p<-read.delim("D:/NatureComms_succK/rardin_protein_quant.txt", stringsAsFactors = F)
head(p)

x<-unique(p$Accession..)[1]
protein_quants<-list()
for(x in unique(p$Accession..)){
  if(nchar(x)>1){
    protein_quants[[x]]<-sum(p[which(p$Accession..==x),"Mean.KO"])
  }
  
}

length(protein_quants)

prot_quant_vec<-c()

for(x in res$UNIPROTKB){
  print(protein_quants[[x]])
  print(length(protein_quants[[x]]))
  if(length(protein_quants[[x]])==1){
    prot_quant_vec<-c(prot_quant_vec, protein_quants[[x]])
  }
  if(length(protein_quants[[x]])==0){
    prot_quant_vec<-c(prot_quant_vec, NA)
  }
}


prot_quant_vec

res$protein_quant<-prot_quant_vec

res$nsites_over_protquant<- res$nsites/res$protein_quant


tca_genes<-c("Cs", "Aco2", "Idh2", "Idh3a", "Idh3g", "Dld", "Dlst", "Suclg1", "Sucla2", 
             "Suclg2", "Sdha", "Fh", "Mdh")



#write.table(res, "D:/NatureComms_succK/R/outputs/rardin_succinyl_norm.txt", 
#            sep="\t", row.names = F)





############ remove the genes from TCA cycle to another table ###########

nrow(res)
length(tca_genes)
wo_tca<-res[-which(res$firstgene %in% tca_genes),]
par(cex=1.4)

## sites norm to protein quant, both distributions
hist(log2(wo_tca$nsites_over_protquant), breaks=15,
     main="Normalized to protein LFQ",
     xlab="# sites/protein abundance",
     xlim=c(-23,-10))

hist(log2(res[res$firstgene %in% tca_genes,"nsites_over_protquant"]), 
     add=T, col="red", breaks=4)

text(-21,20, "p-value = 6E-4")

legend(-15,30, legend=c("all proteins", "TCA proteins"), fill=c("white", "red"))
wilcox.test(log2(wo_tca$nsites_over_protquant), 
       log2(res[res$firstgene %in% tca_genes,"nsites_over_protquant"]))



#### divided by length?
hist(log2(wo_tca$nsites_per_len), breaks=15,
     main="Normalized to protein length",
     xlab="# sites/protein length", xlim=c(-11, 2))

hist(log2(res[res$firstgene %in% tca_genes,"nsites_per_len"]), 
     add=T, col="red", breaks=5)
legend(-4,28, legend=c("all proteins", "TCA proteins"), fill=c("white", "red"))

wilcoxres<-wilcox.test(log2(wo_tca$nsites_per_len), 
            log2(res[res$firstgene %in% tca_genes,"nsites_per_len"]))

text(-0.5,7, paste("p-value =",round(wilcoxres$p.value,4) ))



#### divided by #lysine?
hist(log2(wo_tca$num_k), breaks=15,
     main="Normalized to # lysines",
     xlab="# sites/ # protein lysines", xlim=c(-9, 5))

hist(log2(res[res$firstgene %in% tca_genes,"num_k"]), 
     add=T, col="red", breaks=5)
legend(-0.5,28, legend=c("all proteins", "TCA proteins"), fill=c("white", "red"))

wilcoxres<-wilcox.test(log2(wo_tca$num_k), 
                       log2(res[res$firstgene %in% tca_genes,"num_k"]))

text(-7,20, paste("p-value =",round(wilcoxres$p.value,4) ))


