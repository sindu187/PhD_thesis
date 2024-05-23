#

#complete gene information
Chr_genes<-as.data.frame(read.csv("/home/sinduja/work/rdata/geneData.csv"))

#installing packages
library("DESeq2") 
library(ggplot2)
library(RColorBrewer)
library(stringr)
library(dplyr)
library(ggpubr)
library(pheatmap)

##### peiferdata #########

load("/home/sinduja/work/work/tpm.RData") #tpm values of genes
peiferdata<-read.csv("/home/sinduja/work/work/breakpointdata/peifer_reevaluated.csv") # input peifer information
peiferdata<-peiferdata[peiferdata$CLASS!="whole" & peiferdata$CLASS!="mbp"  & peiferdata$CLASS!="outlier",] # removing outlier samples
peiferdata<-peiferdata[peiferdata$TRANSLOCATION=="YES",] # samples containg the translocation
peifersamp<-peiferdata$SAMPLE



#DESeq analysis 

load("/home/sinduja/work/rdata/peifer15.RData") # rawdata counts
nx<-exprs(normX)
nx=2^nx
dim (nx)

#genecount filtration
nx_primaryfilter<-nx[apply(nx,1,function(x)length(x[as.integer(x)>3])>=2),]# filtration taking more than 3 read count counts in 2 samples 
nx_secondaryfilter<-as.data.frame(nx_primaryfilter[(rowVars(nx_primaryfilter)>1),])
total_gene<-nx_secondaryfilter
total_gene<-total_gene[,sort(colnames(total_gene))]
total_gene<-total_gene[,colnames(total_gene) %in% peiferdata$SAMPLE]
sample<-peiferdata$SAMPLE
total_gene <-data.matrix(total_gene)
mode(total_gene)<-"integer"
coldata<-DataFrame(translocation= peiferdata$TRANSLOCATION,mycn=peiferdata$MYCN_status, type=peiferdata$TYPE,gain=peiferdata$X17QGAIN, loss= peiferdata$X11QLOSS)
rownames(coldata)<-sample
coldata$translocation<-factor(coldata$translocation,levels = c("NO","YES"))
coldata$mycn<-factor(coldata$mycn,levels=c(0,1,2))
coldata$class<-factor(coldata$class,levels=c("Translocation_NMYCN","Indirect_translocation","Translocation_LR","NMYCN_17qgain","MYCN_17qgain","No_translocation_LR","No_translocation_NMYCN","No_translocation_MYCN"))

dds<-DESeqDataSetFromMatrix(countData = total_gene,colData = coldata,design =~ translocation) #DESeq analysis
dds<-DESeq(dds,fitType = 'local')

#PCA plot
vds<-vst(dds)
pca<-plotPCA(vdd, intgroup=c("translocation","type","gain","loss","class"), returnData=F)
ggplot(pca$data, aes(pca$data$PC1, pca$data$PC2, color=class, shape=gain)) +
  geom_point(size=5)+ theme_pubr() +
  labs(x = pca$labels$x, y = pca$labels$y) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "bold.italic", size = 10)) +
  theme(axis.text.y = element_text(face = "bold.italic", size = 10)) +
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 15, b = 0, l = 0))) + 
  theme(axis.title.x = element_text(margin = margin(t = 15, r = 0, b = 0, l = 0))) +
  theme(strip.text = element_text(size=15)) +
  # theme(legend.title = element_blank()) +
  theme(axis.title.x = element_text(size = 12, face = "bold")) +
  theme(axis.title.y = element_text(size = 12, face = "bold"))
dev.off()



###peiferall###
result <- results(dds)
peiferall<-result[!is.na(result$padj),]
peiferall$external_gene_name<-rownames(peiferall)
rownames(peiferall)<-NULL
peiferall<-as.data.frame(peiferall)
peiferall_genes<-inner_join(peiferall,Chr_genes, by="external_gene_name")
pall_remove_chr1<-peiferall_genes[grep("^CHR",peiferall_genes$chromosome_name),]
peiferallgenes<-peiferall_genes[!peiferall_genes$chromosome_name %in% pall_remove_chr1$chromosome_name,]

# removing confounding by mycn amplification
load("/home/sinduja/work/work/rnadata/peifer/DDS.RData")
result <- results(DDS)
mycngene<-result[result$padj<0.05 & !is.na(result$padj),] #signficant genes between samples with and without mycn amplification
peiferall_nomycngene<-peiferallgenes[!peiferallgenes$external_gene_name %in% rownames(mycngene),]

#signficant genes
peifer_sig<-peiferall_nomycngene[peiferall_nomycngene$padj<0.05,] #significant genes between samples with and without translocation


