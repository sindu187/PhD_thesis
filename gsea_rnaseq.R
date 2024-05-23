BiocManager::install("Genomicsranges",force=TRUE)
BiocManager::install("pathview",lib ="/home/sinduja/R/x86_64-pc-linux-gnu-library/4.1")
BiocManager::install("devtools",lib ="/home/sinduja/R/x86_64-pc-linux-gnu-library/4.1")
BiocManager::install("ggraph")
BiocManager::install("clusterProfiler", force = TRUE)
BiocManager::install("wordcloud", force = TRUE)
BiocManager::install("ggplotly")
BiocManager::install("enrichplot")
BiocManager::install("ggpubr")
BiocManager::install("ggridges")
BiocManager::install("tidyverse")
BiocManager::install("genefilter", force=TRUE)
BiocManager::install("freetype2", force=TRUE)

library(dplyr)
library(ggplotly)
library(RCurl)
library(GenomicRanges)
library(ggtree)
library(clusterProfiler)
library(enrichplot)
library(pheatmap)

library(devtools)
library(ggplot2)
library(ggnewscale)
library(ggridges)
library(tidyverse)
library(pathfindR)


# SET THE DESIRED ORGANISM HERE
organism = "org.Hs.eg.db"
#BiocManager::install(organism)

library(organism, character.only = TRUE)
library(clusterProfiler)
library(wordcloud)
library(DOSE)


#Input preparation
# reading in data from deseq2 analysis
peiferall<-read.csv("/home/sinduja/work/work/breakpointdata/peiferallgenes_nomycn.csv")
df<-peiferall

# we want the log2 fold change 
original_gene_list <- df$log2FoldChange
# name the vector
names(original_gene_list) <- df$ensembl_gene_id
# omit any NA values 
gene_list<-na.omit(original_gene_list)

# sort the list in decreasing order (required for clusterProfiler)

gene_list = sort(gene_list, decreasing = TRUE)

sig_genes_df = subset(df, padj < 0.05)

genes <- sig_genes_df$log2FoldChange
names(genes) <- sig_genes_df$ensembl_gene_id
genes <- na.omit(genes)
genes<- names(genes)

#gene ontology  over-representation analysis

go_enrich <- enrichGO(gene = genes,
                      universe = names(gene_list),
                      OrgDb = organism, 
                      keyType = 'ENSEMBL',
                      readable = T,
                      ont = "BP",
                      pvalueCutoff = 0.05, 
                      qvalueCutoff = 0.10)


#keys
head(keys(org.Hs.eg.db, keytype="GENETYPE"))

#plots
upsetplot(go_enrich)
dotplot(go_enrich)
cnetplot(go_enrich, categorySize="pvalue", foldChange=gene_list)

selected_pathways= go_enrich$Description[5]
sample(go_enrich$Description,c(4:7))
cnetplot(go_enrich,showCategory = selected_pathways)




#gene ontology analysis

gse1 <- gseGO(geneList=gene_list, 
              ont ="BP", 
              keyType = "ENSEMBL",
              eps = 0 ,
              OrgDb = org.Hs.eg.db,
              minGSSize = 10, 
              maxGSSize = 800, 
              pvalueCutoff = 0.05,
              verbose = TRUE,
              pAdjustMethod = "BH")

 #disease enrichment analysis 
genes1<-bitr(genes, "ENSEMBL","ENTREZID" , organism, drop = TRUE) #convert to entrez
dotplot(gse1)
de <- names(gene_list)[abs(gene_list) > 1] #fold change >1

genes_list1<-bitr(de, "ENSEMBL","ENTREZID" , organism, drop = TRUE) 


edo <-enrichDGN(
  gene=genes1$ENTREZID,
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  universe = genes_list1$ENTREZID  ,
  minGSSize = 10,
  maxGSSize = 500,
  qvalueCutoff = 0.2,
  readable = FALSE
)

edo <-enrichDGN(genes_list1$ENTREZID)
barplot(edo, showCategory = 20)



#KEGG

# Convert gene IDs for gseKEGG function
# We will lose some genes here because not all IDs will be converted
ids<-bitr(names(original_gene_list), fromType = "ENSEMBL", toType = "ENTREZID", OrgDb=organism)
# remove duplicate IDS (here I use "ENSEMBL", but it should be whatever was selected as keyType)
dedup_ids = ids[!duplicated(ids[c("ENSEMBL")]),]

# Create a new dataframe df2 which has only the genes which were successfully mapped using the bitr function above
df2 = df[df$ensembl_gene_id %in% dedup_ids$ENSEMBL,]

# Create a new column in df2 with the corresponding ENTREZ IDs
df2$Y = dedup_ids$ENTREZID

# Create a vector of the gene unuiverse
kegg_gene_list <- df2$log2FoldChange

# Name vector with ENTREZ ids
names(kegg_gene_list) <- df2$Y

# omit any NA values 
kegg_gene_list<-na.omit(kegg_gene_list)

# sort the list in decreasing order (required for clusterProfiler)
kegg_gene_list = sort(kegg_gene_list, decreasing = TRUE)


kegg_organism = "hsa"
kk2 <- gseKEGG(geneList     = kegg_gene_list,
               organism     = kegg_organism,
               nPerm        = 10000,
               minGSSize    = 3,
               maxGSSize    = 800,
               pvalueCutoff = 0.05,
               pAdjustMethod = "none",
               keyType       = "ncbi-geneid")
#svg("/home/sinduja/work/work/GSEA_1/filteredmycn_kegg.svg", width=12,height=10)
dotplot(kk2, showCategory = c(15), title = "Enriched Pathways" , split=".sign") + facet_grid(.~.sign)
p1 <- gseaplot2(kk2, geneSetID = 1, title = kk2$Description[42])
kk3 <- setReadable(kk2,'org.Hs.eg.db', 'ENTREZID')
cnetplot(kk3, showCategory=kk3$Description[42],Size="pvalue", circular = TRUE, colorEdge = TRUE,foldChange=kegg_gene_list)
p1
#dev.off()


# upregulated and downregulated
upregulated_genes=df[df$log2FoldChange>0.00,]
upregulated_p<-upregulated_genes[order(upregulated_genes$stat,decreasing = TRUE),]
upregulated_b1<-upregulated_genes[order(upregulated_genes$log2FoldChange,decreasing = TRUE),]
downregulated_genes=df[df$log2FoldChange<=0.00,]
downregulated<-downregulated_genes[order(downregulated_genes$stat),]

upregulated_genes11<-upregulated_genes[upregulated_genes$chromosome_name==11,]
upregulated11<-upregulated_genes11[order(upregulated_genes11$log2FoldChange,decreasing = TRUE),]
ggplot(upregulated11,aes(start_position))+geom_bar(fill = "#0073C2FF")+labs(x = "genomeposition",y="Frequency",title="Genomeposition of upregulated chr11 genes")



qplot(downregulated_genes2$start_position,
      geom="histogram",
      binwidth = 10000000,  
      main = "chr2 deletion", 
      xlab = "Age",  
      fill=I("blue"), 
      col=I("red"), 
      alpha=I(.2),
      xlim=c(0,2.5e+08), ylim=c(0,15))

