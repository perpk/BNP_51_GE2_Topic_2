library("DESeq2")

df <- read.table('matrix.txt', header = TRUE, sep=",")
conditions <- as.factor(c("tumor", "tumor", "tumor", "normal", "normal", "tumor", "tumor", "tumor", "normal", "normal"))

install.packages("BiocManager")
BiocManager::install("DESeq2")

#split names and data
gene_names<- df[,1]
gene_data<- df[,2:11]
# εκτέλεση εργαλείου DESeq2
dsq<-DESeqDataSetFromMatrix(gene_data, DataFrame(conditions), ~ conditions)
dsq<-DESeq(dsq)
#εξαγωγή αποτελεσμάτων από αντικείμενο DESeqDataSet
res<-results(dsq, contrast = c("conditions", "normal" ,"tumor"))
#μετατροπή από DESQeqResults σε dataframe
res <- as.data.frame(res)
#Προσθήκη στήλης gene_names
res$genes<- gene_names
#αφαίρεση NA τιμών
filt1 <- !is.na(res$padj)
res2<- res[filt1,]
res4 <- res2[res2$padj < 10^-8,]

nrow(res2[res2$log2FoldChange < 0, ])
nrow(res2[res2$log2FoldChange > 0, ])
nrow(res2)

threshold <- 10^-8
res2$diffexpressed = res2$padj < threshold
library(ggplot2)
ggplot(data=res2, aes(x=log2FoldChange, y=-log10(padj), col=diffexpressed)) +
  geom_point() + theme_minimal()
