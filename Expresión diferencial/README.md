# Análisis de expresión diferencial

## Creamos la matriz

```
module load conda
source activate /home/fabianc.salgado/data/POCtrinity
module load R/3.6.3

abundance_estimates_to_matrix.pl --est_method RSEM \
--gene_trans_map /home/fabianc.salgado/shared/paula_torres/gasteracantha/trinity/trinity_without_2000/Trinity_2000.fasta.gene_trans_map \
--out_prefix Trinity_trans_2000 /home/fabianc.salgado/shared/paula_torres/gasteracantha/mapeo/mapeo_2000/RSEM.isoforms.results/*isoforms.results
```

## EdgeR

```
library(edgeR)
library(limma)
library(Glimma)
library(DESeq2)
library(ggrepel)
library(gplots)
library(RColorBrewer)
x <- read.table("/home/paula/Documents/gasteracantha/edgeR/2000/matriz/Trinity_trans_2000.isoform.counts.matrix", header=T)
x<-as.matrix(x)
##Analysis with edgeR
group<-factor(c("O","O","O","W","W","W","Y","Y","Y"))
y<- DGEList(counts=x,group=group)
levels(y$samples$group)
[1] "O" "W" "Y"
design <- model.matrix(~0+group, data=y$samples)
colnames(design) <- levels(y$samples$group)
design
y<-calcNormFactors(y)
y<-estimateDisp(y,design)
fit <- glmQLFit(y, design)
my.contrasts <- makeContrasts(OvsW=O-W, WvsO=W-O, WvsY= W-Y, YvsW= Y-W, OvsY= O-Y, YvsO=Y-O, levels=design)
qlf.OvsW <- glmQLFTest(fit, contrast=my.contrasts[,"OvsW"])
#solucion
fit <- glmFit(y, design)
> qlf.OvsW <- glmLRT(fit, contrast=my.contrasts[,"OvsW"])
##para saber cuantos transcritos hay
summary(decideTests(qlf.OvsW,p.value=0.05))
topTags(qlf.OvsW)

##GRAFICOS
#PCA
points <- c(0,1,2)
colors <- c("blue", "darkgreen", "red")
#pdf("plotMDS_2.pdf")
plotMDS(y, col=colors[group],pch=points[group])
legend("topleft",legend=levels(group),pch=points,col=colors)
dev.off()

points<-23y
col<-c("orange","white","yellow")
col<-rep(col,each=3)
coge<-c("O","W","Y")
pdf("pca_all_snps.pdf")
plot(pc$x,pch=points,col="black",bg=col,cex=2)
legend("topright",legend=coge,pch=points,col="black",pt.bg= c("orange","white","yellow"))
dev.off()



plotBCV(y)
pdf("plotBCV.pdf")
plotBCV(y)
dev.off()

pdf("glmQLFit.pdf")
plotQLDisp(fit)
dev.off()

pdf("qlf.OvsW.pdf")
plotMD(qlf.OvsW)
dev.off()

pdf("qlf.OvsY.pdf")
plotMD(qlf.OvsY)
dev.off()

pdf("qlf.WvsY.pdf")
plotMD(qlf.WvsY)
dev.off()

#volcano_plot
topTagsOvsW<-topTags(qlf.OvsW)
volcanoOvsW<-cbind(topTagsOvsW$table$logFC, -log10(topTagsOvsW$table$FDR))
colnames(volcanoOvsW) <-c("logFC", "negLogPval")
head (volcanoOvsW)
pdf("volcanoOvsW.pdf")
plot(volcanoOvsW,pch=19)
dev.off()

#Volcano_real
v <- voom(y,design,plot = TRUE)
fit <- lmFit(v)
fit.cont<-contrasts.fit(fit,my.contrasts[,"OvsW"])
fit.cont <- eBayes(fit.cont)
dim(fit.cont)
#pdf("volcanoprueba.pdf")
volcanoplot(fit.cont,coef=1,highlight=100)
dev.off()
volcanoplot(fit.cont,coef=1,highlight=100,names=rownames(v$E))


#heatmap500
logcounts <- cpm(y,log=TRUE)
var_genes <- apply(logcounts, 1, var)
select_var <- names(sort(var_genes, decreasing=TRUE))[1:100]
highly_variable_lcpm <- logcounts[select_var,]
dim(highly_variable_lcpm)
mypalette <- scale_fill_viridis(discrete = TRUE)
morecols <- colorRampPalette(mypalette)
pdf("heatmap.pdf")
heatmap.2(highly_variable_lcpm,col=rev(morecols(50)),trace="none", main="Top 500 most variable genes across samples",key=FALSE)
dev.off()

##heatmap bonito
library(heatmaply)
heat_df<-as.data.frame(highly_variable_lcpm)
heatmaply(heat_df,showticklabels=c(TRUE,FALSE),margins = c(10, 10, 20, 10),file="heatmaply_plot.png",width=600,height=800)
#----------------------------------------------------------------------------------------------------------------------------------------------
##Analysis with DeSeq2
condition<-rep(c("Orange","White","Yellow"),each=3)
type<-rep("paired-end",9)
names<-c()
for (i in seq(1:9)){
        names<-c(names,gsub(" ", "",paste("PTRUR00",i)))
}
coldata<-data.frame(condition,type)
rownames(coldata)<-names
##calculate diffexpressed
dds <- DESeqDataSetFromMatrix(countData = round(x),
                              colData = coldata,
                              design = ~ condition)
dds<-DESeq(dds)
##make contrast
res_OvsY <- results(dds contrast=c("condition","Orange","Yellow"))
res_YvsW <- results(dds, contrast=c("condition","Yellow","White"))
res__OvsW <- results(dds, contrast=c("condition","Orange","White"))
##organizamos dataframe para hacer volcanoplot
YvsW<-data.frame(gene_symbol=rownames(res_YvsW),log2FoldChange=res_YvsW$log2FoldChange,pvalue=res_YvsW$pvalue)
YvsW$diffexpressed <- "NO"
YvsW$diffexpressed[YvsW$log2FoldChange > 0.6 & YvsW$pvalue < 0.05] <- "UP"
qlf_OvsW_table_id$diffexpressed[qlf_OvsW_table_id$log2FoldChange > 0.6 & qlf_OvsW_table_id$PValue < 0.05] <- "UP"
# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
YvsW$diffexpressed[YvsW$log2FoldChange < -0.6 & YvsW$pvalue < 0.05] <- "DOWN"
qlf_OvsW_table_id$diffexpressed[qlf_OvsW_table_id$log2FoldChange < -0.6 & qlf_OvsW_table_id$PValue < 0.05] <- "DOWN"

YvsW$delabel <- NA
YvsW$delabel[YvsW$diffexpressed != "NO"] <- YvsW$gene_symbol[YvsW$diffexpressed != "NO"]
##añadir gene_id al data frame 
#Cargamos los datos de los id
blastx_id<-read.table("/home/paula/Documents/gasteracantha/anotation/trinity_nr_blastx_id.txt",sep="\t",head=F)
colnames(blastx_id)<-c("gene_symbol","gene_id")
#añadimos la columna del id a 
OvsW_with_id<-merge(OvsW,blastx_id, by= "gene_symbol",all.x=TRUE)
OvsW_with_id$gene_id[OvsW_with_id$diffexpressed == "NO"] <- NA
##Extraer transcripts id
##Extraer transcripts id
YvsW_Up<-YvsW %>% filter(diffexpressed=='UP') 
rownames(YvsW_Up)<-YvsW_Up$gene_symbol
gene_symbol<-rownames(YvsW_Up)
write(gene_symbol,"YvsW_transcript_id_up.txt")
#extraer secuencias fasta de los up
for j in $(cat YvsW_transcript_id_up.txt)    ; do grep -A1 "$j" Trinity_2000.fasta >> transcripts_YvsW.fasta ; done

##organizacion de datos para volcano con edgeR results
qlf_YvsW_table_id<- qlf.YvsW$table    
qlf_YvsW_table_id$gene_symbol<-rownames(qlf_YvsW_table_id)
qlf_YvsW_table_id$diffexpressed <- "NO"
qlf_YvsW_table_id$diffexpressed[qlf_YvsW_table_id$logFC > 1 & qlf_YvsW_table_id$PValue < 0.05] <- "UP"
qlf_YvsW_table_id$diffexpressed[qlf_YvsW_table_id$logFC < -1 & qlf_YvsW_table_id$PValue < 0.05] <- "DOWN"
qlf_YvsW_table_id<-merge(qlf_YvsW_table_id,blastx_id, by= "gene_symbol",all.x=TRUE)
qlf_YvsW_table_id$gene_id[qlf_YvsW_table_id$diffexpressed == "NO"] <- NA


##volcano_plot
mynamestheme <- theme(plot.title = element_text(size = (16),hjust=0.5))
volcano_OvsY_with_id<-ggplot(data=OvsY_with_id, aes(x=log2FoldChange, y=-log10(pvalue), col=diffexpressed)) +
        geom_point() + 
        geom_text_repel(label=OvsY_with_id$gene_id) +
        xlim(-40,40)+
        theme(panel.background = element_rect(fill = "white",colour = "black")) +
        #geom_text_repel() +
        scale_color_manual(values=c("#FFCF20FF","gray62","#541352FF")) +
        geom_vline(xintercept=c(-0.6, 0.6),linetype = "dashed", col="black") +
        geom_hline(yintercept=-log10(0.05), linetype = "dashed", col="black") + mynamestheme +
        #annotate("rect", xmin = 3, xmax = 30, ymin = -log10(0.01), ymax = 16, alpha=.2, fill="#541352FF") +
        #annotate("rect", xmin = -5, xmax = -30, ymin = -log10(0.01), ymax = 16, alpha=.2, fill="#FFCF20FF") +
        ggtitle("Volcano plot Orange vs Yellow") 

volcano_YvsW_with_id<-ggplot(data=YvsW_with_id, aes(x=log2FoldChange, y=-log10(pvalue), col=diffexpressed)) +
        geom_point() + 
        geom_text_repel(label=YvsW_with_id$gene_id,box.padding = 0.5) +
        xlim(-40,40)+
        theme(panel.background = element_rect(fill = "white",colour = "black")) +
        #geom_text_repel() +
        scale_color_manual(values=c("#FFCF20FF","gray62","#541352FF")) +
        geom_vline(xintercept=c(-0.6, 0.6),linetype = "dashed", col="black") +
        geom_hline(yintercept=-log10(0.05), linetype = "dashed", col="black") + mynamestheme +
        ggtitle("Volcano plot Yellow vs White") 

volcano_OvsW_with_id<-ggplot(data=OvsW_with_id, aes(x=log2FoldChange, y=-log10(pvalue), col=diffexpressed)) +
        geom_point() + 
        geom_text_repel(label=OvsW_with_id$gene_id,box.padding = 0.5) +
        xlim(-40,40)+
        theme(panel.background = element_rect(fill = "white",colour = "black")) +
        #geom_text_repel() +
        scale_color_manual(values=c("#FFCF20FF","gray62","#541352FF")) +
        geom_vline(xintercept=c(-0.6, 0.6),linetype = "dashed", col="black") +
        geom_hline(yintercept=-log10(0.05), linetype = "dashed", col="black") + mynamestheme +
        ggtitle("Volcano plot Orange vs White") 

##qlf_volcano
qlf_volcano_OvsW_with_id<-ggplot(data=qlf_OvsW_table_id, aes(x=logFC, y=-log10(PValue), col=diffexpressed)) +
        geom_point() + 
        geom_text_repel(label=qlf_OvsW_table_id$gene_id.y) +
        xlim(-40,40)+
        theme(panel.background = element_rect(fill = "white",colour = "black")) +
        #geom_text_repel() +
        scale_color_manual(values=c("#FFCF20FF","gray62","#541352FF")) +
        geom_vline(xintercept=c(-0.6, 0.6),linetype = "dashed", col="black") +
        geom_hline(yintercept=-log10(0.05), linetype = "dashed", col="black") + mynamestheme +
        #annotate("rect", xmin = 3, xmax = 30, ymin = -log10(0.01), ymax = 16, alpha=.2, fill="#541352FF") +
        #annotate("rect", xmin = -5, xmax = -30, ymin = -log10(0.01), ymax = 16, alpha=.2, fill="#FFCF20FF") +
        ggtitle("Volcano plot Orange vs White") 

qlf_volcano_OvsY_with_id<-ggplot(data=qlf_OvsY_table_id, aes(x=logFC, y=-log10(PValue), col=diffexpressed)) +
        geom_point() + 
        geom_text_repel(label=qlf_OvsY_table_id$gene_id) +
        xlim(-40,40)+
        theme(panel.background = element_rect(fill = "white",colour = "black")) +
        #geom_text_repel() +
        scale_color_manual(values=c("#FFCF20FF","gray62","#541352FF")) +
        geom_vline(xintercept=c(-0.6, 0.6),linetype = "dashed", col="black") +
        geom_hline(yintercept=-log10(0.05), linetype = "dashed", col="black") + mynamestheme +
        #annotate("rect", xmin = 3, xmax = 30, ymin = -log10(0.01), ymax = 16, alpha=.2, fill="#541352FF") +
        #annotate("rect", xmin = -5, xmax = -30, ymin = -log10(0.01), ymax = 16, alpha=.2, fill="#FFCF20FF") +
        ggtitle("Volcano plot Orange vs Yellow") 

qlf_volcano_YvsW_with_id<-ggplot(data=qlf_YvsW_table_id, aes(x=logFC, y=-log10(PValue), col=diffexpressed)) +
        geom_point() + 
        geom_text_repel(label=qlf_YvsW_table_id$gene_id) +
        xlim(-40,40)+
        theme(panel.background = element_rect(fill = "white",colour = "black")) +
        #geom_text_repel() +
        scale_color_manual(values=c("#FFCF20FF","gray62","#541352FF")) +
        geom_vline(xintercept=c(-1, 1),linetype = "dashed", col="black") +
        geom_hline(yintercept=-log10(0.05), linetype = "dashed", col="black") + mynamestheme +
        #annotate("rect", xmin = 3, xmax = 30, ymin = -log10(0.01), ymax = 16, alpha=.2, fill="#541352FF") +
        #annotate("rect", xmin = -5, xmax = -30, ymin = -log10(0.01), ymax = 16, alpha=.2, fill="#FFCF20FF") +
        ggtitle("Volcano plot Yellow vs White") 

png("qlf_volvano_YvsW_id_high.png",width=880,height=480)
volcano_YvsW_with_id
dev.off()


pdf("volcano_OvsY.pdf",width=15,height=8)
volcano_OvsY
dev.off()

png("qlf_volcano_YvsW_with_id.png")
volcano_OvsY
dev.off()

write.table(heat_df, file = "heat_df.txt", sep = " ", row.names = T)
```
