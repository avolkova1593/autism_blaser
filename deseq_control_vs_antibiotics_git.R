###This script provides comparison between trasciption of
###antibiotics (STAT + STAT-Birth) and controlconditioins
###Written by Angelina Volkova 01/14/2020

library(DESeq2)
library(ggplot2)
library(dplyr)
library(ComplexHeatmap)
library(circlize)
set.seed(1234)
setwd ("/Users/av1936/Desktop/Google_Drive/collaborative_projects/blaser")

file <- read.table("feature_counts.txt", sep="\t", header = T,row.names=1)
samples_fc <- read.csv("RNA_Seq_analysis/FC_samples/samplesheet_2_f.csv")
samples_hip <- read.csv("RNA_Seq_analysis/A_samples/samplesheet_2_a.csv")

####Get the names of the samples from the samplesheet files
a_control_ids <- gsub(" ","_",as.character(samples_hip[samples_hip$condition%in%c("CONTROL"),]$id))
a_stat_ids <- gsub(" ","_",as.character(samples_hip[samples_hip$condition%in%c("STAT"),]$id))
a_stat_birth_ids <- gsub(" ","_",as.character(samples_hip[samples_hip$condition%in%c("STAT-BIRTH"),]$id))

f_control_ids <- gsub(" ","_",as.character(samples_fc[samples_fc$condition%in%c("CONTROL"),]$id))
f_stat_ids <- gsub(" ","_",as.character(samples_fc[samples_fc$condition%in%c("STAT"),]$id))
f_stat_birth_ids <- gsub(" ","_",as.character(samples_fc[samples_fc$condition%in%c("STAT-BIRTH"),]$id))

a_all <- c(a_control_ids,a_stat_ids,a_stat_birth_ids)
f_all <- c(f_control_ids,f_stat_ids,f_stat_birth_ids)

##Remove numbers after the dots in the feature file
names(file) <- sub("\\..*","",names(file))

###Split the feature file into amygdala and fc
f_file <- file[,names(file)%in%c(f_all),]
f_samples <- c(names(f_file))
names(f_file) <- sub("\\..*","",names(f_file))

a_file <- file[,names(file)%in%c(a_all),]
a_samples <- c(names(a_file))
names(a_file) <- sub("\\..*","",names(a_file))

###Change the names of both data frames and get corresponding conditions
names(a_file)[names(a_file)%in%c(a_control_ids)]<- "cont"
names(a_file)[names(a_file)%in%c(a_stat_ids)]<- "antibiotics"
names(a_file)[names(a_file)%in%c(a_stat_birth_ids)]<- "antibiotics"
names(f_file)[names(f_file)%in%c(f_control_ids)]<- "cont"
names(f_file)[names(f_file)%in%c(f_stat_ids)]<- "antibiotics"
names(f_file)[names(f_file)%in%c(f_stat_birth_ids)]<- "antibiotics"
a_conditions <- c(names(a_file))
f_conditions <- c(names(f_file))
names(a_file) <- a_samples
names(f_file) <- f_samples

###Set samples data frames
a_df_samples <- data.frame(row.names = c(names(a_file)), 
                      vector = as.factor(a_conditions))
f_df_samples <- data.frame(row.names = c(names(f_file)), 
                           vector = as.factor(f_conditions))

###Make data frames deseq ready
a_df_deseq <- DESeqDataSetFromMatrix(countData = a_file, colData = a_df_samples, design = ~ vector)
f_df_deseq <- DESeqDataSetFromMatrix(countData = f_file, colData = f_df_samples, design = ~ vector)

###Remove genes with counts less than 2
a_df_deseq <- a_df_deseq[rowSums(counts(a_df_deseq) >= 2) >= 2]
f_df_deseq <- f_df_deseq[rowSums(counts(f_df_deseq) >= 2) >= 2]

###Run Deseq
a_processed <- DESeq(a_df_deseq)
f_processed <- DESeq(f_df_deseq)

###Get deseq results
a_results <- as.data.frame(results(a_processed, contrast = c("vector", "cont","antibiotics")))
f_results <- as.data.frame(results(f_processed, contrast = c("vector", "cont", "antibiotics")))

###Subset deseq results
a_new_results_sig <- a_results[a_results$padj<0.01,]
a_new_results_sig <- a_new_results_sig[!is.na(a_new_results_sig$padj),]
a_new_results_sig <- a_new_results_sig%>%mutate(direction = ifelse((log2FoldChange >= 0 & padj < 0.01), "upregulated",
                                                                             ifelse((log2FoldChange <= 0 & padj < 0.01),
                                                                                    "downregulated", "not significant"))) 

f_new_results_sig <- f_results[f_results$padj<0.01,]
f_new_results_sig <- f_new_results_sig[!is.na(f_new_results_sig$padj),]
f_new_results_sig <- f_new_results_sig%>%mutate(direction = ifelse((log2FoldChange >= 0 & padj < 0.01), "upregulated",
                                                                   ifelse((log2FoldChange <= 0 & padj < 0.01),
                                                                         "downregulated", "not significant"))) 

###Select the genes with logfc > 1.5
a_results_logfc <- a_new_results_sig[abs(a_results$log2FoldChange)>=1.5,]
f_results_logfc <- f_new_results_sig[abs(f_results$log2FoldChange)>=1.5,]

###Save the gene expression
write.table(a_new_results_sig, file="/Users/av1936/Desktop/a_volcano_logfc_anti.txt", 
            sep="\t", row.names = T, col.names = T,quote=F)
write.table(f_new_results_sig, file="/Users/av1936/Desktop/f_volcano_logfc_anti.txt", 
            sep="\t", row.names = T, col.names = T,quote=F)
write.table(a_results_logfc, file="/Users/av1936/Desktop/a_volcano_logfc1_5_anti.txt", 
            sep="\t", row.names = T, col.names = T,quote=F)
write.table(f_results_logfc, file="/Users/av1936/Desktop/f_volcano_logfc1_5_anti.txt", 
            sep="\t", row.names = T, col.names = T,quote=F)

###Get only significant results
a_results_no_na <- a_results_logfc[!is.na(a_results_logfc$padj),]
f_results_no_na <- f_results_logfc[!is.na(f_results_logfc$padj),]
a_results_sig <- a_results[a_results_no_na$padj<0.05,]
f_results_sig <- f_results[f_results_no_na$padj<0.05,]

###Calculate row means across vectors
a_mean_cont <- rowMeans(counts(a_processed, normalized=TRUE)[, a_processed$vector == "cont"])
a_mean_anti <- rowMeans(counts(a_processed, normalized=TRUE)[, a_processed$vector == "antibiotics"])
f_mean_cont <- rowMeans(counts(f_processed, normalized=TRUE)[, f_processed$vector == "cont"])
f_mean_anti <- rowMeans(counts(f_processed, normalized=TRUE)[, f_processed$vector == "antibiotics"])

###Plots

### Log Fold change
a_logfc <- lfcShrink(a_processed, contrast = c("vector", "cont", "antibiotics")) 
a_logfc <- cbind(as.data.frame(a_logfc), a_mean_cont, a_mean_anti)
f_logfc <- lfcShrink(f_processed, contrast = c("vector", "cont", "antibiotics")) 
f_logfc <- cbind(as.data.frame(f_logfc), f_mean_cont, f_mean_anti)

###Volcano plots
a_logfc_vol <- a_logfc%>%mutate(threshold = ifelse((log2FoldChange >= 0 & padj < 0.01), "up_sig",
                                                   ifelse((log2FoldChange <= 0 & padj < 0.01),
                                                   "down_sig", "not_sig"))) 
a_logfc_vol$threshold[is.na(a_logfc_vol$threshold)] <- "not_sig"

ggplot(a_logfc_vol, aes(x = log2FoldChange, y = -log10(pvalue))) + 
  geom_point(aes(colour = threshold), size=1) +
  ggtitle("Amygdala \nCONTROL vs STAT + STAT-BIRTH") +
  xlab("Shrunken LogFC")+
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_colour_manual(values=c("blue", "grey", "red3"),
                      name="Threshold 0.01", 
                      breaks=c("down_sig", "not_sig", "up_sig"), 
                      labels=c("Down", "None", "Up"))+
  theme_bw()+
  theme(axis.text=element_text(size=14,colour="black"),
        axis.title=element_text(size=14,colour="black"),
        text = element_text(size=14,colour="black"))


f_logfc_vol <- f_logfc%>%mutate(threshold = ifelse((log2FoldChange >= 0 & padj < 0.01), "up_sig",
                                                   ifelse((log2FoldChange <= 0 & padj < 0.01),
                                                          "down_sig", "not_sig"))) 
f_logfc_vol$threshold[is.na(f_logfc_vol$threshold)] <- "not_sig"

ggplot(f_logfc_vol, aes(x = log2FoldChange, y = -log10(pvalue))) + 
  geom_point(aes(colour = threshold), size=1) +
  ggtitle("Frontal Cortex \nCONTROL vs STAT + STAT-BIRTH") +
  xlab("Shrunken LogFC")+
  theme(plot.title = element_text(hjust = 0.5)) +
  #ylim(0,20) + xlim(-2,2) + 
  scale_colour_manual(values=c("blue", "grey", "red3"),
                      name="Threshold 0.01", 
                      breaks=c("down_sig", "not_sig", "up_sig"), 
                      labels=c("Down", "None", "Up"))+
  theme_bw()+
  theme(axis.text=element_text(size=14,colour="black"),
        axis.title=element_text(size=14,colour="black"),
        text = element_text(size=14,colour="black"))

###Heatmaps

###Heatmap of the genes with padj more than 0.01 and log fold change more than 1.5
alpha <- 0.01
a_top_genes <- rownames(a_results)[a_results$padj <= alpha & abs(a_results$log2FoldChange)>=1.5 &!is.na(a_results$padj)]
a_file_top <- a_file[row.names(a_file)%in%c(a_top_genes),]
a_file_log <- a_file_top
a_file_log[a_file_log<0] = -log2(abs(a_file_log[a_file_log<0]))
a_file_log[a_file_log>0] = log2(a_file_log[a_file_log>0])
a_file_log <- as.data.frame(t(scale(t(a_file_log))))
is.na(a_file_log)<-sapply(a_file_log, is.infinite)
a_file_log[is.na(a_file_log)]<-0 
a_conditions <- factor(a_conditions, levels=c("cont","antibiotics"))
a_top_anno <- HeatmapAnnotation(Condition = a_conditions,
                               col = list(Condition = c("cont" = "forestgreen", 
                                                        "antibiotics" = "magenta4")))
a_save_genes <- a_results[a_results$padj <= alpha & abs(a_results$log2FoldChange)>=1.5 &!is.na(a_results$padj),]
write.table(a_save_genes,file="A_genes_p_001_and_logfc.txt",sep="\t", col.names=T,
            row.names = T, quote=F)
col_fun = colorRamp2(c(-1, 0, 3), c("blue", "white", "red"))
Heatmap(a_file_log,column_order = order(a_conditions),col=col_fun, 
        top_annotation = a_top_anno, show_column_names = F,cluster_columns = F,
        heatmap_legend_param = list(title = "Z-score"),
        column_title="Amygdala \np<0.01 & logFC>1.5")

f_top_genes <- rownames(f_results)[f_results$padj <= alpha & abs(f_results$log2FoldChange)>=1.5 &!is.na(f_results$padj)]
f_file_top <- f_file[row.names(f_file)%in%c(f_top_genes),]
f_file_log <- f_file_top
f_file_log[f_file_log<0] = -log2(abs(f_file_log[f_file_log<0]))
f_file_log[f_file_log>0] = log2(f_file_log[f_file_log>0])
f_file_log <- as.data.frame(t(scale(t(f_file_log))))
is.na(f_file_log)<-sapply(f_file_log, is.infinite)
f_file_log[is.na(f_file_log)]<-0 
f_conditions <- factor(f_conditions, levels=c("cont","antibiotics"))
f_top_anno <- HeatmapAnnotation(Condition = f_conditions,
                                col = list(Condition = c("cont" = "forestgreen", 
                                                         "antibiotics" = "magenta4")))
f_save_genes <- f_results[f_results$padj <= alpha & abs(f_results$log2FoldChange)>=1.5 &!is.na(f_results$padj),]
write.table(f_save_genes,file="F_genes_p_001_and_logfc.txt",sep="\t", col.names=T,
            row.names = T, quote=F)
col_fun = colorRamp2(c(-1, 0, 3), c("blue", "white", "red"))
Heatmap(f_file_log,column_order = order(f_conditions),col=col_fun, 
        top_annotation = f_top_anno, show_column_names = F,cluster_columns = F,
        heatmap_legend_param = list(title = "Z-score"),row_names_gp = grid::gpar(fontsize = 8),
        column_title="Frontal Cortex \np<0.01 & logFC>1.5")

###All genes with p value < 0.01
a_top_genes_all <- a_results[a_results$padj <= alpha &!is.na(a_results$padj),]
write.table(a_top_genes_all,file="A_all_genes_threshold_001.txt",sep="\t", col.names=T,
            row.names = T, quote=F)
f_top_genes_all <- f_results[f_results$padj <= alpha & !is.na(f_results$padj),]
write.table(f_top_genes_all,file="F_all_genes_threshold_001.txt",sep="\t", col.names=T,
            row.names = T, quote=F)

###PCA Plots

###Visualize screeplot 
log_a_file <- scale(a_file)
log_a_file[log_a_file<0] = -log2(abs(log_a_file[log_a_file<0]))
log_a_file[log_a_file>0] = log2(log_a_file[log_a_file>0])
a_pca <-prcomp(t(log_a_file))
screeplot(a_pca, type='lines')

###Create PCA object to be plotted with ggplot
a_plot_pca <- as.data.frame(a_pca$x)[,c(1,2)]
a_plot_pca$Condition <- c(as.character(a_conditions))

###Plot PCA
ggplot(a_plot_pca, aes(x = PC1, y = PC2, color = Condition)) +
  geom_point(size = 2) +
  xlab("PC1 (8.6%)") + ylab("PC2 (6.5%)") + 
  scale_color_manual(values = c("magenta4","forestgreen")) + 
  ggtitle("Amygdala \nPCA on All Genes")+
  theme_bw()+
  theme(axis.text=element_text(size=14,colour="black"),
        axis.title=element_text(size=14,colour="black"),
        text = element_text(size=14,colour="black"))
write.table(a_plot_pca,"pca_amygdala_antibiotics.txt",
            sep="\t",quote=F,row.names = F,col.names = T)
log_f_file <- scale(f_file)
log_f_file[log_f_file<0] = -log2(abs(log_f_file[log_f_file<0]))
log_f_file[log_f_file>0] = log2(log_f_file[log_f_file>0])
f_pca <-prcomp(t(log_f_file))
screeplot(f_pca, type='lines')

###Create PCA object to be plotted with ggplot
f_plot_pca <- as.data.frame(f_pca$x)[,c(1,2)]
f_plot_pca$Condition <- c(as.character(f_conditions))
ggplot(f_plot_pca, aes(x = PC1, y = PC2, color = Condition)) +
  geom_point(size = 2) +
  xlab("PC1 (15.7%)") + ylab("PC2 (5.4%)") + 
  scale_color_manual(values = c("magenta4","forestgreen")) + 
  ggtitle("Frontal Cortex \nPCA on All Genes")+
  theme_bw()+
  theme(axis.text=element_text(size=14,colour="black"),
        axis.title=element_text(size=14,colour="black"),
        text = element_text(size=14,colour="black"))
write.table(f_plot_pca,"pca_frontal_antibiotics.txt",
            sep="\t",quote=F,row.names = F,col.names = T)
