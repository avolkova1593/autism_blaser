library(DESeq2)
library(ggplot2)
library(dplyr)
library(ComplexHeatmap)
library(circlize)

set.seed(1234)
setwd ("/Users/av1936/Desktop/Google_Drive/collaborative_projects/blaser_alzheimers")

file <- read.table("feature_counts.txt", sep="\t", header = T,row.names=1)
samples_fc <- read.csv("RNA_Seq_analysis/FC_samples/samplesheet_2.csv")
samples_hip <- read.csv("RNA_Seq_analysis/A_samples/samplesheet_2.csv")
####Get the names of the samples from the samplesheet files
a_control_ids <- sub(" ","_",as.character(samples_hip[samples_hip$condition%in%c("CONTROL"),]$id))
a_stat_ids <- sub(" ","_",as.character(samples_hip[samples_hip$condition%in%c("STAT"),]$id))
a_stat_birth_ids <- sub(" ","_",as.character(samples_hip[samples_hip$condition%in%c("STAT-BIRTH"),]$id))
f_control_ids <- sub(" ","_",as.character(samples_fc[samples_fc$condition%in%c("CONTROL"),]$id))
f_stat_ids <- sub(" ","_",as.character(samples_fc[samples_fc$condition%in%c("STAT"),]$id))
f_stat_birth_ids <- sub(" ","_",as.character(samples_fc[samples_fc$condition%in%c("STAT-BIRTH"),]$id))
a_control_ids <- sub(" ","_",a_control_ids)
a_stat_ids <- sub(" ","_",a_stat_ids)
f_control_ids <- sub(" ","_",f_control_ids)
f_stat_ids <- sub(" ","_",f_stat_ids)
a_stat <- c(a_control_ids,a_stat_ids)
a_stat_birth <- c(a_control_ids,a_stat_birth_ids)
a_all <- c(a_control_ids,a_stat_ids,a_stat_birth_ids)
f_stat <- c(f_control_ids,f_stat_ids)
f_stat_birth <- c(f_control_ids,f_stat_birth_ids)
f_all <- c(f_control_ids,f_stat_ids,f_stat_birth_ids)

###Remove numbers after the dots in the feature file
names(file) <- sub("\\..*","",names(file))
###Split the feature file into amygdala and fc
a_file <- file[,names(file)%in%c(a_all),]
a_file_stat <- file[,names(file)%in%c(a_stat),]
a_file_stat_birth <- file[,names(file)%in%c(a_stat_birth),]
a_samples_stat <- c(names(a_file_stat))
a_samples_stat_birth <- c(names(a_file_stat_birth))
names(a_file) <- sub("\\..*","",names(a_file))
names(a_file_stat) <- sub("\\..*","",names(a_file_stat))
names(a_file_stat_birth) <- sub("\\..*","",names(a_file_stat_birth))

f_file <- file[,names(file)%in%c(f_all),]
f_file_stat <- file[,names(file)%in%c(f_stat),]
f_file_stat_birth <- file[,names(file)%in%c(f_stat_birth),]
f_samples_stat <- c(names(f_file_stat))
f_samples_stat_birth <- c(names(f_file_stat_birth))
names(f_file) <- sub("\\..*","",names(f_file))
names(f_file_stat) <- sub("\\..*","",names(f_file_stat))
names(f_file_stat_birth) <- sub("\\..*","",names(f_file_stat_birth))

###change the names of both data frames and get corresponding conditions
names(a_file)[names(a_file)%in%c(a_control_ids)]<- "cont"
names(a_file)[names(a_file)%in%c(a_stat_ids)]<- "stat"
names(a_file)[names(a_file)%in%c(a_stat_birth_ids)]<- "stat_birth"

names(a_file_stat)[names(a_file_stat)%in%c(a_control_ids)]<- "cont"
names(a_file_stat)[names(a_file_stat)%in%c(a_stat_ids)]<- "stat"
names(a_file_stat_birth)[names(a_file_stat_birth)%in%c(a_control_ids)]<- "cont"
names(a_file_stat_birth)[names(a_file_stat_birth)%in%c(a_stat_birth_ids)]<- "stat_birth"

names(f_file)[names(f_file)%in%c(f_control_ids)]<- "cont"
names(f_file)[names(f_file)%in%c(f_stat_ids)]<- "stat"
names(f_file)[names(f_file)%in%c(f_stat_birth_ids)]<- "stat_birth"

names(f_file_stat)[names(f_file_stat)%in%c(f_control_ids)]<- "cont"
names(f_file_stat)[names(f_file_stat)%in%c(f_stat_ids)]<- "stat"
names(f_file_stat_birth)[names(f_file_stat_birth)%in%c(f_control_ids)]<- "cont"
names(f_file_stat_birth)[names(f_file_stat_birth)%in%c(f_stat_birth_ids)]<- "stat_birth"

a_conditions <- c(names(a_file))
a_conditions_stat <- c(names(a_file_stat))
a_conditions_stat_birth <- c(names(a_file_stat_birth))
f_conditions <- c(names(f_file))
f_conditions_stat <- c(names(f_file_stat))
f_conditions_stat_birth <- c(names(f_file_stat_birth))

names(a_file_stat) <- a_samples_stat
names(a_file_stat_birth) <- a_samples_stat_birth

names(f_file_stat) <- f_samples_stat
names(f_file_stat_birth) <- f_samples_stat_birth

###set samples data frames
a_df_samples_stat <- data.frame(row.names = c(names(a_file_stat)), 
                      vector = as.factor(a_conditions_stat))
a_df_samples_stat_birth <- data.frame(row.names = c(names(a_file_stat_birth)), 
                                vector = as.factor(a_conditions_stat_birth))
f_df_samples_stat <- data.frame(row.names = c(names(f_file_stat)), 
                           vector = as.factor(f_conditions_stat))
f_df_samples_stat_birth <- data.frame(row.names = c(names(f_file_stat_birth)), 
                                vector = as.factor(f_conditions_stat_birth))
###Make data frames deseq ready
a_df_deseq_stat <- DESeqDataSetFromMatrix(countData = a_file_stat, colData = a_df_samples_stat, design = ~ vector)
a_df_deseq_stat_birth <- DESeqDataSetFromMatrix(countData = a_file_stat_birth, colData = a_df_samples_stat_birth, design = ~ vector)

f_df_deseq_stat <- DESeqDataSetFromMatrix(countData = f_file_stat, colData = f_df_samples_stat, design = ~ vector)
f_df_deseq_stat_birth <- DESeqDataSetFromMatrix(countData = f_file_stat_birth, colData = f_df_samples_stat_birth, design = ~ vector)

###Remove genes with counts less than 2
a_df_deseq_stat <- a_df_deseq_stat[rowSums(counts(a_df_deseq_stat) >= 2) >= 2]
a_df_deseq_stat_birth <- a_df_deseq_stat_birth[rowSums(counts(a_df_deseq_stat_birth) >= 2) >= 2]

f_df_deseq_stat <- f_df_deseq_stat[rowSums(counts(f_df_deseq_stat) >= 2) >= 2]
f_df_deseq_stat_birth <- f_df_deseq_stat_birth[rowSums(counts(f_df_deseq_stat_birth) >= 2) >= 2]

###Run Deseq
a_processed_stat <- DESeq(a_df_deseq_stat)
a_processed_stat_birth <- DESeq(a_df_deseq_stat_birth)

f_processed_stat <- DESeq(f_df_deseq_stat)
f_processed_stat_birth <- DESeq(f_df_deseq_stat_birth)

###Get deseq results
a_results_stat <- results(a_processed_stat, contrast = c("vector", "cont","stat"))
a_results_stat_birth <- results(a_processed_stat_birth, contrast = c("vector", "cont","stat_birth"))

f_results_stat <- results(f_processed_stat, contrast = c("vector", "cont", "stat"))
f_results_stat_birth <- results(f_processed_stat_birth, contrast = c("vector", "cont", "stat_birth"))

###Select the genes with logfc > 1.5
a_results_logfc_stat <- a_results_stat[abs(a_results_stat$log2FoldChange)>=1.5,]
a_results_logfc_stat_birth <- a_results_stat_birth[abs(a_results_stat_birth$log2FoldChange)>=1.5,]

f_results_logfc_stat <- f_results_stat[abs(f_results_stat$log2FoldChange)>=1.5,]
f_results_logfc_stat_birth <- f_results_stat_birth[abs(f_results_stat_birth$log2FoldChange)>=1.5,]

###Get only significant results
a_results_no_na_stat <- a_results_logfc_stat[!is.na(a_results_logfc_stat$padj),]
a_results_no_na_stat_birth <- a_results_logfc_stat_birth[!is.na(a_results_logfc_stat_birth$padj),]

f_results_no_na_stat <- f_results_logfc_stat[!is.na(f_results_logfc_stat$padj),]
f_results_no_na_stat_birth <- f_results_logfc_stat_birth[!is.na(f_results_logfc_stat_birth$padj),]

a_results_sig_stat <- a_results_stat[a_results_no_na_stat$padj<0.05,]
a_results_sig_stat_birth <- a_results_stat_birth[a_results_no_na_stat_birth$padj<0.05,]

f_results_sig_stat <- f_results_stat[f_results_no_na_stat$padj<0.05,]
f_results_sig_stat_birth <- f_results_stat_birth[f_results_no_na_stat_birth$padj<0.05,]

###calculate row means across vectors
a_mean_cont_stat <- rowMeans(counts(a_processed_stat, normalized=TRUE)[, a_processed_stat$vector == "cont"])
a_mean_anti_stat <- rowMeans(counts(a_processed_stat, normalized=TRUE)[, a_processed_stat$vector == "stat"])
a_mean_cont_stat_birth <- rowMeans(counts(a_processed_stat_birth, normalized=TRUE)[, a_processed_stat_birth$vector == "cont"])
a_mean_anti_stat_birth <- rowMeans(counts(a_processed_stat_birth, normalized=TRUE)[, a_processed_stat_birth$vector == "stat_birth"])

f_mean_cont_stat <- rowMeans(counts(f_processed_stat, normalized=TRUE)[, f_processed_stat$vector == "cont"])
f_mean_anti_stat <- rowMeans(counts(f_processed_stat, normalized=TRUE)[, f_processed_stat$vector == "stat"])
f_mean_cont_stat_birth <- rowMeans(counts(f_processed_stat_birth, normalized=TRUE)[, f_processed_stat_birth$vector == "cont"])
f_mean_anti_stat_birth <- rowMeans(counts(f_processed_stat_birth, normalized=TRUE)[, f_processed_stat_birth$vector == "stat_birth"])
###Plots

### Log Fold change
a_logfc_stat <- lfcShrink(a_processed_stat, contrast = c("vector", "cont", "stat")) 
a_logfc_stat <- cbind(as.data.frame(a_logfc_stat), a_mean_cont_stat, a_mean_anti_stat)
a_logfc_stat_birth <- lfcShrink(a_processed_stat_birth, contrast = c("vector", "cont", "stat_birth")) 
a_logfc_stat_birth <- cbind(as.data.frame(a_logfc_stat_birth), a_mean_cont_stat_birth, a_mean_anti_stat_birth)

f_logfc_stat <- lfcShrink(f_processed_stat, contrast = c("vector", "cont", "stat")) 
f_logfc_stat <- cbind(as.data.frame(f_logfc_stat), f_mean_cont_stat, f_mean_anti_stat)
f_logfc_stat_birth <- lfcShrink(f_processed_stat_birth, contrast = c("vector", "cont", "stat_birth")) 
f_logfc_stat_birth <- cbind(as.data.frame(f_logfc_stat_birth), f_mean_cont_stat_birth, f_mean_anti_stat_birth)
###Volcano plots
a_logfc_vol_stat <- a_logfc_stat%>%mutate(threshold = ifelse((log2FoldChange >= 0 & padj < 0.01), "up_sig",
                                                   ifelse((log2FoldChange <= 0 & padj < 0.01),
                                                   "down_sig", "not_sig"))) 
a_logfc_vol_stat$threshold[is.na(a_logfc_vol_stat$threshold)] <- "not_sig"

a_logfc_vol_stat_birth <- a_logfc_stat_birth%>%mutate(threshold = ifelse((log2FoldChange >= 0 & padj < 0.01), "up_sig",
                                                             ifelse((log2FoldChange <= 0 & padj < 0.01),
                                                                    "down_sig", "not_sig"))) 
a_logfc_vol_stat_birth$threshold[is.na(a_logfc_vol_stat_birth$threshold)] <- "not_sig"

ggplot(a_logfc_vol_stat, aes(x = log2FoldChange, y = -log10(pvalue))) + 
  geom_point(aes(colour = threshold), size=1) +
  ggtitle("Amygdala \nCONTROL vs STAT") +
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

ggplot(a_logfc_vol_stat_birth, aes(x = log2FoldChange, y = -log10(pvalue))) + 
  geom_point(aes(colour = threshold), size=1) +
  ggtitle("Amygdala \nCONTROL vs STAT-BIRTH") +
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

f_logfc_vol_stat <- f_logfc_stat%>%mutate(threshold = ifelse((log2FoldChange >= 0 & padj < 0.01), "up_sig",
                                                   ifelse((log2FoldChange <= 0 & padj < 0.01),
                                                          "down_sig", "not_sig"))) 
f_logfc_vol_stat$threshold[is.na(f_logfc_vol_stat$threshold)] <- "not_sig"

f_logfc_vol_stat_birth <- f_logfc_stat_birth%>%mutate(threshold = ifelse((log2FoldChange >= 0 & padj < 0.01), "up_sig",
                                                             ifelse((log2FoldChange <= 0 & padj < 0.01),
                                                                    "down_sig", "not_sig"))) 
f_logfc_vol_stat_birth$threshold[is.na(f_logfc_vol_stat_birth$threshold)] <- "not_sig"
ggplot(f_logfc_vol_stat, aes(x = log2FoldChange, y = -log10(pvalue))) + 
  geom_point(aes(colour = threshold), size=1) +
  ggtitle("Frontal Cortex \nCONTROL vs STAT") +
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

ggplot(f_logfc_vol_stat_birth, aes(x = log2FoldChange, y = -log10(pvalue))) + 
  geom_point(aes(colour = threshold), size=1) +
  ggtitle("Frontal Cortex \nCONTROL vs STAT-BIRTH") +
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
a_top_genes_stat <- rownames(a_results_stat)[a_results_stat$padj <= alpha & abs(a_results_stat$log2FoldChange)>=1.5 &!is.na(a_results_stat$padj)]
a_top_genes_stat_birth <- rownames(a_results_stat_birth)[a_results_stat_birth$padj <= alpha & abs(a_results_stat_birth$log2FoldChange)>=1.5 &!is.na(a_results_stat_birth$padj)]
a_top_genes <- unique(c(a_top_genes_stat,a_top_genes_stat_birth))
a_file_top <- a_file[row.names(a_file)%in%c(a_top_genes),]
#write.table(as.data.frame(a_top_genes),"A_top_genes_stat_and_stat_birth.txt", row.names = F,col.names = F,
 #           quote=F, sep="\t")
a_file_log <- a_file_top
a_file_log[a_file_log<0] = -log2(abs(a_file_log[a_file_log<0]))
a_file_log[a_file_log>0] = log2(a_file_log[a_file_log>0])
a_file_log <- as.data.frame(t(scale(t(a_file_log))))
is.na(a_file_log)<-sapply(a_file_log, is.infinite)
a_file_log[is.na(a_file_log)]<-0 
a_conditions <- factor(a_conditions, levels=c( "cont","stat","stat_birth"))
a_top_anno <- HeatmapAnnotation(Condition = a_conditions,
                               col = list(Condition = c("cont" = "forestgreen", 
                                                        "stat" = "deeppink3",
                                                        "stat_birth"="blue")))
a_save_genes_stat <- a_results_stat[a_results_stat$padj <= alpha & abs(a_results_stat$log2FoldChange)>=1.5 &!is.na(a_results_stat$padj),]
a_save_genes_stat_birth <- a_results_stat_birth[a_results_stat_birth$padj <= alpha & abs(a_results_stat_birth$log2FoldChange)>=1.5 &!is.na(a_results_stat_birth$padj),]
#write.table(a_save_genes_stat,file="A_genes_stat_p_001_and_logfc.txt",sep="\t", col.names=T,
 #           row.names = T, quote=F)
#write.table(a_save_genes_stat_birth,file="A_genes_stat_birth_p_001_and_logfc.txt",sep="\t", col.names=T,
 #           row.names = T, quote=F)
col_fun = colorRamp2(c(-1, 0, 3), c("blue", "white", "red"))

Heatmap(a_file_log,column_order = order(a_conditions),col=col_fun, 
        top_annotation = a_top_anno, show_column_names = F,cluster_columns = F,
        heatmap_legend_param = list(title = "Z-score"),
        column_title="Amygdala \np<0.01 & logFC>1.5")

f_top_genes_stat <- rownames(f_results_stat)[f_results_stat$padj <= alpha & abs(f_results_stat$log2FoldChange)>=1.5 &!is.na(f_results_stat$padj)]
f_top_genes_stat_birth <- rownames(f_results_stat_birth)[f_results_stat_birth$padj <= alpha & abs(f_results_stat_birth$log2FoldChange)>=1.5 &!is.na(f_results_stat_birth$padj)]
f_top_genes <- unique(c(f_top_genes_stat,f_top_genes_stat_birth))
f_file_top <- f_file[row.names(f_file)%in%c(f_top_genes),]
#write.table(as.data.frame(a_top_genes),"A_top_genes_stat_and_stat_birth.txt", row.names = F,col.names = F,
#           quote=F, sep="\t")
f_file_log <- f_file_top
f_file_log[f_file_log<0] = -log2(abs(f_file_log[f_file_log<0]))
f_file_log[f_file_log>0] = log2(f_file_log[f_file_log>0])
f_file_log <- as.data.frame(t(scale(t(f_file_log))))
is.na(f_file_log)<-sapply(f_file_log, is.infinite)
f_file_log[is.na(f_file_log)]<-0 
f_conditions <- factor(f_conditions, levels=c( "cont","stat","stat_birth"))
f_top_anno <- HeatmapAnnotation(Condition = f_conditions,
                                col = list(Condition = c("cont" = "forestgreen", 
                                                         "stat" = "deeppink3",
                                                         "stat_birth"="blue")))
f_save_genes_stat <- f_results_stat[f_results_stat$padj <= alpha & abs(f_results_stat$log2FoldChange)>=1.5 &!is.na(f_results_stat$padj),]
f_save_genes_stat_birth <- f_results_stat_birth[f_results_stat_birth$padj <= alpha & abs(f_results_stat_birth$log2FoldChange)>=1.5 &!is.na(f_results_stat_birth$padj),]
#write.table(a_save_genes_stat,file="A_genes_stat_p_001_and_logfc.txt",sep="\t", col.names=T,
#           row.names = T, quote=F)
#write.table(a_save_genes_stat_birth,file="A_genes_stat_birth_p_001_and_logfc.txt",sep="\t", col.names=T,
#           row.names = T, quote=F)
col_fun = colorRamp2(c(-1, 0, 3), c("blue", "white", "red"))

Heatmap(f_file_log,column_order = order(f_conditions),col=col_fun, 
        top_annotation = a_top_anno, show_column_names = F,cluster_columns = F,
        heatmap_legend_param = list(title = "Z-score"),
        column_title="Frontal Cortex \np<0.01 & logFC>1.5")









###Heatmap of all genes with p value < 0.01
a_top_genes_all_stat <- a_results_stat[a_results_stat$padj <= alpha &!is.na(a_results_stat$padj),]
write.table(a_top_genes_all_stat,file="A_stat_genes_threshold_001.txt",sep="\t", col.names=T,
            row.names = T, quote=F)
a_top_genes_all_stat_birth <- a_results_stat_birth[a_results_stat_birth$padj <= alpha &!is.na(a_results_stat_birth$padj),]
write.table(a_top_genes_all_stat_birth,file="A_stat_birth_genes_threshold_001.txt",sep="\t", col.names=T,
            row.names = T, quote=F)
# a_file_top_all <- a_file[row.names(a_file)%in%c(a_top_genes_all),]
# a_file_log_all <- log2(a_file_top_all)
# is.na(a_file_log_all)<-sapply(a_file_log_all, is.infinite)
# a_file_log_all[is.na(a_file_log_all)]<-0 
# a_top_anno_all <- HeatmapAnnotation(Condition = a_conditions,
#                                 col = list(Condition = c("cont" = "forestgreen", "antibiotics" = "goldenrod2")))
# Heatmap(scale(a_file_log_all),
#         top_annotation = a_top_anno_all, show_column_names = F,
#         show_row_names = F,
#         heatmap_legend_param = list(title = "Z score"),
#         column_title="Amygdala \n All Genes with p<0.01")
# 
f_top_genes_all_stat <- f_results_stat[f_results_stat$padj <= alpha &!is.na(f_results_stat$padj),]
write.table(f_top_genes_all_stat,file="FC_stat_genes_threshold_001.txt",sep="\t", col.names=T,
            row.names = T, quote=F)
f_top_genes_all_stat_birth <- f_results_stat_birth[f_results_stat_birth$padj <= alpha &!is.na(f_results_stat_birth$padj),]
write.table(f_top_genes_all_stat_birth,file="FC_stat_birth_genes_threshold_001.txt",sep="\t", col.names=T,
            row.names = T, quote=F)
# f_file_top_all <- f_file[row.names(f_file)%in%c(f_top_genes_all),]
# f_file_log_all <- log2(f_file_top_all)
# is.na(f_file_log_all)<-sapply(f_file_log_all, is.infinite)
# f_file_log_all[is.na(f_file_log_all)]<-0 
# f_top_anno_all <- HeatmapAnnotation(Condition = f_conditions,
#                                 col = list(Condition = c("cont" = "forestgreen", "antibiotics" = "goldenrod2")))
# Heatmap(scale(f_file_log_all),
#         top_annotation = a_top_anno, show_column_names = F,
#         heatmap_legend_param = list(title = "Z score"),
#         column_title="Frontal Cortex \n All Genes with p<0.01")

###PCA Plots
### Compare with PCA
# visualize screeplot 
log_a_file <- scale(a_file)
log_a_file[log_a_file<0] = -log2(abs(log_a_file[log_a_file<0]))
log_a_file[log_a_file>0] = log2(log_a_file[log_a_file>0])
a_pca <-prcomp(t(log_a_file))
screeplot(a_pca, type='lines')
# Create PCA object to be plotted with ggplot
a_plot_pca <- as.data.frame(a_pca$x)[,c(1,2)]
a_plot_pca$Condition <- a_conditions
# Plot PCA
ggplot(a_plot_pca, aes(x = PC1, y = PC2, color = Condition)) +
  geom_point(size = 2) +
  xlab("PC1 (8.6%)") + ylab("PC2 (6.5%)") + 
  scale_color_manual(values = c("forestgreen","deeppink3","blue")) + 
  ggtitle("Amygdala \nPCA on All Genes")+
  theme_bw()+
  theme(axis.text=element_text(size=14,colour="black"),
        axis.title=element_text(size=14,colour="black"),
        text = element_text(size=14,colour="black"))

log_f_file <- scale(f_file)
log_f_file[log_f_file<0] = -log2(abs(log_f_file[log_f_file<0]))
log_f_file[log_f_file>0] = log2(log_f_file[log_f_file>0])
f_pca <-prcomp(t(log_f_file))
screeplot(f_pca, type='lines')
# Create PCA object to be plotted with ggplot
f_plot_pca <- as.data.frame(f_pca$x)[,c(1,2)]
f_plot_pca$Condition <- f_conditions
# Plot PCA
ggplot(f_plot_pca, aes(x = PC1, y = PC2, color = Condition)) +
  geom_point(size = 2) +
  xlab("PC1 (15.7%)") + ylab("PC2 (5.4%)") + 
  scale_color_manual(values = c("forestgreen","deeppink3","blue")) + 
  ggtitle("Frontal Cortex \nPCA on All Genes")+
  theme_bw()+
  theme(axis.text=element_text(size=14,colour="black"),
        axis.title=element_text(size=14,colour="black"),
        text = element_text(size=14,colour="black"))









