library(TCGAbiolinks)
library(dplyr)
library(DT)
library(SummarizedExperiment)


# # =========================================================================
# # clinical
# project_name <- 'TCGA-SKCM' 
# # project_name <- 'TCGA-UVM' 
# 
# clinical_data <- GDCquery_clinic(project = project_name, type = "clinical")
# 
# 
# # gene expression
# Expr_df <- GDCquery(project = project_name,
#                     data.category = "Transcriptome Profiling",
#                     data.type = "Gene Expression Quantification",
#                     workflow.type = "HTSeq - FPKM")
# GDCdownload(Expr_df, method = "api", files.per.chunk = 100)
# expdat <- GDCprepare(query = Expr_df)
# Expr_matrix <- assay(expdat)
# 
# 
# # gene id transformation, only keep the genes that can be transformed into geneid
# library(clusterProfiler)
# gene<-bitr(rownames(Expr_matrix),"ENSEMBL","SYMBOL","org.Hs.eg.db")
# Expr_matrix<-cbind(rownames(Expr_matrix),Expr_matrix)
# colnames(Expr_matrix)[1]<-"ENSEMBL"
# Gene_expr_df<-merge(gene,Expr_matrix,by="ENSEMBL")
# rm(Expr_matrix)
# 
# 
# # delete duplicate genes
# Gene_expr_df<-Gene_expr_df[-which(duplicated(Gene_expr_df[,2])),]
# rownames(Gene_expr_df) <- Gene_expr_df[,2]
# Gene_expr_df <- Gene_expr_df[,c(-1,-2)]
# 
# 
# 
# # the name of TCGA sample splited by '-'，
# # The first three columns are patient numbers，
# # The fourth column is the type: 11: normal; 01: tumor，
# # https://gdc.cancer.gov/resources-tcga-users/tcga-code-tables/sample-type-codes 
# # 06	Metastatic
# # 01	Primary Solid Tumor
# # 11	Solid Tissue Normal
# # 07	Additional Metastatic
# group <- strsplit(colnames(Gene_expr_df)[-1][-1],"[-]")
# class<-sapply(group,function(I){I[4]})
# class[which(grepl("11",class)==TRUE)]<-"Solid Tissue Normal"
# class[which(grepl("06",class)==TRUE)]<-"Metastatic"
# class[which(grepl("07",class)==TRUE)]<-"Additional Metastatic"
# class[which(grepl("01",class)==TRUE)]<-"Primary Solid Tumor"
# clinical_data <- cbind(clinical_data, class)
# # Metastatic: 364
# # Additional Metastatic:1
# # Primary Solid Tumor: 102
# # Solid Tissue Normal: 1
# 
# # merge survival data to expression data according to sample information
# newid<-lapply(strsplit(colnames(Gene_expr_df),'-'),
#               function(i){paste0(i[1:3],collapse = '-')})
# newid<-sapply(1:length(newid),function(i){newid[[i]]})
# 
# Gene_expr_df<-rbind(Gene_expr_df,newid)
# rownames(Gene_expr_df)[dim(Gene_expr_df)[1]]<-'newid'
# 
# colnames(clinical_data)[1]<-'newid'
# df_OS<-merge(clinical_data,t(Gene_expr_df),by="newid")
# 
# df_OS_dropdu<-df_OS[-which(duplicated(df_OS[,1])),]
# rownames(df_OS_dropdu)<-df_OS_dropdu[,1]
# df_OS_dropdu<-df_OS_dropdu[,-1]








#########################################################################
#########################################################################
# sample screen
df_tumor <- df_OS_dropdu[df_OS_dropdu[69]!='Solid Tissue Normal',]
# df_tumor <- df_OS_dropdu[df_OS_dropdu[69]=='Solid Tissue Normal',]
df_tumor <- cbind(rep(0, dim(df_tumor)[1]), df_tumor)
  
  
# Dividing the sample according to genes
gene <- 'XPNPEP3'
gexp <- as.numeric(df_tumor[, gene])
sample_div <- quantile(gexp, 0.5)
survival_var_choose <- gexp
survival_var_choose[gexp<=sample_div] <- 'low expr'
survival_var_choose[gexp>sample_div] <- 'high expr'
print(paste0('propotion:', 
             sum(survival_var_choose=='high expr')/length(gexp)))
colnames(df_tumor)[1] <- gene
df_tumor[,1] <- survival_var_choose


# status
df_tumor[df_tumor$vital_status=='Dead',]$vital_status <- 2
df_tumor[df_tumor$vital_status=='Alive',]$vital_status <- 1
df_tumor$vital_status <- as.numeric(df_tumor$vital_status)


# ==============================================================================
library("survminer")
require("survival")
# head(lung)

# tow group
fit<- survfit(Surv(days_to_death, vital_status) ~ XPNPEP3, data = df_tumor)
# save: 650, 510
ggsurvplot(fit,  size = 1,  
           linetype = "strata", 
           break.time.by = 2000, 
           risk.table.col = "strata",
           legend = c(0.8,0.75),
           legend.title = "", 
           palette = "Dark2",
           conf.int = TRUE, 
           pval = TRUE,
           pval.coord = c(7000,0.5),
           pval.size = 4,
           pval.method = TRUE,
           pval.method.size = 4,
           pval.method.coord = c(7000,0.6),
           risk.table = TRUE,
           surv.median.line = "hv",
           legend.labs = c("high expression", "low expression"),
           risk.table.y.text.col = TRUE,
           xlab = "Follow up time(d)",
           title = gene)
           # font.title = 'Helvetica',
           # font.subtitle = 'Helvetica',
           # font.caption = 'Helvetica',
           # font.x = 'Helvetica',
           # font.y = 'Helvetica',
           # font.tickslab = 'Helvetica',
           # font.legend = 'Helvetica')

