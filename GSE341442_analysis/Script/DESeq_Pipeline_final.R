#installing packages needed
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::version()
BiocManager::install("DESeq2")
BiocManager::install("S4Arrays") #required with DESeq2
BiocManager::install("pheatmap")
BiocManager::install("clusterProfiler") 

#Reading the data

library(readr)
library(dplyr)
library(tidyr)
counts <- read.csv("GSE341442_Counts_allSamples.csv", row.names = 1)
exp <- counts[rowSums(counts) != 0, ]

#Normalisation
library(DESeq2)
library(S4Vectors) 
library(S4Arrays)
library(SummarizedExperiment) 

sample_info <- data.frame(row.names = colnames(exp), condition = c("2% glucose/1 mM K⁺","2% glucose/1 mM K⁺","2% glucose/1 mM K⁺","2% glucose/7 mM K⁺","2% glucose/7 mM K⁺", "2% glucose/7 mM K⁺","0.5% glucose/1 mM K⁺","0.5% glucose/1 mM K⁺","0.5% glucose/1 mM K⁺","0.5% glucose/7 mM K⁺" ,"0.5% glucose/7 mM K⁺","0.5% glucose/7 mM K⁺"))

rep_count <- table(sample_info$condition)
print(rep_count)

dds <- DESeqDataSetFromMatrix(countData = exp, colData = sample_info, design = ~ condition)
dds <- estimateSizeFactors(dds)
keep <- rowMeans(counts(dds)) >= 10
dds <- dds[keep,]
normalised_counts <- counts(dds, normalized = TRUE)
norm <- as.data.frame(normalised_counts)

#PCA visualization
library(ggplot2)
vsd <- vst(dds)
pca_data <- prcomp(t(assay(vsd)))
pca_df <- as.data.frame(pca_data$x) 
pca_df$condition <- sample_info$condition
pca <- ggplot(pca_df,aes(x = PC1, y = PC2, color = condition)) + geom_point(size=3) + theme_minimal() + ggtitle("PCA of samples")
print(pca)
#Heatmap

library(pheatmap)
normalised_counts_log <- log2(normalised_counts + 1)

# Use top 20 most variable genes to see if samples cluster together

pheatmap(normalised_counts_log[1:20, ], annotation_col = sample_info , cluster_cols = TRUE, main = "Do replicates look similar?", fontsize= 6, fontsize_row = 4, cellheight = 8, width = 8, height= 8)


library(EnhancedVolcano)
library(clusterProfiler)
library(org.Sc.sgd.db)
library(AnnotationDbi)
library(stringr)

#performing differential gene expression
dds <- DESeq(dds)

#Differential gene expression results for six conditions, top 10 genes, volcano plot and Gene Ontology

#six different conditions where differential expression is seen
contrasts_list <- list(
  "Condition 1" = c("condition", "2% glucose/7 mM K⁺", "2% glucose/1 mM K⁺"),
  "Condition 2" = c("condition", "2% glucose/7 mM K⁺", "0.5% glucose/1 mM K⁺"),
  "Condition 3" = c("condition", "2% glucose/7 mM K⁺", "0.5% glucose/7 mM K⁺"),
  "Condition 4" = c("condition", "2% glucose/1 mM K⁺", "0.5% glucose/1 mM K⁺"),
  "Condition 5" = c("condition", "2% glucose/1 mM K⁺", "0.5% glucose/7 mM K⁺"),
  "Condition 6" = c("condition", "0.5% glucose/7 mM K⁺", "0.5% glucose/1 mM K⁺")
)

all_results <- list()
all_sig <- list()
all_ego <- list()

dir.create("DEG_results", showWarnings = FALSE) # folder for CSVs
dir.create("GO_results", showWarnings = FALSE) # folder for plots

for(nm in names(contrasts_list)){
  cat("\n=== Running:", nm, "===\n")
  
  # 1. DESeq2 results
  res <- results(dds, contrast = contrasts_list[[nm]])
  res$sgd <- rownames(res)
  res$systematic <- mapIds(org.Sc.sgd.db, keys = rownames(res),
                           column = "SGD", keytype = "SGD", multiVals = "first")
  res$genename <- mapIds(org.Sc.sgd.db, keys = rownames(res),
                         column = "GENENAME", keytype = "SGD", multiVals = "first")
  
  # 2. Filter sig genes
  sig <- res[!is.na(res$padj) & res$padj < 0.05 & abs(res$log2FoldChange) > 1, ]
  
  # 3. Volcano plot
  p <- EnhancedVolcano(res,
                       lab = res$genename,
                       x = 'log2FoldChange',
                       y = 'padj',
                       title = nm,
                       pCutoff = 0.05,
                       FCcutoff = 1.0)
  print(p)
  ggsave(paste0("DEG_results/volcano_", nm, ".png"), plot = p, width = 8, height = 8)
  
  # 4. Save CSV
  write.csv(as.data.frame(sig), paste0("DEG_results/DEG_", nm, ".csv"), row.names = TRUE)
  
  # 5. GO enrichment - only if we have sig genes
  if(nrow(sig) > 0){
    ego <- enrichGO(gene = sig$systematic, # use systematic, not rownames
                    OrgDb = org.Sc.sgd.db,
                    keyType = "SGD",
                    ont = "BP",
                    pAdjustMethod = "BH",
                    qvalueCutoff = 0.05)
    
    all_ego[[nm]] <- ego
    print(dotplot(ego, showCategory = 20) + ggtitle(nm))
    ggsave(paste0("GO_results/dotplot_", nm, ".png"),
           dotplot(ego, showCategory = 20) + ggtitle(nm), width = 8, height = 10)
  } else {
    cat("No significant genes for", nm, "\n")
  }
  
  all_results[[nm]] <- res
  all_sig[[nm]] <- sig
}


