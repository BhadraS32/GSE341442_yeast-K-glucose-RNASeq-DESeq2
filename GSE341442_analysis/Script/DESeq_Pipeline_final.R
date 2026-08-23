#installing packages needed
#if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
#BiocManager::version()
#BiocManager::install("DESeq2")
#BiocManager::install("S4Arrays") #required with DESeq2
#BiocManager::install("pheatmap")
#BiocManager::install("clusterProfiler") 


library(readr) ;library(dplyr) ;library(tidyr); library(ggplot2); library(pheatmap)
library(DESeq2); library(EnhancedVolcano); library(clusterProfiler); library(org.Sc.sgd.db); 
library(AnnotationDbi); 
#1.Reading the data
counts <- read.csv("~/GSE341442_analysis/Data/GSE341442_Counts_allSamples.csv", row.names = 1)
#To remove genes from C.elegans (WB genes) that were contaminated in the sample
counts <- counts[!grepl("^WB", rownames(counts), ignore.case=TRUE), ]
#TO remove genes which have 0 expression in all samples
exp <- counts[rowSums(counts) != 0, ]

#2.Normalisation
#Metadata information
sample_info <- data.frame(row.names = colnames(exp), condition = c("2% glucose/1 mM K⁺","2% glucose/1 mM K⁺","2% glucose/1 mM K⁺","2% glucose/7 mM K⁺","2% glucose/7 mM K⁺", "2% glucose/7 mM K⁺","0.5% glucose/1 mM K⁺","0.5% glucose/1 mM K⁺","0.5% glucose/1 mM K⁺","0.5% glucose/7 mM K⁺" ,"0.5% glucose/7 mM K⁺","0.5% glucose/7 mM K⁺"))
rep_count <- table(sample_info$condition)
print(rep_count)

dds <- DESeqDataSetFromMatrix(countData = exp, colData = sample_info, design = ~ condition)
dds <- estimateSizeFactors(dds)
keep <- rowMeans(counts(dds)) >= 10
dds <- dds[keep,]

normalised_counts <- counts(dds, normalized = TRUE)
norm <- as.data.frame(normalised_counts)

#3PCA visualization
dir.create("~/GSE341442_analysis/results/PCA_and_Heatmap", recursive = TRUE, showWarnings = FALSE)
vsd <- vst(dds)
pca_data <- prcomp(t(assay(vsd)))
pca_df <- as.data.frame(pca_data$x) 
pca_df$condition <- sample_info$condition
pca <- ggplot(pca_df,aes(x = PC1, y = PC2, color = condition)) + geom_point(size=3) + theme_minimal() + ggtitle("PCA of samples")
print(pca)
ggsave(paste0("~/GSE341442_analysis/results/PCA_and_Heatmap/PCA", ".png"), plot = pca, width = 8, height = 8)

#Heatmap
normalised_counts_log <- log2(normalised_counts + 1)
# Use top 20 most variable genes to see if samples cluster together
heatmap <- pheatmap(normalised_counts_log[1:20, ], annotation_col = sample_info , cluster_cols = TRUE, main = "Do replicates look similar?", fontsize= 6, fontsize_row = 4, cellheight = 8, width = 8, height= 15)
print(heatmap)
ggsave(paste0("~/GSE341442_analysis/results/PCA_and_Heatmap/Heatmap", ".png"), plot = heatmap, width = 8, height = 15)

3.#performing differential gene expression
dds <- DESeq(dds)

#Differential gene expression results for six conditions, top 10 genes, volcano plot and Gene Ontology
#six different conditions where differential expression is seen
contrasts_list <- list(
  "Condition_1" = c("condition", "2% glucose/7 mM K⁺", "2% glucose/1 mM K⁺"),
  "Condition_2" = c("condition", "2% glucose/7 mM K⁺", "0.5% glucose/1 mM K⁺"),
  "Condition_3" = c("condition", "2% glucose/7 mM K⁺", "0.5% glucose/7 mM K⁺"),
  "Condition_4" = c("condition", "2% glucose/1 mM K⁺", "0.5% glucose/1 mM K⁺"),
  "Condition_5" = c("condition", "2% glucose/1 mM K⁺", "0.5% glucose/7 mM K⁺"),
  "Condition_6" = c("condition", "0.5% glucose/7 mM K⁺", "0.5% glucose/1 mM K⁺")
)

all_results <- list()
all_sig <- list()
all_ego <- list()

dir.create("~/GSE341442_analysis/results/DEG_results", recursive = TRUE, showWarnings = FALSE) # folder for CSVs
dir.create("~/GSE341442_analysis/results/GO_results", recursive = TRUE, showWarnings = FALSE) # folder for plots

#4. Extracting DEGs from conditions, filtering based on padj and log2fc values, Volcano plot and Gene Ontology Enrichment
for(nm in names(contrasts_list)){
  cat("\n=== Running:", nm, "===\n")
  
  # 1. DESeq2 results
  res <- results(dds, contrast = contrasts_list[[nm]])
  
  res$sgd <- rownames(res)
  res$systematic <- mapIds(org.Sc.sgd.db, keys = rownames(res),
                           column = "SGD", keytype = "SGD", multiVals = "first")
  res$genename <- mapIds(org.Sc.sgd.db, keys = rownames(res),
                         column = "GENENAME", keytype = "SGD", multiVals = "first")
  
  # 2. Filter significant genes and find upregualted and downregulated genes
  sig <- res[!is.na(res$padj) & res$padj < 0.05 & abs(res$log2FoldChange) > 1, ]
  up <- sum(sig$log2FoldChange > 0, na.rm=TRUE)
  down <- sum(sig$log2FoldChange < 0, na.rm=TRUE)
  #Top 10 upregulated genes

  print(head(sig[order(sig$log2FoldChange, decreasing = TRUE), ], 10))
  #Top 10 downregulated genes
  print(head(sig[order(sig$log2FoldChange, decreasing = FALSE), ], 10))
  #Number of up and downregulated genes
  print(up)
  print(down)
  
  # 3. Volcano plot
  p <- EnhancedVolcano(res,
                       lab = res$genename,
                       x = 'log2FoldChange',
                       y = 'padj',
                       title = nm,
                       pCutoff = 0.05,
                       FCcutoff = 1.0)
  print(p)
  ggsave(paste0("~/GSE341442_analysis/results/DEG_results/volcano_", nm, ".png"), plot = p, width = 8, height = 8)
  
  # 4. Save CSV
  write.csv(as.data.frame(sig), paste0("~/GSE341442_analysis/results/DEG_results/DEG_", nm, ".csv"), row.names = TRUE)
  
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
    ggsave(paste0("~/GSE341442_analysis/results/GO_results/dotplot_", nm, ".png"),
           dotplot(ego, showCategory = 20) + ggtitle(nm), width = 8, height = 10)
  } else {
    cat("No significant genes for", nm, "\n")
  }
  
  all_results[[nm]] <- res
  all_sig[[nm]] <- sig
}


