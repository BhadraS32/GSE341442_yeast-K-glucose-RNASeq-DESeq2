# GSE341442_yeast-K-glucose-RNASeq-DESeq2
# Yeast Transcriptomics under Low K⁺ and Glucose Deprivation (GSE341442)

RNA-seq analysis of *S. cerevisiae* using DESeq2. 6 contrasts comparing potassium & glucose conditions.

## Dataset
- **GEO:** GSE341442
- **Samples:** 12 (4 conditions x 3 replicates)
- **Organism:** Saccharomyces cerevisiae S288C

## Experimental Design

| Contrast | Comparison | Question |
| :--- | :--- | :--- |
| Condition 1 | 2% glucose/7 mM K⁺ vs 2% glucose/1 mM K⁺ | Effect of K⁺ in standard glucose media |
| Condition 2 | 2% glucose/7 mM K⁺ vs 0.5% glucose/1 mM K⁺ | Comparing standard glucose media and high K⁺ to limited glucose media and low K⁺ |
| Condition 3 | 2% glucose/7 mM K⁺ vs 0.5% glucose/7 mM K⁺ | Glucose starvation with high potassium- comparing standard glucose media to low glucose media |
| Condition 4 | 2% glucose/1 mM K⁺ vs 0.5% glucose/1 mM K⁺ | Glucose starvation with low potassium- comparing standard glucose media to low glucose media|
| Condition 5 | 2% glucose/1 mM K⁺ vs 0.5% glucose/7 mM K⁺ | Comparing standard glucose, low potassium with low glucose, high potassium media |
| Condition 6 | 0.5% glucose/7 mM K⁺ vs 0.5% glucose/1 mM K⁺| Comparing media with high and low potassium in glucose starvation |

## Results Summary

| Contrast | Upregulated | Downregulated | Top Gene |
| :--- | :--- | :--- | :--- |
| Condition 1 | 415 | 421 | IRT1 |   
| Condition 2 | 727 | 806 | YRO2 |   
| Condition 3 | 618 | 433 | TPO2 |   
| Condition 4 | 418 | 576 | ENA2 |   
| Condition 5 | 804 | 608 | CIS3 |   
| Condition 6 | 661 | 949 | ENA2 |   

*Cutoff: padj < 0.05, |log2FC| > 1*

## Pipeline

| Step | Tool | Package |
| :--- | :--- | :--- |
| QC & Filtering | R | DESeq2 |
| Normalization | VST | DESeq2 |
| DEG | Wald test | DESeq2 |
| Visualization | Volcano, PCA, Heatmap | EnhancedVolcano, pheatmap |
| Enrichment | GO BP | clusterProfiler |

## Folder Structure
GSE341442_analysis/
├── Data/ - GSE341442_Counts_allSamples.csv
├── results/DEG_results/ - 6 DEG CSVs 6 volcano PNGs 
├── results/GO_results/ - GO dotplots 
├── results/PCA_and_Heatmap/ - PCA and heatmap
└── script/ - DESeq2_pipeline_final.R

