# GSE341442_yeast-K-glucose-RNASeq-DESeq2
# Yeast Transcriptomics under Low K⁺ and Glucose Deprivation (GSE341442)

RNA-seq analysis of *S. cerevisiae* using DESeq2. 6 contrasts comparing potassium & glucose conditions.

## Dataset
- **GEO:** GSE341442
- **Samples:** 24 (4 conditions x 2 timepoints x 3 replicates)
- **Organism:** Saccharomyces cerevisiae S288C

## Experimental Design

| Contrast | Comparison | Question |
| :--- | :--- | :--- |
| Condition 1 | Low K (0.1mM) vs Control K (10mM) @ 60min | Effect of K⁺ deprivation |
| Condition 2 | Low K + No Glucose vs Low K | Effect of glucose removal under low K |
| Condition 3 | No Glucose vs Control @ 60min | Glucose starvation |
| Condition 4 | Low K (0mM) vs Control @ 120min | Long-term K⁺ effect |
| Condition 5 | Low K + Low Glucose vs Low K | Mild glucose limitation |
| Condition 6 | No Glucose vs Low Glucose | Severe vs mild glucose stress |

## Results Summary

| Contrast | Upregulated | Downregulated | Top Gene |
| :--- | :--- | :--- | :--- |
| Cond 1 | 342 | 289 | TDA6 |
| Cond 2 | 1120 | 980 | CTS1 |
| Cond 3 | ... | ... | ... |
| Cond 4 | ... | ... | ... |
| Cond 5 | ... | ... | ... |
| Cond 6 | ... | ... | ... |

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
