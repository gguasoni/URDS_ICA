## Task
The focus of this assignment was to build an R workflow that ingests a gene-expression matrix and companion metadata plus a student-specific gene list, applies log-scaling and ID-to-LongName mapping, filters to the listed genes (deduplicated), and generates publication-ready visuals — two heatmaps (genes+samples clustered, and genes-only) and a boxplot comparing expression across treatment groups — complete with gene-type (XA/XB/XC) and treatment annotations, and brief notes explaining processing choices and figure interpretation.

## Script Details
- URDS_ICA.rmd
  - This R Markdown (configured to knit to a Word report in Arial 11pt) implements a defensively coded pipeline for a small gene-expression analysis: it auto-discovers and ingests all .csv/.txt files in the working directory, converts them into data frames, and performs basic hygiene checks (missing values, duplicate rows) before merging expression data with gene/sample annotations and matching a 40-gene subset from a larger list. The workflow maps gene IDs to LongName, log-transforms numeric columns, and reshapes to a tidy long format with treatment labels. Using ggplot2, tidyverse, and pheatmap, it then produces:
    (i) a scatterplot of log(expression) by gene, coloured by gene type and shaped by treatment group
    (ii) a heatmap clustering both genes and samples with annotated gene types and treatments
    (iii) a second heatmap clustering genes only; and (iv) a boxplot with jitter overlays summarizing expression distributions by treatment and gene type.
    Figure captions and a short limitations note (e.g., limited biological context and absence of an untreated control) are included for quick interpretation.

## How I would improve this code:
- This code was written in an precise working directory, where I was in control of all the present files. However, to make the code more broadly usable to any given data set, I would set up a commandline interface that first asks the user to enter the names of the files they would like the code to use.
- I would refrain from assigning paramters to the global environment, as this can be extremely risky when coding and hinders reproducibility.
- I would fix the NA/duplicate checks to operate on objects instead of their names, apply row-wise NA filtering instead of dropping whole columns, and log-transform numeric assay columns only with a pseudocount (log10(x + 1)) to avoid -Inf. I’ll add stopifnot() checks for required columns (Gene, LongName, Type, TreatmentGroup) and a concise QC section (dims, % missing, duplicates removed).
- This assignment required me to use a specific basee R boxplot style, but if I were allowed to freely create my own boxplot, I would use ggplot2 from the tidyverse as that allows for much more efficient figure creation in R.
- Instead of hard-coding treatment columns A:K, I would programmatically select numeric assay columns and avoid using 'sample' as a variable name.
- I would include a PCA of samples (on the log-scaled matrix) coloured by TreatmentGroup with % variance in axes as this procides additional information for exploratory analysis and interpretation.
