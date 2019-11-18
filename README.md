# Resistome Analysis

### Install dependencies in R:
* dplyr
* magrittr
* purrr
* reshape2
* ggplot2
* RColorBrewer
* grid
* vegan
* cluster
* ComplexHeatmap
* circlize
* dendextend
* DESeq2
* ggrepel
* stringr
* metafor

```R
install.packages(c("dplyr", "magrittr", "purrr", "reshape2", "ggplot2", "RColorBrewer", "vegan", "cluster", "circlize", "dendextend", "ggrepel", "stringr", "metafor"))

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
    
BiocManager::install(c("DESeq2", "ComplexHeatmap"))
```

### Install resistomeAnalysis package in R:

```R
install.packages("devtools")
library(devtools)
install_github("blue-moon22/resistomeAnalysis")
```
