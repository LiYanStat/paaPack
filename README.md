The package **paaPack** provides R implementation of the Principal Amalgamation Analysis for microbiome data

## Overview

The **paaPack** is designed to perform hierarchical Principal Amalgamation Analysis (HPAA) with or without the 
guidance of taxonomic tree structure, and provide several useful graphical tools for visualizing the results of HPAA, 
including 

1) hierarchical dendrograms to visualize the full path of amalgamations, 
2) the scree plot showing the percentage change 
in the diversity loss along with the changes of number of compositions, and 
3) the ordination plot showing the changes in the between-sample distance patterns before and after HPAA with any given number of principal compositions (PC).

## Installation

One may install the latest version under development as follows:

```R
install.packages(c("devtools"))
devtools::install_github("LiYanStat/paaPack")
```

## Get Started

The **paaPack** provides a main function `hPAA()` to perform the hierarchical Principal Amalgamation Analysis (HPAA). 
To use the function, the analyst should provide:

+ The compositional data with row representing compositions for each subject and 
  column recording the original components/taxa.
+ The taxonomic tree structure as a vector with dimension same as the number of taxa. The format of the taxonomic 
  vector is the same as output format of the commonly used bioinformatics data processing software **Mothur**. 
  That is, each element of the vector denotes the full taxonomic ranks from kingdom to genus, species level of the taxon, 
  with ranks separated by semicolon. For example a typical element of the taxonomic vector could be
  `k_Bacteria;p_Actinobacteria;c_Actinobacteria;o_Bifidobacteriales;f_Bifidobacteriaceae;g_Bifidobacterium;s_longum`.
  The taxonomic structure is optional. If taxonomy is not provided, the unconstrained HPAA without tree guidance is performed.
  
+ The diversity measures used in the HPAA analysis, and indicate whether strong or weak taxonomic hierarchy is applied in the analysis.

Then a set of plotting methods are provided taking the object of class "hPAA" from the `hPAA()` as input, including `plotHPAA()` 
for the dendrogram showing the full path of hierarchical amalgamation, `plotLine()` for the scree plot showing the percentage change in 
the diversity loss with the changing in the number of principal compositions, and `plotMDS()` for the ordination plot showing the changes 
in the between-sample distance patterns before and after HPAA. In each function, a group of graphical arguments for shaping the figures could be 
specified. For details, the analyst could refer to the documentations of the functions using `help()`.

```R
library(paaPack)
#### functions in paaPack
help(hPAA)  #### fit HPAA models
help(plotHPAA)  #### dendrogram showing the hierarchical amalgamation
help(plotLine)  #### scree plot showing the percentage
help(plotMDS)  #### ordination plot
```

## Demo

A detailed demo with numerical examples is provided as Rmarkdown file **demo/paaPack.rmd** for reference.
