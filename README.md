# LncPath
Identifying the Pathways Regulated by LncRNA Sets of Interest

> Identifies pathways synergisticly regulated by the interested lncRNA(long non-coding RNA) sets based on a lncRNA-mRNA(messenger RNA) interaction network. 1) The lncRNA-mRNA interaction network was built from the protein-protein interactions and the lncRNA-mRNA co-expression relationships in 28 RNA-Seq data sets. 2) The interested lncRNAs can be mapped into networks as seed nodes and a random walk strategy will be performed to evaluate the rate of each coding genes influenced by the seed lncRNAs. 3) Pathways regulated by the lncRNA set will be evaluated by a weighted Kolmogorov-Smirnov statistic as an ES Score. 4) The p value and false discovery rate value will also be calculated through a permutation analysis. 5) The running score of each pathway can be plotted and the heat map of each pathway can also be plotted if an expression profile is provided. 6) The rank and scores of the gene list of each pathway can be printed.

### how to install

Install the **development version** from Github:

```R
Installation method：

1. library(devtools); 
   install_github("hanjunwei-lab/LncPath")
2. install.packages("LncPath")

Use：
library(LncPath)
```

Please cite the following article when using `LncPath`:

Han, J., S. Liu, Z. Sun, Y. Zhang, F. Zhang, C. Zhang, D. Shang, H. Yang, F. Su, Y. Xu, C. Li, H. Ren, and X. Li, LncRNAs2Pathways: Identifying the pathways influenced by a set of lncRNAs of interest based on a global network propagation method. Sci Rep, 2017. 7: p. 46566.