# Ensemble of Gene Set Enrichment Analyses

<div align="center">
<img src="inst/logo/EGSEA_logo.png" align="middle" width=250 />
<br />
<sub>Credit: Roberto Bonelli </sub>
</div>

This package is part of the **Bioconductor** project and implements the Ensemble of Gene Set Enrichment Analyses (EGSEA) method for gene set testing.

**Author:** Monther Alhamdoosh, Luyi Tian, Milica Ng and Matthew Ritchie

**Maintainer:** Monther Alhamdoosh <m.hamdoosh at gmail.com>

Citation (from within R, enter ```citation("EGSEA")```):

Alhamdoosh M, Ng M, Wilson N, Sheridan J, Huynh H, Wilson M and Ritchie M (2017). “Combining multiple tools outperforms individual methods in gene set enrichment analyses.” *Bioinformatics*, 33(3). doi: 10.1093/bioinformatics/btw623.

Alhamdoosh M, Law CW, Tian L et al. Easy and efficient ensemble gene set testing with EGSEA [version 1; peer review: 1 approved, 3 approved with reservations]. F1000Research 2017, 6:2010. doi: 10.12688/f1000research.12544.1

# Installation

To install the *stable release* of this package, start R and enter:
```{r}
## try http:// if https:// URLs are not supported
if (!requireNamespace("BiocManager", quietly=TRUE))
    install.packages("BiocManager")
BiocManager::install("EGSEA")
```

To install the *development version* of this package, start R and enter:
```{r}
install.packages("devtools") # if you have not installed "devtools" package
devtools::install_github("malhamdoosh/EGSEA")
```

# Documentation

To view documentation for the version of this package installed in your system, start R and enter:
```{r}
browseVignettes("EGSEA")
```


