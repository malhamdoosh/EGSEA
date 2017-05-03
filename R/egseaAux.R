#Ensemble of Gene Set Enrichment Analyses 
# 
# Author: Monther Alhamdoosh, E:m.hamdoosh@gmail.com

#' @title EGSEA auxiliary functions
#' @description It lists the accepted sorting methods for analysis results
#' @return It returns a character vector of the accepted values for the sort.by 
#' argument in egsea
#' @export
#' 
#' @name egsea.sort
#' @aliases egsea.sort,egsea-aux
#' @rdname egsea-aux
#' 
#' @examples 
#' egsea.sort()
#' 

egsea.sort <-function(){
    return(c(c("p.value", "p.adj", "avg.rank", "med.rank", "min.rank", 
                            "min.pvalue", "vote.rank", "avg.logfc.dir", "avg.logfc", "direction",
                            "significance"),egsea.base()))
}


#' @title EGSEA P-value Combining Options
#' @description It lists the p-value combining methods
#' @return It returns a character vector of available methods for the 
#' combineMethod argument in egsea
#' @export
#' 
#' @name egsea.combine
#' @aliases egsea.combine,egsea-aux
#' @rdname egsea-aux
#' 
#' @examples 
#' egsea.combine()
#' 

egsea.combine <- function(){
    return(c("fisher", "wilkinson", "average", "logitp", 
                    "sump", "sumz", "median"))
}


#' @title EGSEA Base GSE Methods
#' @description It lists the supported GSEA methods. Since EGSEA base methods 
#' are implemented in the Bioconductor project, the most recent version of 
#' each individual method is always used. 
#' @details These methods include: 
#' \pkg{ora}[1], \pkg{globaltest}[2], \pkg{plage}[3], \pkg{safe}[4], \pkg{zscore}[5],
#'  \pkg{gage}[6], \pkg{ssgsea}[7], \pkg{roast}[8], \pkg{fry}[8], \pkg{padog}[9],
#'  \pkg{camera}[10] and \pkg{gsva}[11]. 
#' The \code{ora, gage, camera} and \code{gsva} methods depend on a competitive null 
#' hypothesis while the remaining seven methods are based on a self-contained hypothesis. 
#' Conveniently, EGSEA is not limited to these twelve
#' GSE methods and new GSE tests can be easily integrated into the framework. \cr
#' \cr
#' Note: the execution time of base methods can vary depending on the size of gene set collections, 
#' number of samples, number of genes and number of contrasts.  When a gene set collection of 
#' around 200 gene sets was tested on a dataset of 17,500 genes, 8 samples and 2 contrasts, the
#' execution  time of base methods in ascending order was as follows:
#' \code{globaltest; safe; gage; gsva; zscore; plage; fry; camera; ora; ssgsea; padog}. When the
#' same dataset was tested on a large gene set collection of 3,700 gene sets, the execution
#' time of base methods in ascending order was as follows: 
#' \code{globaltest; camera; fry; zscore; plage; safe; gsva; ora; gage; padog; ssgsea}. Apparently, the
#' size of gene set collection plays a key role in the execution time of most of the base
#' methods. The reduction rate of execution time between the large and small gene set
#' collections varied between 18\% and 88\%. \code{camera, fry, plage, zscore} and \code{ora} showed the least
#' reduction rate of execution time. As a result, there is no guarantee that a single combination 
#' of base methods would run faster than other combinations. It is worth mentioning that 
#' our simulation results showed that the increasing number of base methods in the EGSEA 
#' analysis is desirable to achieve high performance.  
#'
#' @references 
#' [1] Tavazoie, S. et al. (1999). Systematic determination of genetic network architecture.
#' Nature Genetics, 22(3), 281-5.\cr 
#' [2] Goeman, J. J. et al. (2004). A global test for groups of genes: testing association with a
#' clinical outcome. Bioinformatics, 20(1), 93-9.\cr 
#' [3] Tomfohr, J. et al. (2005). Pathway level analysis of gene expression using singular
#' value decomposition. BMC Bioinformatics, 6, 225.\cr 
#' [4] Barry, W. T. et al. (2005). Significance analysis of functional categories in gene
#' expression studies: a structured permutation approach. Bioinformatics, 21(9), 1943-9.\cr 
#' [5] Lee, E. et al. (2008). Inferring pathway activity toward precise disease classification.
#' PLoS Computational Biology, 4(11), e1000217.\cr 
#' [6] Luo, W. et al. (2009). GAGE: generally applicable gene set enrichment for pathway
#' analysis. BMC Bioinformatics, 10, 161.\cr 
#' [7] Barbie, D. A. et al. (2009). Systematic RNA interference reveals that oncogenic KRASdriven
#' cancers require TBK1. Nature, 462(7269), 108-12.\cr 
#' [8] Wu, D. et al. (2010). ROAST: rotation gene set tests for complex microarray
#' experiments. Bioinformatics, 26(17), 2176-82.\cr 
#' [9] Tarca, A. L. et al. (2009). A novel signaling pathway impact analysis. Bioinformatics,
#' 25(1), 75-82.\cr 
#' [10] Wu, D. and Smyth, G. K. (2012). Camera: a competitive gene set test accounting for
#' inter-gene correlation. Nucleic Acids Research, 40(17), e133.\cr 
#' [11] Hanzelmann, S. et al. (2013). GSVA: gene set variation analysis for microarray and
#' RNA-seq data. BMC Bioinformatics, 14, 7.
#' 
#' 
#' @return It returns a character vector of supported GSE methods. 
#' @export
#' 
#' @name egsea.base
#' @aliases egsea.base,egsea-aux
#' @rdname egsea-aux
#' 
#' @examples
#' egsea.base()
#' 

egsea.base <- function(){
    return(c("camera", "roast", "safe", "gage", "padog", "plage", "zscore", 
                    "gsva", "ssgsea", "globaltest", "ora", "fry")) 
    # , "SPIA", "gsea", "samgs"
}

#' @title EGSEA logo
#' @description This function writes out the official EGSEA package logo
#' @details This function generates a PNG file of the EGSEA logo, which can
#' be used to acknowledge EGSEA in presentations/reports. 
#' The logo was designed by Roberto Bonelli from The Walter and Eliza Hall Institute
#' of Medical Research. 
#' @param out.dir character, the target directory to which the logo will be written. 
#' @return a PNG file. 
#' 
#' @export 
#' 
#' @name egsea.logo
#' @aliases egsea.logo,egsea-aux
#' @rdname egsea-aux

egsea.logo <- function(out.dir = "./"){
    logo.file = system.file("logo", "EGSEA_logo.png", package="EGSEA")
    if (file.exists(logo.file)){
        file.copy(logo.file, out.dir)
    }else{
        cat("EGSEA logo was not installed\n")
    }
}