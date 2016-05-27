#Ensemble of Gene Set Enrichment Analyses 
# 
# Author: Monther Alhamdoosh, E:m.hamdoosh@gmail.com
#R (>= 3.0.1),
#gage (>= 2.14.4),
#PADOG (>= 1.6.0),
#GSVA (>= 1.12.0),
#globaltest (>= 5.18.0),
#limma (>= 3.20.9),
#edgeR (>= 3.6.8),
#pathview (>= 1.4.2),
#HTMLUtils (>= 0.1.5),
#hwriter (>= 1.3.2),
#gplots (>= 2.14.2),
#ggplot2 (>= 1.0.0),
#safe (>= 3.4.0),
#topGO (>= 2.16.0),
#stringi (>= 0.5.0),
#XML (>= 3.98.1.3),
#GO.db (>= 3.1.2),    
#EGSEAdata (>= 0.99.0)
###############################################################################
#' @name EGSEA-package
#' @aliases EGSEA
#' @docType package
#' @title Ensemble of Gene Enrichment Analysis (EGSEA)
#' @author Monther Alhamdoosh, Milica Ng and Matthew Ritchie 
#' @description This packages provides the implementatino of the EGSEA 
#' algorithm 
#' and addition functions to help perform GSE analysis
NULL



#' @title Ensemble of Gene Set Enrichment Analyses Function
#' 
#' @description This is the main function to carry out gene set enrichment 
#' analysis using the
#'  EGSEA algorithm. This function is aimed to extend the limma-voom pipeline 
#' of RNA-seq analysis. 
#'
#' @details EGSEA, an acronym for \emph{Ensemble of Gene Set Enrichment 
#' Analyses}, utilizes the 
#' analysis results of eleven prominent GSE algorithms from the literature to 
#' calculate 
#' collective significance scores for gene sets. These methods include: 
#' \pkg{ora}, 
#' \pkg{globaltest}, \pkg{plage}, \pkg{safe}, \pkg{zscore}, \pkg{gage}, 
#' \pkg{ssgsea},
#' \pkg{roast}, \pkg{fry}, \pkg{padog}, \pkg{camera} and \pkg{gsva}. 
#' The ora, gage, camera and gsva methods depend on a competitive null 
#' hypothesis while the 
#' remaining seven methods are based on a self-contained hypothesis. 
#' Conveniently, the 
#' algorithm proposed here is not limited to these eleven GSE methods and new 
#' GSE tests 
#' can be easily integrated into the framework. This function takes the voom 
#' object and 
#' the contrast matrix as parameters. 
#' The results of EGSEA can be seen using the \code{\link{topSets}} function. 
#'
#' @param voom.results list, an EList object generated using the  
#' \code{\link[limma]{voom}} function. 
#' Entrez Gene IDs should be used as row names.  
#' @param contrasts double, an N x L matrix indicates the contrast of the 
#' linear model coefficients for
#' which the test is required. N is number of experimental conditions and L is 
#' number of contrasts.
#' @param logFC double, an K x L matrix indicates the log2 fold change of each 
#' gene for each contrast. 
#' K is the number of genes included in the analysis. If logFC=NULL, the logFC 
#' values are 
#' estimated using the  \code{\link[limma]{ebayes}} for each contrast. 
#' @param gs.annots list, indexed collections of gene sets. It is generated 
#' using one of these functions: 
#'  \code{\link{buildIdxEZID}}, \code{\link{buildMSigDBIdxEZID}}, 
#' \code{\link{buildKEGGIdxEZID}},
#' \code{\link{buildGeneSetDBIdxEZID}}, and \code{\link{buildCustomIdxEZID}}. 
#' @param symbolsMap dataframe, an K x 2 matrix stores the gene symbol of each 
#' Entrez Gene ID. It 
#' is used for the heatmap visualization. The order of rows should match that 
#' of the
#'  \strong{voom.results}. Default symbolsMap=NULL. 
#' @param baseGSEAs character, a vector of the gene set tests that should be 
#' included in the 
#' ensemble. Type  \code{\link{egsea.base}} to see the supported GSE methods. 
#' By default, all
#' supported methods are used.
#' @param minSize integer, the minimum size of a gene set to be included in the 
#' analysis. 
#' Default minSize= 2. 
#' @param display.top integer, the number of top gene sets to be displayed in 
#' the EGSEA report.
#' You can always access the list of all tested gene sets using the returned 
#' gsa list. 
#' Default is 20.
#' @param combineMethod character, determines how to combine p-values from 
#' different 
#' GSEA method. Type  \code{\link{egsea.combine}}() to see supported methods.
#' @param combineWeights double, a vector determines how different GSEA methods 
#' will be weighted. 
#' Its values should range between 0 and 1. This option is not supported 
#' currently.
#' @param sort.by character, determines how to order the analysis results in 
#' the stats table. Type  
#' \code{\link{egsea.sort}}() to see all available options.
#' @param egsea.dir character, directory into which the analysis results are 
#' written out. 
#' @param kegg.dir character, the directory of KEGG pathway data file (.xml) 
#' and image file (.png). 
#' Default kegg.dir=paste0(egsea.dir, "/kegg-dir/").
#' @param logFC.cutoff numeric, cut-off threshold of logFC and is used for 
#' Sginificance Score 
#' and Regulation Direction Calculations. Default logFC.cutoff=0.
#' @param sum.plot.axis character, the x-axis of the summary plot. All the 
#' values accepted by the 
#' \strong{sort.by} parameter can be used. Default sum.plot.axis="p.value".
#' @param sum.plot.cutoff numeric, cut-off threshold to filter the gene sets of 
#' the summary plots
#' based on the values of the \strong{sum.plot.axis}. Default 
#' sum.plot.cutoff=NULL.
#' @param vote.bin.width numeric, the bin width of the vote ranking. Default 
#' vote.bin.width=5.
#' @param num.threads numeric, number of CPU threads to be used. Default 
#' num.threads=2.
#' @param report logical, whether to generate the EGSEA interactive report. It 
#' takes longer time
#' to run. Default is True.
#' @param print.base logical, whether to write out the results of the 
#' individual GSE methods.
#'  Default FALSE.  
#' @param verbose logical, whether to print out progress messages and warnings. 
#' 
#' @return A list of elements, each with two/three elements that store the top 
#' gene sets and the detailed analysis
#' results for each contrast and the comparative analysis results.
#' 
#' @seealso \code{\link{topSets}}, \code{\link{egsea.base}}, 
#' \code{\link{egsea.sort}}, 
#' \code{\link{buildIdxEZID}}, \code{\link{buildMSigDBIdxEZID}}, 
#' \code{\link{buildKEGGIdxEZID}},
#' \code{\link{buildGeneSetDBIdxEZID}}, and \code{\link{buildCustomIdxEZID}}
#'
#' @importFrom limma camera mroast fry lmFit  contrasts.fit eBayes topTable plotMDS
#' @importFrom globaltest gt gt.options result
#' @importFrom PADOG padog
#' @importFrom GSVA gsva
#' @importFrom safe safe safe.toptable getCmatrix
#' @importFrom parallel mclapply
#' @importFrom AnnotationDbi Ontology
#' @importFrom stats p.adjust pchisq phyper
#' @importFrom grDevices pdf dev.off png 
#' @importFrom graphics par
#' @importFrom utils write.table browseURL write.csv data capture.output
#' @importFrom gage gage kegg.gsets
#' @importFrom metap logitp meanp sumlog sump sumz wilkinsonp
#' @import hwriter HTMLUtils stringi ggplot2 pathview gplots  Biobase topGO
#' @export
#' @references 
#' Monther Alhamdoosh, Milica Ng, Nicholas J. Wilson, Julie M. Sheridan, Huy 
#' Huynh, Michael J. Wilson
#' and Matthew E. Ritchie. Combining multiple tools outperforms individual 
#' methods in gene set enrichment
#' analyses. 
#' 
#' @examples
#' library(EGSEAdata)
#' data(il13.data)
#' v = il13.data$voom
#' contrasts = il13.data$contra
#' gs.annots = buildIdxEZID(entrezIDs=rownames(v$E), species="human", 
#' msigdb.gsets="none", 
#'          kegg.updated=FALSE, kegg.exclude = c("Metabolism"))
#' # set report = TRUE to generate the EGSEA interactive report
#' gsa = egsea(voom.results=v, contrasts=contrasts,  gs.annots=gs.annots, 
#'          symbolsMap=v$genes, baseGSEAs=egsea.base()[-c(2,5,6,9,12)], 
#' 			display.top = 5, sort.by="avg.rank", 
#' 			egsea.dir="./il13-egsea-report", 
#'          num.threads = 2, report = FALSE)
#' topSets(gsa) 
#' 


egsea <- function(voom.results, contrasts, logFC=NULL,
        gs.annots, symbolsMap=NULL, 
        baseGSEAs=egsea.base(),
        minSize=2, display.top=20, 
        combineMethod="fisher", combineWeights = NULL,      
        sort.by="p.adj", 
        egsea.dir="./", kegg.dir=NULL, 
        logFC.cutoff=0, sum.plot.axis="p.adj", sum.plot.cutoff=NULL, 
        vote.bin.width=5,
        num.threads=4, report = TRUE,
        print.base = FALSE, 
        verbose=FALSE){
    if (is.null(sum.plot.cutoff)){
        if (sum.plot.axis %in% c("p.value", "p.adj"))
            sum.plot.cutoff = 1     
        else
            sum.plot.cutoff = 10000
    }
    if (length(baseGSEAs) == 1){
        sort.by = "p.adj"
    }
    set.seed(581986) # to gurantee reproducibility of results
    if (verbose)    
        return(egsea.main(voom.results, contrasts, gs.annots, baseGSEAs, 
combineMethod, 
                        combineWeights,sort.by,  egsea.dir, 
                        kegg.dir, logFC, symbolsMap, minSize, display.top, 
logFC.cutoff, sum.plot.cutoff, sum.plot.axis, 
                        vote.bin.width, print.base, verbose, num.threads, 
report))
    else
        suppressWarnings(
                return(egsea.main(voom.results, contrasts, gs.annots, 
baseGSEAs, combineMethod, 
                                combineWeights,sort.by,  egsea.dir, 
                                kegg.dir, logFC, symbolsMap, minSize, 
display.top, logFC.cutoff, sum.plot.cutoff, 
                                sum.plot.axis, vote.bin.width,print.base, 
verbose, num.threads, report)))
    
}


#' @title Ensemble of Gene Set Enrichment Analyses Function
#' 
#' @description This is the main function to carry out gene set enrichment 
#' analysis using the
#'  EGSEA algorithm. This function is aimed to use the raw count matrix to 
#' perform the EGSEA analysis. 
#'
#' @details EGSEA, an acronym for \emph{Ensemble of Gene Set Enrichment 
#' Analyses}, utilizes the 
#' analysis results of eleven prominent GSE algorithms from the literature to 
#' calculate 
#' collective significance scores for gene sets. These methods include: 
#' \pkg{ora}, 
#' \pkg{globaltest}, \pkg{plage}, \pkg{safe}, \pkg{zscore}, \pkg{gage}, 
#' \pkg{ssgsea},
#' \pkg{roast}, \pkg{fry}, \pkg{padog}, \pkg{camera} and \pkg{gsva}. 
#' The ora, gage, camera and gsva methods depend on a competitive null 
#' hypothesis while the 
#' remaining seven methods are based on a self-contained hypothesis. 
#' Conveniently, the 
#' algorithm proposed here is not limited to these eleven GSE methods and new 
#' GSE tests 
#' can be easily integrated into the framework. This function takes the raw 
#' count matrix, 
#' the experimental group of each sample, the design matrix and the contrast 
#' matrix as parameters. 
#' It performs TMM normalization and then applies \link[limma]{voom} to 
#' calculate the logCPM and weighting factors. 
#' The results of EGSEA can be seen using the \code{\link{topSets}} function. 
#'
#' @param counts double, numeric matrix of read counts where genes are the rows 
#' and samples are
#' the columns.
#' @param group character, vector or factor giving the experimental 
#' group/condition for each sample/library
#' @param design double, numeric matrix giving the design matrix of the linear 
#' model fitting.
#' @param contrasts double, an N x L matrix indicates the contrast of the 
#' linear model coefficients for
#' which the test is required. N is number of experimental conditions and L is 
#' number of contrasts.
#' @param logFC double, an K x L matrix indicates the log2 fold change of each 
#' gene for each contrast. 
#' K is the number of genes included in the analysis. If logFC=NULL, the logFC 
#' values are 
#' estimated using the  \code{\link[limma]{eBayes}} for each contrast. 
#' @param gs.annots list, indexed collections of gene sets. It is generated 
#' using one of these functions: 
#'  \code{\link{buildIdxEZID}}, \code{\link{buildMSigDBIdxEZID}}, 
#' \code{\link{buildKEGGIdxEZID}},
#' \code{\link{buildGeneSetDBIdxEZID}}, and \code{\link{buildCustomIdxEZID}}. 
#' @param symbolsMap dataframe, an K x 2 matrix stores the gene symbol of each 
#' Entrez Gene ID. It 
#' is used for the heatmap visualization. The order of rows should match that 
#' of the 
#' \strong{counts}. Default symbolsMap=NULL. 
#' @param baseGSEAs character, a vector of the gene set tests that should be 
#' included in the 
#' ensemble. Type  \code{\link{egsea.base}} to see the supported GSE methods. 
#' By default, all
#' supported methods are used.
#' @param minSize integer, the minimum size of a gene set to be included in the 
#' analysis. 
#' Default minSize= 2. 
#' @param display.top integer, the number of top gene sets to be displayed in 
#' the EGSEA report.
#' You can always access the list of all tested gene sets using the returned 
#' gsa list. 
#' Default is 20.
#' @param combineMethod character, determines how to combine p-values from 
#' different 
#' GSEA method. Type  \code{\link{egsea.combine}}() to see supported methods.
#' @param combineWeights double, a vector determines how different GSEA methods 
#' will be weighted. 
#' Its values should range between 0 and 1. This option is not supported 
#' currently.
#' @param sort.by character, determines how to order the analysis results in 
#' the stats table. Type  
#' \code{\link{egsea.sort}}() to see all available options.
#' @param egsea.dir character, directory into which the analysis results are 
#' written out. 
#' @param kegg.dir character, the directory of KEGG pathway data file (.xml) 
#' and image file (.png). 
#' Default kegg.dir=paste0(egsea.dir, "/kegg-dir/").
#' @param logFC.cutoff numeric, cut-off threshold of logFC and is used for 
#' Sginificance Score 
#' and Regulation Direction Calculations. Default logFC.cutoff=0.
#' @param sum.plot.axis character, the x-axis of the summary plot. All the 
#' values accepted by the 
#' \strong{sort.by} parameter can be used. Default sum.plot.axis="p.value".
#' @param sum.plot.cutoff numeric, cut-off threshold to filter the gene sets of 
#' the summary plots
#' based on the values of the \strong{sum.plot.axis}. Default 
#' sum.plot.cutoff=NULL.
#' @param vote.bin.width numeric, the bin width of the vote ranking. Default 
#' vote.bin.width=5.
#' @param num.threads numeric, number of CPU threads to be used. Default 
#' num.threads=2.
#' @param report logical, whether to generate the EGSEA interactive report. It 
#' takes longer time
#' to run. Default is True.
#' @param print.base logical, whether to write out the results of the 
#' individual GSE methods.
#'  Default FALSE.  
#' @param verbose logical, whether to print out progress messages and warnings. 
#' 
#' @return A list of elements, each with two/three elements that store the top 
#' gene sets and the detailed analysis
#' results for each contrast and the comparative analysis results.
#' 
#' @seealso \code{\link{topSets}}, \code{\link{egsea.base}}, 
#' \code{\link{egsea.sort}}, 
#' \code{\link{buildIdxEZID}}, \code{\link{buildMSigDBIdxEZID}}, 
#' \code{\link{buildKEGGIdxEZID}},
#' \code{\link{buildGeneSetDBIdxEZID}}, and \code{\link{buildCustomIdxEZID}}
#'
#' @importFrom  edgeR calcNormFactors DGEList
#' @importFrom limma voom
#' @importFrom stats model.matrix
#' @export
#' @references 
#' Monther Alhamdoosh, Milica Ng, Nicholas J. Wilson, Julie M. Sheridan, Huy 
#' Huynh, Michael J. Wilson
#' and Matthew E. Ritchie. Combining multiple tools outperforms individual 
#' methods in gene set enrichment
#' analyses. 
#' 
#' @examples
#' library(EGSEAdata)
#' data(il13.data.cnt)
#' cnt = il13.data.cnt$counts
#' group = il13.data.cnt$group
#' design = il13.data.cnt$design
#' contrasts = il13.data.cnt$contra
#' genes = il13.data.cnt$genes
#' gs.annots = buildIdxEZID(entrezIDs=rownames(cnt), species="human", 
#' msigdb.gsets="none",
#'          kegg.updated=FALSE, kegg.exclude = c("Metabolism"))
#' # set report = TRUE to generate the EGSEA interactive report
#' gsa = egsea.cnt(counts=cnt, group=group, design=design, contrasts=contrasts, 
#'          gs.annots=gs.annots, 
#'          symbolsMap=genes, baseGSEAs=egsea.base()[-c(2,5,6,9,12)], 
#' display.top = 5,
#'           sort.by="avg.rank", 
#' egsea.dir="./il13-egsea-cnt-report", 
#'          num.threads = 2, report = FALSE)
#' topSets(gsa) 
#' 

egsea.cnt <- function(counts, group, design = NULL, contrasts, logFC=NULL,
        gs.annots, symbolsMap=NULL, 
        baseGSEAs=egsea.base(),
        minSize=2, display.top=20, 
        combineMethod="fisher", combineWeights = NULL,      
        sort.by="p.adj", 
        egsea.dir="./", kegg.dir=NULL, 
        logFC.cutoff=0, sum.plot.axis="p.adj", sum.plot.cutoff=NULL, 
        vote.bin.width=5,
        num.threads=4, report = TRUE,
        print.base = FALSE, 
        verbose=FALSE){
    d = DGEList(counts, group=group)
    d = calcNormFactors(d, method="TMM")    
    if (is.null(design)){
        design = model.matrix(~0+group)
        colnames(design) = levels(factor(group))
    }
    voom.results = voom(d, design=design, plot=FALSE)
    return (egsea(voom.results, contrasts, logFC,
                    gs.annots, symbolsMap, 
                    baseGSEAs,
                    minSize, display.top, 
                    combineMethod, combineWeights,      
                    sort.by, 
                    egsea.dir, kegg.dir, 
                    logFC.cutoff, sum.plot.axis, sum.plot.cutoff, 
                    vote.bin.width,
                    num.threads, report,
                    print.base , 
                    verbose))
    
}

#' @title Over-representation Analysis with EGSEA Reporting Capabilities 
#' 
#' @description This is the main function to carry out gene set enrichment 
#' analysis using the
#'  over-representation analysis (ORA) only. 
#'
#' @details This function takes a list of Entrez gene IDs and uses the gene set 
#' collections 
#' from \pkg{EGSEAdata} or a custom-built collection to find over-represented 
#' gene sets in 
#' this list. It takes the advantage of the existing EGSEA reporting 
#' capabilities and generate 
#' an interative report for the ORA analysis. 
#' The results can be explored using the \code{\link{topSets}} function. 
#'
#' @param entrezIDs character, a vector of Entrez Gene IDs to be tested for ORA.
#' @param universe character, a vector of Enterz IDs to be used as a background 
#' list. If 
#' universe=NULL, the background list is created from the \pkg{AnnotationDbi} 
#' package.
#' @param logFC double, is a matrix or vector of log fold changes of the same 
#' length of entrezIDs.
#' If logFC=NULL, 1 is used as a default value. Then, the regulation direction 
#' in heatmaps
#' and pathway maps is not indicative to the gene regulation direction.
#' @param title character, a short description of the experimental contrast. 
#' @param gs.annots list, indexed collections of gene sets. It is generated 
#' using one of these functions: 
#'  \code{\link{buildIdxEZID}}, \code{\link{buildMSigDBIdxEZID}}, 
#' \code{\link{buildKEGGIdxEZID}},
#' \code{\link{buildGeneSetDBIdxEZID}}, and \code{\link{buildCustomIdxEZID}}. 
#' @param symbolsMap dataframe, an K x 2 matrix stores the gene symbol of each 
#' Entrez Gene ID. It 
#' is used for the heatmap visualization. The order of rows should match that 
#' of the 
#' \strong{entrezIDs}. Default symbolsMap=NULL. 
#' @param minSize integer, the minimum size of a gene set to be included in the 
#' analysis. 
#' Default minSize= 2. 
#' @param display.top integer, the number of top gene sets to be displayed in 
#' the EGSEA report.
#' You can always access the list of all tested gene sets using the returned 
#' gsa list. 
#' Default is 20.
#' @param sort.by character, determines how to order the analysis results in 
#' the stats table. 
#' It takes "p.value", "p.adj" or "Significance". 
#' @param egsea.dir character, directory into which the analysis results are 
#' written out. 
#' @param kegg.dir character, the directory of KEGG pathway data file (.xml) 
#' and image file (.png). 
#' Default kegg.dir=paste0(egsea.dir, "/kegg-dir/").
#' @param logFC.cutoff numeric, cut-off threshold of logFC and is used for 
#' Sginificance Score 
#' and Regulation Direction Calculations. Default logFC.cutoff=0.
#' @param sum.plot.axis character, the x-axis of the summary plot. All the 
#' values accepted by the 
#' \strong{sort.by} parameter can be used. Default sum.plot.axis="p.adj".
#' @param sum.plot.cutoff numeric, cut-off threshold to filter the gene sets of 
#' the summary plots
#' based on the values of the \strong{sum.plot.axis}. Default 
#' sum.plot.cutoff=NULL.
#' @param vote.bin.width numeric, the bin width of the vote ranking. Default 
#' vote.bin.width=5.
#' @param num.threads numeric, number of CPU threads to be used. Default 
#' num.threads=2.
#' @param report logical, whether to generate the EGSEA interactive report. It 
#' takes longer time
#' to run. Default is True.
#' @param print.base logical, whether to write out the results of the 
#' individual GSE methods.
#'  Default FALSE.  
#' @param verbose logical, whether to print out progress messages and warnings. 
#' 
#' @return A list of elements, each with two/three elements that store the top 
#' gene sets and the detailed analysis
#' results for each contrast and the comparative analysis results.
#' 
#' @seealso \code{\link{topSets}},
#' \code{\link{buildIdxEZID}}, \code{\link{buildMSigDBIdxEZID}}, 
#' \code{\link{buildKEGGIdxEZID}},
#' \code{\link{buildGeneSetDBIdxEZID}}, and \code{\link{buildCustomIdxEZID}}
#'
#' @importFrom org.Hs.eg.db org.Hs.egGO2ALLEGS
#' @importFrom org.Mm.eg.db org.Mm.egGO2ALLEGS
#' @importFrom org.Rn.eg.db org.Rn.egGO2ALLEGS
#' @export
#' @references 
#' Monther Alhamdoosh, Milica Ng, Nicholas J. Wilson, Julie M. Sheridan, Huy 
#' Huynh, Michael J. Wilson
#' and Matthew E. Ritchie. Combining multiple tools outperforms individual 
#' methods in gene set enrichment
#' analyses. 
#' 
#' @examples
#' library(EGSEAdata)
#' data(il13.data)
#' voom.results = il13.data$voom
#' contrast = il13.data$contra
#' library(limma)
#' vfit = lmFit(voom.results, voom.results$design)
#' vfit = contrasts.fit(vfit, contrast)
#' vfit = eBayes(vfit)
#' top.Table = topTable(vfit, coef=1, number=Inf, p.value=0.05, lfc=1)
#' deGenes = as.character(top.Table$FeatureID)
#' logFC =  top.Table$logFC
#' names(logFC) = deGenes
#' gs.annots = buildIdxEZID(entrezIDs=deGenes, species="human", 
#' msigdb.gsets="none",
#'          kegg.updated=FALSE, kegg.exclude = c("Metabolism"))
#' # set report = TRUE to generate the EGSEA interactive report
#' gsa = egsea.ora(entrezIDs=deGenes, universe= 
#' as.character(voom.results$genes[,1]),
#'              logFC =logFC, title="X24IL13-X24",  
#' gs.annots=gs.annots, 
#'              symbolsMap=top.Table[, c(1,2)], display.top = 5,
#'               egsea.dir="./il13-egsea-ora-report", num.threads = 2, 
#' 				report = FALSE)
#' topSets(gsa) 

egsea.ora <- function(entrezIDs, universe=NULL, logFC=NULL, title=NULL, 
        gs.annots, symbolsMap=NULL, 
        minSize=2, display.top=20, sort.by = "p.adj",   
        egsea.dir="./", kegg.dir=NULL, 
        logFC.cutoff=0, sum.plot.axis="p.adj", sum.plot.cutoff=NULL, 
        vote.bin.width=5,
        num.threads=4, report = TRUE,
        print.base = FALSE, 
        verbose=FALSE){
    voom.results = list(ids=as.character(entrezIDs))
    if (!is.null(universe))
        voom.results$featureIDs = universe
    if (is.null(logFC)){
        logFC = matrix(1, length(entrezIDs), 1)
        rownames(logFC) = entrezIDs     
    } else if (!is.matrix(logFC)){
        tmp = logFC
        logFC = matrix(1, length(entrezIDs), 1)
        logFC[, 1] = tmp
        if (!is.null(names(tmp)))
            rownames(logFC) = as.character(names(tmp))      
        else
            rownames(logFC) = voom.results$ids
    }
    if (is.null(title))
        colnames(logFC) = c("Experimental contrast")
    else
        colnames(logFC) = title
    
    dummyContrasts = matrix(0, 2, 1)
    dummyContrasts[, 1] = c(1,-1)
    colnames(dummyContrasts) = colnames(logFC)
    
    return (egsea(voom.results, contrasts = dummyContrasts, logFC,
                    gs.annots, symbolsMap, 
                    baseGSEAs = c("ora"),
                    minSize, display.top, 
                    combineMethod = "average", combineWeights = NULL,       
                    sort.by, 
                    egsea.dir, kegg.dir, 
                    logFC.cutoff, sum.plot.axis, sum.plot.cutoff, 
                    vote.bin.width,
                    num.threads, report,
                    print.base , 
                    verbose))

}

#' @title Table of Top Gene Sets from an EGSEA Analysis
#' @description Extract a table of the top-ranked gene sets from an EGSEA 
#' analysis.
#' 
#' @param object EGSEAResults object, the analysis result object from  \code{\link{egsea}}, 
#' \code{\link{egsea.cnt}}
#' or  \code{\link{egsea.ora}}. 
#' @param contrast contrast column number or column name specifying which 
#' contrast is of interest.
#' if contrast = 0 or "comparison" and the number of contrasts greater than 1, 
#' the comparative gene sets are 
#' retruned. 
#' @param gs.label the number or label of the gene set collection of interest.
#' @param sort.by character, determines how to order the analysis results in 
#' the stats table. 
#' The accepted values depend on the function used to generate the EGSEA 
#' results.
#' @param number integer, maximum number of gene sets to list
#' @param names.only logical, whether to display the EGSEA statistics or not. 
#' @param verbose logical, whether to print out progress messages and warnings. 
#' 
#' @export
#' @return 
#' A dataframe of top gene sets with the calculated statistics for each if 
#' names.only = FALSE.
#' 
#' @name topSets
#' @aliases topSets,EGSEAResults-method
#' 
#' @examples
#' library(EGSEAdata)
#' data(il13.data)
#' v = il13.data$voom
#' contrasts = il13.data$contra
#' gs.annots = buildIdxEZID(entrezIDs=rownames(v$E), species="human", 
#' msigdb.gsets="none", 
#'          kegg.updated=FALSE, kegg.exclude = c("Metabolism"))
#' gsa = egsea(voom.results=v, contrasts=contrasts,  gs.annots=gs.annots, 
#'          symbolsMap=v$genes, 
#' baseGSEAs=egsea.base()[-c(2,5,6,9,12)], display.top = 5,
#'           sort.by="avg.rank",  
#'          num.threads = 2, report=FALSE)
#' topSets(gsa, contrast=1, gs.label="kegg", number = 10)
#' topSets(gsa, contrast=1, gs.label=1, sort.by="ora", number = 10, 
#' names.only=FALSE)
#' topSets(gsa, contrast=0, gs.label="kegg", number = 10)
#' 

setGeneric(name="topSets",
        def = function(object, contrast=1, gs.label=1, sort.by=NULL, 
                number = 10, names.only=TRUE, verbose = TRUE){
            standardGeneric("topSets")
        }
)

setMethod(f = "topSets",
        signature(object = "EGSEAResults"),
        definition = function(object, contrast=1, gs.label=1, sort.by=NULL, 
                number = 10, names.only=TRUE, verbose = TRUE){
            tryCatch({              
                        if (is.numeric(contrast))
                            if (contrast != 0)
                                contrast = object@contrasts[contrast]
                            else
                                contrast = "comparison"
                        if (verbose)
                            cat(paste0("Extracting the top gene sets of the collection \n",
                                object@gs.annots[[gs.label]]$name, " for the contrast ",
                                contrast, "\n Sorted by ", 
                                ifelse(is.null(sort.by), object@sort.by, sort.by)
                                , "\n"))
                        if (tolower(contrast) == "comparison") 
                            top.gs = object@results[[gs.label]][["comparison"]][["test.results"]]
                        else
                            top.gs = object@results[[gs.label]][["test.results"]][[contrast]]
                        if (! is.null(sort.by)){
                            top.gs = top.gs[order(top.gs[,sort.by],
                                            decreasing=(sort.by == "Significance")), ]
                        }   
                        top.gs = cbind(Rank=seq(1, nrow(top.gs)), top.gs)
                        
                        number = ifelse(number <= nrow(top.gs), number, 
                                nrow(top.gs))
                        #print(names(top.gs))
                        top.gs = top.gs[1:number, ]
                        if (names.only)
                            return(rownames(top.gs))
                        else{
                            top.gs = as.data.frame(top.gs)
                            top.gs[, "Direction"] = 
                                    as.character(lapply(as.numeric(top.gs[, "Direction"]), 
                                                    function(x) if (x > 0) "Up" else if (x < 0) 
                                                            "Down" else "No Change"))
                            return(top.gs)
                        }
                    }, 
                    error = function(e){
                        cat(paste0("ERROR: topSets(...) encountered an error:\n", e ))
                    })
            return(NULL)
        }

)


#' @title Show the content of the EGSEAResults object
#' @description Show the content of the EGSEAResults object
#' 
#' @param object EGSEAResults object, the analysis result object from  \code{\link{egsea}}, 
#' \code{\link{egsea.cnt}}
#' or  \code{\link{egsea.ora}}. 
#' 
#' @export
#' @return 
#' NULL
#' 
#' @name show
#' @aliases show,EGSEAResults-method
#' 
#' @examples
#' library(EGSEAdata)
#' data(il13.data)
#' v = il13.data$voom
#' contrasts = il13.data$contra
#' gs.annots = buildIdxEZID(entrezIDs=rownames(v$E), species="human", 
#' msigdb.gsets="none", 
#'          kegg.updated=FALSE, kegg.exclude = c("Metabolism"))
#' gsa = egsea(voom.results=v, contrasts=contrasts,  gs.annots=gs.annots, 
#'          symbolsMap=v$genes, 
#' baseGSEAs=egsea.base()[-c(2,5,6,9,12)], display.top = 5,
#'           sort.by="avg.rank", num.threads = 2, report=FALSE)
#' show(gsa)

#setGeneric(name="show",
#        def = function(object){
#            standardGeneric("show")
#        }
#)

setMethod(f = "show",
        signature(object = "EGSEAResults"),
        definition = function(object){
            cat("An object of class \"EGSEAResults\"\n")
            cat(paste0("\tTotal number of genes: ", 
                            length(object@gs.annots[[1]]$featureIDs), "\n"))
            cat(paste0("\tTotal number of samples: ", object@sampleSize, "\n"))
            cat(paste0("\tContrasts: ", paste(object@contrasts, collapse=", "), "\n"))
            cat(paste0("\tBase GSE methods: ", 
                            paste(object@baseMethods, collapse=", "), "\n"))
            cat(paste0("\tP-values combining method: ", object@combineMethod, "\n"))
            cat(paste0("\tSorting statistic: ", object@sort.by, "\n"))
            cat(paste0("\tOrganism: ", object@gs.annots[[1]]$species, "\n"))
            cat(paste0("\tHTML report generated: ", ifelse(object@report, "Yes", "No"), "\n"))
            if (object@report)
                cat(paste0("\tHTML report directory: ", object@report.dir, "\n"))
            cat("\tTested gene set collections: \n")
            for (gs.annot in object@gs.annots){
                cat(paste0("\t\t", gs.annot$name, " (", gs.annot$label , "): ",
                                length(gs.annot$idx), " gene sets\n"))
            }
            cat("Use summary(object) and topSets(object, ...) to explore this object.\n")
        }

)


#' @title Summary of the EGSEAResults object
#' @description Summary of the EGSEAResults object
#' 
#' @param object EGSEAResults object, the analysis result object from  \code{\link{egsea}}, 
#' \code{\link{egsea.cnt}}
#' or  \code{\link{egsea.ora}}. 
#' 
#' @export
#' @return 
#' NULL
#' 
#' @name summary
#' @aliases summary,EGSEAResults-method
#' 
#' @examples
#' library(EGSEAdata)
#' data(il13.data)
#' v = il13.data$voom
#' contrasts = il13.data$contra
#' gs.annots = buildIdxEZID(entrezIDs=rownames(v$E), species="human", 
#' msigdb.gsets="none", 
#'          kegg.updated=FALSE, kegg.exclude = c("Metabolism"))
#' gsa = egsea(voom.results=v, contrasts=contrasts,  gs.annots=gs.annots, 
#'          symbolsMap=v$genes, 
#' baseGSEAs=egsea.base()[-c(2,5,6,9,12)], display.top = 5,
#'           sort.by="avg.rank", num.threads = 2, report=FALSE)
#' summary(gsa)

#setGeneric(name="summary",
#        def = function(object){
#            standardGeneric("summary")
#        }
#)

setMethod(f = "summary",
        signature(object = "EGSEAResults"),
        definition = function(object){
            # for each collection and contrast, print top 10 + comparison if exists            
            for (label in names(object@results)){
                cat(paste0("*Top 10 gene sets in the ", object@gs.annots[[label]]$name, 
                                " collection\n"))
                for (contrast in 1:length(object@contrasts)){     
                    cat(paste0("**Contrast ", object@contrasts[contrast], "\n"))
                    cat(paste(topSets(object, contrast, label, verbose=FALSE), collapse=" | "))
                    cat("\n")
                }
                if ("comparison" %in% names(object@results[[label]])){
                    cat(paste0("**Comparison analysis \n"))
                    cat(paste(topSets(object, "comparison", label, verbose=FALSE), collapse=" | "))
                    cat("\n")
                }
            }
        }

)

#' @title Plot a heatmap for a given gene set
#' @description Summary of the EGSEAResults object
#' 
#' @param object EGSEAResults object, the analysis result object from  \code{\link{egsea}}, 
#' \code{\link{egsea.cnt}}
#' or  \code{\link{egsea.ora}}. 
#' @param gene.set character, the name of the gene set. 
#' See the output of \code{\link{topSets}}.
#' @param gs.label character or numeric, the index/name of the gene set collection.
#' See names(object@@gs.annots) for valid values.
#' @param contrast character or numeric, the index/name of the contrast or 0/"comparison". 
#' See object@@contrasts for valid values.
#' @param file.name character, the name of the heatmap file without an extension.
#' @param format character, takes "pdf" or "png".
#' @param verbose logical, whether to print out progress messages and warnings.
#' 
#' @export
#' @return 
#' NULL
#' 
#' @name plotHeatmap
#' @aliases plotHeatmap,EGSEAResults-method
#' 
#' @examples
#' library(EGSEAdata)
#' data(il13.data)
#' v = il13.data$voom
#' contrasts = il13.data$contra
#' gs.annots = buildIdxEZID(entrezIDs=rownames(v$E), species="human", 
#' msigdb.gsets="none", 
#'          kegg.updated=FALSE, kegg.exclude = c("Metabolism"))
#' gsa = egsea(voom.results=v, contrasts=contrasts,  gs.annots=gs.annots, 
#'          symbolsMap=v$genes, 
#' baseGSEAs=egsea.base()[-c(2,5,6,9,12)], display.top = 5,
#'           sort.by="avg.rank", num.threads = 2, report=FALSE)
#' plotHeatmap(gsa, "Asthma")
#' plotHeatmap(gsa, "Asthma", contrast = "comparison", file.name = "asthma.hm.cmp")

setGeneric(name="plotHeatmap",
        def = function(object, gene.set, gs.label=1, contrast=1, file.name="heatmap", 
                format = "pdf", verbose=TRUE){
            standardGeneric("plotHeatmap")
        }
)

setMethod(f = "plotHeatmap",
        signature="EGSEAResults",
        definition = function(object, gene.set, gs.label=1, contrast=1, 
                file.name="heatmap", format = "pdf", verbose=TRUE){            
            tryCatch({            
                if (is.numeric(contrast))
                    if (contrast != 0)
                        contrast = object@contrasts[contrast]
                    else
                        contrast = "comparison"
                if (verbose)
                    cat(paste0("Plotting heatmap for ", gene.set, 
                        " from the collection \n",
                        object@gs.annots[[gs.label]]$name, " and for the contrast ",
                        contrast, "\n"))
                if (contrast %in% object@contrasts){
                    suppressWarnings(generateHeatMap(gene.set, object@gs.annots[[gs.label]], 
                            object@logFC[, contrast], 
                            object@symbolsMap, file.name,
                            format, print.csv = FALSE))
                }else if (tolower(contrast) == "comparison"){
                    suppressWarnings(generateHeatMap(gene.set, object@gs.annots[[gs.label]], 
                            object@logFC, 
                            object@symbolsMap, file.name,
                            format, print.csv = FALSE))
                }else{
                    stop("Unrecognized contrast value. 
						Use one of the object@contrasts or a numeric value.")
                }
            }, 
            error = function(e){
                cat(paste0("ERROR: plotHeatmap(...) encountered an error:\n", e ))
            })            
        }
)

#' @title Plot a pathway map for a given KEGG pathway
#' @description Summary of the EGSEAResults object
#' 
#' @param object EGSEAResults object, the analysis result object from  \code{\link{egsea}}, 
#' \code{\link{egsea.cnt}}
#' or  \code{\link{egsea.ora}}. 
#' @param gene.set character, the name of the pathway. 
#' See the output of \code{\link{topSets}}.
#' @param gs.label character or numeric, the index/name of the KEGG pathways collection.
#' See names(object@@gs.annots)[grep("^kegg", names(object@@gs.annots))] for valid values.
#' @param contrast character or numeric, the index/name of the contrast or 0/"comparison". 
#' See object@@contrasts for valid values.
#' @param file.name character, the name of the heatmap file without an extension.
#' @param verbose logical, whether to print out progress messages and warnings.
#' 
#' 
#' @export
#' @return 
#' NULL
#' 
#' @name plotPathway
#' @aliases plotPathway,EGSEAResults-method
#' 
#' @examples
#' library(EGSEAdata)
#' data(il13.data)
#' v = il13.data$voom
#' contrasts = il13.data$contra
#' gs.annots = buildIdxEZID(entrezIDs=rownames(v$E), species="human", 
#' msigdb.gsets="none", 
#'          kegg.updated=FALSE, kegg.exclude = c("Metabolism"))
#' gsa = egsea(voom.results=v, contrasts=contrasts,  gs.annots=gs.annots, 
#'          symbolsMap=v$genes, 
#' baseGSEAs=egsea.base()[-c(2,5,6,9,12)], display.top = 5,
#'           sort.by="avg.rank", num.threads = 2, report=FALSE)
#' plotPathway(gsa, "Asthma")
#' plotPathway(gsa, "Asthma", contrast="comparison", file.name = "asthma.map.cmp")

setGeneric(name="plotPathway",
        def = function(object, gene.set, gs.label=1, contrast=1, file.name="pathway"
                , verbose=TRUE){
            standardGeneric("plotPathway")
        }
)

setMethod(f = "plotPathway",
        signature="EGSEAResults",
        definition = function(object, gene.set, gs.label=1, contrast=1, 
                file.name="pathway", verbose=TRUE){            
            tryCatch({               
                if (is.numeric(gs.label))
                    gs.label = names(object@gs.annots)[gs.label]
                stopifnot(length(grep("^kegg", tolower(gs.label))) == 1)
                if (is.numeric(contrast))
                    if (contrast != 0)
                        contrast = object@contrasts[contrast]
                    else
                        contrast = "comparison"
                if (verbose)
                    cat(paste0("Plotting pathway map for ", gene.set, 
                        " from the collection \n",
                        object@gs.annots[[gs.label]]$name, " and for the contrast ",
                        contrast, "\n"))
                if (contrast %in% object@contrasts){
                    suppressWarnings(generatePathway(gene.set, object@gs.annots[[gs.label]], 
                                    object@logFC[, contrast], 
                                    file.name = file.name))
                }else if (tolower(contrast) == "comparison"){
                    suppressWarnings(generateComparisonPathway(gene.set, 
                                    object@gs.annots[[gs.label]], 
                                    object@logFC, 
                                    file.name = file.name))
                }else{
                    stop("Unrecognized contrast value. 
                                    Use one of the object@contrasts or a numeric value.")
                }
            }, 
            error = function(e){
                cat(paste0("ERROR: plotPathway(...) encountered an error:\n", e ))
            })            
        }
)

#' @title Plot a multi-dimensional scaling (MDS) plot for the gene set rankings
#' @description Plot a multi-dimensional scaling (MDS) plot for the gene set rankings
#' of the base GSE methods
#' 
#' @param object EGSEAResults object, the analysis result object from  \code{\link{egsea}}, 
#' \code{\link{egsea.cnt}}
#' or  \code{\link{egsea.ora}}. 
#' @param gs.label character or numeric, the index/name of the KEGG pathways collection.
#' See names(object@@gs.annots)[grep("^kegg", names(object@@gs.annots))] for valid values.
#' @param contrast character or numeric, the index/name of the contrast or 0/"comparison". 
#' See object@@contrasts for valid values.
#' @param file.name character, the name of the heatmap file without an extension.
#' @param format character, takes "pdf" or "png".
#' @param verbose logical, whether to print out progress messages and warnings. 
#' 
#' @export
#' @return 
#' NULL
#' 
#' @name plotMDS
#' @aliases plotMDS,EGSEAResults-method
#' 
#' @examples
#' library(EGSEAdata)
#' data(il13.data)
#' v = il13.data$voom
#' contrasts = il13.data$contra
#' gs.annots = buildIdxEZID(entrezIDs=rownames(v$E), species="human", 
#' msigdb.gsets="none", 
#'          kegg.updated=FALSE, kegg.exclude = c("Metabolism"))
#' gsa = egsea(voom.results=v, contrasts=contrasts,  gs.annots=gs.annots, 
#'          symbolsMap=v$genes, 
#' baseGSEAs=egsea.base()[-c(2,5,6,9,12)], display.top = 5,
#'           sort.by="avg.rank", num.threads = 2, report=FALSE)
#' plotMDS(gsa)

setGeneric(name="plotMDS",
            def = function(object, gs.label=1, contrast=1, 
                file.name="methods.mds", format = "pdf",
                verbose = TRUE){
            standardGeneric("plotMDS")
        }
)

setMethod(f = "plotMDS",
        signature="EGSEAResults",
        definition = function(object, gs.label=1, contrast=1, 
                file.name="methods.mds", format = "pdf",
                verbose = TRUE){                
            tryCatch({         
                    if (is.numeric(contrast))
                        if (contrast != 0)
                            contrast = object@contrasts[contrast]
                        else
                            contrast = "comparison"
                    if (tolower(contrast) == "comparison")
                        stop("plotMDS(...) is not supported for the comparison analysis")
                    if (length(object@baseMethods) < 2){
                        stop("plotMDS(...) requires at least two base methods.")
                    }
                    if (verbose)
                        cat(paste0("Generating MDS plot for the collection \n",
                            object@gs.annots[[gs.label]]$name, " and for the contrast ",
                            contrast, "\n"))
                    capture.output(generateMDSMethodsPlot(object@results[[gs.label]]
                            [["test.results"]][[contrast]], 
                            object@baseMethods, 
                            file.name, format))
                }, 
                error = function(e){
                    cat(paste0("ERROR: plotMDS(...) encountered an error:\n", e ))
                }
            )
            
        }
)


#' @title Generate a Summary plot for EGSEA analysis
#' @description Generate a Summary plot for EGSEA analysis
#' 
#' @param object EGSEAResults object, the analysis result object from  \code{\link{egsea}}, 
#' \code{\link{egsea.cnt}}
#' or  \code{\link{egsea.ora}}. 
#' @param gs.label character or numeric, the index/name of the KEGG pathways collection.
#' See names(object@@gs.annots)[grep("^kegg", names(object@@gs.annots))] for valid values.
#' @param contrast character or numeric, the index/name of the contrast. To compare two
#' contrasts, pass their indexes/names.
#' See object@@contrasts for valid values.
#' @param file.name character, the name of the heatmap file without an extension.
#' @param format character, takes "pdf" or "png".
#' @param x.axis character, the x-axis of the summary plot. All the 
#' values accepted by the 
#' \strong{sort.by} parameter can be used. Default x.axis="p.value".
#' @param x.cutoff numeric, cut-off threshold to filter the gene sets of 
#' the summary plots
#' based on the values of the \strong{x.axis}. Default 
#' x.cutoff=NULL.
#' @param verbose logical, whether to print out progress messages and warnings. 
#' 
#' @export
#' @return 
#' NULL
#' 
#' @name plotSummary
#' @aliases plotSummary,EGSEAResults-method
#' 
#' @examples
#' library(EGSEAdata)
#' data(il13.data)
#' v = il13.data$voom
#' contrasts = il13.data$contra
#' gs.annots = buildIdxEZID(entrezIDs=rownames(v$E), species="human", 
#' msigdb.gsets="none", 
#'          kegg.updated=FALSE, kegg.exclude = c("Metabolism"))
#' gsa = egsea(voom.results=v, contrasts=contrasts,  gs.annots=gs.annots, 
#'          symbolsMap=v$genes, 
#' baseGSEAs=egsea.base()[-c(2,5,6,9,12)], display.top = 5,
#'           sort.by="avg.rank", num.threads = 2, report=FALSE)
#' plotSummary(gsa)
#' plotSummary(gsa, contrast=c(1,2), file.name = "summary.cmp")
#' 
setGeneric(name="plotSummary",
            def = function(object, gs.label=1, contrast=1, 
                    file.name="summary", format = "pdf",
                    x.axis = "p.adj", x.cutoff = NULL,
                    verbose = TRUE){
                standardGeneric("plotSummary")
            }
)

setMethod(f = "plotSummary",
        signature="EGSEAResults",
        definition = function(object, gs.label=1, contrast=1, 
                file.name="summary", format = "pdf",
                x.axis = "p.adj", x.cutoff = NULL,
                verbose = TRUE){                
            tryCatch({         
                    if (is.numeric(contrast))                            
                        contrast = object@contrasts[contrast]                                  
                    if (is.null(x.cutoff)){
                        if (x.axis %in% c("p.value", "p.adj"))
                            x.cutoff = 1     
                        else
                            x.cutoff = 10000
                    }                        
                    if (length(contrast) == 2){
                        if (verbose)
                            cat(paste0("Generating Summary plots for the collection \n",
                            object@gs.annots[[gs.label]]$name, " and for the comparison ",
                            paste(contrast, collapse=" vs "), "\n"))                        
                        plot.data = generatePlotData.comparison(
                            object@results[[gs.label]][["test.results"]][contrast], 
                            object@results[[gs.label]][["comparison"]][["test.results"]], 
                            object@gs.annots[[gs.label]], 
                            x.axis)
                        if (x.axis %in% c("p.value", "p.adj")){
                            generateSummaryPlots(plot.data, file.name,
                                paste0("-log10(p-value) for ", contrast[1]),
                                paste0("-log10(p-value) for ", contrast[2]),
                                format = format)
                        }else{
                            generateSummaryPlots(plot.data, file.name,
                                paste0("1 / ", x.axis, " for ", contrast[1]),
                                paste0("1 / ", x.axis, " for ", contrast[2]),
                                format = format)
                        }
                    }else if (length(contrast) == 1){
                        if (verbose)
                            cat(paste0("Generating Summary plots for the collection \n",
                            object@gs.annots[[gs.label]]$name, " and for the contrast ",
                            contrast, "\n"))
                        plot.data = generatePlotData(
                                object@results[[gs.label]]
                                        [["test.results"]][[contrast]], 
                                object@gs.annots[[gs.label]], 
                                x.cutoff, x.axis)
                        if (x.axis %in% c("p.value", "p.adj"))
                            generateSummaryPlots(plot.data, file.name, 
                                    format = format)
                        else
                            generateSummaryPlots(plot.data, file.name, 
                                    Xlab = paste0("1 / ", x.axis),
                                    format = format)
                    }else{
                        stop("Wrong number of contrasts. Max is 2.")
                    }                 
                }, 
                error = function(e){
                    cat(paste0("ERROR: plotMDS(...) encountered an error:\n", e ))
                }
            )
            
        }
)

#' @title EGSEA Result Sorting Options
#' @description It lists the accepted sorting methods for analysis results
#' @return It returns a character vector of the accepted values for the sort.by 
#' argument in egsea
#' @export
#' @examples 
#' egsea.sort()

egsea.sort <-function(){
    return(c(c("p.value", "p.adj", "avg.rank", "med.rank", "min.rank", 
"min.pvalue", "vote.rank",
                            "Significance"),egsea.base()))
}


#' @title EGSEA P-value Combining Options
#' @description It lists the p-value combining methods
#' @return It returns a character vector of available methods for the 
#' combineMethod argument in egsea
#' @export
#' @examples 
#' egsea.combine()

egsea.combine <- function(){
    return(c("fisher", "wilkinson", "average", "logitp", 
            "sump", "sumz"))
}


#' @title EGSEA Base GSE Methods
#' @description It lists the supported GSEA methods
#' @return It returns a character vector of supported GSE methods. 
#' @export
#' @examples
#' egsea.base()

egsea.base <- function(){
    return(c("camera", "roast", "safe", "gage", "padog", "plage", "zscore", 
            "gsva", "ssgsea", "globaltest", "ora", "fry")) 
        # , "SPIA", "gsea", "samgs"
}






#TODO: use BiocParallel instead of parallel

#TODO: add to topSets() anno = T/F whether to add annotation info, e.g., ID, etc 

#TODO: showSet(gsa, gene.set, contrast = NULL or value, gs.label) display genes size, anno, etc.
# if contrast = NULL, no analysis results. Otherwise, display analysis results and
# return list of genes with fold changes
#TODO: 


# R CMD build --resave-data EGSEA 
# R CMD check EGSEA_0.99.0.tar.gz 
# R CMD BiocCheck EGSEA_0.99.0.tar.gz 
# R CMD INSTALL EGSEA_0.99.0.tar.gz 


# R CMD build --resave-data EGSEAdata
# R CMD check EGSEAdata_0.99.0.tar.gz 
# R CMD BiocCheck EGSEAdata_0.99.0.tar.gz 


## Create test units, type inside the EGSEA directory
# devtools::use_testthat()