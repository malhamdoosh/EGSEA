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
#' @description This packages provides the implementation of the EGSEA 
#' algorithm and addition functions to help perform GSE analysis. 
#' @references 
#' Monther Alhamdoosh, Milica Ng, Nicholas J. Wilson, Julie M. Sheridan, Huy Huynh, 
#' Michael J. Wilson, Matthew E. Ritchie; Combining multiple tools outperforms 
#' individual methods in gene set enrichment analyses. Bioinformatics 2017; 
#' 33 (3): 414-424. doi: 10.1093/bioinformatics/btw623
#' 
#' @rdname EGSEA-package
NULL

#' @title Core functions to perform ensemble of gene set enrichment analysis (EGSEA)
#' 
#' @description \code{egsea} is the main function to carry out gene set enrichment 
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
#' algorithm proposed here is not limited to these twelve GSE methods and new 
#' GSE tests 
#' can be easily integrated into the framework. This function takes the voom 
#' object and 
#' the contrast matrix as parameters. 
#' The results of EGSEA can be seen using the \code{\link{topSets}} function. \cr 
#' \cr
#' EGSEA report is an interactive HTML report that is generated if \code{report=TRUE} to
#' enable a swift navigation through the results of an EGSEA analysis. The following pages 
#' are generated for each gene set collection and contrast/comparison: \cr
#' 1. Stats Table page shows the detailed statistics of the EGSEA analysis for the 
#' \code{display.top} gene sets. It shows the EGSEA scores, individual rankings and 
#' additional annotation for each gene set. Hyperlinks to the source of each gene set
#' can be seen in this table when they are available. The "Direction" column shows the regulation
#' direction of a gene set which is calculated based on the \code{logFC}, which is
#' either calculated from the limma differential expression analysis or provided by the user. 
#' The \code{logFC.cutoff} and \code{fdr.cutoff} are applied for this calculation. 
#' The calculations of the EGSEA 
#' scores can be seen in the references section. The method \code{topSets} can be used to
#' generate custom Stats Table. \cr
#' 2. Heatmaps page shows the heatmaps of the gene fold changes for the gene sets that are
#' presented in the Stats Table page. Red indicates up-regulation
#' while blue indicates down-regulation. Only genes that appear in the input expression/count 
#' matrix are visualized in the heat map. Gene names are coloured based on their 
#' statistical significance in the \code{limma} differential expression analysis. 
#' The "Interpret Results" link below each heat map allows the user to download the 
#' original heat map values along with additional statistics from \code{limma} DE analysis (
#' if available) so that they can be used to perform further analysis in R, e.g., customizing 
#' the heat map visualization. Additional heat maps can be generated and customized
#'  using the method \code{plotHeatmap}. \cr
#' 3. Summary Plots page  shows the methods ranking plot along with the summary plots of 
#' EGSEA analysis. The method plot uses multidimensional scaling (MDS) to visualize the 
#' ranking of individual methods on a given gene set collection. The summary plots are 
#' bubble plots that visualize the distribution of gene sets based on the EGSEA
#' Significance Score and another EGSEA score (default, p-value). 
#' Two summary plots are generated: ranking and directional plots. Each gene set is 
#' reprersented with a bubble which is coloured based on the EGSEA ranking (in ranking
#' plots ) or gene set regulation direction (in directional plots) and sized based on the 
#' gene set cardinality (in ranking plots) or EGSEA Significance score (in directional plots).
#' Since the EGSEA "Significance Score" is proportional to the p-value and the 
#' absolute fold changes, it could be useful to highlight gene sets that
#' have high Significance scores. The blue labels on the summary plot indicate 
#' gene sets that do not appear in the top 10 list of gene sets based on the "sort.by" 
#' argument (black labels) yet they appear in the top 5 list of gene sets based on 
#' the EGSEA "Significance Score". If two contrasts are provided, the rank is calculated 
#' based on the "comparison" analysis results and the "Significance Score" is calculated 
#' as the mean. If \code{sort.by = NULL}, the slot \code{sort.by} of the \code{object} 
#' is used to order gene sets.
#' The method \code{plotSummary} can be used to customize the Summary plots by 
#' changing the x-axis score
#' and filtering bubbles based on the values of the x-axis. The method \code{plotMethods} can be
#' used to generate Methods plots. \cr
#' 4. Pathways page shows the KEGG pathways for the gene sets that are presented in the
#' Stats Table of a KEGG gene set collection. The gene fold changes are overlaid on the 
#' pathway maps and coloured based on the gene regulation direction: blue for down-regulation
#' and red for up-regulation. The method \code{plotPathway} can be used to generate
#' additional pathway maps. Note that this page only appears if a KEGG gene set collection
#' is used in the EGSEA analysis. \cr
#' 5. Go Graphs page shows the Gene Ontology graphs for top 5 GO terms in each of 
#' three GO categories: Biological Processes (BP), Molecular Functions (MF), 
#' and Cellular Components (CC). Nodes are coloured based on the default \code{sort.by}
#' score where red indicates high significance and yellow indicates low significance. 
#' The method \code{plotGOGraph} can be used to customize GO graphs by 
#' changing the default sorting score and the number of significance nodes that can be
#' visualized. It is recommended that a small number of nodes is selected. Note that
#' this page only appears if a Gene Ontology gene set collection is used, i.e., for
#' the c5 collection from MSigDB or the gsdbgo collection from GeneSetDB. \cr
#' \cr
#' Finally, the "Interpret Results" hyperlink in the EGSEA report allows the user to download
#' the fold changes and limma analysis results and thus improve the interpretation of the results.
#' \cr 
#' Note that the running time of this function significantly increseas when 
#' \code{report = TRUE}. For example, the analysis in the example section below
#' was conducted on the $203$ signaling and disease KEGG pathways using a MacBook Pro 
#' machine that had a 2.8 GHz Intel Core i7 CPU and 16 GB of RAM. The execution time 
#' varied between 23.1 seconds (single thread) to 7.9 seconds (16 threads) when the HTML 
#' report generation was disabled. The execution time took 145.5 seconds when the report 
#' generation was enabled using 16 threads.
#' 
#' @param voom.results list, an EList object generated using the  
#' \code{\link[limma]{voom}} function, which it has at least 
#' three elements: \code{E} log2 normalized expression values,
#' \code{design} design matrix, \code{targets} sample information,
#' which must have a column named \strong{group}. If a matrix 
#' \code{weights} is provided, it will be used in \code{limma}-based
#' methods. Entrez Gene IDs should be used as row names.  
#' @param contrasts double, an N x L matrix indicates the contrasts of the 
#' linear model coefficients for which the test is required if the 
#' design matrix does not have an intercept. N is number 
#' of columns of the design matrix and L is number of contrasts. This
#' matrix should be based on the primary factor of the design matrix. 
#' If the design matrix includes an intercept, this parameter can be 
#' a vector of integers that specify the columns of the design matrix. 
#' If this parameter is \code{NULL}, all pairwise comparisons based on 
#' \code{group}  or \code{voom.results$targets$group} are created, 
#' assuming that \code{group} is the primary factor in the \code{design}
#' matrix. Likewise, all the coefficients of the primary factor are used
#' if the \code{design} matrix has an intercept. 
#' @param logFC double, an K x L matrix indicates the log2 fold change of each 
#' gene for each contrast. K is the number of genes included in the analysis. 
#' If logFC=NULL, the logFC values are 
#' estimated using the  \code{\link[limma]{ebayes}} for each contrast. 
#' For \code{egsea.ora}, it can be a matrix or vector of the same 
#' length of \code{entrezIDs}. If logFC=NULL, 1 is used as a default value. Then, the 
#' regulation direction in heatmaps and pathway maps is not indicative of the 
#' gene regulation direction.
#' @param gs.annots list, list of objects of class GSCollectionIndex. It is generated 
#' using one of these functions: 
#'  \code{\link{buildIdx}}, \code{\link{buildMSigDBIdx}}, 
#' \code{\link{buildKEGGIdx}},
#' \code{\link{buildGeneSetDBIdx}}, and \code{\link{buildCustomIdx}}. 
#' @param symbolsMap dataframe, an K x 2 matrix stores the gene symbol of each 
#' Entrez Gene ID. The first column must be the Entrez Gene IDs and the 
#' second column must be the Gene Symbols. It is used for the heatmap visualization. 
#' In \code{egsea} and \code{egsea.cnt}, the number of rows should match that 
#' of the  \strong{voom.results} and \strong{counts}, respectively. 
#' Default symbolsMap=NULL. 
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
#' For \code{egsea.ora}, it takes "p.value", "p.adj" or "Significance". 
#' @param report.dir character, directory into which the analysis results are 
#' written out. 
#' @param kegg.dir character, the directory of KEGG pathway data file (.xml) 
#' and image file (.png). 
#' Default kegg.dir=paste0(report.dir, "/kegg-dir/").
#' @param logFC.cutoff numeric, cut-off threshold of logFC and is used for 
#' the calculation of Sginificance Score 
#' and Regulation Direction. Default logFC.cutoff=0.
#' @param fdr.cutoff numeric, cut-off threshold of DE genes and is used
#' for the calculation of Significance Score and Regulation Direction. 
#' Default fdr.cutoff = 0.05. 
#' @param sum.plot.axis character, the x-axis of the summary plot. All the 
#' values accepted by the 
#' \strong{sort.by} parameter can be used. Default sum.plot.axis="p.value".
#' @param sum.plot.cutoff numeric, cut-off threshold to filter the gene sets of 
#' the summary plots
#' based on the values of the \strong{sum.plot.axis}. Default 
#' sum.plot.cutoff=NULL.
#' @param vote.bin.width numeric, the bin width of the vote ranking. Default 
#' vote.bin.width=5.
#' @param num.threads numeric, number of CPU cores to be used. Default 
#' num.threads=2.
#' @param report logical, whether to generate the EGSEA interactive report. It 
#' takes longer time
#' to run. Default is True.
#' @param interactive logical, whether to generate interactive tables and plots. 
#' Note this might dramatically increase the size of the EGSEA report.
#' @param keep.base logical, whether to write out the results of the 
#' individual GSE methods.
#'  Default FALSE.  
#' @param verbose logical, whether to print out progress messages and warnings. 
#' @param keep.limma logical, whether to store the results of the limma analysis
#' in the EGSEAResults object.
#' @param keep.set.scores logical, whether to calculate the gene set enrichment scores
#' per sample for the methods that support this option, i.e., "ssgsea".

#' 
#' @return \code{egsea} returns an object of the class EGSEAResults, which
#' stores the top gene sets and the detailed analysis
#' results for each contrast and the comparative analysis results.
#' 
#' @seealso \code{\link{topSets}}, \code{\link{egsea.base}}, 
#' \code{\link{egsea.sort}}, 
#' \code{\link{buildIdx}}, \code{\link{buildMSigDBIdx}}, 
#' \code{\link{buildKEGGIdx}},
#' \code{\link{buildGeneSetDBIdx}}, and \code{\link{buildCustomIdx}}
#'
#' @importFrom limma camera mroast fry lmFit  contrasts.fit eBayes topTable plotMDS
#' @importFrom globaltest gt gt.options result
#' @importFrom PADOG padog
#' @importFrom GSVA gsva
#' @importFrom safe safe safe.toptable getCmatrix
#' @importFrom parallel mclapply detectCores
#' @importFrom AnnotationDbi Ontology
#' @importFrom stats p.adjust pchisq phyper quantile
#' @importFrom grDevices pdf dev.off png colorRampPalette 
#' @importFrom graphics par legend
#' @importFrom utils write.table browseURL write.csv data capture.output packageVersion timestamp
#' @importFrom gage gage kegg.gsets
#' @importFrom metap logitp meanp sumlog sump sumz wilkinsonp
#' @importFrom methods new slot
#' @import hwriter HTMLUtils stringi ggplot2 pathview gplots  Biobase topGO
#' @export
#' @references 
#' Monther Alhamdoosh, Milica Ng, Nicholas J. Wilson, Julie M. Sheridan, Huy Huynh, 
#' Michael J. Wilson, Matthew E. Ritchie; Combining multiple tools outperforms 
#' individual methods in gene set enrichment analyses. Bioinformatics 2017; 
#' 33 (3): 414-424. doi: 10.1093/bioinformatics/btw623
#' 
#' @name egsea
#' @aliases egsea,egsea-main
#' @rdname egsea-main
#' 
#' @examples
#' # Example of egsea
#' library(EGSEAdata)
#' data(il13.data)
#' v = il13.data$voom
#' contrasts = il13.data$contra
#' gs.annots = buildIdx(entrezIDs=rownames(v$E), species="human", 
#' msigdb.gsets="none", 
#'          kegg.updated=FALSE, kegg.exclude = c("Metabolism"))
#' # set report = TRUE to generate the EGSEA interactive report
#' gsa = egsea(voom.results=v, contrasts=contrasts,  gs.annots=gs.annots, 
#'          symbolsMap=v$genes, baseGSEAs=egsea.base()[-c(2,5,6,9,12)], 
#' 			display.top = 5, sort.by="avg.rank", 
#' 			report.dir="./il13-egsea-report", 
#'          num.threads = 2, report = FALSE)
#' topSets(gsa) 
#' 


egsea <- function(voom.results, contrasts = NULL, logFC=NULL,
        gs.annots, symbolsMap=NULL, 
        baseGSEAs=egsea.base(),
        minSize=2, display.top=20, 
        combineMethod="wilkinson", combineWeights = NULL,      
        sort.by="p.adj", 
        report.dir=NULL, kegg.dir=NULL, 
        logFC.cutoff = 0, fdr.cutoff = 0.05, 
        sum.plot.axis="p.adj", sum.plot.cutoff=NULL, 
        vote.bin.width=5,
        num.threads=4, report = TRUE,
        interactive = FALSE,
        keep.base = FALSE, 
        verbose=FALSE,
        keep.limma=TRUE, keep.set.scores=FALSE
        ){

    set.seed(581986) # to guarantee reproducibility of results
    if (verbose)    
        return(egsea.main(voom.results, contrasts, gs.annots, baseGSEAs, 
            combineMethod, combineWeights,sort.by,  report.dir, 
            kegg.dir, logFC, symbolsMap, minSize, display.top, 
            logFC.cutoff, fdr.cutoff, sum.plot.cutoff, sum.plot.axis, 
            vote.bin.width, keep.base, verbose, num.threads, 
            report, keep.limma, keep.set.scores, interactive))
    else
        suppressWarnings(
            return(egsea.main(voom.results, contrasts, gs.annots, 
                baseGSEAs, combineMethod, combineWeights, sort.by, 
                report.dir, kegg.dir, logFC, symbolsMap, minSize, 
                display.top, logFC.cutoff, fdr.cutoff,
                sum.plot.cutoff, 
                sum.plot.axis, vote.bin.width,keep.base, 
                verbose, num.threads, report, keep.limma, 
                keep.set.scores, interactive)))    
}

#' @title Ensemble of Gene Set Enrichment Analyses Function
#' 
#' @description \code{egsea.cnt} is the main function to carry out gene set 
#' enrichment analysis using the EGSEA algorithm. This function is aimed to 
#' use raw RNASeq count matrix to perform the EGSEA analysis. 
#'
#' @param counts double, an K x M numeric matrix of read counts where genes are
#' the rows and samples are the columns. In this case, TMM normalization is used
#' to calculate the normalization factors. \code{counts} can be also a 
#' \code{DGEList} object. In this case, the normalization factors can be 
#' calculated by the user prior to invoking this method. This is REQUIRED. 
#' @param group character, vector or factor giving the experimental 
#' group/condition for each sample/library.  This is REQUIRED.
#' @param design double, an M x N numeric matrix giving the design matrix of 
#' the linear model fitting. N is the number of coefficients in model. 
#' If this parameter is 
#' \code{NULL}, \code{model.matrix(~0+group)} is used to create a deisgn matrix.
#' @inheritParams contrasts 
#' @inheritParams logFC 
#' @inheritParams gs.annots 
#' @inheritParams symbolsMap 
#' @inheritParams baseGSEAs 
#' @inheritParams minSize 
#' @inheritParams display.top 
#' @inheritParams combineMethod 
#' @inheritParams combineWeights 
#' @inheritParams sort.by 
#' @inheritParams report.dir 
#' @inheritParams kegg.dir 
#' @inheritParams logFC.cutoff 
#' @inheritParams fdr.cutoff 
#' @inheritParams sum.plot.axis 
#' @inheritParams sum.plot.cutoff 
#' @inheritParams vote.bin.width 
#' @inheritParams num.threads 
#' @inheritParams report 
#' @inheritParams interactive
#' @inheritParams keep.base 
#' @inheritParams verbose 
#' @inheritParams keep.limma 
#' @inheritParams keep.set.scores 
#' 
#' @return \code{egsea.cnt} returns an object of the class EGSEAResults, which
#' stores the top gene sets and the detailed analysis
#' results for each contrast and the comparative analysis results.
#' 
#'
#' @importFrom  edgeR calcNormFactors DGEList
#' @importFrom limma voom
#' @importFrom stats model.matrix
#' @export
#' 
#' @name egsea.cnt
#' @aliases egsea.cnt,egsea-main
#' @rdname egsea-main
#' 
#' @examples
#' # Example of egsea.cnt
#' library(EGSEAdata)
#' data(il13.data.cnt)
#' cnt = il13.data.cnt$counts
#' group = il13.data.cnt$group
#' design = il13.data.cnt$design
#' contrasts = il13.data.cnt$contra
#' genes = il13.data.cnt$genes
#' gs.annots = buildIdx(entrezIDs=rownames(cnt), species="human", 
#' msigdb.gsets="none",
#'          kegg.updated=FALSE, kegg.exclude = c("Metabolism"))
#' # set report = TRUE to generate the EGSEA interactive report
#' gsa = egsea.cnt(counts=cnt, group=group, design=design, contrasts=contrasts, 
#'          gs.annots=gs.annots, 
#'          symbolsMap=genes, baseGSEAs=egsea.base()[-c(2,5,6,9,12)], 
#' display.top = 5,
#'           sort.by="avg.rank", 
#' report.dir="./il13-egsea-cnt-report", 
#'          num.threads = 2, report = FALSE)
#' topSets(gsa) 
#' 

egsea.cnt <- function(counts, group, design = NULL, contrasts = NULL, logFC=NULL,
        gs.annots, symbolsMap=NULL, 
        baseGSEAs=egsea.base(),
        minSize=2, display.top=20, 
        combineMethod="wilkinson", combineWeights = NULL,      
        sort.by="p.adj", 
        report.dir=NULL, kegg.dir=NULL, 
        logFC.cutoff=0, fdr.cutoff = 0.05, sum.plot.axis="p.adj", sum.plot.cutoff=NULL, 
        vote.bin.width=5,
        num.threads=4, report = TRUE,
        interactive = FALSE,
        keep.base = FALSE, 
        verbose=FALSE, keep.limma=TRUE, keep.set.scores=FALSE){
    stopifnot(!is.null(counts))
    if (is.matrix(counts)){
        d = DGEList(counts, group=group)
        d = calcNormFactors(d, method="TMM")  
    }else{
        stopifnot(class(counts) == "DGEList")
        d = counts
    }
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
                    report.dir, kegg.dir, 
                    logFC.cutoff, fdr.cutoff, 
                    sum.plot.axis, sum.plot.cutoff, 
                    vote.bin.width,
                    num.threads, report, interactive,
                    keep.base , 
                    verbose, keep.limma, keep.set.scores
                    ))
    
}

#' @title Over-representation Analysis with EGSEA Reporting Capabilities 
#' 
#' @description \code{egsea.ora} is the main function to carry out 
#' over-representation analysis (ORA) on gene set collections using
#' a list of genes. 
#'
#' @details \code{egsea.ora} takes a list of gene IDs and uses the 
#' gene set collections from \pkg{EGSEAdata} or a custom-built collection to 
#' find over-represented gene sets in this list. It takes the advantage of 
#' the existing EGSEA reporting capabilities and generate an interative report 
#' for the ORA analysis. 
#' The results can be explored using the \code{\link{topSets}} function. 
#'
#' @param geneIDs character, a vector of Gene IDs to be tested for ORA. They 
#' must be Entrez IDs if \pkg{EGSEAdata} collections are used. 
#' @param universe character, a vector of Enterz IDs to be used as a background 
#' list. If universe=NULL, the background list is created from the 
#' \pkg{AnnotationDbi} package.
#' @inheritParams logFC 
#' @param title character, a short description of the experimental contrast. 
#' @inheritParams gs.annots 
#' @inheritParams symbolsMap 
#' @inheritParams minSize 
#' @inheritParams display.top 
#' @inheritParams sort.by 
#' @inheritParams report.dir 
#' @inheritParams kegg.dir 
#' @inheritParams sum.plot.axis 
#' @inheritParams sum.plot.cutoff 
#' @inheritParams num.threads
#' @inheritParams report 
#' @inheritParams interactive
#' @inheritParams verbose 
#' 
#' @return \code{egsea.ora} returns an object of the class EGSEAResults, which
#' stores the top gene sets and the detailed analysis
#' results.
#' 
#'
#' @importFrom org.Hs.eg.db org.Hs.egGO2ALLEGS
#' @importFrom org.Mm.eg.db org.Mm.egGO2ALLEGS
#' @importFrom org.Rn.eg.db org.Rn.egGO2ALLEGS
#' @export
#' 
#' @name egsea.ora
#' @aliases egsea.ora,egsea-main
#' @rdname egsea-main
#' 
#' @examples
#' # Example of egsea.ora
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
#' gs.annots = buildIdx(entrezIDs=deGenes, species="human", 
#' msigdb.gsets="none",
#'          kegg.updated=FALSE, kegg.exclude = c("Metabolism"))
#' # set report = TRUE to generate the EGSEA interactive report
#' gsa = egsea.ora(geneIDs=deGenes, universe= 
#' as.character(voom.results$genes[,1]),
#'              logFC =logFC, title="X24IL13-X24",  
#' gs.annots=gs.annots, 
#'              symbolsMap=top.Table[, c(1,2)], display.top = 5,
#'               report.dir="./il13-egsea-ora-report", num.threads = 2, 
#' 				report = FALSE)
#' topSets(gsa) 
#' 

egsea.ora <- function(geneIDs, universe=NULL, logFC=NULL, title=NULL, 
        gs.annots, symbolsMap=NULL, 
        minSize=2, display.top=20, sort.by = "p.adj",   
        report.dir=NULL, kegg.dir=NULL, 
        sum.plot.axis="p.adj", sum.plot.cutoff=NULL,
        num.threads=4, report = TRUE,    
        interactive = FALSE,
        verbose=FALSE){
    voom.results = list(ids=as.character(geneIDs))
    if (!is.null(universe))
        voom.results$featureIDs = universe
    if (is.null(logFC)){
        logFC = matrix(1, length(geneIDs), 1)
        rownames(logFC) = geneIDs     
    } else if (!is.matrix(logFC)){
        tmp = logFC
        logFC = matrix(1, length(geneIDs), 1)
        logFC[, 1] = tmp
        if (!is.null(names(tmp)))
            rownames(logFC) = as.character(names(tmp))      
        else
            rownames(logFC) = voom.results$ids
    }
    if (is.null(title))
        title = c("ExperimentalContrast")
    title = gsub("[[:punct:]]+", "", title)
    title = gsub("[[:space:]]+", "", title)
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
                    report.dir, kegg.dir, 
                    logFC.cutoff = 0, fdr.cutoff = 1, 
                    sum.plot.axis, sum.plot.cutoff, 
                    vote.bin.width = 5,
                    num.threads, report, interactive,
                    keep.base = FALSE, 
                    verbose))

}


#' @title Ensemble of Gene Set Enrichment Analyses Function
#' 
#' @description \code{egsea.ma} is the main function to carry out gene set 
#' enrichment analysis using the EGSEA algorithm. This function is aimed to 
#' use a microarray expression matrix to perform the EGSEA analysis. 
#'
#' @param expr double, an K x M numeric matrix of intensities where genes
#'  are the rows and samples are the columns. In this case, it is assumed
#' that the expression values are already normalized log2 intensities
#' and filtered, i.e. rows corresponding to control and low-quality 
#' probes removed. Row names of \code{expr} must match the first column in 
#' \code{probe.annot}.  This is REQUIRED.
#' @inheritParams group 
#' @param probe.annot double, an K x n numeric matrix where rows are 
#' probes, the first column 
#' contains probe IDs, the second column contains Entrez Gene IDs, 
#' and the third column (optional) contains the Gene Symbols. 
#' Symbols are used for the heatmap visualization
#' Additional annotation columns can be added for  
#' the probes, e.g., Chromosome, QualityScore, etc. These
#' additional columns are not used.  This is REQUIRED.
#' @param probeMap.method character, the method to be used in mapping the 
#' Probe IDs into Enterz Gene IDs. Accepted methods include: 
#' "avg", "med", "var", "sum" and "iqr".
#' EGSEA selects the probe with the highest
#' average, median, variance, sum or IQR of expression, respectively,  
#' as a representative for each expressed gene. 
#' @inheritParams design 
#' @inheritParams contrasts 
#' @inheritParams logFC 
#' @inheritParams gs.annots  
#' @inheritParams baseGSEAs 
#' @inheritParams minSize 
#' @inheritParams display.top 
#' @inheritParams combineMethod 
#' @inheritParams combineWeights 
#' @inheritParams sort.by 
#' @inheritParams report.dir 
#' @inheritParams kegg.dir 
#' @inheritParams logFC.cutoff 
#' @inheritParams fdr.cutoff 
#' @inheritParams sum.plot.axis 
#' @inheritParams sum.plot.cutoff 
#' @inheritParams vote.bin.width 
#' @inheritParams num.threads 
#' @inheritParams report 
#' @inheritParams interactive
#' @inheritParams keep.base 
#' @inheritParams verbose 
#' @inheritParams keep.limma 
#' @inheritParams keep.set.scores 
#' 
#' @return \code{egsea.ma} returns an object of the class EGSEAResults, which
#' stores the top gene sets and the detailed analysis
#' results for each contrast and the comparative analysis results.
#' 
#'
#' @importFrom limma neqc 
#' @importFrom stats model.matrix median var IQR
#' @importClassesFrom limma EListRaw EList
#' @import Glimma
#' @export
#' 
#' @name egsea.cnt
#' @aliases egsea.cnt,egsea-main
#' @rdname egsea-main
#' 
#' @examples
#' # Example of egsea.ma
#' library(Glimma)
#' data(arraydata)
#' expr = arraydata$arrays$E
#' group = as.factor(arraydata$targets$Condition)
#' levels(group) = c("DPcreEzh2", "DPev", "LumcreEzh2", "Lumev")
#' probe.annot = arraydata$arrays$genes[-1]
#' design <- model.matrix(~0+ group + as.factor(arraydata$targets$Experiment))
#' colnames(design)[1:4] = levels(group)
#' colnames(design)[5:6] = c("Exp2", "Exp3")
#' contr = makeContrasts("DPEzh2KOvsWT" = DPcreEzh2-DPev, 
#'                         "LumEzh2KOvsWT" = LumcreEzh2-Lumev,
#'                       levels=colnames(design))
#' 
#' gs.annots = buildIdx(entrezIDs=unique(probe.annot[, 2]),
#' 		 	species="mouse", 
#' 			msigdb.gsets="none",
#'          kegg.updated=FALSE, kegg.exclude = c("Metabolism")) 
#' 
#' # set report = TRUE to generate the EGSEA interactive report
#' gsa = egsea.ma(expr=expr, group=group, 
#' 			probe.annot = probe.annot,
#' 			design=design, contrasts=contr, 
#'          gs.annots=gs.annots, 
#'          baseGSEAs=egsea.base()[-c(2,5,6,9,12)], 
#' 			display.top = 5,
#'           sort.by="avg.rank", 
#' 			report.dir="./ezh2-egsea-ma-report", 
#'          num.threads = 2, report = FALSE)
#' topSets(gsa) 
#' 

egsea.ma <- function(expr, group, probe.annot, probeMap.method = "avg", 
        design = NULL, contrasts = NULL, logFC=NULL,
        gs.annots, 
        baseGSEAs=egsea.base(),
        minSize=2, display.top=20, 
        combineMethod="wilkinson", combineWeights = NULL,      
        sort.by="p.adj", 
        report.dir=NULL, kegg.dir=NULL, 
        logFC.cutoff=0, fdr.cutoff = 0.05, sum.plot.axis="p.adj", sum.plot.cutoff=NULL, 
        vote.bin.width=5,
        num.threads=4, report = TRUE,
        interactive = FALSE,
        keep.base = FALSE, 
        verbose=FALSE, keep.limma=TRUE, keep.set.scores=FALSE){
    stopifnot(!is.null(group))
    stopifnot(!is.null(probe.annot))  
    stopifnot(identical(rownames(expr), 
                    as.character(probe.annot[, 1])))
    probeMap.method = tolower(probeMap.method)
    stopifnot(probeMap.method %in% c("avg", "med", "var", "sum", "iqr"))
    if (is.null(design)){
        design = model.matrix(~0+group)
        colnames(design) = levels(factor(group))
    }
    if (is.matrix(expr)){
        elist = new("EList")
        elist$E = expr
    }else{
        stop("The 'expr' parameter must be a matrix")        
    }
    targets = data.frame(SampleID=colnames(expr), group=group)
    elist$design = design
    elist$targets = targets    
    elist$genes = probe.annot
    elist.filtered  = mapProbesIntoEntrezIDs(elist, probe.annot, probeMap.method)
    if (ncol(probe.annot) > 2){
        ezs = unique(probe.annot[, 2])
        symbolsMap = probe.annot[match(ezs, probe.annot[, 2]), c(2,3)]
    }else
        symbolsMap = NULL
    set.seed(581986) # to guarantee reproducibility of results
    if (verbose)    
        return(egsea.main(elist.filtered, contrasts, gs.annots, baseGSEAs, 
                        combineMethod, combineWeights,sort.by,  report.dir, 
                        kegg.dir, logFC, symbolsMap, minSize, display.top, 
                        logFC.cutoff, fdr.cutoff, sum.plot.cutoff, sum.plot.axis, 
                        vote.bin.width, keep.base, verbose, num.threads, 
                        report, keep.limma, keep.set.scores, interactive))
    else
        suppressWarnings(
                return(egsea.main(elist.filtered, contrasts, gs.annots, 
                                baseGSEAs, combineMethod, combineWeights, sort.by, 
                                report.dir, kegg.dir, logFC, symbolsMap, minSize, 
                                display.top, logFC.cutoff, fdr.cutoff,
                                sum.plot.cutoff, 
                                sum.plot.axis, vote.bin.width,keep.base, 
                                verbose, num.threads, report, keep.limma, 
                                keep.set.scores, interactive)))       
}

mapProbesIntoEntrezIDs <- function(elist, probe.annot, probeMap.method){
    ezs = unique(probe.annot[, 2])
    sel = c()
    for (ez in ezs){
        reps = which(probe.annot[, 2] == ez)
        if (length(reps) > 1){
            if (probeMap.method == "avg")
            # APPROACH 1: mean 
                exprs = rowMeans(elist$E[reps, ])
            else if (probeMap.method == "med")
            # APPROACH 2: median
                exprs = sapply(1:length(reps), 
                            function(x) median(elist$E[reps[x], ]))
            else if  (probeMap.method == "var")
            # APPROACH 3: variance
                exprs = sapply(1:length(reps),
                            function(x) var(elist$E[reps[x], ]))
            else if  (probeMap.method == "sum")
            # APPROACH 4: sum
                exprs = rowSums(elist$E[reps, ])
            else if (probeMap.method == "iqr")
            # APPROACH 5: IQR (used in Agilent's GeneSpring)
                exprs = sapply(1:length(reps), 
                    function(x) IQR(elist$E[reps[x], ]))            
            sel = c(sel, reps[which(max(exprs) == exprs)])
        }else
            sel = c(sel, reps)
    }
    elist.filtered = elist[sel, ]
    rownames(elist.filtered) = ezs
    return(elist.filtered)
    
}

#TODO: use BiocParallel instead of parallel

#TODO: add to topSets() anno = T/F whether to add annotation info, e.g., ID, etc 

#TODO: add support to use linear models in camera, roast and fry

#TODO: explore the use of ReportingTools for html generation
#TODO: Build a universal I/O interface for the wrappers


# rsync -av --exclude ".*" ../EGSEA .
# R-dev CMD build --resave-data EGSEA / R CMD build --resave-data EGSEA
# R-dev CMD check EGSEA_*.tar.gz   / R CMD check EGSEA_*.tar.gz 
# R-dev CMD BiocCheck EGSEA_*.tar.gz / R CMD BiocCheck EGSEA_*.tar.gz  
# R-dev CMD INSTALL EGSEA_*.tar.gz / R CMD INSTALL EGSEA_*.tar.gz  


# R-dev CMD build --resave-data EGSEAdata
# R-dev CMD check EGSEAdata_0.99.0.tar.gz 
# R-dev CMD BiocCheck EGSEAdata_0.99.0.tar.gz 


## Create test units, type inside the EGSEA directory
# devtools::use_testthat()