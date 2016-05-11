#Ensemble of Gene Set Enrichment Analyses 
# 
# Author: Monther Alhamdoosh, E:m.hamdoosh@gmail.com
###############################################################################

writeResultsToHTML <- function(contrast.name, gsa.results, gs.annot, method, 
file.name){
    title = paste0("Gene Set Enrichment Analysis on ", gs.annot$name, " 
using ", method, " (", 
            contrast.name, ")")
    path = strsplit(file.name, "/")[[1]]    
    file = path[length(path)]
    file = gsub(".txt", ".html", file)
    path = paste0(paste(path[1:length(path) -1], collapse = "/"), "/")
    if (!is.data.frame(gsa.results)){
        rows = row.names(gsa.results)
        gsa.results = data.frame(gsa.results)   
        gsa.results = lapply(gsa.results, as.numeric)
        gsa.results = as.data.frame(gsa.results)
        row.names(gsa.results) = rows
    }   
    annot.frame = gs.annot$anno[match(rownames(gsa.results), 
gs.annot$anno[,"GeneSet"]),]
    if ("BroadUrl" %in% colnames(annot.frame)){
        annot.frame[, "BroadUrl"] = hmakeTag("a", "Go to Broad", 
href=annot.frame[, "BroadUrl"]) 
    }
    if ("SourceURL" %in% colnames(annot.frame)){
        annot.frame[, "SourceDB"] = hmakeTag("a", annot.frame[, 
"SourceDB"], href=annot.frame[, "SourceURL"])
        annot.frame = annot.frame[, "SourceURL" != 
colnames(annot.frame)]
    }
    table.data = data.frame(annot.frame, gsa.results)   

    
    capture.output(HTMLsortedTable(table.data, title, title, file=file, path=path))
    
    annot.frame = gs.annot$anno[match(rownames(gsa.results), 
gs.annot$anno[,"GeneSet"]),]
    table.data = data.frame(Rank=as.numeric(seq(1, nrow(gsa.results))), 
            annot.frame,
            gsa.results)
    write.table(table.data,
            file=file.name, 
            sep="\t", 
            quote=FALSE, 
            row.names=FALSE)
}

writeEGSEAResultsToHTML <- function(contrast.name, gsa.results, gs.annot, 
file.name, comparison=FALSE){   
    
    title = paste0("Gene Set Enrichment Analysis on ", gs.annot$name, " 
using EGSEA (",
            contrast.name, ")") 
    path = strsplit(file.name, "/")[[1]]    
    file = path[length(path)]
    file = gsub(".txt", ".html", file)
    path = paste0(paste(path[1:length(path) -1], collapse = "/"), "/")
    if (!is.data.frame(gsa.results)){
        rows = row.names(gsa.results)
        gsa.results = data.frame(gsa.results)   
        gsa.results = lapply(gsa.results, as.numeric)
        gsa.results = as.data.frame(gsa.results)
        row.names(gsa.results) = rows
    }   
    ids = gs.annot$anno[match(rownames(gsa.results), gs.annot$anno[,2]), 
"ID"]
    if (!comparison){
        heatmaps = paste0("../hm-top-gs-", gs.annot$label, "/", 
                sub(" - ", "-", contrast.name), "/", ids, 
".heatmap.pdf")
        csvfiles = paste0("../hm-top-gs-", gs.annot$label, "/", 
                sub(" - ", "-", contrast.name), "/", ids, 
".heatmap.csv")
        
    }
    else{
        heatmaps = paste0("../hm-top-gs-", gs.annot$label, "/", ids, 
".heatmap.multi.pdf")       
        csvfiles = paste0("../hm-top-gs-", gs.annot$label, "/", ids, 
".heatmap.multi.csv")       
    }
    
    heatmaps = hmakeTag("a", "Show Map", href=heatmaps)     
    csvfiles = hmakeTag("a", "Download values", href=csvfiles)
    heatmaps = paste(heatmaps, csvfiles)
    
    if (length(grep("^kegg", gs.annot$label)) == 0){
        table.data = data.frame(Rank=as.numeric(seq(1, 
nrow(gsa.results))), 
                Heatmaps=heatmaps, 
                gs.annot$anno[match(rownames(gsa.results), 
gs.annot$anno[,"GeneSet"]),],
                gsa.results)
        
        write.table(table.data[, -c(which(colnames(table.data) == 
"Heatmaps"))],
                file=file.name, 
                sep="\t", 
                quote=FALSE, 
                row.names=FALSE)
        
    }
    else{
        if (!comparison){
            pathways = paste0("../pv-top-gs-", gs.annot$label, "/", 
                    sub(" - ", "-", contrast.name), "/", 
ids, ".pathview.png")           
        } else{
            pathways = paste0("../pv-top-gs-", gs.annot$label, "/",
                    ids, ".pathview.multi.png")
        }
        pathways = hmakeTag("a", "Show Pathway", href=pathways)
        pathways = paste(pathways, csvfiles)
        table.data = data.frame(Rank=as.numeric(seq(1, 
nrow(gsa.results))),
                Heatmaps=heatmaps,
                Pathways=pathways,
                gs.annot$anno[match(rownames(gsa.results), 
gs.annot$anno[,"GeneSet"]),],
                gsa.results)
        write.table(table.data[, -c(which(colnames(table.data) == 
"Heatmaps"), which(colnames(table.data) == "Pathways"))],
                file=file.name, 
                sep="\t", 
                quote=FALSE, 
                row.names=FALSE)
    }
    for (i in 1:ncol(table.data)){
        if (is.numeric(table.data[, i])){
            if (!anyOccur(colnames(table.data)[i], c("p.value", 
"p.adj", "min.pvalue")))
                table.data[, i] = round(table.data[,i], 2)
            else
                table.data[, i] = round(table.data[,i], 6)
        }
    }

    if ("BroadUrl" %in% colnames(table.data)){
        table.data[, "BroadUrl"] = hmakeTag("a", "Go to Broad", 
                href=table.data[, "BroadUrl"], 
'target'='_blank')  
    }
    if ("SourceURL" %in% colnames(table.data)){
        table.data[, "SourceDB"] = hmakeTag("a", table.data[, 
"SourceDB"], href=table.data[, "SourceURL"])
        table.data = table.data[, "SourceURL" != colnames(table.data)]
    }
    if ("GeneSet" %in% colnames(table.data)){
        table.data[, "GeneSet"] = gsub("_", " ", table.data[, 
"GeneSet"]) 
    }
    if ("Description" %in% colnames(table.data)){
        
        desc = c()
        for ( t in table.data[, "Description"]){
            l = stri_length(t)
            if (is.na(l) || l <= 50)
                desc = c(desc, t)
            else
                desc = c(desc, paste0("<span title='", t,"'>", 
substr(t, 1, 50), " ... </span>"))
        }
        table.data[, "Description"] = desc
    }
    if ("Direction" %in% colnames(table.data)){
   
        table.data[, "Direction"] = 
as.character(lapply(as.numeric(table.data[, "Direction"]), 
                function(x) if (x > 0) "Up" else if (x < 0) 
"Down" else "No Change"))

    }
    capture.output(HTMLsortedTable(table.data, title, title, file=file, path=path))
    
}

anyOccur <- function(string, list){
    occur = FALSE
    for (x in list){
        if (length(grep(x, string)) > 0){
            occur = TRUE
            break
        }           
    }
    return(occur)
}

generateSummaryPage.comparison <- function(contrast.names, gs.annot, file.name){
    path = strsplit(file.name, "/")[[1]]    
    file = path[length(path)]
    file = gsub(".txt", ".html", file)
    path = paste0(paste(path[1:length(path) -1], collapse = "/"), "/")
    
    title = paste0("Gene Set Enrichment Analysis on ", gs.annot$name, " 
using EGSEA (Comparison Analysis)")
    if (length(contrast.names) == 2 ){
        img.files = paste0("./", gs.annot$label, c("-summary.rank", 
"-summary.dir"),".png")
        plot.titles = c("Summary plot based on gene set rank and size", 
        "Summary plot based on regulation direction and significance")
    }
    else if (length(contrast.names) > 2){
        img.files = c()
        plot.titles = c()
        anchor.names = c()
        anchors.table = matrix("-", length(contrast.names), 
length(contrast.names))
        colnames(anchors.table) = contrast.names
        rownames(anchors.table) = contrast.names
        for (i in 1:(length(contrast.names)-1)){
            for (j in (i+1):length(contrast.names)){
                anchor.names = c(anchor.names, c(paste0('anch', 
i,j),paste0('anch', i,j)))
                anchors.table[i, j] = hmakeTag("a", "View", 
style="text-decoration: none",
                        href=paste0("#", 'anch', i,j))
                anchors.table[j, i] = hmakeTag("a", "View", 
style="text-decoration: none",
                        href=paste0("#", 'anch', i,j))
                img.files = c(img.files, paste0("./", 
gs.annot$label, 
                                paste0('-', 
i,j), c("-summary.rank", "-summary.dir"),
                                ".png"))
                titles = c("Summary plot based on gene set rank 
and size <br/>", 
                        "Summary plot based on 
regulation direction and significance <br/>")
                titles = paste0(titles, contrast.names[i]," | 
", contrast.names[j])
                plot.titles = c(plot.titles, titles)
                
            }
        }
    }
    else{
        img.files = c()
        plot.titles = c()
    }
    
    pdf.files = sub(".png", ".pdf", img.files)  
    images = hmakeTag("a", hmakeTag("img", src=img.files, width=500), 
href=pdf.files)
    pdfs = hmakeTag("a", "Download pdf file", href=pdf.files)
    content = paste(images, plot.titles, pdfs, sep="<br/>")
    if (length(contrast.names) > 2){
        anchors = hmakeTag("a", "", name=anchor.names)
        content = paste(anchors, content, sep="")
    }   
    if (length(content) %% 2 != 0){     
        content = c(content, "")# to make sure there are multiple of 2s
    }   
    content = matrix(content, ncol=2, byrow = TRUE) 

    p = openPage(file, dirname=path, title=title)
    hwrite(title, heading=1, br=TRUE, page=p)
    if (length(contrast.names) > 2){
        hwrite(anchors.table, align="center", border=1, width="1000px", 
                cellspacing=0, row.bgcolor='#ffffaa',
                row.style=list('font-weight:bold'),
                col.bgcolor='#ffffaa',
                col.style=list('font-weight:bold'),
                page=p)
    }
    hwrite(content, align="center", border=0, page=p)
    closePage(p)    
}

generateAllGOgraphsPage.comparison <- function(contrast.names, gs.annot, 
file.name){
    title = paste0("Gene Set Enrichment Analysis on ", gs.annot$name, " 
using EGSEA (Comparison Analysis)") 
    img.files = character(0)
    plot.titles = character(0)
    for (contrast.name in contrast.names){
        img.files = c(img.files, paste0("./", sub(" - ", "-", 
contrast.name), "-", gs.annot$label, "-top-", 
            c("BP", "MF", "CC"), ".png"))
        plot.titles = c(plot.titles, paste0(c("Top GO Biological 
Process (", "Top GO Molecular Functions (", 
        "Top GO Cellular Components ("), contrast.name, ")."))
    }
    pdf.files = sub(".png", ".pdf", img.files)
    
    images = hmakeTag("a", hmakeTag("img", src=img.files, width=600), 
href=pdf.files)
    pdfs = hmakeTag("a", "Download pdf file", href=pdf.files)
        
    content = paste(images, plot.titles, pdfs, sep="<br/>")
    if (length(content) %% 2 != 0){     
        content = c(content, "")# to make sure there are multiple of 2s
    }   
    content = matrix(content, ncol=2, byrow = TRUE)
    
    path = strsplit(file.name, "/")[[1]]    
    file = path[length(path)]
    file = gsub(".txt", ".html", file)
    path = paste0(paste(path[1:length(path) -1], collapse = "/"), "/")
    
    p = openPage(file, dirname=path, title=title)
    hwrite(title, heading=1, br=TRUE, page=p)
    hwrite(content, align="center", border=0, page=p)
    closePage(p)
}

generateAllGOgraphsPage <- function(contrast.name, gs.annot, file.name){
    title = paste0("Gene Set Enrichment Analysis on ", gs.annot$name, " 
using EGSEA (",
            contrast.name, ")") 
    
    path = strsplit(file.name, "/")[[1]]    
    file = path[length(path)]
    file = gsub(".txt", ".html", file)
    path = paste0(normalizePath(paste(path[1:length(path) -1], collapse = 
"/")), "/")
    
    tmp = c("BP", "MF", "CC")
    img.files = c()
    plot.titles = c()
    cat.names = c("Top GO Biological Process.", "Top GO Molecular 
Functions.", 
            "Top GO Cellular Components.")
    names(cat.names) = tmp
    for (cat in tmp){
        f= paste0(sub(" - ", "-", contrast.name), "-", gs.annot$label, 
"-top-", 
                cat, ".png")
        if (file.exists(paste0(path, f))){
            img.files = c(img.files, f)
            plot.titles = c(plot.titles, cat.names[cat])
        }
    }
    if (length(img.files) > 0){
        img.files = paste0("./", img.files)
        pdf.files = sub(".png", ".pdf", img.files)
        
        images = hmakeTag("a", hmakeTag("img", src=img.files, 
width=600), href=pdf.files)
        pdfs = hmakeTag("a", "Download pdf file", href=pdf.files)
        
        content = paste(images, plot.titles, pdfs, sep="<br/>")
        if (length(content) %% 2 != 0){     
            content = c(content, "")# to make sure there are 
# multiple of 2s
        }   
        content = matrix(content, ncol=2, byrow = TRUE)
    }
    else
        content = hmakeTag("h3", "GO graphs could not be generated!")
    p = openPage(file, dirname=path, title=title)
    hwrite(title, heading=1, br=TRUE, page=p)
    hwrite(content, align="center", border=0, page=p)
    closePage(p)
}

generateSummaryPage <- function(contrast.name, gs.annot, file.name){
    title = paste0("Gene Set Enrichment Analysis on ", gs.annot$name, " 
using EGSEA (",
            contrast.name, ")") 
    img.files = paste0("./", sub(" - ", "-", contrast.name), "-", 
gs.annot$label, "-summary.rank.png")
    img.files = c(img.files, paste0("./", sub(" - ", "-", contrast.name), 
"-", gs.annot$label, "-summary.dir.png"))
    path = strsplit(file.name, "/")[[1]]    
    path = paste0(paste(path[1:length(path) -1], collapse = "/"), "/")
    mds.file = paste0(sub(" - ", "-", contrast.name), "-", gs.annot$label, 
"-methods.png")
    if (file.exists(paste0(path, mds.file)))
        img.files = c(img.files, paste0("./", mds.file))
    pdf.files = sub(".png", ".pdf", img.files)
    
    images = hmakeTag("a", hmakeTag("img", src=img.files, width=600), 
href=pdf.files)
    pdfs = hmakeTag("a", "Download pdf file", href=pdf.files)
    plot.titles = c("Summary plot based on gene set rank and size", 
            "Summary plot based on regulation direction and significance")
    if (file.exists(paste0(path, mds.file)))
        plot.titles = c(plot.titles,
                "MDS plot for the gene set ranking in different base methods.")
    
    content = paste(images, plot.titles, pdfs, sep="<br/>")
    if (length(content) %% 2 != 0){     
        content = c(content, "")# to make sure there are multiple of 2s
    }   
    content = matrix(content, ncol=2, byrow = TRUE)
    
    path = strsplit(file.name, "/")[[1]]    
    file = path[length(path)]
    file = gsub(".txt", ".html", file)
    path = paste0(paste(path[1:length(path) -1], collapse = "/"), "/")
    
    p = openPage(file, dirname=path, title=title)
    hwrite(title, heading=1, br=TRUE, page=p)
    hwrite(content, align="center", border=0, page=p)
    closePage(p)    
}

generateAllHeatmapsPage <- function(contrast.name, gsa.results, gs.annot, 
file.name, comparison=FALSE){
    title = paste0("Gene Set Enrichment Analysis on ", gs.annot$name, " 
using EGSEA (",
            contrast.name, ")") 
    ids = gs.annot$anno[match(rownames(gsa.results), gs.annot$anno[,2]), 
"ID"]
    names = gs.annot$anno[match(rownames(gsa.results), gs.annot$anno[,2]), 
"GeneSet"]
    if(!comparison){
        heatmaps.img = paste0("./", sub(" - ", "-", contrast.name), 
"/", ids, ".heatmap.png")
        heatmaps = paste0("./", sub(" - ", "-", contrast.name), "/", 
ids, ".heatmap.pdf")
        csvfiles = paste0("./", sub(" - ", "-", contrast.name), "/", 
ids, ".heatmap.csv")
    } else{
        heatmaps.img = paste0("./", ids, ".heatmap.multi.png")
        heatmaps = paste0("./", ids, ".heatmap.multi.pdf")
        csvfiles = paste0("./", ids, ".heatmap.multi.csv")
    }
    heatmaps.img = hmakeTag("a", hmakeTag("img", src=heatmaps.img, 
width=300, center=TRUE), href=heatmaps) 
    heatmaps = hmakeTag("a", "Large Map", href=heatmaps)    
    csvfiles = hmakeTag("a", "Download values", href=csvfiles)
    heatmaps = paste(heatmaps, csvfiles)
    heatmaps = paste(heatmaps.img, ids, substr(names, 1, 30), heatmaps, 
sep="<br />")   
#   head(heatmaps)
    if (length(heatmaps) %% 5 != 0)
        heatmaps = c(heatmaps, rep("", 5 - (length(heatmaps) %% 5 )))
    hm.table = matrix(heatmaps, ncol=5, byrow = TRUE)
            
            
    path = strsplit(file.name, "/")[[1]]    
    file = path[length(path)]
    file = gsub(".txt", ".html", file)
    path = paste0(paste(path[1:length(path) -1], collapse = "/"), "/")
    
    
    p = openPage(file, dirname=path, title=title)
    hwrite(title, heading=1, br=TRUE, page=p)
    hwrite(hm.table, col.width = 200, page=p)
    closePage(p)        
    
}

generateAllPathwaysPage <- function(contrast.name, gsa.results, gs.annot, 
file.name, comparison=FALSE){
    title = paste0("Gene Set Enrichment Analysis on ", gs.annot$name, " 
using EGSEA (",
            contrast.name, ")") 
    ids = gs.annot$anno[match(rownames(gsa.results), gs.annot$anno[,2]), 
"ID"]
    names = gs.annot$anno[match(rownames(gsa.results), gs.annot$anno[,2]), 
"GeneSet"]
    if (!comparison){
        pathways.img = paste0("./", sub(" - ", "-", contrast.name), 
"/", ids, ".pathview.png")
        csvfiles = paste0("../hm-top-gs-",gs.annot$label, "/", sub(" - 
", "-", contrast.name), "/", ids, ".heatmap.csv")
        pathways = paste0("./", sub(" - ", "-", contrast.name), "/", 
ids, ".pathview.png")
    } else{
        pathways.img = paste0("./", ids, ".pathview.multi.png")
        csvfiles = paste0("../hm-top-gs-",gs.annot$label, "/", ids, 
".heatmap.multi.csv")
        pathways = paste0("./", ids, ".pathview.multi.png")
    }   
    pathways.img = hmakeTag("a", hmakeTag("img", src=pathways.img, 
width=300, center=TRUE), href=pathways)     
    csvfiles = hmakeTag("a", "Download values", href=csvfiles)  
    pathways = hmakeTag("a", "Large Pathway", href=pathways)
    pathways = paste(pathways, csvfiles)
    pathways = paste(pathways.img,ids, substr(names, 1, 30), pathways, 
sep="<br />")   
#   head(heatmaps)
    if (length(pathways) %% 5 != 0)
        pathways = c(pathways, rep("", 5 - (length(pathways) %% 5 )))
    hm.table = matrix(pathways, ncol=5, byrow = TRUE)
    
    
    path = strsplit(file.name, "/")[[1]]    
    file = path[length(path)]
    file = gsub(".txt", ".html", file)
    path = paste0(paste(path[1:length(path) -1], collapse = "/"), "/")
    
    
    p = openPage(file, dirname=path, title=title)
    hwrite(title, heading=1, br=TRUE, page=p)
    hwrite(hm.table, col.width = 200, page=p)
    closePage(p)
}

getNumberofSamples <- function(voom.results, contrast){
    if (is.null(voom.results$design)){
        return("N/A")
    }
    samples = c()
    sam.idx = colnames(voom.results$E)
    for(i in 1:ncol(contrast)){
        d = voom.results$design[, contrast[,i] != 0]
        if (is.null(ncol(d))){
            samples = sam.idx[ d == 1]
        }else if (ncol(d) > 1){         
            for (j in 1:ncol(d))
                samples = c(samples, sam.idx[ d[,j] == 1])
        }
        else
            stop("Invalid contrasts selected.")
    }
    return(length(unique(samples)))
}

generateEGSEAReport <- function(voom.results, contrast, gs.annots, baseGSEAs, 
combineMethod, 
        sort.by,  egsea.dir, 
        kegg.dir, logFC.cal, symbolsMap){
    gs.labels = sapply(gs.annots, function(x) x$label)
    gs.names =  sapply(gs.annots, function(x) x$name)
    ranked.gs.dir = paste0("./ranked-gene-sets-", combineMethod)
    top.gs.dir = paste0("./top-gene-sets-", combineMethod)
    pv.dir = paste0("./pv-top-gs-", gs.labels[grep("^kegg", gs.labels)], 
"/")
    hm.dir = paste0("./hm-top-gs-", gs.labels, "/")
    summary.dir = "./summary/"
    go.dir = "./go-graphs/"
    
    
    p = openPage("index.html", dirname=egsea.dir, title="Ensemble of Gene 
Set Enrichment Analyses (EGSEA) - Report")
    hwrite("Ensemble of Gene Set Enrichment Analyses (EGSEA) - Report", 
heading=1, br=TRUE, page=p) 
    hwrite("Analysis Parameters", heading=2, br=TRUE, page=p)
    #### write analysis parameters
    hwrite(paste0(hmakeTag("b", "Total number of genes: ") ,
                    length(gs.annots[[1]]$featureIDs)), 
br=TRUE, page=p)
    hwrite(paste0(hmakeTag("b", "Total number of samples: " ),
                    getNumberofSamples(voom.results, 
contrast)), br=TRUE, page=p)
    hwrite(paste0(hmakeTag("b", "Number of contrasts: ") ,ncol(contrast)), 
br=TRUE, page=p)
    hwrite(paste0(hmakeTag("b","Base GSEA methods: ") ,paste(baseGSEAs, 
collapse=",")), br=TRUE, page=p)
    hwrite(paste0(hmakeTag("b","P-value combine method: " ), 
combineMethod), br=TRUE, page=p)
    hwrite(paste0(hmakeTag("b","Sorting statistic: " ),sort.by), br=TRUE, 
page=p)
    hwrite(paste0(hmakeTag("b","Fold changes calculated: " ),logFC.cal), 
br=TRUE, page=p)    
    hwrite(paste0(hmakeTag("b","Gen IDs - Symbols mapping used: " ), 
                    ifelse(!is.null(symbolsMap), "Yes", 
"No")), br=TRUE, page=p)
    hwrite(paste0(hmakeTag("b","Organism: " ), gs.annots[[1]]$species), 
br=TRUE, page=p)
    
    
    hwrite("Analysis Results", heading=2, br=TRUE, page=p)
    main = ""
    ## Main Page Generation 
    for (i in 1:ncol(contrast)){        
        
        file.name = paste0(ranked.gs.dir, "/ranked-", gs.labels, 
"-gene-sets-", 
                sub(" - ", "-", colnames(contrast)[i]), 
'.html')        
        temp = paste0(gs.names, " (", hmakeTag("a", "Stats Table", 
href=file.name))
            
        file.name = paste0(hm.dir,  
                sub(" - ", "-", colnames(contrast)[i]), 
'-allHeatmaps.html')
        temp = paste0(temp, ", ", hmakeTag("a", "Heatmaps" , 
href=file.name))
                    
        kegg.idx =  grep("^kegg", gs.labels)
        if (length(kegg.idx) != 0){
            file.name = paste0(pv.dir, 
                    sub(" - ", "-", colnames(contrast)[i]), 
'-allPathways.html')    
            temp[kegg.idx] = paste0(temp[kegg.idx], ", ", 
hmakeTag("a", "Pathways" , href=file.name))
        }
        
        go.idx = which(gs.labels == "c5")
        if (length(go.idx) !=0 ){
            file.name = paste0(go.dir, 
                    sub(" - ", "-", colnames(contrast)[i]), 
'-allGOgraphs.html')    
            temp[go.idx] = paste0(temp[go.idx], ", ", hmakeTag("a", 
"GO Graphs" , href=file.name))
        }
        
        file.name = paste0(summary.dir, sub(" - ", "-", 
colnames(contrast)[i]), "-", gs.labels, 
                "-summary.html")
        temp = paste0(temp, ", ", hmakeTag("a", "Summary Plots", 
href=file.name))
        
        file.name = paste0(ranked.gs.dir, "/ranked-", gs.labels, 
"-gene-sets-", 
                sub(" - ", "-", colnames(contrast)[i]), '.txt')
        temp = paste0(temp, ", ", hmakeTag("a", "Download 
Stats",target="_blank", href=file.name))
            
        
        temp = paste0(temp, ")")
        
        temp = hmakeTag("il", paste(hmakeTag("b", 
colnames(contrast)[i]), hmakeTag("ul", paste(hmakeTag("li", temp), 
collapse="\n")), sep="\n"))
        main = paste(main, temp, sep="\n")
    }
    
    if (ncol(contrast) > 1){
        file.name = paste0(ranked.gs.dir, "/ranked-", gs.labels, 
"-gene-sets-compare.html")      
        
        temp = paste0(gs.names, " (", hmakeTag("a", "Stats Table", 
href=file.name))
        
        file.name = paste0(hm.dir, 'allHeatmaps.html')
        temp = paste0(temp, ", ", hmakeTag("a", "Heatmaps" , 
href=file.name))
        
        kegg.idx =  grep("^kegg", gs.labels)
        if (length(kegg.idx) != 0){
            file.name = paste0(pv.dir, 'allPathways.html')
            temp[kegg.idx] = paste0(temp[kegg.idx], ", ", 
hmakeTag("a", "Pathways" , href=file.name))
        }
        
        go.idx = which(gs.labels == "c5")
        if (length(go.idx) !=0 ){
            file.name = paste0(go.dir, 'allGOgraphs.html')  
            temp[go.idx] = paste0(temp[go.idx], ", ", hmakeTag("a", 
"GO Graphs" , href=file.name))
        }
        
        file.name = paste0(summary.dir, gs.labels, "-summary.html")
        temp = paste0(temp, ", ", hmakeTag("a", "Summary Plots", 
href=file.name))
        
        file.name = paste0(ranked.gs.dir, "/ranked-", gs.labels, 
"-gene-sets-compare.txt")
        temp = paste0(temp, ", ", hmakeTag("a", "Download Stats" , 
target="_blank", href=file.name))
        
        temp = paste0(temp, ")")
        
        temp = hmakeTag("il", paste(hmakeTag("b", "Comparison 
Analysis"), hmakeTag("ul", paste(hmakeTag("li", temp), collapse="\n")), 
sep="\n"))
        main = paste(main, temp, sep="\n")
    }

    hwrite(hmakeTag("ul", main, style="list-style-type:square"), page=p)
    hwrite(hmakeTag("footer", "Report generated by EGSEA package. For any 
inquiry, please contact the authors."))
    closePage(p)
    
}