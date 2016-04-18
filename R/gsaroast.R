# Wrapper function to run ROAST on different contrasts

runroast <- function(voom.results, contrast, gs.annot,  ranked.gs.dir="", 
output = TRUE,
        num.workers=4){     
    # run ROAST and write out ranked 'gene sets' for each 'contrast'    

    file.name = paste0(ranked.gs.dir, "/roast-ranked-", gs.annot$label, 
"-gene-sets-", 
            sub(" - ", "-", colnames(contrast)), '.txt')        
    roast.results = vector("list", ncol(contrast))  
    for(i in 1:ncol(contrast)){
        print(paste0("   Running ROAST for ", colnames(contrast)[i]))
        roast.results[[i]] = mroast(y=voom.results, 
                index=gs.annot$idx[as.numeric(gs.annot$anno[, 
"NumGenes"]) < 2000], 
                design=voom.results$design, 
                contrast=contrast[,i], nrot=9999)
        # returns PropDown/PropUp ==> proportion of genes that are 
# down/up-regulated
        roast.results[[i]] = roast.results[[i]][order(roast.results[[i]][, 
                                "PValue"]),]
        roast.results[[i]] = cbind(Rank=seq(1, 
nrow(roast.results[[i]])), roast.results[[i]])
        if (!output)
            next
        writeResultsToHTML(colnames(contrast)[i], roast.results[[i]], 
gs.annot, "ROAST", file.name[i])
    }
    names(roast.results) = colnames(contrast)
    return(roast.results)
}