# Wrapper function to run ORA on different contrasts


runora <- function(voom.results, contrast, gs.annot,  
        ranked.gs.dir="", output=TRUE,
        num.workers=4, verbose = TRUE){     
    # run hypergeomteric test and write out ranked 'gene sets' for each 
'contrast'
    # The p-value you want is the probability of getting 100 white balls in 
    # a sample of size 400 from an urn with 3000 white balls and 12000 
# black balls. 
    # phyper(100, 3000, 12000, 400)
    #

    file.name = paste0(ranked.gs.dir, "/ora-ranked-", gs.annot$label, 
"-gene-sets-", 
            sub(" - ", "-", colnames(contrast)), '.txt')    
    if (!is.null(voom.results$E)){
        vfit = lmFit(voom.results, voom.results$design)
        vfit = contrasts.fit(vfit, contrast)
        vfit = eBayes(vfit)
        pvalue.cut=0.05
        logfc.cut=1     
        universe = gs.annot$featureIDs
    }
    else if (!is.null(voom.results$featureIDs)){
        universe = voom.results$featureIDs
    }
    else{
        # following topGO approach
        if (gs.annot$species == "Homo sapiens")
            EG.GO = AnnotationDbi::toTable(org.Hs.egGO2ALLEGS)
        else if (gs.annot$species == "Mus musculus")
            EG.GO = AnnotationDbi::toTable(org.Mm.egGO2ALLEGS)
        else if (gs.annot$species == "Rattus norvegicus")
            EG.GO = AnnotationDbi::toTable(org.Rn.egGO2ALLEGS)
        else
            stop("Unsupported species for ORA universe.")
        d = duplicated(EG.GO[,c("gene_id", "go_id", "Ontology")])
        EG.GO = EG.GO[!d, ]     
        universe = as.character(unique(EG.GO$gene_id))# fix the univese
    }
    
    ora.results = vector("list", ncol(contrast))    
    for(i in 1:ncol(contrast)){     
        if (verbose)
            print(paste0("   Running ORA for ", colnames(contrast)[i]))
        else
            cat(".")
        if (!is.null(voom.results$E)){
            deGenes = rownames(topTable(vfit, coef=i, number=Inf, 
                        p.value=pvalue.cut, 
                        lfc=logfc.cut))     
            if (length(deGenes) == 0){
                deGenes = rownames(topTable(vfit, coef=i, number=Inf, 
                                p.value=pvalue.cut, 
                                lfc=0))
                if (verbose)
                    print("ORA used a cut-off logFC = 0")
            }
                
        }else{
            deGenes = voom.results$ids          
        }
        deGenes = deGenes[which(deGenes %in% universe)]
        ora.stats = runora.contra(deGenes, gs.annot, universe)
        ora.results[[i]] = ora.stats
#       print(head(ora.results[[i]]))
        ora.results[[i]] = 
ora.results[[i]][order(ora.results[[i]][,"p.adj"]),]        
        ora.results[[i]] = cbind(Rank=seq(1, nrow(ora.results[[i]])), 
ora.results[[i]])
        if (!output)
            next
        writeResultsToHTML(colnames(contrast)[i], ora.results[[i]], 
gs.annot, "ORA", file.name[i])
    }
    names(ora.results) = colnames(contrast)
    return(ora.results)
}


runora.contra <- function(deGenes, gs.annot, universe){
    tmp = rep(NA, length(gs.annot$original))
    ora.stats = data.frame(p.value=tmp, p.adj = tmp) # , NumGenes=tmp
    totalDE = length(deGenes)
    n = length(universe) - totalDE
    for (j in 1:length(gs.annot$original)){
        gset = gs.annot$original[[j]]
        totalDEinS = length(intersect(gset, deGenes)) 
        totalSinUniverse = length(intersect(gset, universe))        
        ora.stats[j, "p.value"] = phyper(q = totalDEinS- 0.5, m= 
totalDE, n = n, 
                k = totalSinUniverse, lower.tail = FALSE)
        #ora.stats[j, "NumGenes"] = totalDEinS
    }
    ora.stats[, "p.adj"] = p.adjust(ora.stats[, "p.value"], method = "BH")
    row.names(ora.stats) = names(gs.annot$original)
    return ( ora.stats)
}




