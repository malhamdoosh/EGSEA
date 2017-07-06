# Wrapper function to run ORA on different contrasts


runora <- function(voom.results, contrast, gs.annot, 
        num.workers=4, verbose = TRUE){     
    # run hypergeomteric test and write out ranked 'gene sets' for each 
# 'contrast'
    # The p-value you want is the probability of getting 100 white balls in 
    # a sample of size 400 from an urn with 3000 white balls and 12000 
# black balls. 
    # phyper(100, 3000, 12000, 400)
    #
    if (is.matrix(contrast)){
        contr.names = colnames(contrast)
        contr.num = ncol(contrast)
    }else{
        contr.names = names(contrast)
        contr.num = length(contrast)
    }
      
    if (!is.null(voom.results$E)){
        vfit = lmFit(voom.results, voom.results$design)
        if (is.matrix(contrast)){
            vfit = contrasts.fit(vfit, contrast)
            coefs = 1:contr.num
        }else{
            coefs = contrast
        }
        vfit = eBayes(vfit)
        pvalue.cut=0.05
        logfc.cut=1     
        universe = gs.annot@featureIDs
    }
    else if (!is.null(voom.results$featureIDs)){ # ORA Analysis 
        universe = voom.results$featureIDs
    }
    else{
        # following topGO approach
        if (gs.annot@species == "Homo sapiens")
            EG.GO = AnnotationDbi::toTable(org.Hs.egGO2ALLEGS)
        else if (gs.annot@species == "Mus musculus")
            EG.GO = AnnotationDbi::toTable(org.Mm.egGO2ALLEGS)
        else if (gs.annot@species == "Rattus norvegicus")
            EG.GO = AnnotationDbi::toTable(org.Rn.egGO2ALLEGS)
        else
            stop("Unsupported species for ORA universe.")
        d = duplicated(EG.GO[,c("gene_id", "go_id", "Ontology")])
        EG.GO = EG.GO[!d, ]     
        universe = as.character(unique(EG.GO$gene_id))# fix the univese
    }
    
    ora.results = vector("list", contr.num)    
    for(i in 1:contr.num){     
        if (verbose)
            message("   Running ORA for ", contr.names[i])
        else
            message(".", appendLF = FALSE)
        if (!is.null(voom.results$E)){
            deGenes = rownames(topTable(vfit, coef=coefs[i], number=Inf, 
                        p.value=pvalue.cut, 
                        lfc=logfc.cut))     
            if (length(deGenes) == 0){
                deGenes = rownames(topTable(vfit, coef=i, number=Inf, 
                                p.value=pvalue.cut, 
                                lfc=0))
                if (verbose)
                    message("ORA used a cut-off logFC = 0")
            }
                
        }else{
            deGenes = voom.results$ids          
        }
        deGenes = deGenes[which(deGenes %in% universe)]
        ora.stats = runora.collection(deGenes, gs.annot, universe)
        ora.results[[i]] = ora.stats
#       print(head(ora.results[[i]]))
        ora.results[[i]] = 
ora.results[[i]][order(ora.results[[i]][,"p.value"]),]   # order by p.value     
        ora.results[[i]] = cbind(Rank=seq(1, nrow(ora.results[[i]])), 
ora.results[[i]])
        
    }
    names(ora.results) = colnames(contrast)
    return(ora.results)
}


runora.collection <- function(deGenes, gs.annot, universe){
    tmp = rep(NA, length(gs.annot@original))
    ora.stats = data.frame(p.value=tmp, p.adj = tmp) # , NumGenes=tmp
    totalDE = length(deGenes)
    n = length(universe) - totalDE
    for (j in 1:length(gs.annot@original)){
        gset = gs.annot@original[[j]]
        totalDEinS = length(intersect(gset, deGenes)) 
        totalSinUniverse = length(intersect(gset, universe))        
        ora.stats[j, "p.value"] = phyper(q = totalDEinS- 0.5, m= 
totalDE, n = n, 
                k = totalSinUniverse, lower.tail = FALSE)
        #ora.stats[j, "NumGenes"] = totalDEinS
    }
    ora.stats[, "p.adj"] = p.adjust(ora.stats[, "p.value"], method = "BH")
    row.names(ora.stats) = names(gs.annot@original)
    return ( ora.stats)
}




