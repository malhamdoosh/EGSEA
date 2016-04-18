# Wrapper function to run CAMERA on different contrasts

runcamera <- function(voom.results, contrast, gs.annot,  ranked.gs.dir="", 
output=TRUE,
        num.workers=4){     
    # run CAMERA and write out ranked 'gene sets' for each 'contrast'   
    file.name = paste0(ranked.gs.dir, "/camera-ranked-", gs.annot$label, 
"-gene-sets-", 
            sub(" - ", "-", colnames(contrast)), '.txt')        
    camera.results = vector("list", ncol(contrast)) 
    for(i in 1:ncol(contrast)){     
        print(paste0("   Running CAMERA for ", colnames(contrast)[i]))
        camera.results[[i]] = camera(y=voom.results, 
index=gs.annot$idx, 
                design=voom.results$design, 
                contrast=contrast[,i])
#       print(head(camera.results[[i]]))
        camera.results[[i]] = 
camera.results[[i]][order(camera.results[[i]][,"PValue"]),]     
        camera.results[[i]] = cbind(Rank=seq(1, 
nrow(camera.results[[i]])), camera.results[[i]])
        if (!output)
            next
        writeResultsToHTML(colnames(contrast)[i], camera.results[[i]], 
gs.annot, "CAMERA", file.name[i])
    }
    names(camera.results) = colnames(contrast)
    return(camera.results)
}
