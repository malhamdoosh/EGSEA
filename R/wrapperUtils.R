# TODO: Add comment
# 
# Author: monther
###############################################################################

# function to extract logCPM matrices for methods that work with two groups only 

prepareTwoGroupsData <- function(voom.results, contrast, gs.annot,
        min.samples = NULL, verbose = FALSE){    
    if (!is.matrix(contrast)){
        # find the reference samples 
        if (is.null(voom.results$targets) || 
                is.null(voom.results$targets$group))
            stop(paste0("The data frame 'targets' of the object 'voom.results' ",
                         "must have a column named 'group'."))
        ref.group = levels(factor(voom.results$targets$group))[1]  
        if (verbose)
            message("   Reference group is identified as ", ref.group)
        cnt.sam.indx = which(voom.results$targets$group == ref.group)
    }
    # use gene indexes instead of gene IDs
    groupData = list()    
    gsets = list()
    for (j in 1:length(gs.annot@idx)){
        gsets[[j]] = as.character(gs.annot@idx[[j]])
    }
    names(gsets) = names(gs.annot@idx)
    groupData[["gsets"]] = gsets
    # Extract a logCPM matrix for each contrast
    data.log = voom.results$E
    rownames(data.log) = as.character(seq(1, nrow(data.log)))
    design = voom.results$design       
    sam.idx = 1:ncol(data.log)
    groupData[["data"]] = list()
    #set.seed(05081986)
    contr.num = ifelse(is.matrix(contrast), ncol(contrast), length(contrast))
    for(i in 1:contr.num){
        if (is.matrix(contrast)){
            # find the indexes of the treatment group samples 
            d = design[, contrast[,i] > 0]
            if (is.null(ncol(d))){
                tre.sam.indx = sam.idx[ d == 1]
            }else if (ncol(d) > 1){
                tre.sam.indx = c()
                for (j in 1:ncol(d))
                    tre.sam.indx = c(tre.sam.indx, sam.idx[ d[,j] 
                                            == 1])
            }
            else
                stop("Invalid contrasts selected.")                   
            # find the indexes of the control group samples 
            d = design[, contrast[,i] < 0]
            if (is.null(ncol(d))){
                cnt.sam.indx = sam.idx[ d == 1]
            }else if (ncol(d) > 1){
                cnt.sam.indx = c()
                for (j in 1:ncol(d))
                    cnt.sam.indx = c(cnt.sam.indx, sam.idx[ d[,j] 
                                            == 1])
            }
            else
                stop("Invalid contrasts selected.")     
        }else{
            tre.sam.indx = sam.idx[ design[, contrast[i]] == 1]
        }
        # Check if a minimum number of samples is required
        if (! is.null(min.samples)){
            if (length(tre.sam.indx) == 1)
                tre.sam.indx = rep(tre.sam.indx, min.samples)
            else if (length(tre.sam.indx) < min.samples)
                tre.sam.indx = c(tre.sam.indx, 
                        sample(tre.sam.indx, min.samples - length(tre.sam.indx)))
            
            if (length(cnt.sam.indx) == 1)
                cnt.sam.indx = rep(cnt.sam.indx, min.samples)
            else if (length(cnt.sam.indx) < min.samples)
                cnt.sam.indx = c(cnt.sam.indx, 
                        sample(cnt.sam.indx, min.samples - length(cnt.sam.indx)))
        }
        # logCPM matrix has control samples then treatment samples
        data.log.sel = data.log[, c(cnt.sam.indx, tre.sam.indx)]
        groupData$data[[i]] = list()
        groupData$data[[i]][["logCPM"]] = data.log.sel        
        # group1 is control / reference
        groupData$data[[i]][["group1"]] = seq(1,length(cnt.sam.indx))
        groupData$data[[i]][["group2"]] = seq(length(cnt.sam.indx) + 
                        1,ncol(data.log.sel))
    }
    names(groupData$data) = colnames(contrast)
    return(groupData)
}


getNumberofSamples <- function(voom.results, contrast){
    if (is.null(voom.results$design)){
        return(0)
    }
    if (is.matrix(contrast)){
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
    }else{
        return(nrow(voom.results$design))
    }    
}


