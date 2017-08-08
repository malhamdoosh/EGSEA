#Ensemble of Gene Set Enrichment Analyses 
# 
# Author: Monther Alhamdoosh, E:m.hamdoosh@gmail.com


## test egsea.sort

sort.opts = c(c("p.value", "p.adj", "vote.rank", "avg.rank", "med.rank", 
                "min.pvalue", "min.rank",  "avg.logfc", "avg.logfc.dir", "direction",
                "significance"),egsea.base())
expect_identical(egsea.sort(), sort.opts)


# test egsea.base
expect_identical(egsea.base(), c("camera", "roast", "safe", "gage", "padog", 
                "plage", "zscore", "gsva", "ssgsea", 
                    "globaltest", "ora", "fry"))


# test egsea.combine
expect_identical(egsea.combine(), c("fisher", "wilkinson", "average", "logitp", 
                                    "sump", "sumz", "votep", "median"))


# test buildCustomIdxEZID
entrezIDs = c("2180", "2181",  "2182",  "2194",  "23205", "23305", "27349", "31",    "32",   
 "51703", "54995", "55301", "81616")
gs.anno = buildCustomIdx(geneIDs=entrezIDs, gsets=list("set1"=entrezIDs))
slots = slotNames(gs.anno)
expect_true("original" %in% slots)
expect_true("idx" %in% slots)
expect_equal(length(gs.anno@original), length(gs.anno@idx))
expect_equal(length(gs.anno@original[[1]]), length(gs.anno@idx[[1]]))
expect_identical(entrezIDs, gs.anno@featureIDs)