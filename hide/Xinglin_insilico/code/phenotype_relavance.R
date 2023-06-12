library(GSEABenchmarkeR)
geo2kegg <- loadEData("geo2kegg")
names(geo2kegg)
geo2kegg[[1]]
geo2kegg <- maPreproc(geo2kegg[1:5])
geo2kegg[[1]]
se <- geo2kegg[[1]]
table(se$GROUP)


geo2kegg <- runDE(geo2kegg, de.method="limma", padj.method="flexible")
rowData(geo2kegg[[1]], use.names=TRUE)

library(EnrichmentBrowser)
kegg.gs <- getGenesets(org="hsa", db="go")
kegg.ora.res <- runEA(geo2kegg[[1]], method="ora", gs=kegg.gs, perm=0)
kegg.ora.res

kegg.ora.res$ora$GSE1297
mala.kegg <- readRDS("/Library/Frameworks/R.framework/Versions/4.1/Resources/library/GSEABenchmarkeR/extdata/malacards/GO_BP.rds")
sapply(mala.kegg, nrow)

obs.score <- evalRelevance(kegg.ora.res$ora$GSE1297$ranking, mala.kegg$ALZ)
gs.names <- kegg.ora.res$ora$GSE1297$ranking$GENE.SET
gs.ids <- substring(gs.names, 1, 10)
rand.scores <- compRand(mala.kegg$ALZ, gs.ids, perm=1000)
summary(rand.scores)
(sum(rand.scores >= obs.score) + 1) / 1001
compOpt(mala.kegg$ALZ, gs.ids)
###note: GSE1297+limma+ORA+GO, obs phenotype relavence =  1814, opt PR = 3064, p = 0.000999

##CRFE
fnlist <- function(x, fil){z <- deparse(substitute(x))
    cat(z, "\t", file=fil)
    nams=names(x) 
    for (i in seq_along(x) ){cat(nams[i], "\t", x[[i]], "\n", 
                              file=fil, append=TRUE) }
}
fnlist(kegg.gs,"/Users/JasonJia/Desktop/CRFE/Xinglin_insilico/phenotype_relevance/anno_GO.txt")


kegg.gs.df = data.frame(names(kegg.gs), )

genes = as.data.frame(rowData(geo2kegg[[1]], use.names=TRUE))
genes = genes[order(genes$ADJ.PVAL),]
rownames(genes)
fileConn<-file("/Users/JasonJia/Desktop/CRFE/Xinglin_insilico/phenotype_relevance/GSE1297genes.txt")
writeLines(rownames(genes), fileConn)
close(fileConn)

out = data.frame(row.names(genes), genes$ADJ.PVAL)
write.table(out, "/Users/JasonJia/Desktop/CRFE/Xinglin_insilico/phenotype_relevance/GSE1297genes_copy.txt",sep ="\t", row.names = F, quote = T, col.names = F)
