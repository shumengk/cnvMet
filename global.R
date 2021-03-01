geneAUC <- readRDS("gene_auc.rds")
geneROCPlot <- readRDS("gene_rocPlot.rds")
coadOnlyAUC <- readRDS("coadOnly_gene_auc.rds")
coadOnlyROC <- readRDS("coadOnly_gene_rocPlot.rds")
cnvROCPlot <- readRDS("cnv_DMR.rds")
cnvAUC <- readRDS("cnv_auc.rds")
cnvIndROC <- readRDS("cnv_rocPlot.rds")


geneList <- unique(geneAUC$gene)
cnvList <- sort(c(unique(cnvAUC$chrLoc),"DMR-overlapped"))
