###########################
## Update h2o
if ("package:h2o" %in% search()) { detach("package:h2o", unload=TRUE) }
if ("h2o" %in% rownames(installed.packages())) { remove.packages("h2o") }
pkgs <- c("RCurl","jsonlite")
for (pkg in pkgs) {
  if (! (pkg %in% rownames(installed.packages()))) { install.packages(pkg) }
}
install.packages("h2o", type="source", 
                 repos="http://h2o-release.s3.amazonaws.com/h2o/rel-zermelo/4/R")
###########################
library("h2o")
library("bigrquery")
library("tidyr")
library("dplyr")


## Set up Google Bigquery project name
project<- "tcga-masked"

## Downlaod clinicalData for masked CNV data
clinicalQuery <- paste(
  "SELECT
     cd.case_gdc_id, cd.project_short_name, cd.age_at_diagnosis, cd.race, cd.gender, cd.ethnicity
   FROM
   `isb-cgc.TCGA_bioclin_v0.Clinical` cd")
# clinicalData_masked is used to extract metadata, i.e. cancer_type, for masked CNV dataset
clinicalData_masked <- query_exec(clinicalQuery, project,max_pages=Inf, use_legacy_sql=FALSE)

# clinicalData is used to extract gender and cancer type data for unmasked CNVs
#clinicalData<-readRDS("clinicalData.rds")

##### Upstream/Downstream analysis
## Make the dataframe to store data
cnv_auc <- data.frame(matrix(ncol=6, nrow=20)) # Store auc for each CNV
rocPlot <- data.frame() # Store tpr/fpr values for ROC Plots
metCount <- data.frame() # Store methylation probe count data

# Declare DNA segment (from )
geneName <- "6_149661_254283:1 CNV Downstream" 
chr <- paste0('chr',sapply(strsplit(geneName,"_"),`[`,1))
start <- as.numeric(sapply(strsplit(geneName,"_"),`[`,2))
end <- as.numeric(sapply(strsplit(geneName,"_"),`[`,3))
plotName <- paste0(chr,":",start,"-",end)

# (1) Random sampling of the same size up/downstream + within CNV
# (2) percentage of CNV length instead of 1MB
cnv_len <- end-start
upstream_mid <- end+(end-start)
downstream_mid <- start-(end-start)
lb <- start-1e6
lb1 <- lb+(end-start)

ub <- end+1e6
ub1 <- ub-(end-start)

## Plot methylation probe distribution
# Extract probe IDs within 2MB of CNVs
billing <- "tcga-masked"
cpg_query <- paste("SELECT
                 cn.CpG_probe_id AS probe_ID,
                 cn.chromosome AS Chr,
                 cn.position AS Start_Pos,
                 FROM
                 `isb-cgc.platform_reference.GDC_hg38_methylation_annotation` cn
                 WHERE   		         	
                 (cn.chromosome = 'chr6' AND			
                   cn.position>",lb, "AND 	
                   cn.position<",ub,")")

tb <- bq_project_query(billing, cpg_query)
cpg_spread <- bq_table_download(tb)
cpg_spread$cpg_name <- geneName
metCount <- rbind(metCount,cpg_spread)

# Plot density of cpg probes
library(ggplot2)
ggplot(cpg_spread) + 
  geom_histogram(aes(x=Start_Pos),binwidth=cnv_len)+theme_bw()+
  annotate("rect", xmin = start, xmax = end, ymin = 0, ymax = Inf,
           alpha = .2)+annotate(geom="text", x=start, y=50, label=plotName)+
  xlab('Chromosome Location')+ggtitle(paste("Methylation Probe Count within 2Mb Window of",plotName))

## Extract methylation data for a specific gene or for a specific CNV segment
query <- paste("WITH PID as (SELECT
                 cn.CpG_probe_id AS probe_ID,
                 cn.chromosome AS Chr,
                 cn.position AS Start_Pos,
                 FROM
                 `isb-cgc.platform_reference.GDC_hg38_methylation_annotation` cn
                 WHERE   		         	
                 (cn.chromosome = 'chr6' AND			
                   cn.position>",start-cnv_len, "AND 	
                   cn.position<",start,"))
               
               SELECT
               met.probe_ID, met.beta_value, met.case_gdc_id, met.sample_barcode
               FROM `isb-cgc.TCGA_hg38_data_v0.DNA_Methylation_chr6` met
               JOIN PID ON (PID.probe_ID = met.probe_ID)")

# Get data from Google Bigquery
tb <- bq_project_query(billing, query)
gene <- bq_table_download(tb,page_size=1e6)

# Transform data from long to wide format (one probe/column)
gene <- gene[!gene$probe_ID=="",]
gene_wide <- gene%>% 
  pivot_wider(names_from=probe_ID,values_from=beta_value,
                                  values_fn=list(beta_value=mean))

# Get aliquot id to get clinical data
gene_wide$Tissue <- clinicalData_masked[match(gene_wide$case_gdc_id,
                                       clinicalData_masked$case_gdc_id),"project_short_name"]
met_dat <- as.data.frame(gene_wide)
rm(gene_wide)

# Extract sample type information from barcode: 11 = solid normal tissue normal, 
# 10 = blood derived normal, 01 = Solid tumor
met_dat$sampleType <- sapply(strsplit(sapply(strsplit(met_dat$sample_barcode,"-"), 
                                             `[`, 4),"[aA-zZ]+"),`[`,1)
met_dat$Tissue <- ifelse(met_dat$Tissue=="TCGA-COAD","COAD","Normal")

## Differentiate COAD from other cancer types
# If too many predictors: Pick n random beta values for model building
if (ncol(met_dat)>=50) {
  set.seed(123)
  random_bv <- sample(3:ncol(met_dat)-3,4)
  met_dat1 <- met_dat[,c(random_bv,which(colnames(met_dat)=="Tissue"))]
} else {
    met_dat1 <- met_dat
}


response <- "Tissue"
predictors <- setdiff(names(met_dat1), c(response,"case_gdc_id",
                      "sample_barcode","aliquot_id","sampleType"))
met_dat1[[response]] <- as.factor(met_dat1[[response]])  

## Differentiate COAD tumor tissue vs normal tissue
# coad_dat <- filter(met_dat,Tissue=="COAD")
# coad_dat$Tissue <- ifelse(coad_dat$sampleType==11,"Normal",coad_dat$Tissue)
# coad_dat1 <- coad_dat[,c(random_bv,which(colnames(coad_dat)=="Tissue"))]
# #coad_dat1 <- coad_dat
# response <- "Tissue"
# predictors_coad <- setdiff(names(met_dat1), c(response,"case_gdc_id",
#                            "sample_barcode","aliquot_id","sampleType"))
# coad_dat1[[response]] <- as.factor(coad_dat1[[response]])  

# Start h2o
h2o.init(nthreads=6,min_mem_size = "4g")

# CHANGE FOR EVERY CNV
j=7
cnv_auc[j,1] <- geneName
# Pick random number to decide which iteration to save tpr,fpr
plotNum <- sample(2:6,1)
for (i in 2:6){
  
  ## COAD vs. other cancer types
  # Load data into h2o
  train1.hex<- as.h2o(met_dat1, destination_frame = "train1.hex")
  model1 <- h2o.gbm(x=predictors,
                   y=response,
                   training_frame = train1.hex,
                   nfolds=10,
                   ntrees=100)
  # Record the AUC
  h2o.auc(model1, train=FALSE, xval=TRUE)
  cnv_auc[j,i] <- h2o.auc(model1, train=FALSE, xval=TRUE)
 
  ## COAD normal vs. tumor tissue
  # train_coad.hex<- as.h2o(coad_dat1, destination_frame = "train_coad.hex")        
  # model_coad <- h2o.gbm(x=predictors_coad,
  #                   y=response,
  #                   training_frame = train_coad.hex,
  #                   nfolds=10,
  #                   ntrees=100)
  # # Record the AUC
  # coad_auc[j,i] <- h2o.auc(model_coad, train=FALSE, xval=TRUE)
  
  if (i==plotNum){
    tpr <- h2o.performance(model1,xval=T)@metrics$thresholds_and_metric_scores["tpr"]
    fpr <- h2o.performance(model1,xval=T)@metrics$thresholds_and_metric_scores["fpr"]
    tmp <- cbind(tpr,fpr)
    tmp$chrLoc <- geneName
    rocPlot <- rbind(rocPlot,tmp)

    # tpr_coad <- h2o.performance(model_coad,xval=T)@metrics$thresholds_and_metric_scores["tpr"]
    # fpr_coad <- h2o.performance(model_coad,xval=T)@metrics$thresholds_and_metric_scores["fpr"]
    # tmp1 <- cbind(tpr_coad,fpr_coad)
    # tmp1$chrLoc <-  geneName
    # coad_rocPlot <- rbind(coad_rocPlot,tmp1)
  }
  gc()
  h2o:::.h2o.garbageCollect()
  h2o:::.h2o.garbageCollect()
  h2o:::.h2o.garbageCollect()
}

h2o.shutdown()
# Plot upstream/downstream regions of a CNV
ggplot(cnvROC1,aes(fpr,tpr,colour=chrLoc))+geom_line() + 
 geom_segment(aes(x=0,y=0,xend = 1, yend = 1),linetype = 2,col='grey')+
  theme_bw()


saveRDS(gene_auc,file="cnv_up_downstream_auc.rds")
saveRDS(cnvROC,file="cnv_up_downstream_rocPlot.rds")
saveRDS(coad_auc,file="coadOnly_cnv_up_downstream_auc.rds")
saveRDS(coad_rocPlot,file="coad_Only_cnv_up_downstream_rocPlot.rds")
h2o.shutdown(prompt=FALSE)

#############
## Processing for DMR-overlapped CNV
# cnv_query <- paste("#standardSQL
#      with topCNVs as
#                      (SELECT
#                      Chr,
#                      Start_Pos,
#                      End_Pos,
#                      COUNT(*) AS total
#                      FROM
#                      `tcga-unmasked.cnvMet.cnv_DMR` cn
#                      GROUP BY
#                      Chr,
#                      Start_Pos,
#                      End_Pos
#                      ORDER BY
#                      total DESC
#                      limit 50)
#         SELECT
#                      Aliquot_ID AS Aliquot_ID,
#                      cn.Chr AS Chr,
#                      cn.Start_Pos AS Start_Pos,
#                      cn.End_Pos AS End_Pos,
#                      cn.Num_Probes AS Num_Probes,
#                      avg(Segment_Mean) AS Segment_Mean,
#                      cn.Cancer_Type AS Cancer_Type
#                      FROM
#                      `tcga-unmasked.TCGA_unmasked_copyNum.TCGA_CNV_unmasked` cn
#                      JOIN
#                      topCNVs top ON (top.Chr=cn.Chr and top.Start_Pos=cn.Start_Pos and top.End_Pos=cn.End_Pos)
#                      group by Aliquot_ID,cn.Chr,cn.Start_Pos,cn.End_Pos,Num_Probes,Cancer_Type")
# 
# cnv_gene <- query_exec(cnv_query,project,max_pages=Inf, use_legacy_sql=FALSE)
# cnv_gene1 <- unite(cnv_gene, "chromosome_loc", c("Chr","Start_Pos","End_Pos"))
# cnv_gene_wide <- cnv_gene1 %>% pivot_wider(names_from = chromosome_loc,values_from=Segment_Mean)
# 
# 
# # Isolate COAD samples
# cnv_gene_wide$Cancer_Type <- ifelse(cnv_gene_wide$Cancer_Type=="TCGA-COAD","COAD","Normal")
# response <- "Cancer_Type"
# predictors <- setdiff(names(cnv_gene_wide), c(response, "Aliquot_ID","Num_Probes"))
# cnv_gene_wide[[response]] <- as.factor(cnv_gene_wide[[response]])  
# 
# # Start h2o
# h2o.init(nthreads=4,min_mem_size = "4g")
# 
# # Load data into h2o
# train.hex<- as.h2o(cnv_gene_wide, destination_frame = "train.hex")        
# 
# model <- h2o.gbm(x=predictors,
#                  y=response,
#                  training_frame = train.hex,
#                  nfolds=10,
#                  ntrees=100)
# h2o.auc(model,xval=TRUE)
# tpr <- h2o.performance(model,xval=T)@metrics$thresholds_and_metric_scores["tpr"]
# fpr <- h2o.performance(model,xval=T)@metrics$thresholds_and_metric_scores["fpr"]
# tmp <- cbind(tpr,fpr)
# ggplot(tmp,aes(fpr,tpr))+geom_line() + 
#   geom_segment(aes(x=0,y=0,xend = 1, yend = 1),
#                linetype = 2,col='grey')+theme_bw()+
#   ggtitle(geneIn)
# h2o.shutdown(prompt=FALSE)
