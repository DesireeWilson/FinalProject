#loading libraries:
library(tidyverse)
library(ggforce)
library(GEOquery)
library(gt)
library(gtsummary)
library(pheatmap)
library(limma)
library(illuminaHumanv4.db)
library(survminer)
library(caret)
library(randomForest)
library(gbm)
library(kknn)
library(glmnet)
library(mice)
library(prophet)
library(reshape2)


#setting working directory:
thepath <- "C:/Users/wilso/OneDrive/Documents/GitHub/FinalProject/GSE70769/"
setwd(thepath)


#downloading the data:
gse70769 <- getGEO("GSE70769", GSEMatrix=TRUE)


#dimensions of the pheno data
dim(pData(phenoData(gse70769[[1]]))) # 94 rows and 51 columns; so 94 samples


#selecting columns:
pheno <- pData(phenoData(gse70769[[1]]))[, c(41, 42, 45, 47, 48, 51)]
rename_cols <- c("bcrRelapse", "clinicalStage", "pathologicalStage", "PSAdx", "timeToRecurrence", "tumorGleason")
colnames(pheno) <- rename_cols
#pheno_rename <- pheno %>% gt()

#selecting columns:
pheno_bcr <- pheno %>% dplyr::select(bcrRelapse, clinicalStage, tumorGleason)

# summarize the data with our package
pheno_summary <- tbl_summary(pheno_bcr)

#saving to a file:
gt::gtsave(as_gt(pheno_summary), file = "Graphs/SummaryTable.png")

#taking a peak at meta-data accompanying this dataset:
head(gse70769)
# $GSE70769_series_matrix.txt.gz
# ExpressionSet (storageMode: lockedEnvironment)
# assayData: 47323 features, 94 samples 
# element names: exprs 
# protocolData: none
# phenoData
# sampleNames: GSM1817906 GSM1817907 ... GSM1817999 (94
#                                                    total)
# varLabels: title geo_accession ... tumour gleason:ch1 (51
#                                                        total)
# varMetadata: labelDescription
# featureData
# featureNames: ILMN_1343291 ILMN_1343295 ... ILMN_3311190
# (47323 total)
# fvarLabels: ID Species ... GB_ACC (30 total)
# fvarMetadata: Column Description labelDescription
# experimentData: use 'experimentData(object)'
# pubMedIds: 26501111 
# Annotation: GPL10558 

#okay, I need to figure out the data type for gse70769:
typeof(gse70769) #okay, it's a list

#now I need to determine how many elements are in the list:
length(gse70769) #okay, so one

#now, I want to determine the data type of this element:
typeof(gse70769[[1]]) #S4, not helpful; let's try class()
class(gse70769[[1]]) 
# [1] "ExpressionSet"
# attr(,"package")
# [1] "Biobase"
#SWEET!

#obtaining gene expression matrix from gse70769:
exprsdat <- Biobase::exprs(gse70769[[1]])


#I suspect this data is already normalized (which is great)
#I will perform a density plot just in case:
#To do this, first I need to transform the data from a
#wide format to a long format:
pdf(file = "Graphs/DensityPlot.pdf"
    , width = 11
    , height = 8.5)
last_idx <- length(colnames(exprsdat))
exprsdat_long <- gather(data = data.frame(exprsdat)
                    , key = "sampleName"
                    , value = "log2exprs"
                    , colnames(exprsdat)[1]:colnames(exprsdat)[last_idx]
                    , factor_key = TRUE)
ggplot(exprsdat_long, aes(x=log2exprs, color=sampleName)) + geom_density()
dev.off()
#yes, the data is already normalized :)

#okay, so now I am going to do a correlation matrix
#followed by a heatmap:
pdf(file = "Graphs/CorrelationMatrixHeatmap.pdf"
    , width = 11
    , height = 8.5)
exprsdat_cor <- cor(x = exprsdat
                    , use = "pairwise.complete.obs")
pheatmap(as.matrix(exprsdat_cor))
dev.off()

#So, it looks like there are two samples that may be outliers.
#I'm going to :
exprsdat_cormtest <- cor.mtest(exprsdat)
all((rowSums(exprsdat_cormtest$p) > 0) == FALSE) # all true so no outliers
#so from here, I will perform differential gene expression



#before I perform the differential gene expression analysis, I want to see
#if there is a difference in biochemical recurrent rate between Gleason
#7=3+4 vs 7=4+3 prostate cancer cases:
pheno_table_gleason7 <- pheno %>% dplyr::select(bcrRelapse, tumorGleason) %>% 
  filter((tumorGleason == "7=3+4" | tumorGleason == "7=4+3") & (bcrRelapse != "N/A")) %>%
  group_by(bcrRelapse, tumorGleason) %>%
  summarize(n = n())

pheno_pivot_table_gleason7 <- dcast(data = pheno_table_gleason7
                           ,formula = bcrRelapse ~ tumorGleason
                           ,fun.aggregate=sum
                           ,value.var="n")

#performing chi-square on the values:
chisq.test(x = pheno_pivot_table_gleason7[,-1])
#Pearson's Chi-squared test with Yates' continuity correction
#
#data:  pheno_pivot_table_gleason7[, -1]
#X-squared = 1.0486, df = 1, p-value = 0.3058
#WooHoo! It's not significant!!!! So Gleason is probably NOT a confounding factor


#now I want to create a pheno table with just the gleason7 samples:
pheno_gleason7 <- pheno %>% dplyr::select(bcrRelapse, tumorGleason) %>% 
  filter((tumorGleason == "7=3+4" | tumorGleason == "7=4+3") & (bcrRelapse != "N/A"))

#selecting only the samples with Gleason 7:
selectedSamples <- rownames(pheno_gleason7)
selectedSamples_idx <- which(colnames(exprsdat) %in% selectedSamples)
exprsdat_gleason7 <- exprsdat[,selectedSamples_idx]
selectedSamples == colnames(exprsdat_gleason7)

#now perform differential expression analysis:
gleason7 <- factor(pheno_gleason7$tumorGleason, levels = c("7=3+4", "7=4+3"))
design <- model.matrix(~0+gleason7)
colnames(design) <- c("gleason3plus4", "gleason4plus3")
fit <- lmFit(exprsdat_gleason7, design)
contrast.matrix <- makeContrasts("gleason3plus4 - gleason4plus3", levels=design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)
topTable(fit2, coef=1, adjust="BH")
#there are two probes that are significant, ILMN_2337740  and ILMN_1720965 
#both probes measure the expression of TULP4 gene
#I need to write the code to prove this but I just did a Google search
#TULP4 is involved in protein ubiquitination. I dn't recall seeing any
#studies involved in protein ubiquitination AND prostate cancer. Interesting
#angle.

#I need to do a boxplot of this gene in 3+4 vs 4+3 prostate cancer cases.
#I need to transform data from wide to long
twoProbes <- c("ILMN_2337740", "ILMN_1720965")
twoProbes_idx <- which(rownames(exprsdat_gleason7) %in% twoProbes)
exprsdat_gleason7_twoProbes <- exprsdat_gleason7[twoProbes_idx,]
exprsdat_gleason7_twoProbes[,1:3]
exprsdat_gleason7_twoProbes_long <- gather(data.frame(exprsdat_gleason7_twoProbes),
                                           sampleID, logExprs)
gleason7cat <- pheno_gleason7$tumorGleason
gleason7cat <- gsub("7=3+4", "g3plus4", gleason7cat, fixed = TRUE)
gleason7cat <- gsub("7=4+3", "g4plus3", gleason7cat, fixed = TRUE)

probeID <- rep(rownames(exprsdat_gleason7_twoProbes), times = ncol(exprsdat_gleason7_twoProbes))
condition <- rep(gleason7cat, times = ncol(exprsdat_gleason7_twoProbes))
exprsdat_gleason7_twoProbes_long_probeID <- data.frame(probeID, 
                                                       exprsdat_gleason7_twoProbes_long$sampleID,
                                                       gleason7cat,
                                                       exprsdat_gleason7_twoProbes_long$logExprs)
colnames(exprsdat_gleason7_twoProbes_long_probeID) <- c("probeID", "sampleID", "gleason7cat", "logExprs")
head(exprsdat_gleason7_twoProbes_long_probeID)
pdf(file = "Graphs/ViolinPlotTwoTULP4probes.pdf"
    , width = 11
    , height = 8.5)
ggplot(exprsdat_gleason7_twoProbes_long_probeID,
       aes(probeID, logExprs, fill = gleason7cat, color = gleason7cat)) +
  geom_violin() +
  ggtitle("Log Expression of Two Probes Associated with TULP4 Gene")
dev.off()
#I don't think the probes are worth invesigating for survival analysis
#now, I want to start looking at the genes from my dissertation. 
#first, I have to link the probes to the gene symbols:
probe_ids <- rownames(exprsdat_gleason7)
# The below function call will return a datafram with probe_id, gene symbol
# and ŕefgene_id for your data
probeToSymbol_df <- select(illuminaHumanv4.db, 
       keys = probe_ids, 
       columns=c("SYMBOL"), 
       keytype="PROBEID")
genesOfInterest <- c("MYBPC1", "SLC7A4", "CDC42EP5", "BCAS1", "PAK1IP1", 
                     "ANPEP", "SLC23A1", "HS.270778", "GLB1L2", "GCNT2",
                     "TNFRSF19", "CUX2", "POTEG", "NCAPD3", "LOC728606", 
                     "SLC22A3", "C9ORF61", "SRD5A2", "FLJ31568", "SERPINF2", 
                     "PGM5", "LOC645993", "ZDHHC8P", "CNTNAP2", "F12", "CXCL14",
                     "CTHRC1", "ZNF467", "CDC20")
genes_idx <- which(probeToSymbol_df$SYMBOL %in% genesOfInterest)
exprsdat_gleason7_idx <- which(rownames(exprsdat_gleason7) %in% probeToSymbol_df$PROBEID[genes_idx])
exprsdat_gleason7_probesOfInterest <- exprsdat_gleason7[exprsdat_gleason7_idx,]
rownames(exprsdat_gleason7_probesOfInterest) == probeToSymbol_df$PROBEID[genes_idx]
#all are true! :)
#I'm going to rename the rownames and merge probe ID to gene symbol together:
probeID_geneSymbol <- paste0(rownames(exprsdat_gleason7_probesOfInterest), 
                             ":", probeToSymbol_df$SYMBOL[genes_idx])
exprsdat_gleason7_probesOfInterest_probeIDSymbol_merge <- exprsdat_gleason7_probesOfInterest
rownames(exprsdat_gleason7_probesOfInterest_probeIDSymbol_merge) <- probeID_geneSymbol

#going to perform differential gene analysis but focusing on the p-value, not adjusted p-value
fit <- lmFit(exprsdat_gleason7_probesOfInterest_probeIDSymbol_merge, design)
contrast.matrix <- makeContrasts("gleason3plus4 - gleason4plus3", levels=design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)
topTable(fit2, coef=1, adjust="BH")
#okay so only SRD5A2 and PGM5 are differentially expressed between 4+3 vs 3+4
#the unfortunate part is the log fold change is like 0.26 and 0.17 respectively
#I'll take a look at these probes in the survival analysis.
sample_idx <- which(rownames(pheno) %in% rownames(pheno_gleason7))
pheno_survival <- pheno[sample_idx,]
pheno_survival <- data.frame(sampleID, pheno_survival)
#now, I need to add the probe expression to the pheno data as well:
pheno_survival_withProbes <- data.frame(pheno_survival, t(exprsdat_gleason7_twoProbes), t(exprsdat_gleason7_probesOfInterest_probeIDSymbol_merge))
colnames(pheno_survival_withProbes)[8:9] <- c("ILMN_1720965.TULP4", "ILMN_2337740.TULP4")
pheno_survival_withProbes_tibble <- as_tibble(pheno_survival_withProbes)
pheno_survival_withProbes_tibble$tumorGleason <-
  case_match(
  pheno_survival_withProbes_tibble$tumorGleason,
  "7=3+4" ~ "g3+4",
  "7=4+3" ~ "g4+3",
)
pheno_survival_withProbes_tibble$bcrRelapse <-
  case_match(
    pheno_survival_withProbes_tibble$bcrRelapse,
    "Y" ~ TRUE,
    "N" ~ FALSE,
  )
#before I do the survival analysis, let me first create violin plots:
exprsdat_gleason7_allProbes_long <- gather(data.frame(t(pheno_survival_withProbes[ ,c(8:53)])),
                                           sampleID, logExprs)
gleason7cat <- pheno_survival_withProbes$tumorGleason
gleason7cat <- gsub("7=3+4", "g3plus4", gleason7cat, fixed = TRUE)
gleason7cat <- gsub("7=4+3", "g4plus3", gleason7cat, fixed = TRUE)
gleason7cat <- rep(gleason7cat, times = ncol(pheno_survival_withProbes[ ,c(8:53)]))
probeID <- rep(colnames(pheno_survival_withProbes)[8:53], times = nrow(pheno_survival_withProbes))
exprsdat_gleason7_allProbes_long_probeID <- data.frame(probeID, 
                                                       exprsdat_gleason7_allProbes_long$sampleID,
                                                       gleason7cat,
                                                       exprsdat_gleason7_allProbes_long$logExprs)
colnames(exprsdat_gleason7_allProbes_long_probeID) <- c("probeID", "sampleID", "gleason7cat", "logExprs")
head(exprsdat_gleason7_allProbes_long_probeID)
pdf(file = "Graphs/ViolinPlotsAllprobes.pdf"
    ,onefile = TRUE
    , width = 11
    , height = 8.5)
numOfPages <- round( (ncol(pheno_survival_withProbes[ ,c(8:53)]) / 3) + 0.499999)
for (currentPage in 1:numOfPages){
  print(ggplot(exprsdat_gleason7_allProbes_long_probeID,
         aes(x = probeID, y = logExprs, fill = gleason7cat, color = gleason7cat)) +
    geom_violin() +
    facet_wrap_paginate( ~ probeID, scales = "free_x", ncol = 3, nrow = 1, page = currentPage) +
    ggtitle("Log Expression of All Probes Associated with Key Genes"))
}
dev.off()
#I cannot perform a survival analysis based on the expression of the probes.
#I would need to establish a cut off for expression high and expression low for 
#each of the probes and I cannto do this since the expression profiles are
#so close. Instead, I will perform survival analysis based on biochemical
#recurrence and Gleason score.
pheno_survival_withProbes_tibble$timeToRecurrence <- as.numeric(pheno_survival_withProbes_tibble$timeToRecurrence)
s <- survfit(Surv(timeToRecurrence, bcrRelapse) ~ tumorGleason, data = pheno_survival_withProbes_tibble)
plot(s)
pdf(file = "Graphs/KaplanMeierPlots.pdf"
    ,onefile = TRUE
    , width = 11
    , height = 8.5)
ggsurvplot(s, data = pheno_survival_withProbes_tibble,
           pval = TRUE, risk.table=TRUE,
           legend.labs=c("Gleason3+4","Gleason4+3"),
           title="Kaplan-Meier Curve for Prostate Cancer Cases with Biochemical Recurrence",
           risk.table.height=0.15)
dev.off()
#There is no difference between the two groups. It doesn't matter if the 
#patient has Gleason 3+4 or Gleason 4+3. 























