#loading libraries:
library(tidyverse)
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
pheno_bcr <- pheno %>% select(bcrRelapse, clinicalStage, tumorGleason)

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

pheno_pivot_table_gleason7 <- dcast(data = pheno_table
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
