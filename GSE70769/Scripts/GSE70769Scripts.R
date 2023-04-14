#loading libraries:
library(tidyverse)
library(GEOquery)
library(gt)
library(gtsummary)
library(pheatmap)


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
exprsdat <- exprs(gse70769[[1]])


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

