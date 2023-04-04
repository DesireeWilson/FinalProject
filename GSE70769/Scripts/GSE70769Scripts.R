#loading libraries:
library(tidyverse)
library(GEOquery)
library(gt)
library(gtsummary)


#setting working directory:
thepath <- "C:/Users/wilso/Documents/FinalProject/GSE70769/"
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

