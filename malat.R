#BiocManager::install("TCGAbiolinks")
#BiocManager::install("TCGAutils")
#BiocManager::install("curatedTCGAData")
#BiocManager::install("TCGAutils")

library(curatedTCGAData)
library(TCGAutils)
library(MultiAssayExperiment)
library(patchwork)

head(
  curatedTCGAData(
    diseaseCode = "COAD", assays = "*", version = "1.1.38"
  )
)


#library(TCGAbiolinks)
library(dplyr)

browseVignettes("curatedTCGAData")

# Get all codes on TCGA
data('diseaseCodes', package = "TCGAutils")
diseaseCodes[,c(1,4)]

# Get sample types on TCGA
data(sampleTypes, package = "TCGAutils")
head(sampleTypes)

# Get 5 Cancers data from TCGA 
#readData = curatedTCGAData(c("CHOL","LIHC","CESC","HNSC","COAD"), assays =  "RNASeq2GeneNorm", 
#                           version = '2.0.1',
#                           dry.run = FALSE)



################ 1 ###############
# $`CESC_RNASeq2GeneNorm-20160128`
# 01  06  11 
# 304   2   3 
################ 2 ###############
# $`CHOL_RNASeq2GeneNorm-20160128`
# 01 11 
# 36  9 
################# 3 ############### 
# $`COAD_RNASeq2GeneNorm_illuminaga-20160128`
# 01 
# 191 
################# 4 ###############
# $`COAD_RNASeq2GeneNorm_illuminahiseq-20160128`
# 01  02  06  11 
# 283   1   1  41 
################# 5 ###############
# $`HNSC_RNASeq2GeneNorm-20160128`
# 01  06  11 
# 520   2  44 
################# 6 ###############
# $`LIHC_RNASeq2GeneNorm-20160128`
# 01  02  11 
# 371   2  50 

# 1-CESC 2-CHOL 3-COAD 4-COAD 5-HNSC 6-LIHC




COADData = curatedTCGAData(c("COAD"), assays =  c("RNASeq2Gene", "RNASeq2GeneNorm"),
                            version = '2.0.1',
                            dry.run = FALSE)

# Preform simple exploration on the dataset
typeof(COADData)
rownames(COADData)
colnames(COADData)
colData(COADData)
getClinicalNames("COAD")
table(table(sampleMap(COADData)$primary)) 

# Get how many sample types in each experiment/cancer 
sampleTables(COADData)
sampleTables(COADData,vial = TRUE)

# extract assay data

coadraw.data <- (assay(COADData[[1]]))

as.integer(coadraw.data)
round(coadraw.data)

coad.pritumor.data <- coadraw.data[,which(substr(colnames(coadraw.data),14,16) == "01A")]

coad.normal.data <- coadraw.data[,which(substr(colnames(coadraw.data),14,16) == "11A")]


# Data pre-processing and clean up
# Assay data mapping to clinical data and create sample info object
# preform DESeq2 analysis
# Select EMT gene markers , cancer stem cells , metastatic and find correlation


summarise(test, avg = mean(test$`TCGA-2W-A8YY-01A-11R-A37O-07`))
test <- as.data.frame(CESC.cancer)

test["MALAT1",]

# Get all samples IDs from the first cancer
colnames(COADData[[1]])
# Get the RNA assays of normal sample/patients from first cancer CESC
CESC.normal <- assay(COADData[[1]])[,which(substr(colnames(COADData[[1]]),14,16) == "11A")]

dim(CESC.normal)

CESC.normal <- log2(CESC.normal + 1)

# Get the counts of MALAT-1 gene from all normal samples/patients in CESC
CESC.normal.malat <- as.data.frame(CESC.normal["MALAT1",])

# Get the RNA assays of diseased sample/patients from first cancer CESC
CESC.cancer <- assay(COADData[[1]])[,which(substr(colnames(COADData[[1]]),14,16) == "01A")]

dim(CESC.cancer)

CESC.cancer <- log2(CESC.cancer + 1)

# Get the counts of MALAT-1 gene from all diseased samples/patients in CHOL
CESC.cancer.malat <- as.data.frame(CESC.cancer["MALAT1",])

# Conduct normalization and transformation
# Preform analysis using Limma
# Visualize 

# rafalib::mypar(1,2)
boxplot(c(CESC.normal.malat,CESC.cancer.malat),ylab = "expression",xlab = "MALAT-1 between cancer and normal in Cervical squamous cell carcinoma and endocervical adenocarcinoma")

colnames(CESC.normal.malat) <- "Normal MALAT-1"
colnames(CESC.cancer.malat) <- "Cancer MALAT-1"


library(ggplot2)
# create a data frame with some sample data
data <- data.frame(
  group = rep(c("A", "B"), each = 50),
  value = c(rnorm(50, mean = 0, sd = 1), rnorm(50, mean = 2, sd = 1))
)


dim(CESC.cancer.malat)

df.data <- data.frame(c(CESC.normal.malat,CESC.cancer.malat))

ggplot(data, aes(x = group, y = value)) +
  geom_boxplot()


colnames(df.data) <- c("group1","group2")

library(ggplot2)

# create a data frame with the two sets of data
df <- data.frame(
  Group = c(rep("Normal", length(malat.normal)), rep("Cancer", length(malat.cancer))),
  MALAT = c(malat.normal, malat.cancer)
)

# create the boxplot using ggplot
ggplot(df, aes(x = Group, y = MALAT)) +
  geom_boxplot() +
  labs(x = "Group", y = "MALAT") +
  ggtitle("Distribution of MALAT by Group")


t.test(CESC.normal.malat,CESC.cancer.malat)