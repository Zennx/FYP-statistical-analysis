library(edgeR)
library(dplyr)
library(pheatmap)

################### Data Merging ###################
# Function to read and prepare each file
setwd("/Users/Zen/Library/CloudStorage/OneDrive-InternationalMedicalUniversity/Biomedical Science BM122/FYP research/Transcriptomics Analysis - Zen/Data/FeatureCounts")
read_and_prepare <- function(file) {
  df <- read.table(file, sep = "\t", header = TRUE)
  
  # Extract sample name from file name
  sample_name <- gsub(".tabular", "", file)
  
  # Rename the second column (counts) to the sample name
  colnames(df)[2] <- sample_name
  
  return(df)
}
# List all files in the directory ending with .txt
files <- list.files(pattern = "*.tabular")

# Read the first file to initialize the merged dataframe
fc <- read_and_prepare(files[1])

# Loop over the remaining files and perform left_join
for (file in files[-1]) {
  df <- read_and_prepare(file)
  fc <- left_join(fc, df, by = colnames(df)[1])  # Join by the first column, typically the gene identifier
}

# Optional: Handle any missing values
fc <- fc %>%
  replace(is.na(.), 0)

# View the first few rows of the combined data
head(fc)

######################### EdgeR ############################

group <- factor(c(1,1,2,2,3,3,4,4))
design <- model.matrix(~group)
bcv <- 0.6
counts <- matrix( rnbinom(40,size=1/bcv^2,mu=10), 20,2)
dgelist <- DGEList(counts=fc, group=group)
dgelist <- estimateCommonDisp(dgelist)
dgelist <- estimateGLMTagwiseDisp(dgelist, design)
et <- exactTest(dgelist, dispersion=bcv^2)
fit <- glmQLFit(dgelist, design)
qlf.3vs1 <- glmQLFTest(fit, coef=3)
tr <- glmTreat(fit, coef=1.5, lfc=1)
topTags(tr, n=20)
logcpm <- cpm(dgelist, log=TRUE)
heatmap <- pheatmap(logcpm, cluster_rows=TRUE, cluster_cols=TRUE, show_rownames=FALSE, show_colnames=TRUE, main="Heatmap of logCPM values")

print(heatmap)
