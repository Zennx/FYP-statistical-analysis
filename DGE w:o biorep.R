library(edgeR)
library(dplyr)
library(ggplot2)
library(EnhancedVolcano)
library(UpSetR)
library(clusterProfiler)

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
#specify elements for design matrix
targets <- data.frame(samples=c("MR220CL2_NS", "MR220CL2_S", "NMR152_NS", "NMR152_S", "ML82_2_NS", "ML82_2_S", "ML125_2_NS", "ML125_2_S"),
                      Lines=c("MR220CL2", "MR220CL2", "NMR152", "NMR152", "ML82_2", "ML82_2", "ML125_2", "ML125_2"),
                      Treat=c("NS", "S", "NS", "S", "NS", "S", "NS", "S"))

Group <- factor(paste(targets$Lines,targets$Treat,sep="."))
targets <- cbind(targets,Group=Group)

group <- targets$Group
Line <- targets$Line
Treat <- targets$Treat

#generate designmatrix. The concept of the matrix is to compare how each line respond to the given treatment aka stress (S) vs non-stress (NS)
design <- model.matrix(~0+Line + Line:Treat, data=fc)

#counts <- matrix( rnbinom(40,size=1/bcv^2,mu=10), 20,2)
#generate DGEList object
dgelist <- DGEList(counts=fc, group=group)
#filter out lowly expressed genes
keep <- filterByExpr(dgelist, min.count=10, min.total.count=15)
dgelist <- dgelist[keep, , keep.lib.sizes=FALSE]

#normalise the counts
dgelist<- calcNormFactors(dgelist)

################### PCA/MDS plot ###################
#plotMDS(dgelist, top = 500, pch = 16, cex = 0.5, main = "EdgeR DEG MDS plot", gene.selection = "common")

#get log counts per million (logCPM) values
cpm <- cpm(dgelist,normalized.lib.sizes = TRUE, log = TRUE)

#convert the logcpm values to eigenvectors
pca_res <- prcomp(t(cpm), center = TRUE, scale. = TRUE)

#generate the pca dataframe
pca <-as.data.frame(pca_res$x) %>% 
  tibble::rownames_to_column("Sample") %>% 
  mutate(Lines = targets$Lines) %>% 
  mutate(Treat = targets$Treat) %>%
  select(Sample,Lines,Treat, PC1:PC8)

#plot the pca plot using ggplot
pca_plot <- ggplot(pca, aes(PC1, PC2,color = Lines, label = Sample) ) + 
  geom_point(size = 3) +
  theme_bw() + geom_text_repel(aes(label = Sample)) + ggtitle("PCA plot")


################### Differential Expression ###################
#given that this experiments has no bioreplicates, we provide an estimate of the biological coefficient of variation
bcv <- 0.3 # biological coefficient of variation, sqrt dispersion

#test the exactTest function (simplest way to compare between the conditions) without complex model fitting (unapplicatble in this case)
et <- exactTest(dgelist, dispersion=bcv^2)
#view exact test results 
topTags(et)
#> topTags(et)
#Comparison of groups:  ML125-2.S-ML125-2.NS 
#Geneid    logFC    logCPM       PValue          FDR
#34231 Os10g0542100 16.95273  7.353160 2.805980e-14 1.033835e-09
#27955 Os08g0104400 12.66004  9.518500 1.163951e-13 1.114090e-09


#dgelist$samples$group
#[1] MR220CL2.NS MR220CL2.S  NMR152.NS   NMR152.S    ML82-2.NS   ML82-2.S    ML125-2.NS  ML125-2.S  
#Levels: ML125-2.NS ML125-2.S ML82-2.NS ML82-2.S MR220CL2.NS MR220CL2.S NMR152.NS NMR152.S
####the above only does the default comparison of pair 1:2 . to get all the required comparisons for the purpose of this *f-ed* up experiments state all the contrast required ie. comparisons then do all the required pairwise comparisons

#specify the comparisosn
#A vs B， A= treatment/condition，B=control/Wt
#Compare for each line, S vs NS
#MR220CL2.SvsNS = 'MR220CL2.S'-'MR220CL2.NS', #parent (P0)
#NMR152.SvsNS = 'NMR152.S'-'NMR152.NS', #generation 1 (F1)
#ML82_2.SvsNS = 'ML82_2.S'-'ML82_2.NS', #generation2-line1 (F2L1)
#ML125_2.SvsNS = 'ML125_2.S'-'ML125_2.NS', #generation2-line2 (F2L2)


#Contrasts for each stressed samples VS parental MR220CL2 (S)
#NMR152vsMR220 = NMR152.S - MR220CL2.S,
#ML82vsMR220 = "ML82-2.S" - MR220CL2.S,
#ML125vsMR220 = "ML125-2.S" - MR220CL2.S,

################### Pairwise comparison ###################
###test the pairwise comparison to see if the contrast is correct
#Compare for each line, S vs NS
ML125_2.SvsNS_et <- exactTest(dgelist, dispersion=bcv^2, pair = 1:2 )
ML82_2.SvsNS_et <- exactTest(dgelist, dispersion=bcv^2, pair = 3:4 )
MR220CL2.SvsNS_et <- exactTest(dgelist, dispersion=bcv^2, pair = 5:6 )
NMR152.SvsNS_et <- exactTest(dgelist, dispersion=bcv^2, pair = 7:8 )

#Contrasts for each stressed samples VS parental MR220CL2 (S)
NMR152vsMR220_et <- exactTest(dgelist, dispersion=bcv^2, pair = c(6,8))
ML82vsMR220_et <- exactTest(dgelist, dispersion=bcv^2, pair = c(6,4))
ML125vsMR220_et <- exactTest(dgelist, dispersion=bcv^2, pair = c(6,2))

#view exact test results
topTags(ML125_2.SvsNS_et)
topTags(ML82_2.SvsNS_et)
topTags(MR220CL2.SvsNS_et)
topTags(NMR152.SvsNS_et)
topTags(NMR152vsMR220_et)
topTags(ML82vsMR220_et)
topTags(ML125vsMR220_et)

#> topTags(ML125_2.SvsNS_et)
#Comparison of groups:  ML125_2.S-ML125_2.NS
################### Volcano plot ###################

#Duplicate the geneid, logFC and PValue columns from the table of the 3vs1 comparison to new DF
#volcanodata <- data.frame(Geneid = c(qlf.3vs1$genes$Geneid), logFC = c(qlf.3vs1$table$logFC), PValue = c(qlf.3vs1$table$PValue))
ML125_2.SvsNS_et_df <- ML125_2.SvsNS_et %>% data.frame()
ML82_2.SvsNS_et_df <- ML82_2.SvsNS_et %>% data.frame()
MR220CL2.SvsNS_et_df <- MR220CL2.SvsNS_et %>% data.frame()
NMR152.SvsNS_et_df <- NMR152.SvsNS_et %>% data.frame()
NMR152vsMR220_et_df <- NMR152vsMR220_et %>% data.frame()
ML82vsMR220_et_df <- ML82vsMR220_et %>% data.frame()
ML125vsMR220_et_df <- ML125vsMR220_et %>% data.frame()

#plot the volcano plot for ML125_2.SvsNS_et 
{
volcaplt <-ML125_2.SvsNS_et_df

keyvals <- ifelse(
  volcaplt$logFC < -2 & volcaplt$PValue < 0.05, 'lightgreen',
  ifelse(volcaplt$logFC > 2 & volcaplt$PValue < 0.05, 'red1',
         'grey'))
keyvals[is.na(keyvals)] <- 'black'
names(keyvals)[keyvals == 'red1'] <- 'Overexpressed'
names(keyvals)[keyvals == 'grey'] <- 'N.S. OR |logFC|>=2'
names(keyvals)[keyvals == 'lightgreen'] <- 'Underexpressed'

volcano <- EnhancedVolcano(
  volcaplt,
  lab = volcaplt$Geneid,
  x = "logFC",
  y = "PValue",
  xlim = c(-15, 15),
  ylim = c(0, 7.5),
  pCutoff = 0.05,
  FCcutoff = 2,
  labSize = 2,
  title = 'EdgeR DEG Volcano Plot',
  titleLabSize = 15,
  subtitle = 'ML125-2 Non-Stressed vs ML125-2 Stressed',
  subtitleLabSize = 12,
  legendPosition = 'right',
  drawConnectors = TRUE,
  max.overlaps = 20,
  colCustom = keyvals,
  legendLabSize = 12,
  legendIconSize = 3.0,
  axisLabSize = 12
)
print(volcano)
}
#plot the volcano plot for ML82_2.SvsNS_et
{
volcaplt <- ML82_2.SvsNS_et

keyvals <- ifelse(
  volcaplt$logFC < -2 & volcaplt$PValue < 0.05, 'lightgreen',
  ifelse(volcaplt$logFC > 2 & volcaplt$PValue < 0.05, 'red1',
         'grey'))
keyvals[is.na(keyvals)] <- 'black'
names(keyvals)[keyvals == 'red1'] <- 'Overexpressed'
names(keyvals)[keyvals == 'grey'] <- 'N.S. OR |logFC|>=2'
names(keyvals)[keyvals == 'lightgreen'] <- 'Underexpressed'

volcano <- EnhancedVolcano(
  volcaplt,
  lab = volcaplt$Geneid,
  x = "logFC",
  y = "PValue",
  xlim = c(-15, 15),
  ylim = c(0, 7.5),
  pCutoff = 0.05,
  FCcutoff = 2,
  labSize = 2,
  title = 'EdgeR DEG Volcano Plot',
  titleLabSize = 15,
  subtitle = 'ML82-2 Non-Stressed vs ML82-2 Stressed',
  subtitleLabSize = 12,
  legendPosition = 'right',
  drawConnectors = TRUE,
  max.overlaps = 20,
  colCustom = keyvals,
  legendLabSize = 12,
  legendIconSize = 3.0,
  axisLabSize = 12
)
print(volcano)
}
#plot the volcano plot for MR220CL2.SvsNS_et
{
  volcaplt <- MR220CL2.SvsNS_et_df
    
    keyvals <- ifelse(
      volcaplt$table$logFC < -2 & volcaplt$PValue < 0.05, 'lightgreen',
      ifelse(volcaplt$table$logFC > 2 & volcaplt$PValue < 0.05, 'red1',
             'grey'))
  keyvals[is.na(keyvals)] <- 'black'
  names(keyvals)[keyvals == 'red1'] <- 'Overexpressed'
  names(keyvals)[keyvals == 'grey'] <- 'N.S. OR |logFC|>=2'
  names(keyvals)[keyvals == 'lightgreen'] <- 'Underexpressed'
  
  volcano <- EnhancedVolcano(
    volcaplt,
    lab = volcaplt$Geneid,
    x = "logFC",
    y = "PValue",
    xlim = c(-15, 15),
    ylim = c(0, 7.5),
    pCutoff = 0.05,
    FCcutoff = 2,
    labSize = 2,
    title = 'EdgeR DEG Volcano Plot',
    titleLabSize = 15,
    subtitle = 'MR220CL2 Non-Stressed vs MR220CL2 Stressed',
    subtitleLabSize = 12,
    legendPosition = 'right',
    drawConnectors = TRUE,
    max.overlaps = 20,
    colCustom = keyvals,
    legendLabSize = 12,
    legendIconSize = 3.0,
    axisLabSize = 12
  )
  print(volcano)
}
#plot the volcano plot for NMR152.SvsNS_et
{
  volcaplt <- NMR152.SvsNS_et
    
    keyvals <- ifelse(
      volcaplt$logFC < -2 & volcaplt$PValue < 0.05, 'lightgreen',
      ifelse(volcaplt$logFC > 2 & volcaplt$PValue < 0.05, 'red1',
             'grey'))
  keyvals[is.na(keyvals)] <- 'black'
  names(keyvals)[keyvals == 'red1'] <- 'Overexpressed'
  names(keyvals)[keyvals == 'grey'] <- 'N.S. OR |logFC|>=2'
  names(keyvals)[keyvals == 'lightgreen'] <- 'Underexpressed'
  
  volcano <- EnhancedVolcano(
    volcaplt,
    lab = volcaplt$Geneid,
    x = "logFC",
    y = "PValue",
    xlim = c(-15, 15),
    ylim = c(0, 7.5),
    pCutoff = 0.05,
    FCcutoff = 2,
    labSize = 0,
    title = 'EdgeR DEG Volcano Plot',
    titleLabSize = 15,
    subtitle = 'NMR152 Non-Stressed vs NMR152 Stressed',
    subtitleLabSize = 12,
    legendPosition = 'right',
    drawConnectors = FALSE,
    max.overlaps = 20,
    colCustom = keyvals,
    legendLabSize = 12,
    legendIconSize = 3.0,
    axisLabSize = 12
  )
  print(volcano)
}
#plot the volcano plot for NMR152vsMR220_et
{
  volcaplt <- NMR152vsMR220_et
    
    keyvals <- ifelse(
      volcaplt$logFC < -2 & volcaplt$PValue < 0.05, 'lightgreen',
      ifelse(volcaplt$logFC > 2 & volcaplt$PValue < 0.05, 'red1',
             'grey'))
  keyvals[is.na(keyvals)] <- 'black'
  names(keyvals)[keyvals == 'red1'] <- 'Overexpressed'
  names(keyvals)[keyvals == 'grey'] <- 'N.S. OR |logFC|>=2'
  names(keyvals)[keyvals == 'lightgreen'] <- 'Underexpressed'
  
  volcano <- EnhancedVolcano(
    volcaplt,
    lab = volcaplt$Geneid,
    x = "logFC",
    y = "PValue",
    xlim = c(-15, 15),
    ylim = c(0, 7.5),
    pCutoff = 0.05,
    FCcutoff = 2,
    labSize = 0,
    title = 'EdgeR DEG Volcano Plot',
    titleLabSize = 15,
    subtitle = 'NMR152 Stressed vs MR220 Stressed',
    subtitleLabSize = 12,
    legendPosition = 'right',
    drawConnectors = FALSE,
    max.overlaps = 20,
    colCustom = keyvals,
    legendLabSize = 12,
    legendIconSize = 3.0,
    axisLabSize = 12
  )
  print(volcano)
}
#plot the volcano plot for ML82vsMR220_et
{
  volcaplt <- ML82vsMR220_et
    
    keyvals <- ifelse(
      volcaplt$logFC < -2 & volcaplt$PValue < 0.05, 'lightgreen',
      ifelse(volcaplt$logFC > 2 & volcaplt$PValue < 0.05, 'red1',
             'grey'))
  keyvals[is.na(keyvals)] <- 'black'
  names(keyvals)[keyvals == 'red1'] <- 'Overexpressed'
  names(keyvals)[keyvals == 'grey'] <- 'N.S. OR |logFC|>=2'
  names(keyvals)[keyvals == 'lightgreen'] <- 'Underexpressed'
  
  volcano <- EnhancedVolcano(
    volcaplt,
    lab = volcaplt$Geneid,
    x = "logFC",
    y = "PValue",
    xlim = c(-15, 15),
    ylim = c(0, 7.5),
    pCutoff = 0.05,
    FCcutoff = 2,
    labSize = 0,
    title = 'EdgeR DEG Volcano Plot',
    titleLabSize = 15,
    subtitle = 'ML82 Stressed vs MR220 Stressed',
    subtitleLabSize = 12,
    legendPosition = 'right',
    drawConnectors = FALSE,
    max.overlaps = 20,
    colCustom = keyvals,
    legendLabSize = 12,
    legendIconSize = 3.0,
    axisLabSize = 12
  )
  print(volcano)
}
#plot the volcano plot for ML125vsMR220_et
{
  volcaplt <- ML125vsMR220_et
    
    keyvals <- ifelse(
      volcaplt$logFC < -2 & volcaplt$PValue < 0.05, 'lightgreen',
      ifelse(volcaplt$logFC > 2 & volcaplt$PValue < 0.05, 'red1',
             'grey'))
  keyvals[is.na(keyvals)] <- 'black'
  names(keyvals)[keyvals == 'red1'] <- 'Overexpressed'
  names(keyvals)[keyvals == 'grey'] <- 'N.S. OR |logFC|>=2'
  names(keyvals)[keyvals == 'lightgreen'] <- 'Underexpressed'
  
  volcano <- EnhancedVolcano(
    volcaplt,
    lab = volcaplt$Geneid,
    x = "logFC",
    y = "PValue",
    xlim = c(-15, 15),
    ylim = c(0, 7.5),
    pCutoff = 0.05,
    FCcutoff = 2,
    labSize = 0,
    title = 'EdgeR DEG Volcano Plot',
    titleLabSize = 15,
    subtitle = 'ML125 Stressed vs MR220 Stressed',
    subtitleLabSize = 12,
    legendPosition = 'right',
    drawConnectors = FALSE,
    max.overlaps = 20,
    colCustom = keyvals,
    legendLabSize = 12,
    legendIconSize = 3.0,
    axisLabSize = 12
  )
  print(volcano)
}

####Zen - plot the rest of the comparison's volcano plot

#####do not just change the variable and run because you risk taking the wrong data!!!!!!!!

################### Heatmap merged ###################
SvsNS <- data.frame(ML82_2 = ML82_2.SvsNS_et_df$logFC, 
                    NMR152 = NMR152.SvsNS_et_df$logFC,
                    MR220CL2 = MR220CL2.SvsNS_et_df$logFC,
                    ML125_2 = ML125_2.SvsNS_et_df$logFC) %>%
         select(MR220CL2, NMR152, ML125_2, ML82_2)

SvsS <- data.frame(ML82_2 = ML82vsMR220_et_df$logFC, 
                    NMR152 = NMR152vsMR220_et_df$logFC,
                    ML125_2 = ML125vsMR220_et_df$logFC)

SvsS <-filter_if(SvsS, is.numeric, all_vars((.) != 0))

#heatmap(SvsNS, cluster_rows=TRUE, cluster_cols=TRUE, show_rownames=FALSE, show_colnames=TRUE, main="Heatmap of logFC values")
#heatmap to visualise expression across sample
#include column annotation for the heatmap
my_sample_col <- data.frame(
  sample = c("MR220CL2", "NMR152", "ML125_2", "ML82_2"),
  group = c("P0", "F1", rep("F2", 2))
) %>% tibble::column_to_rownames("sample")
#for this use logcpm generated
# Define the sample group names and corresponding colors
group_names <- c("P0", "F1", "F2")
colours <- c("pink", "magenta4", "firebrick3")


# Create a named list for annotation colors
anno_color <- list(group = setNames(colours, group_names))
DEGmapSvsNS <- pheatmap(SvsNS, 
                        cluster_rows=TRUE, 
                        cluster_cols=FALSE, 
                        show_rownames=FALSE, 
                        show_colnames=TRUE, 
                        main="Heatmap of logFC values between lines for Stressed vs NonStressed", 
                        col = viridis(50), 
                        annotation = my_sample_col,
                        annotation_colors = anno_color,# add the group information
                        scale="row")


#2nd heatmap on SvsS
my_sample_col2 <- data.frame(
  sample = c("NMR152", "ML125_2", "ML82_2"),
  group = c("F1", rep("F2", 2))
) %>% tibble::column_to_rownames("sample")
#for this use logcpm generated
# Define the sample group names and corresponding colors
group_names <- c("F1", "F2")
colours <- c("magenta4", "firebrick3")

DEGmapSvsS <- pheatmap(SvsS, 
                        clustering_method = "complete",
                        cluster_rows=TRUE, 
                        cluster_cols=FALSE, 
                        show_rownames=FALSE,
                        show_colnames=TRUE, 
                        main="LogFC values between Stressed lines vs Stressed MR220CL2", 
                        col = viridis(50), 
                        annotation = my_sample_col2,
                        annotation_colors = anno_color,# add the group information
                        scale="row")



#DEGmapSvsNS <- heatmap(SvsNS, main="Heatmap of logFC values for Stressed vs NonStressed", col = colorRampPalette(c("blue", "white", "red"))(100), scale="none")
#heatmap of each comparison

##combine the log2FC of each comparison into a single dataframe, then plot heatmap with the comparison as x axis (don't need to show gene name because your R will crash)
###this allows you to compare and see for the same set of genes, how is it changing across different comparisons
#two heatmaps - one by lines, one by the stress of the G2 and G3 lines to the stressed parent

################### Top 20 DEGs ###################
MR220CL2.SvsNS_tt <- topTags(MR220CL2.SvsNS_et, sort.by = c("PValue") , n=20)

NMR152.SvsNS_tt <- topTags(NMR152.SvsNS_et, sort.by = "PValue" , n=20)
NMR152vsMR220_tt <- topTags(NMR152vsMR220_et, sort.by = "PValue" , n=20)
ML82vsMR220_tt <- topTags(ML82vsMR220_et, sort.by = "PValue" , n=20)
ML125vsMR220_tt <- topTags(ML125vsMR220_et, sort.by = "PValue" , n=20)

merged_SvsNS <- data_frame(Geneid = MR220CL2.SvsNS_tt$table$Geneid)
merged_SvsNS <- inner_join(merged_SvsNS,
                           MR220CL2.SvsNS_et_df,
                           by = "Geneid") %>% select(`Geneid`, `logFC`) %>% rename(MR220CL2 = logFC)
merged_SvsNS <- inner_join(merged_SvsNS,
                           NMR152.SvsNS_et_df,
                          by = "Geneid") %>% select(`Geneid`, `MR220CL2`, `logFC`) %>% rename(NMR152 = logFC)
merged_SvsNS <- inner_join(merged_SvsNS,
                           ML125_2.SvsNS_et_df,
                           by = "Geneid") %>% select(`Geneid`,`MR220CL2`, `NMR152`,`logFC`) %>% 
                            rename(ML125_2 = logFC)
merged_SvsNS <- inner_join(merged_SvsNS,
                           ML82_2.SvsNS_et_df,
                           by = "Geneid") %>% select(`Geneid`,`MR220CL2`, `NMR152`, `ML125_2`, `logFC`) %>% 
                rename(ML82_2 = logFC) %>% 
                tibble::column_to_rownames("Geneid")


my_sample_col <- data.frame(
  sample = c("MR220CL2", "NMR152", "ML125_2", "ML82_2"),
  group = c("P0", "F1", rep("F2", 2))
) %>% tibble::column_to_rownames("sample")
#for this use logcpm generated
# Define the sample group names and corresponding colors
group_names <- c("P0", "F1", "F2")
colours <- c("pink", "magenta4", "firebrick3")


# Create a named list for annotation colors
anno_color <- list(group = setNames(colours, group_names))
DEGmapSvsNSTT <- pheatmap(merged_SvsNS, 
                        cluster_rows=TRUE, 
                        cluster_cols=FALSE, 
                        show_rownames=TRUE, 
                        show_colnames=TRUE, 
                        main="Heatmap of top 20 DEGs logFC values, Stressed vs NonStressed", 
                        col = colorRampPalette(c("orange", "white", "purple"))(100), 
                        annotation = my_sample_col,
                        annotation_colors = anno_color,# add the group information
                        scale="row", cellheight = 10, cellwidth = 10)

#converts the geneid to gene symbol
#convert_id <- riceidconverter::RiceIDConvert("Os01g0102800",'RAP', toType = 'SYMBOL')

#make another heatmap for the other data
NMR152vsMR220_tt <- topTags(NMR152vsMR220_et, sort.by = "PValue" , n=20)

merged_SvsS <- data_frame(Geneid = NMR152vsMR220_tt$table$Geneid)
merged_SvsS <- inner_join(merged_SvsS,
                           ML125vsMR220_et_df,
                           by = "Geneid") %>% select(`Geneid`, `logFC`) %>% rename(ML125_2 = logFC)
merged_SvsS <- inner_join(merged_SvsS,
                          NMR152vsMR220_et_df,
                          by = "Geneid") %>% select(`Geneid`, `ML125_2`, `logFC`) %>% rename(NMR152 = logFC)
merged_SvsS <- inner_join(merged_SvsS,
                           ML82vsMR220_et_df,
                           by = "Geneid") %>% select(`Geneid`, `ML125_2`, `NMR152`, `logFC`) %>% 
  rename(ML82_2 = logFC) %>% 
  tibble::column_to_rownames("Geneid")

#2nd heatmap on SvsS
my_sample_col2 <- data.frame(
  sample = c("NMR152", "ML125_2", "ML82_2"),
  group = c("F1", rep("F2", 2))
) %>% tibble::column_to_rownames("sample")
#for this use logcpm generated
# Define the sample group names and corresponding colors
group_names <- c("F1", "F2")
colours <- c("magenta4", "firebrick3")

DEGmapSvsS <- pheatmap(merged_SvsS, 
                       clustering_method = "complete",
                       cluster_rows=TRUE, 
                       cluster_cols=FALSE, 
                       show_rownames=TRUE,
                       show_colnames=TRUE, 
                       main="Top 20 DEGs (LogFC) Stressed mutants vs Stressed MR220CL2", 
                       col = colorRampPalette(c("orange", "white", "purple"))(100), 
                       annotation = my_sample_col2,
                       annotation_colors = anno_color,# add the group information
                       scale="row", cellheight = 10, cellwidth = 10)

################### Drought Stress genes ###################
DSgenes <- read.delim("/Users/Zen/Library/CloudStorage/OneDrive-InternationalMedicalUniversity/Biomedical Science BM122/FYP research/Transcriptomics Analysis - Zen/Data/DroughtStress.tsv")
merged_SvsNSDS <- data_frame(Geneid = DSgenes$Locus.ID)
merged_SvsNSDS <- inner_join(merged_SvsNSDS,
                           MR220CL2.SvsNS_et_df,
                           by = "Geneid") %>% select(`Geneid`, `logFC`) %>% rename(MR220CL2 = logFC)
merged_SvsNSDS <- inner_join(merged_SvsNSDS,
                           NMR152.SvsNS_et_df,
                           by = "Geneid") %>% select(`Geneid`, `MR220CL2`, `logFC`) %>% rename(NMR152 = logFC)
merged_SvsNSDS <- inner_join(merged_SvsNSDS,
                           ML125_2.SvsNS_et_df,
                           by = "Geneid") %>% select(`Geneid`,`MR220CL2`, `NMR152`,`logFC`) %>% 
  rename(ML125_2 = logFC)
merged_SvsNSDS <- inner_join(merged_SvsNSDS,
                           ML82_2.SvsNS_et_df,
                           by = "Geneid") %>% select(`Geneid`,`MR220CL2`, `NMR152`, `ML125_2`, `logFC`) %>% 
  rename(ML82_2 = logFC) %>% 
  tibble::column_to_rownames("Geneid")

#heatmap for drought stress genes
my_sample_col <- data.frame(
  sample = c("MR220CL2", "NMR152", "ML125_2", "ML82_2"),
  group = c("P0", "F1", rep("F2", 2))
) %>% tibble::column_to_rownames("sample")
#for this use logcpm generated
# Define the sample group names and corresponding colors
group_names <- c("P0", "F1", "F2")
colours <- c("pink", "magenta4", "firebrick3")


# Create a named list for annotation colors
anno_color <- list(group = setNames(colours, group_names))
DEGmapSvsNSTT <- pheatmap(merged_SvsNSDS, 
                          clustermethod = "ward.D2",
                          cluster_rows=TRUE, 
                          cluster_cols=FALSE, 
                          show_rownames=FALSE, 
                          show_colnames=TRUE, 
                          main="Heatmap of stress-related DEGs (logFC), Stressed vs NonStressed", 
                          col = colorRampPalette(c("orange", "white", "purple"))(100), 
                          annotation = my_sample_col,
                          annotation_colors = anno_color,# add the group information
                          scale="row")
##??upsetplotRvrdgreregreg
#2nd map for the other data
merged_SvsSDS <- data_frame(Geneid = DSgenes$Locus.ID)

merged_SvsSDS <- inner_join(merged_SvsSDS,
                          ML125vsMR220_et_df,
                          by = "Geneid") %>% select(`Geneid`, `logFC`) %>% rename(ML125_2 = logFC)
merged_SvsSDS <- inner_join(merged_SvsSDS,
                          NMR152vsMR220_et_df,
                          by = "Geneid") %>% select(`Geneid`, `ML125_2`, `logFC`) %>% rename(NMR152 = logFC)
merged_SvsSDS <- inner_join(merged_SvsSDS,
                          ML82vsMR220_et_df,
                          by = "Geneid") %>% select(`Geneid`, `NMR152`, `ML125_2`, `logFC`) %>% 
  rename(ML82_2 = logFC) %>% 
  tibble::column_to_rownames("Geneid")

#2nd heatmap on SvsS
my_sample_col2 <- data.frame(
  sample = c("NMR152", "ML125_2", "ML82_2"),
  group = c("F1", rep("F2", 2))
) %>% tibble::column_to_rownames("sample")
#for this use logcpm generated
# Define the sample group names and corresponding colors
group_names <- c("F1", "F2")
colours <- c("magenta4", "firebrick3")

DEGmapSvsS <- pheatmap(merged_SvsSDS, 
                       clustering_method = "ward.D2",
                       cluster_rows=TRUE, 
                       cluster_cols=FALSE, 
                       show_rownames=FALSE,
                       show_colnames=TRUE, 
                       main="Heatmap of stress-related DEGs (logFC), Stressed mutants vs Stressed parent", 
                       col = colorRampPalette(c("orange", "white", "purple"))(100), 
                       annotation = my_sample_col2,
                       annotation_colors = anno_color,# add the group information
                       scale="row")
################# upset plot ###################
#use upset plot to visualise similar and different genes expressed across each line
upsetSvsNS <- limma::decideTests(MR220CL2.SvsNS_et, method = "separate", adjust.method = "BH", PValue = 0.05, lfc = 2)
upsetSvsNS <- upset(fromList(upsetSvsNS), order.by = "freq")

################### GO enrichment ####################DO NOT USE THE GO DEFAULT CUZ IT SUCKS
#use clusterProfiler to do GO enrichment analysis
#use the top 100 DEGs to do the GO analysis
#use the function enrichGO to do the analysis - have to mannualy specify it's "MF", BP" or "CC"
#goto bioconductor and find the Rice Genome ~ Oryza sativa (I hope they have if not..call jie jie)
###Then you can compile report and show it to you supervisor
```
organism = "org.Osativa.eg.db"
BiocManager::install(organism, force = TRUE)
library(organism)
organism <- org.Osativa.eg.db

#No1 - MR220CL2.SvsNS
{
#gene list from DE of SvsNS 
diffexpressed <- 
  ifelse(MR220CL2.SvsNS_et_df$logFC > 2 & MR220CL2.SvsNS_et_df$PValue < 0.05, "UP",
         ifelse(MR220CL2.SvsNS_et_df$logFC < -2 & MR220CL2.SvsNS_et_df$PValue <0.05, "DOWN", 
                "NS"))

diff_gene_names_SvsNS <- MR220CL2.SvsNS_et_df %>% 
  dplyr::filter(diffexpressed == "UP" | diffexpressed == "DOWN") %>% 
  unique()

map_id <- RiceIDConvert(diff_gene_names_SvsNS$Geneid, fromType = "RAP", toType = "SYMBOL")

ego <- enrichGO(map_id$SYMBOL,
                OrgDb = organism,
                keyType = "SYMBOL",
                ont="BP") #BP, MF, or CC
#number of GO terms for
require(DOSE)
#can adjust it using scalebar using the ggplot setting 
barplot(ego) + ggtitle("GO enrichment for MR220CL2 Stressed vs Non-Stressed (BP)")
}
#No2 - NMR152.SvsNS
{
  #gene list from DE of SvsNS 
  diffexpressed <- 
    ifelse(NMR152.SvsNS_et_df$logFC > 2 & NMR152.SvsNS_et_df$PValue < 0.05, "UP",
           ifelse(NMR152.SvsNS_et_df$logFC < -2 & NMR152.SvsNS_et_df$PValue <0.05, "DOWN", 
                  "NS"))
  
  diff_gene_names_SvsNS <- NMR152.SvsNS_et_df %>% 
    dplyr::filter(diffexpressed == "UP" | diffexpressed == "DOWN") %>% 
    unique()
  
  map_id <- RiceIDConvert(diff_gene_names_SvsNS$Geneid, fromType = "RAP", toType = "SYMBOL")
  
  ego <- enrichGO(map_id$SYMBOL,
                  OrgDb = organism,
                  keyType = "SYMBOL",
                  ont="BP") #BP, MF, or CC
  #number of GO terms for
  require(DOSE)
  #can adjust it using scalebar using the ggplot setting 
  barplot(ego) + ggtitle("GO enrichment for MR220CL2 Stressed vs Non-Stressed (BP)")
#NOTE: NO GO TERM FOUND FOR NMR152.SvsNS for BP MF and CC,  BUT IT'S OKAY
  }
#No3 - ML82-2.SvsNS
{
  #gene list from DE of SvsNS 
  diffexpressed <- 
    ifelse(ML82_2.SvsNS_et_df$logFC > 2 & ML82_2.SvsNS_et_df$PValue < 0.05, "UP",
           ifelse(ML82_2.SvsNS_et_df$logFC < -2 & ML82_2.SvsNS_et_df$PValue <0.05, "DOWN", 
                  "NS"))
  
  diff_gene_names_SvsNS <- ML82_2.SvsNS_et_df %>% 
    dplyr::filter(diffexpressed == "UP" | diffexpressed == "DOWN") %>% 
    unique()
  
  map_id <- RiceIDConvert(diff_gene_names_SvsNS$Geneid, fromType = "RAP", toType = "SYMBOL")
  
  ego <- enrichGO(map_id$SYMBOL,
                  OrgDb = organism,
                  keyType = "SYMBOL",
                  ont="BP") #BP, MF, or CC
  #number of GO terms for
  require(DOSE)
  #can adjust it using scalebar using the ggplot setting 
  barplot(ego) + ggtitle("GO enrichment for MR220CL2 Stressed vs Non-Stressed (BP)")
#NOTE: NO GO TERM FOUND FOR ML82-2.SvsNS for BP MF and CC too,  BUT IT'S OKAY
  }
#No4 - ML125-2.SvsNS
{
  #gene list from DE of SvsNS 
  diffexpressed <- 
    ifelse(ML125_2.SvsNS_et_df$logFC > 2 & ML125_2.SvsNS_et_df$PValue < 0.05, "UP",
           ifelse(ML125_2.SvsNS_et_df$logFC < -2 & ML125_2.SvsNS_et_df$PValue <0.05, "DOWN", 
                  "NS"))
  
  diff_gene_names_SvsNS <- ML125_2.SvsNS_et_df %>% 
    dplyr::filter(diffexpressed == "UP" | diffexpressed == "DOWN") %>% 
    unique()
  
  map_id <- RiceIDConvert(diff_gene_names_SvsNS$Geneid, fromType = "RAP", toType = "SYMBOL")
  
  ego <- enrichGO(map_id$SYMBOL,
                  OrgDb = organism,
                  keyType = "SYMBOL",
                  ont="BP") #BP, MF, or CC
  #number of GO terms for
  require(DOSE)
  #can adjust it using scalebar using the ggplot setting 
  barplot(ego) + ggtitle("GO enrichment for ML125-2 Stressed vs Non-Stressed (BP)")
# No results for MF and CC
  }
#No5 - ML125 vs MR220
{
  #gene list from DE of SvsNS 
  diffexpressed <- 
    ifelse(ML125vsMR220_et_df$logFC > 2 & ML125vsMR220_et_df$PValue < 0.05, "UP",
           ifelse(ML125vsMR220_et_df$logFC < -2 & ML125vsMR220_et_df$PValue <0.05, "DOWN", 
                  "NS"))
  
  diff_gene_names_SvsNS <- ML125vsMR220_et_df %>% 
    dplyr::filter(diffexpressed == "UP" | diffexpressed == "DOWN") %>% 
    unique()
  
  map_id <- RiceIDConvert(diff_gene_names_SvsNS$Geneid, fromType = "RAP", toType = "SYMBOL")
  
  ego <- enrichGO(map_id$SYMBOL,
                  OrgDb = organism,
                  keyType = "SYMBOL",
                  ont="CC") #BP, MF, or CC
  #number of GO terms for
  require(DOSE)
  #can adjust it using scalebar using the ggplot setting 
  barplot(ego) + ggtitle("GO enrichment for MR220CL2 Stressed vs Non-Stressed (BP)")
# also no results
  }
#No6 - NMR152 vs MR220
{
  #gene list from DE of SvsNS 
  diffexpressed <- 
    ifelse(NMR152vsMR220_et_df$logFC > 2 & NMR152vsMR220_et_df$PValue < 0.05, "UP",
           ifelse(NMR152vsMR220_et_df$logFC < -2 & NMR152vsMR220_et_df$PValue <0.05, "DOWN", 
                  "NS"))
  
  diff_gene_names_SvsNS <- NMR152vsMR220_et_df %>% 
    dplyr::filter(diffexpressed == "UP" | diffexpressed == "DOWN") %>% 
    unique()
  
  map_id <- RiceIDConvert(diff_gene_names_SvsNS$Geneid, fromType = "RAP", toType = "SYMBOL")
  
  ego <- enrichGO(map_id$SYMBOL,
                  OrgDb = organism,
                  keyType = "SYMBOL",
                  ont="BP") #BP, MF, or CC
  #number of GO terms for
  require(DOSE)
  #can adjust it using scalebar using the ggplot setting 
  barplot(ego) + ggtitle("GO enrichment for Stressed NMR152 VS MR220CL2 (BP)")
#Note: no results for MF and CC
  }
#No7 - ML82 vs MR220
{
  #gene list from DE of SvsNS 
  diffexpressed <- 
    ifelse(ML82vsMR220_et_df$logFC > 2 & ML82vsMR220_et_df$PValue < 0.05, "UP",
           ifelse(ML82vsMR220_et_df$logFC < -2 & ML82vsMR220_et_df$PValue <0.05, "DOWN", 
                  "NS"))
  
  diff_gene_names_SvsNS <- ML82vsMR220_et_df %>% 
    dplyr::filter(diffexpressed == "UP" | diffexpressed == "DOWN") %>% 
    unique()
  
  map_id <- RiceIDConvert(diff_gene_names_SvsNS$Geneid, fromType = "RAP", toType = "SYMBOL")
  
  ego <- enrichGO(map_id$SYMBOL,
                  OrgDb = organism,
                  keyType = "SYMBOL",
                  ont="BP") #BP, MF, or CC
  #number of GO terms for
  require(DOSE)
  #can adjust it using scalebar using the ggplot setting 
  barplot(ego) + ggtitle("GO enrichment for MR220CL2 Stressed vs Non-Stressed (BP)")
# No results
  }
#No8 - EXPERIMENTAL
{
  #gene list from DE of SvsNS 
  diffexpressed <- 
    ifelse(NMR152vsMR220_et_df$logFC > 2 & NMR152vsMR220_et_df$PValue < 0.05, "UP",
           ifelse(NMR152vsMR220_et_df$logFC < -2 & NMR152vsMR220_et_df$PValue <0.05, "DOWN", 
                  "NS"))
  
  diff_gene_names_SvsNS <- NMR152vsMR220_et_df %>% 
    dplyr::filter(diffexpressed == "UP" | diffexpressed == "DOWN") %>% 
    unique()
  
  map_id <- RiceIDConvert(NMR152vsMR220_et_df$Geneid, fromType = "RAP", toType = "SYMBOL")
  
  ego <- enrichGO(map_id$SYMBOL,
                  OrgDb = organism,
                  keyType = "SYMBOL",
                  ont="BP") #BP, MF, or CC
  #number of GO terms for
  require(DOSE)
  #can adjust it using scalebar using the ggplot setting 
  barplot(ego) + ggtitle("GO enrichment for Stressed NMR152 VS MR220CL2 (BP)")
  #Note: no results for MF and CC
}

