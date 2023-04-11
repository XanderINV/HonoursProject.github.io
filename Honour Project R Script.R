#SCRIPTS FOR HONOURS PROJECT

########Downloading Data from TCGA##########  
library('DESeq2')
library('TCGAbiolinks')
library("TCGAbiolinks")
library("limma")
library("edgeR")
library("glmnet")
library("factoextra")
library("FactoMineR")
library("caret")
library("SummarizedExperiment")
library("gplots")
library("survival")
library("survminer")
library("RColorBrewer")
library("gProfileR")
library("genefilter")
library('airway')
library('tidyverse')
library('apeglm')

#Make a Query to Download the TCGA transcriptome data
query.expression <- GDCquery(project = "TCGA-BRCA",
                             data.category = "Transcriptome Profiling",
                             data.type = "Gene Expression Quantification",
                             barcode = c("TCGA-......................."))
GDCdownload(query.expression)

BRCA.Rnaseq.SE <- GDCprepare(query.expression)

#Create an RDS file that is easier to use in the future
saveRDS(object = BRCA.Rnaseq.SE,
        file = "BRCA_data.RDS",
        compress = FALSE)

setwd("R")
BRCA_data = readRDS(file = "BRCA_data.RDS")

#Obtaining the count matrix through assay function and the clinical data form the RSE BRCA_data
Matrix_Raw_Counts <- assay(BRCA_data, "unstranded")
clinical_data = colData(BRCA_data)


### FULL CLINICAL DATA SCRIPT ####

#Import the counts and clinical data (see Survival Analysis for changes made to this data before importing)
setwd("Lymph nodes predictor")
install.packages('readxl')
library(readxl)
setwd("Lymph nodes predictor/")

Clinical_data <- read_excel("Clinical_lymph.xlsx")
Clinical_data <- as.data.frame(Clinical_lymph)
View(Clinical_data) # n = 1044 (patients with <30 day follow up or time to death have been removed)

#Removing rows with no data or NX for lymph nodes from clinical and counts
Rows_to_remove <- c(which(Clinical_data == "NX", arr.ind=TRUE))
Rows_to_remove <- head(Rows_to_remove, -17)
Clinical_data_1 <- Clinical_data[!(row.names(Clinical_data) %in% Rows_to_remove),]

table(is.na(Clinical_data$ajcc_pathologic_n))
which(is.na(Clinical_data$ajcc_pathologic_n))
Clinical_data_2 <- Clinical_data_1[complete.cases(Clinical_data_1[,"ajcc_pathologic_n"]),]
table(is.na(Clinical_data_2$ajcc_pathologic_n)) # NA removed

#convert lymph node class to binary (1 = )
Clinical_data_2$lymph_class <- ifelse(Clinical_data_2$lymph_class == "N1", 1, 0)  #1 = spread to lymph nodes

#Quick observation of the data
library(DataExplorer)
plot_intro(Clinical_data_2)
plot_bar(Clinical_data_2)
plot_correlation(Clinical_data_2)

col_na <- colSums(is.na(Clinical_data_2))/nrow((Clinical_data_2))
col_na <- col_na[col_na > 0]
barplot(col_na,
        main = "Percentage of Missing Values by Column",
        xlab = "Columns",
        ylab = "Percentage of Missing Values",
        col = "red", las=2, cex.names = 0.7
)

col_only_na <- col_na[col_na == 1]
names(col_only_na) #These columns contain only NA values and should likely be removed (treatments, primary_site, disease type)

#Removing male samples
remove_male <- c(which(Clinical_data_2$gender == "male", arr.ind = TRUE))
Clinical_data_3 <- Clinical_data_2[!(row.names(Clinical_data_2) %in% remove_male),]
table(Clinical_data_3$gender)

#Removing near 0 variance variables
library(caret)
nzv <- nearZeroVar(Clinical_data_2, saveMetrics= TRUE)
nzv_index <- nearZeroVar(Clinical_data_2)
table(nzv$nzv)
Clinical_data_4 <- Clinical_data_3[,-nzv_index]

#Removing columns with useless ID data:
columns_to_remove <- c("sample_submitter_id","sample_id","pathology_report_uuid","diagnosis_id","exposure_id","demographic_id")
Clinical_data_5 <- Clinical_data_4[,!(colnames(Clinical_data_4) %in% columns_to_remove)]

#Remove days to last follow up and days to death since these values are represented together in the survival time column.
remove_daystodeath_followup <- c("days_to_death","days_to_last_follow_up")
Clinical_data_6 <- Clinical_data_5[,!(colnames(Clinical_data_5) %in% remove_daystodeath_followup)]

#Additional columns without relevant clinical data
Colums_irrelevant <- c("sample","submitter_id","bcr_patient_barcode","paper_patient","synchronous_malignancy","barcode","days_to_birth","age_at_index","paper_days_to_death",
                       "paper_days_to_last_followup","paper_days_to_birth","days_to_collection",
                       "year_of_diagnosis","ajcc_staging_system_edition","ajcc_pathologic_n",
                       "year_of_death","year_of_birth","paper_Included_in_previous_marker_papers","paper_vital_status")# synch malig has only no and not reported values
Clinical_data_7 <- Clinical_data_6[,!(colnames(Clinical_data_6) %in% Colums_irrelevant)]

Clinical_data_7[Clinical_data_7 == "not reported"] <- NA
Clinical_data_7[Clinical_data_7 == "NA"] <- NA
col_na_7 <- colSums(is.na(Clinical_data_7))
col_na_7
#converting categorical data into binary where possible
#Removing variables that are implicated in response or future information

Clinical_data_7$vital_status <- ifelse(Clinical_data_7$vital_status == "Dead", 1, 0)  # Dead = 1,  Alive = 0

library(writexl)
write_xlsx(Clinical_data_7,"..\\temporary_clinical.xlsx")




#REMOVAL OF OTHER VARIABLES AND PRE_PROCESSING
#Drop all the rows containing NA, leaves 833 patients
drop_colum <- c("initial_weight","oct_embedded","ajcc_pathologic_stage", "age_at_diagnosis","ajcc_pathologic_stage","primary_diagnosis","ajcc_pathologic_m", "vital_status", "paper_pathologic_stage", "Survival Time")  
Clinical_data_7 <- Clinical_data_7[,!colnames(Clinical_data_7) %in% drop_colum]




##################    TWO OPTIONS HERE      

#1 Is to remove all NA patients with NA values in any column and work with remainder
#2 Is to remove specific clinical variables such as "paper_lncRNA.Clusters","paper_BRCA_Pathology"
# which contain many NA values and have more samples.



# METHOD 1, REMOVING PATIENTS WITH ANY NA VALUES AND KEEPING ALL CLINICAL VARIABLES
Clinical_data_7_NAomit <- na.omit(Clinical_data_7)

#going through each columns and changing accordingly
table(Clinical_data_7_NAomit$ajcc_pathologic_t)
Clinical_data_7_NAomit$Tumor_Size <- ifelse(grepl("^T1", Clinical_data_7_NAomit$ajcc_pathologic_t), "T1",
                                            ifelse(grepl("^T2", Clinical_data_7_NAomit$ajcc_pathologic_t), "T2",
                                                   ifelse(grepl("^T3", Clinical_data_7_NAomit$ajcc_pathologic_t), "T3",
                                                          ifelse(grepl("^T4", Clinical_data_7_NAomit$ajcc_pathologic_t), "T4", NA))))
drop<- c("ajcc_pathologic_t")
Clinical_data_7_NAomit <- Clinical_data_7_NAomit[,!colnames(Clinical_data_7_NAomit) %in% drop]

Clinical_data_7_NAomit$prior_malignancy <- ifelse(Clinical_data_7_NAomit$prior_malignancy == "yes", 1, 0)  # yes = 1,  no = 0

#clean morphology create dummy variables
#Morphology done within EXCEL ,created dummy variables for each morphology too
write_xlsx(Clinical_data_7_NAomit,"..\\Metastasis Clinical Method 1.xlsx")  #File moved to deysnced files clinical folder
Clinical_data_7_NAomit <- read_excel("C:/Users/Xande/OneDrive - National University of Singapore/Project, unsynced files from R/survival/Clinical/5-year Survival Clinical Method 1.xlsx")
Clinical_data_7_NAomit <- Clinical_data_7_NAomit[,-3] #Remove Morphology column

#Dummy variables
Clinical_data_7_NAomit$White <- ifelse(Clinical_data_7_NAomit$race == "White",1,0)
Clinical_data_7_NAomit$Asian  <-  ifelse(Clinical_data_7_NAomit$race == "Asian",1,0)
Clinical_data_7_NAomit$Black <- ifelse(Clinical_data_7_NAomit$race == "black_or_African_American",1,0)
which(Clinical_data_7_NAomit$race == "american indian or alaska native") # need to remove row 416 since theres only one of these
Clinical_data_7_NAomit <- Clinical_data_7_NAomit[-270,] #removed row with alaska native
Clinical_data_7_NAomit <- Clinical_data_7_NAomit[,-3] #removed race column

Clinical_data_7_NAomit$ethnicity <- ifelse(Clinical_data_7_NAomit$ethnicity == "not hispanic or latino",1,0)

Clinical_data_7_NAomit$LumA <- ifelse(Clinical_data_7_NAomit$paper_BRCA_Subtype_PAM50 == "LumA",1,0)
Clinical_data_7_NAomit$LumB <- ifelse(Clinical_data_7_NAomit$paper_BRCA_Subtype_PAM50 == "LumB",1,0)
Clinical_data_7_NAomit$Normal_subtype <- ifelse(Clinical_data_7_NAomit$paper_BRCA_Subtype_PAM50 == "Normal",1,0)
Clinical_data_7_NAomit$Her2 <- ifelse(Clinical_data_7_NAomit$paper_BRCA_Subtype_PAM50 == "Her2",1,0)
Clinical_data_7_NAomit$Basal <- ifelse(Clinical_data_7_NAomit$paper_BRCA_Subtype_PAM50 == "Basal",1,0)
Clinical_data_7_NAomit <- Clinical_data_7_NAomit[,-6] #remove original subtype column

Clinical_data_7_NAomit$T1 <- ifelse(Clinical_data_7_NAomit$Tumor_Size == "T1", 1,0)
Clinical_data_7_NAomit$T2 <- ifelse(Clinical_data_7_NAomit$Tumor_Size == "T2", 1,0)
Clinical_data_7_NAomit$T3 <- ifelse(Clinical_data_7_NAomit$Tumor_Size == "T3", 1,0)
Clinical_data_7_NAomit$T4 <- ifelse(Clinical_data_7_NAomit$Tumor_Size == "T4", 1,0)
Clinical_data_7_NAomit <- Clinical_data_7_NAomit[,-16] #Remove Tumor Size orignial column

Clinical_data_7_NAomit$IDC <- ifelse(Clinical_data_7_NAomit$paper_BRCA_Pathology == "IDC",1,0)
Clinical_data_7_NAomit$ILC <- ifelse(Clinical_data_7_NAomit$paper_BRCA_Pathology == "ILC",1,0)
Clinical_data_7_NAomit$Mixed_P <- ifelse(Clinical_data_7_NAomit$paper_BRCA_Pathology == "Mixed",1,0)
Clinical_data_7_NAomit$Other_P <- ifelse(Clinical_data_7_NAomit$paper_BRCA_Pathology == "Other",1,0)
Clinical_data_7_NAomit <- Clinical_data_7_NAomit[,-5] #Remove Tumor Size orignial column

Clinical_data_7_NAomit$CNV_C1 <- ifelse(Clinical_data_7_NAomit$paper_CNV.Clusters == "C1",1,0)
Clinical_data_7_NAomit$CNV_C2 <- ifelse(Clinical_data_7_NAomit$paper_CNV.Clusters == "C2",1,0)
Clinical_data_7_NAomit$CNV_C3 <- ifelse(Clinical_data_7_NAomit$paper_CNV.Clusters == "C3",1,0)
Clinical_data_7_NAomit$CNV_C4 <- ifelse(Clinical_data_7_NAomit$paper_CNV.Clusters == "C4",1,0)
Clinical_data_7_NAomit$CNV_C5 <- ifelse(Clinical_data_7_NAomit$paper_CNV.Clusters == "C5",1,0)
Clinical_data_7_NAomit$CNV_C6 <- ifelse(Clinical_data_7_NAomit$paper_CNV.Clusters == "C6",1,0)
Clinical_data_7_NAomit <- Clinical_data_7_NAomit[,-5] #Remove Tumor Size orignial column

Clinical_data_7_NAomit$Mutation_C1 <- ifelse(Clinical_data_7_NAomit$paper_Mutation.Clusters == "C1",1,0)
Clinical_data_7_NAomit$Mutation_C2 <- ifelse(Clinical_data_7_NAomit$paper_Mutation.Clusters == "C2",1,0)
Clinical_data_7_NAomit$Mutation_C3 <- ifelse(Clinical_data_7_NAomit$paper_Mutation.Clusters == "C3",1,0)
Clinical_data_7_NAomit$Mutation_C4 <- ifelse(Clinical_data_7_NAomit$paper_Mutation.Clusters == "C4",1,0)
Clinical_data_7_NAomit$Mutation_C5 <- ifelse(Clinical_data_7_NAomit$paper_Mutation.Clusters == "C5",1,0)
Clinical_data_7_NAomit$Mutation_C6 <- ifelse(Clinical_data_7_NAomit$paper_Mutation.Clusters == "C6",1,0)
Clinical_data_7_NAomit$Mutation_C7 <- ifelse(Clinical_data_7_NAomit$paper_Mutation.Clusters == "C7",1,0)
Clinical_data_7_NAomit$Mutation_C8 <- ifelse(Clinical_data_7_NAomit$paper_Mutation.Clusters == "C8",1,0)
Clinical_data_7_NAomit$Mutation_C9 <- ifelse(Clinical_data_7_NAomit$paper_Mutation.Clusters == "C9",1,0)
Clinical_data_7_NAomit$Mutation_C10 <- ifelse(Clinical_data_7_NAomit$paper_Mutation.Clusters == "C10",1,0)
Clinical_data_7_NAomit <- Clinical_data_7_NAomit[,-5] #Remove Mutation column


#Start using loops for this
# convert "paper_DNA.Methylation.Clusters" to character vector
clusters_column <- Clinical_data_7_NAomit$paper_DNA.Methylation.Clusters

for (i in 1:6) {
  # create a new column with the appropriate label
  new_column_name <- paste0("Methylation_C", i)
  Clinical_data_7_NAomit[[new_column_name]] <- 0
  
  # set the appropriate values to 1
  Clinical_data_7_NAomit[[new_column_name]][clusters_column == paste0("C", i)] <- 1
}
Clinical_data_7_NAomit <- Clinical_data_7_NAomit[,-5] #Remove Methylation orignial column

#mRNA Clusters
cluster_labels <- c("C1", "C2", "C3", "C4", "C7")

for (i in seq_along(cluster_labels)) {
  # create a new column with the appropriate label
  new_column_name <- paste0("mRNA_", cluster_labels[i])
  Clinical_data_7_NAomit[[new_column_name]] <- 0
  
  # set the appropriate values to 1
  Clinical_data_7_NAomit[[new_column_name]] <- ifelse(Clinical_data_7_NAomit$paper_mRNA.Clusters == cluster_labels[i], 1, 0)
}
Clinical_data_7_NAomit <- Clinical_data_7_NAomit[,-5] #Remove mRNA orignial column

#miRNA Clusters
cluster_labels <- c("C2", "C3", "C5","C6", "C7")

for (i in seq_along(cluster_labels)) {
  # create a new column with the appropriate label
  new_column_name <- paste0("miRNA_", cluster_labels[i])
  Clinical_data_7_NAomit[[new_column_name]] <- 0
  
  # set the appropriate values to 1
  Clinical_data_7_NAomit[[new_column_name]] <- ifelse(Clinical_data_7_NAomit$paper_miRNA.Clusters == cluster_labels[i], 1, 0)
}
Clinical_data_7_NAomit <- Clinical_data_7_NAomit[,-5] #Remove miRNA orignal column


#lncRNA Clusters
cluster_labels <- c("C1", "C2", "C3", "C4", "C5","C6")
for (i in seq_along(cluster_labels)) {
  # create a new column with the appropriate label
  new_column_name <- paste0("lncRNA_", cluster_labels[i])
  Clinical_data_7_NAomit[[new_column_name]] <- 0
  
  # set the appropriate values to 1
  Clinical_data_7_NAomit[[new_column_name]] <- ifelse(Clinical_data_7_NAomit$paper_lncRNA.Clusters == cluster_labels[i], 1, 0)
}
Clinical_data_7_NAomit <- Clinical_data_7_NAomit[,-5] #Remove lncRNA original column

#Protein Clusters
cluster_labels <- c("C1", "C2", "C4", "C5")
for (i in seq_along(cluster_labels)) {
  # create a new column with the appropriate label
  new_column_name <- paste0("Protein_cluster_", cluster_labels[i])
  Clinical_data_7_NAomit[[new_column_name]] <- 0
  
  # set the appropriate values to 1
  Clinical_data_7_NAomit[[new_column_name]] <- ifelse(Clinical_data_7_NAomit$paper_Protein.Clusters == cluster_labels[i], 1, 0)
}
Clinical_data_7_NAomit <- Clinical_data_7_NAomit[,-5] #Remove protein original column

#Paradigm Clusters
cluster_labels <- c("C2", "C4", "C5", "C6", "C8")
for (i in seq_along(cluster_labels)) {
  # create a new column with the appropriate label
  new_column_name <- paste0("Paradigm_cluster_", cluster_labels[i])
  Clinical_data_7_NAomit[[new_column_name]] <- 0
  
  # set the appropriate values to 1
  Clinical_data_7_NAomit[[new_column_name]] <- ifelse(Clinical_data_7_NAomit$paper_PARADIGM.Clusters == cluster_labels[i], 1, 0)
}
Clinical_data_7_NAomit <- Clinical_data_7_NAomit[,-5] #Remove Paradigm Cluster original column

#Pan_Gyn Clusters
cluster_labels <- c("C1", "C2", "C3", "C4", "C5")
for (i in seq_along(cluster_labels)) {
  # create a new column with the appropriate label
  new_column_name <- paste0("Pan_Gyn_Cluster_", cluster_labels[i])
  Clinical_data_7_NAomit[[new_column_name]] <- 0
  
  # set the appropriate values to 1
  Clinical_data_7_NAomit[[new_column_name]] <- ifelse(Clinical_data_7_NAomit$paper_Pan.Gyn.Clusters== cluster_labels[i], 1, 0)
}
Clinical_data_7_NAomit <- Clinical_data_7_NAomit[,-5] #Remove Pan_Gyn original column


#rownames into patients IDs
Clinical_data_7_NAomit <- as.data.frame(Clinical_data_7_NAomit)
rownames(Clinical_data_7_NAomit) <- Clinical_data_7_NAomit[,1]
Clinical_data_7_NAomit <- Clinical_data_7_NAomit[,-1]

#Convert Factors variables into factors (except age)
Clinical_data_7_NAomit[, -3] <- lapply(Clinical_data_7_NAomit[, -3], as.factor)
sapply(Clinical_data_7_NAomit, class)

#SAVE FINAL PATIENT INFORMATION Data frame IN EXCEL BEFORE PREDICTION
Clinical_data_7_NAomit$patients_IDs <- rownames(Clinical_data_7_NAomit)
write_xlsx(Clinical_data_7_NAomit,"..\\Clinical Info after pre-processing method 1.xlsx")  #File moved to deysnced files clinical folder
Clinical_data_7_NAomit <- Clinical_data_7_NAomit[,-79]



###########################                         SPLIT -> RANk -> TRAIN               
#Split the data
library(caret)
set.seed(42)
trainIndex <- createDataPartition(Clinical_data_7_NAomit$lymph_class, p = .7,
                                  list = FALSE,
                                  times = 1)
Train_m1 <- Clinical_data_7_NAomit[ trainIndex,]
Test_m1 <- Clinical_data_7_NAomit[-trainIndex,]



##############                 SPLIT CLINICAL SAMPLES INTO SUBSETS AND RANKING                  


#Robust RFE
#Creating a table for subsets of samples which ranks feature importance scores. Then we can get the total feature importance score for each feature across all the subsets and generate ROC curves.
library(caret)
set.seed(42)

# create empty list to store results
results_importance <- list()

# create vector of percentages
percentages <- seq(0.8, 0.1, by = -0.1)


lymph_0 <- Train_m1[Train_m1$lymph_class == 0,]
lymph_1 <- Train_m1[Train_m1$lymph_class == 1,]

library(progress)
pb <-progress_bar$new(total = length(percentages))


# loop over each percentage
for (p in percentages) {
  set.seed(42)
  print(p)
  
  n_0 <- round(nrow(lymph_0)*0.7)
  n_1 <- round(nrow(lymph_1)*0.7)
  
  
  # randomly sample p% of patients
  sample_rows_0 <- sample(nrow(lymph_0), n_0)
  sampled_data_0 <- lymph_0[sample_rows_0,]
  
  sample_rows_1 <- sample(nrow(lymph_1), n_1)
  sampled_data_1 <- lymph_1[sample_rows_1,]
  
  sampled_data <- rbind(sampled_data_0,sampled_data_1)
  
  # set up control function for rfe
  ctrl <- rfeControl(
    method = "repeatedcv",
    repeats = 5,
    verbose = FALSE,
    returnResamp = "final",
    functions = rfFuncs,
  )
  
  # perform rfe
  rfe_result <- rfe(
    x = sampled_data[, -1],
    y = sampled_data[, 1],
    sizes = c(1),
    rfeControl = ctrl
  )
  
  # get feature importance scores
  importance <- varImp(rfe_result, scale = FALSE)
  
  # create data frame with feature importance scores
  importance_df <- data.frame(importance[, 1], row.names = rownames(importance))
  colnames(importance_df) <- c("Feature Importance")
  
  # store importance data frame in list of data frames
  results_importance[[as.character(0.1*100)]] <- importance_df
}

# print results
print(results_importance)

#GENERATING THE TABLE WITH IMPORTANCE SCORES FOR EACH Variable ACROSS ALL SUBSETS
# create empty data frame with one column for each subset
install.packages("tidyverse")
library(tidyverse)

# Set rownames as a regular column in each data frame
results_importance <- lapply(results_importance, function(x) {
  rownames_to_column(x, "Variable")  
})


install.packages("purrr")
library(purrr)
library(dplyr)

# Add a new column to each data frame in the list with the same name as the data frame
names(results_importance) <- paste0("df_", seq_along(results_importance))

# Combine all data frames into one
combined_df <- reduce(results_importance, full_join, by = "Variable")

new_colnames <- c("Variable", "80%", "70%", "60%", "50%", "40%", "30%", "20%", "10%")

colnames(combined_df) <- new_colnames

library(dplyr)

combined_df$average <- rowMeans(combined_df[, 2:9], na.rm = TRUE)

combined_df <- combined_df %>% arrange(desc(average))





### CREATE FEATURE IMPORTANCE PLOT

library(ggplot2)

# create a new data frame with the desired ordering
combined_ordered <- combined_df[order(combined_df$average, decreasing = TRUE), ]

# plot the data using ggplot2
ggplot(combined_ordered, aes(x = average, y = reorder(Variable, average))) +
  geom_bar(stat = "identity", fill = "steelblue") +
  xlab("Average Score") +
  ylab("Variable") +
  theme(axis.text.y = element_text(size = 10, hjust = 1)) +
  ggtitle("Variable Importance Scores")

###### REDO USING CROSS VALIDATION AND (ONLY RF FOR NOW) PREDICTING IN TRAIN SET TO OBTAIN TRAIN AUC CURVE



####                              RF MODEL TRAINING
library(randomForest)
library(pROC)
library(caret)

# Define number of variables to use
num_var <- 77

# Initialize empty data frames to store results
results_df <- data.frame(matrix(ncol=4, nrow=num_var))
colnames(results_df) <- c("Variables", "AUC_train", "Sensitivity_train", "Specificity_train")

results_RF_Folds <- data.frame(matrix(ncol=16, nrow=num_var))
colnames(results_RF_Folds) <- c("Variables", 
                                paste(rep(c("AUC", "Sensitivity", "Specificity"), each=5), 1:5, sep="_"))

# Loop over number of variables to include
set.seed(42)
for (i in 1:num_var) {
  # Select top i genes
  selected_var <- combined_df$Variable[1:i]
  
  # Create Train data frame with response variable and selected genes
  Training <- Train_m1[, c("lymph_class", selected_var)]
  
  # Set up K-fold cross-validation
  set.seed(42)
  K <- 5
  folds <- createFolds(Training$lymph_class, K, list = TRUE)
  
  # Initialize empty vectors to store predictions and labels for each fold
  train_predictions <- numeric(nrow(Training))
  train_labels <- numeric(nrow(Training))
  auc_vec <- numeric(K)
  sens_vec <- numeric(K)
  spec_vec <- numeric(K)
  
  # Train and predict using K-fold cross-validation
  for (k in 1:K) {
    # Split data into training and validation sets
    train_indices <- unlist(folds[-k])
    valid_indices <- unlist(folds[k])
    
    train_set <- Training[train_indices, ]
    valid_set <- Training[valid_indices, ]
    
    # Train random forest model on training set
    rf_model <- randomForest(lymph_class ~ ., data=train_set, ntree=500, importance=TRUE)
    
    
    # Make predictions on training set for AUC, sensitivity, and specificity
    train_prediction_probs <- predict(rf_model, valid_set, type = "prob")
    train_prediction_df <- data.frame(prob = train_prediction_probs[,2], class = valid_set$lymph_class)
    
    # Calculate AUC, sensitivity, and specificity of ROC curve for validation and training sets
    auc_vec[k] <- roc(train_prediction_df$class, train_prediction_df$prob)$auc
    sens_vec[k] <- sum(train_prediction_probs[valid_set$lymph_class==1,2] >= 0.5)/sum(valid_set$lymph_class==1)
    spec_vec[k] <- sum(train_prediction_probs[valid_set$lymph_class==0,2] < 0.5)/sum(valid_set$lymph_class==0)
  }
  
  # Calculate average AUC, sensitivity, and specificity across all folds
  avg_auc_train <- mean(auc_vec)
  avg_sens_train <- mean(sens_vec)
  avg_spec_train <- mean(spec_vec)
  
  # Store results in results_df
  results_df[i, "Variables"] <- i
  results_df[i, "AUC_train"] <- avg_auc_train
  results_df[i, "Sensitivity_train"] <- avg_sens_train
  results_df[i, "Specificity_train"] <- avg_spec_train
  
  # Store AUC, sensitivity, and specificity for each fold in results_RF_Folds
  results_RF_Folds[i, c(paste0("AUC_", 1:5))] <- auc_vec
  results_RF_Folds[i, c(paste0("Sensitivity_", 1:5))] <- sens_vec
  results_RF_Folds[i, c(paste0("Specificity_", 1:5))] <- spec_vec
}

#plot
results_RF_Folds$Variables <- seq(1, 77)

library(ggplot2)

results_df$average = (results_df$AUC_train + results_df$Sensitivity_train + results_df$Specificity_train)/3
#6 clinical variable is best but 23 is more balanced

# Reshape the results_RF_Folds dataframe to long format for easier plotting
results_RF_Folds_long <- reshape2::melt(results_RF_Folds, id.vars=c("Variables"), variable.name="Fold", value.name="Value")

# Plot the results using ggplot
ggplot() +
  # Add the AUC values for each fold as points
  geom_point(data=results_RF_Folds_long, aes(x=Variables, y=Value, color=Fold)) +
  # Add the average AUC for each number of variables as a line
  geom_line(data=results_df, aes(x=Variables, y=AUC_train)) +
  # Add axis labels and title
  labs(x="Number of clinical variables", y="AUC", title="Random Forest AUC vs Number of Variables")


#Individual plots

#AUC
library(ggplot2)

# create a long format of results_RF_Folds
results_RF_Folds_long_AUC <- reshape2::melt(results_RF_Folds, id.vars = "Variables", 
                                            variable.name = "Fold", value.name = "AUC")
# filter only the AUC columns
results_RF_Folds_long_AUC <- results_RF_Folds_long_AUC %>%
  filter(str_detect(Fold, "AUC"))
# convert Variables column to numeric
results_RF_Folds_long_AUC$Variables <- as.numeric(results_RF_Folds_long_AUC$Variables)
# calculate the mean AUC for each Variables value
results_RF_Folds_mean <- aggregate(AUC ~ Variables, data = results_RF_Folds_long_AUC, mean)
# plot the AUC curve
ggplot(results_RF_Folds_long_AUC, aes(x = Variables, y = AUC)) +
  geom_point(aes(color = Fold), alpha = 0.5) +
  geom_line(data = results_RF_Folds_mean, aes(x = Variables, y = AUC), size = 1.2) +
  labs(x = "Number of Clinical Variables", y = "AUC") +
  theme_classic() +
  theme(plot.title = element_blank(),
        axis.text.x = element_text(size = 24, color = "black", face = "bold"),
        axis.text.y = element_text(size = 24, color = "black", face = "bold"),
        axis.title.x = element_text(size = 14, color = "black", face = "bold"),
        axis.title.y = element_text(size = 14, color = "black", face = "bold"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(size = 1)) +
  guides(color = FALSE)+
  scale_x_continuous(expand = c(0, 0))
geom_vline(xintercept = 14, linetype = "dotted", color = "blue") +
  geom_segment(aes(x = 14, y = 0.4, xend = 14, yend = 0.8), linetype = "dotted", color = "blue") +
  annotate("text", x = 22, y = 0.45, label = "X = 14", color = "blue", vjust = -1, hjust = 1)




#Sensitivity
# create a long format of results_RF_Folds
results_RF_Folds_long_Sens <- reshape2::melt(results_RF_Folds, id.vars = "Variables", 
                                             variable.name = "Fold", value.name = "Sensitivity")
# filter only the Sensitivity columns
results_RF_Folds_long_Sens <- results_RF_Folds_long_Sens %>%
  filter(str_detect(Fold, "Sensitivity"))
# convert Variables column to numeric
results_RF_Folds_long_Sens$Variables <- as.numeric(results_RF_Folds_long_Sens$Variables)
# calculate the mean Sensitivity for each Variables value
results_RF_Folds_mean_Sens <- aggregate(Sensitivity ~ Variables, data = results_RF_Folds_long_Sens, mean)
# plot the Sensitivity curve
ggplot(results_RF_Folds_long_Sens, aes(x = Variables, y = Sensitivity)) +
  geom_point(aes(color = Fold), alpha = 0.5) +
  geom_line(data = results_RF_Folds_mean_Sens, aes(x = Variables, y = Sensitivity), size = 1.2) +
  labs(x = "Number of Clinical Variables", y = "Sensitivity") +
  theme_classic() +
  theme(plot.title = element_blank(),
        axis.text.x = element_text(size = 24, color = "black", face = "bold"),
        axis.text.y = element_text(size = 24, color = "black", face = "bold"),
        axis.title.x = element_text(size = 14, color = "black", face = "bold"),
        axis.title.y = element_text(size = 14, color = "black", face = "bold"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(size = 1)) +
  guides(color = FALSE)+
  scale_x_continuous(expand = c(0, 0))
geom_vline(xintercept = 14, linetype = "dotted", color = "blue") +
  geom_segment(aes(x = 14, y = 0.4, xend = 14, yend = 0.8), linetype = "dotted", color = "blue")+
  annotate("text", x = 22, y = 0.25, label = "X = 14", color = "blue", vjust = -1, hjust = 1) 


#Specificity
# create a long format of results_RF_Folds
results_RF_Folds_long_Spec <- reshape2::melt(results_RF_Folds, id.vars = "Variables", 
                                             variable.name = "Fold", value.name = "Specificity")
# filter only the Specificity columns
results_RF_Folds_long_Spec <- results_RF_Folds_long_Spec %>%
  filter(str_detect(Fold, "Specificity"))
# convert Variables column to numeric
results_RF_Folds_long_Spec$Variables <- as.numeric(results_RF_Folds_long_Spec$Variables)
# calculate the mean Specificity for each Variables value
results_RF_Folds_mean_Spec <- aggregate(Specificity ~ Variables, data = results_RF_Folds_long_Spec, mean)
# plot the Specificity curve
ggplot(results_RF_Folds_long_Spec, aes(x = Variables, y = Specificity)) +
  geom_point(aes(color = Fold), alpha = 0.5) +
  geom_line(data = results_RF_Folds_mean_Spec, aes(x = Variables, y = Specificity), size = 1.2) +
  labs(x = "Number of Clinical Variables", y = "Specificity") +
  theme_classic() +
  theme(plot.title = element_blank(),
        axis.text.x = element_text(size = 24, color = "black", face = "bold"),
        axis.text.y = element_text(size = 24, color = "black", face = "bold"),
        axis.title.x = element_text(size = 14, color = "black", face = "bold"),
        axis.title.y = element_text(size = 14, color = "black", face = "bold"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(size = 1)) +
  guides(color = FALSE)+
  scale_x_continuous(expand = c(0, 0))







############################                   SVM TRAINING MODEL               ##################################33

#SVM DOESNT SEEM TO WORK WITH PREDICTORS THAT HAVE 0 IMPORTANCE SCORE, WHAT TO DO? REMOVE OR DIFFERENT MODEL?
library(e1071)
library(pROC)
library(caret)

#CREATE DATA FRAME FROM COMBINED DF THAT HAS ONLY FEATURES WITH IMPORTANCE SCORES (NON 0)
combined_df_SVM <- combined_df[combined_df$average != 0, ]


# Define number of variables to use
num_var <- 65

# Initialize empty data frames to store results
results_SVM_df <- data.frame(matrix(ncol=4, nrow=num_var))
colnames(results_SVM_df) <- c("Variables", "AUC_train", "Sensitivity_train", "Specificity_train")

results_SVM_Folds_scaled <- data.frame(matrix(ncol=16, nrow=num_var))
colnames(results_SVM_Folds_scaled) <- c("Variables", 
                                        paste(rep(c("AUC", "Sensitivity", "Specificity"), each=5), 1:5, sep="_"))

# Loop over number of variables to include
set.seed(42)
for (i in 1:num_var) {
  # Select top i genes
  selected_var_scaled <- combined_df_SVM$Variable[1:i]
  
  # Create Train data frame with response variable and selected genes
  Training <- Train_m1[, c("lymph_class", selected_var_scaled)]
  
  # Set up K-fold cross-validation
  set.seed(42)
  K <- 5
  folds <- createFolds(Training$lymph_class, K, list = TRUE)
  
  # Initialize empty vectors to store predictions and labels for each fold
  train_predictions <- numeric(nrow(Training))
  train_labels <- numeric(nrow(Training))
  auc_vec <- numeric(K)
  sens_vec <- numeric(K)
  spec_vec <- numeric(K)
  
  # Train and predict using K-fold cross-validation
  for (k in 1:K) {
    # Split data into training and validation sets
    train_indices <- unlist(folds[-k])
    valid_indices <- unlist(folds[k])
    
    train_set <- Training[train_indices, ]
    valid_set <- Training[valid_indices, ]
    
    
    # Train SVM model on training set
    svm_model <- svm(lymph_class ~ ., data = train_set, kernel = "sigmoid", 
                     type = "C-classification", probability = TRUE,
                     maxiter = 10000)
    
    # Make predictions on training set for AUC, sensitivity, and specificity
    train_prediction_probs <- predict(svm_model, valid_set, type = "prob")
    train_prediction_probs <- as.data.frame(train_prediction_probs)
    train_prediction_df <- data.frame(prob = train_prediction_probs[,1], class = valid_set$lymph_class)
    
    # Calculate AUC, sensitivity, and specificity of ROC curve for validation and training sets
    auc_vec[k] <- roc(as.numeric(train_prediction_df$class),as.numeric(train_prediction_df$prob))$auc
    sens_vec[k] <- sum(train_prediction_df$prob == "1" & train_prediction_df$class == "1")/sum(train_prediction_df$class==1)
    spec_vec[k] <- sum(train_prediction_df$prob == "0" & train_prediction_df$class == "0")/sum(train_prediction_df$class==0)
  }
  
  # Calculate average AUC, sensitivity, and specificity across all folds
  avg_auc_train <- mean(auc_vec)
  avg_sens_train <- mean(sens_vec)
  avg_spec_train <- mean(spec_vec)
  
  # Store results in results_df
  results_SVM_df[i, "Variables"] <- i
  results_SVM_df[i, "AUC_train"] <- avg_auc_train
  results_SVM_df[i, "Sensitivity_train"] <- avg_sens_train
  results_SVM_df[i, "Specificity_train"] <- avg_spec_train
  
  # Store AUC, sensitivity, and specificity for each fold in results_RF_Folds
  results_SVM_Folds_scaled[i, c(paste0("AUC_", 1:5))] <- auc_vec
  results_SVM_Folds_scaled[i, c(paste0("Sensitivity_", 1:5))] <- sens_vec
  results_SVM_Folds_scaled[i, c(paste0("Specificity_", 1:5))] <- spec_vec
}


#plot
results_SVM_Folds_scaled$Variables <- seq(1, num_var)


# Reshape the results_RF_Folds dataframe to long format for easier plotting
results_SVM_Folds_expanded <- reshape2::melt(results_SVM_Folds_scaled, id.vars=c("Variables"), variable.name="Fold", value.name="Value")

# Plot the results using ggplot
ggplot() +
  # Add the AUC values for each fold as points
  geom_point(data=results_SVM_Folds_expanded, aes(x=Variables, y=Value, color=Fold), alpha = 0.5) +
  # Add the average AUC for each number of variables as a line
  geom_line(data=results_SVM_df, aes(x=Variables, y=AUC_train), cex = 1.5) +
  # Add axis labels and title
  labs(x="Number of variables", y="AUC", title="SVM AUC vs Number of Variables")


#AUC
# create a long format of results_RF_Folds
results_SVM_Folds_expanded_AUC <- reshape2::melt(results_SVM_Folds_scaled, id.vars = "Variables", 
                                                 variable.name = "Fold", value.name = "AUC")
# filter only the AUC columns
results_SVM_Folds_expanded_AUC <- results_SVM_Folds_expanded_AUC %>%
  filter(str_detect(Fold, "AUC"))
# convert Variables column to numeric
results_SVM_Folds_expanded_AUC$Variables <- as.numeric(results_SVM_Folds_expanded_AUC$Variables)
# calculate the mean AUC for each Variables value
results_SVM_Folds_mean_scaled <- aggregate(AUC ~ Variables, data = results_SVM_Folds_expanded_AUC, mean)
# plot the AUC curve
ggplot(results_SVM_Folds_expanded_AUC, aes(x = Variables, y = AUC)) +
  geom_point(aes(color = Fold), alpha = 0.5) +
  geom_line(data = results_SVM_Folds_mean_scaled, aes(x = Variables, y = AUC), size = 1.2) +
  labs(x = "Number of Variables", y = "AUC") +
  theme_classic() +
  theme(plot.title = element_blank(),
        axis.text.x = element_text(size = 24, color = "black", face = "bold"),
        axis.text.y = element_text(size = 24, color = "black", face = "bold"),
        axis.title.x = element_text(size = 14, color = "black", face = "bold"),
        axis.title.y = element_text(size = 14, color = "black", face = "bold"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(size = 1)) +
  guides(color = FALSE) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(limits = c(0.4, 0.8), expand = c(0, 0))
#leave out for now
geom_vline(xintercept = 54, linetype = "dotted", color = "red") +
  geom_segment(aes(x = 54, y = 0.5, xend = 54, yend = 0.75), linetype = "dotted", color = "red")+
  annotate("text", x = 69, y = 0.5, label = "X = 54", color = "red", vjust = -1, hjust = 1)


# Specificity
# create a long format of results_SVM_Folds
results_SVM_Folds_expanded_spec <- reshape2::melt(results_SVM_Folds_scaled, id.vars = "Variables", 
                                                  variable.name = "Fold", value.name = "Specificity")
# filter only the Specificity columns
results_SVM_Folds_expanded_spec <- results_SVM_Folds_expanded_spec %>%
  filter(str_detect(Fold, "Specificity"))
# convert Variables column to numeric
results_SVM_Folds_expanded_spec$Variables <- as.numeric(results_SVM_Folds_expanded_spec$Variables)
# calculate the mean Specificity for each Variables value
results_SVM_Folds_mean_spec <- aggregate(Specificity ~ Variables, data = results_SVM_Folds_expanded_spec, mean)
# plot the Specificity curve
ggplot(results_SVM_Folds_expanded_spec, aes(x = Variables, y = Specificity)) +
  geom_point(aes(color = Fold), alpha = 0.5) +
  geom_line(data = results_SVM_Folds_mean_spec, aes(x = Variables, y = Specificity), size = 1.2) +
  labs(x = "Number of Variables", y = "Specificity") +
  theme_classic() +
  theme(plot.title = element_blank(),
        axis.text.x = element_text(size = 24, color = "black", face = "bold"),
        axis.text.y = element_text(size = 24, color = "black", face = "bold"),
        axis.title.x = element_text(size = 14, color = "black", face = "bold"),
        axis.title.y = element_text(size = 14, color = "black", face = "bold"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(size = 1)) +
  guides(color = FALSE)+
  scale_x_continuous(expand = c(0, 0))
#leave out for now
geom_vline(xintercept = 54, linetype = "dotted", color = "red") +
  geom_segment(aes(x = 54, y = 0.5, xend = 54, yend = 0.7), linetype = "dotted", color = "red") +
  annotate("text", x = 69, y = 0.5, label = "X = 54", color = "red", vjust = -1, hjust = 1)



# Sensitivity
# create a long format of results_SVM_Folds
results_SVM_Folds_expanded_sens <- reshape2::melt(results_SVM_Folds_scaled, id.vars = "Variables", 
                                                  variable.name = "Fold", value.name = "Sensitivity")
# filter only the Sensitivity columns
results_SVM_Folds_expanded_sens <- results_SVM_Folds_expanded_sens %>%
  filter(str_detect(Fold, "Sensitivity"))
# convert Variables column to numeric
results_SVM_Folds_expanded_sens$Variables <- as.numeric(results_SVM_Folds_expanded_sens$Variables)
# calculate the mean Sensitivity for each Variables value
results_SVM_Folds_mean_sens <- aggregate(Sensitivity ~ Variables, data = results_SVM_Folds_expanded_sens, mean)
# plot the Sensitivity curve
ggplot(results_SVM_Folds_expanded_sens, aes(x = Variables, y = Sensitivity)) +
  geom_point(aes(color = Fold), alpha = 0.5) +
  geom_line(data = results_SVM_Folds_mean_sens, aes(x = Variables, y = Sensitivity), size = 1.2) +
  labs(x = "Number of Variables", y = "Sensitivity") +
  theme_classic() +
  theme(plot.title = element_blank(),
        axis.text.x = element_text(size = 24, color = "black", face = "bold"),
        axis.text.y = element_text(size = 24, color = "black", face = "bold"),
        axis.title.x = element_text(size = 14, color = "black", face = "bold"),
        axis.title.y = element_text(size = 14, color = "black", face = "bold"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(size = 1)) +
  guides(color = FALSE)+
  scale_x_continuous(expand = c(0, 0))




#Leave out for now
geom_vline(xintercept = 54, linetype = "dotted", color = "red") +
  geom_segment(aes(x = 54, y = 0.5, xend = 54, yend = 0.7), linetype = "dotted", color = "red") +
  annotate("text", x = 69, y = 0.5, label = "X = 54", color = "red", vjust = -1, hjust = 1)






#DECIDE BEST SUBSET
results_SVM_df$Average <- (results_SVM_df$AUC_train + results_SVM_df$Sensitivity_train + results_SVM_df$Specificity_train)/3







########################33                         XGBOOST TRAINING MODEL             ###############################3

#Load required packages
library(xgboost)
library(pROC)

# Define number of variables to use
num_var <- 77

# Initialize empty data frames to store results
results_df_XGB <- data.frame(matrix(ncol=4, nrow=num_var))
colnames(results_df_XGB) <- c("Variables", "AUC_train", "Sensitivity_train", "Specificity_train")

results_XGB_Folds <- data.frame(matrix(ncol=16, nrow=num_var))
colnames(results_XGB_Folds) <- c("Variables", 
                                 paste(rep(c("AUC", "Sensitivity", "Specificity"), each=5), 1:5, sep="_"))

# Loop over number of variables to include
set.seed(42)
for (i in 1:num_var) {
  # Select top i genes
  selected_var_scaled <- combined_df$Variable[1:i]
  
  # Create Train data frame with response variable and selected genes
  Training <- Train_m1[, c("lymph_class", selected_var_scaled)]
  
  # Set up K-fold cross-validation
  set.seed(42)
  K <- 5
  folds <- createFolds(Training$lymph_class, K, list = TRUE)
  
  # Initialize empty vectors to store predictions and labels for each fold
  train_predictions <- numeric(nrow(Training))
  train_labels <- numeric(nrow(Training))
  auc_vec <- numeric(K)
  sens_vec <- numeric(K)
  spec_vec <- numeric(K)
  
  # Train and predict using K-fold cross-validation
  for (k in 1:K) {
    # Split data into training and validation sets
    train_indices <- unlist(folds[-k])
    valid_indices <- unlist(folds[k])
    
    train_set <- Training[train_indices, ]
    valid_set <- Training[valid_indices, ]
    
    # Convert factor variables to numerical values (1 or 0)
    train_set[,-1] <- as.numeric(ifelse(train_set[,-1] == 1, 1, 0))
    valid_set[,-1] <- as.numeric(ifelse(valid_set[,-1] == 1, 1, 0))
    
    train_set$lymph_class <- ifelse(train_set$lymph_class == 1, 1, 0)
    valid_set$lymph_class <- ifelse(valid_set$lymph_class == 1, 1, 0)
    
    # Train XGBoost model on training set
    xgb_model <- xgboost(data = as.matrix(train_set[,-1]), label = train_set$lymph_class, nrounds = 100, objective = "binary:logistic")
    
    # Make predictions on validation set for AUC, sensitivity, and specificity
    
    train_prediction_probs <- predict(xgb_model, as.matrix(valid_set[,-1]))
    train_prediction_df <- data.frame(prob = train_prediction_probs, class = valid_set$lymph_class)
    
    # Calculate AUC, sensitivity, and specificity of ROC curve for validation and training sets
    auc_vec[k] <- roc(train_prediction_df$class,train_prediction_df$prob)$auc
    sens_vec[k] <- sum(train_prediction_probs[valid_set$lymph_class==1] >= 0.5)/sum(valid_set$lymph_class==1)
    spec_vec[k] <- sum(train_prediction_probs[valid_set$lymph_class==0] < 0.5)/sum(valid_set$lymph_class==0)
  }
  
  # Calculate average AUC, sensitivity, and specificity across all folds
  avg_auc_train <- mean(auc_vec)
  avg_sens_train <- mean(sens_vec)
  avg_spec_train <- mean(spec_vec)
  
  # Store results in results_df
  results_df_XGB[i, "Variables"] <- i
  results_df_XGB[i, "AUC_train"] <- avg_auc_train
  results_df_XGB[i, "Sensitivity_train"] <- avg_sens_train
  results_df_XGB[i, "Specificity_train"] <- avg_spec_train
  
  # Store AUC, sensitivity, and specificity for each fold in results_RF_Folds
  results_XGB_Folds[i, c(paste0("AUC_", 1:5))] <- auc_vec
  results_XGB_Folds[i, c(paste0("Sensitivity_", 1:5))] <- sens_vec
  results_XGB_Folds[i, c(paste0("Specificity_", 1:5))] <- spec_vec
}



#plot
results_XGB_Folds$Variables <- seq(1, num_var)


# Reshape the results_RF_Folds dataframe to long format for easier plotting
results_XGB_Folds_expanded <- reshape2::melt(results_XGB_Folds, id.vars=c("Variables"), variable.name="Fold", value.name="Value")

# Plot the results using ggplot
ggplot() +
  # Add the AUC values for each fold as points
  geom_point(data=results_XGB_Folds_expanded, aes(x=Variables, y=Value, color=Fold), alpha = 0.5) +
  # Add the average AUC for each number of variables as a line
  geom_line(data=results_df_XGB, aes(x=Variables, y=AUC_train), cex = 1.5) +
  # Add axis labels and title
  labs(x="Number of variables", y="AUC", title="XGB AUC vs Number of Variables")






#AUC
# create a long format of results_RF_Folds
results_XGB_Folds_expanded_AUC <- reshape2::melt(results_XGB_Folds, id.vars = "Variables", 
                                                 variable.name = "Fold", value.name = "AUC")
# filter only the AUC columns
results_XGB_Folds_expanded_AUC <- results_XGB_Folds_expanded_AUC %>%
  filter(str_detect(Fold, "AUC"))
# convert Variables column to numeric
results_XGB_Folds_expanded_AUC$Variables <- as.numeric(results_XGB_Folds_expanded_AUC$Variables)
# calculate the mean AUC for each Variables value
results_XGB_Folds_mean_scaled <- aggregate(AUC ~ Variables, data = results_XGB_Folds_expanded_AUC, mean)
# plot the AUC curve
ggplot(results_XGB_Folds_expanded_AUC, aes(x = Variables, y = AUC)) +
  geom_point(aes(color = Fold), alpha = 0.5) +
  geom_line(data = results_XGB_Folds_mean_scaled, aes(x = Variables, y = AUC), size = 1.2) +
  labs(x = "Number of Clinical Variables", y = "AUC") +
  theme_classic() +
  theme(plot.title = element_blank(),
        axis.text.x = element_text(size = 24, color = "black", face = "bold"),
        axis.text.y = element_text(size = 24, color = "black", face = "bold"),
        axis.title.x = element_text(size = 14, color = "black", face = "bold"),
        axis.title.y = element_text(size = 14, color = "black", face = "bold"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(size = 1)) +
  guides(color = FALSE)+
  scale_x_continuous(expand = c(0, 0))+
  scale_y_continuous(limits = c(0.4, 0.8), expand = c(0, 0))
#leave out for now
geom_vline(xintercept = 54, linetype = "dotted", color = "red") +
  geom_segment(aes(x = 54, y = 0.5, xend = 54, yend = 0.75), linetype = "dotted", color = "red")+
  annotate("text", x = 69, y = 0.5, label = "X = 54", color = "red", vjust = -1, hjust = 1)




#Sensitivity
# create a long format of results_RF_Folds
results_XGB_Folds_expanded_sens <- reshape2::melt(results_XGB_Folds, id.vars = "Variables", 
                                                  variable.name = "Fold", value.name = "Sensitivity")
# filter only the Sensitivity columns
results_XGB_Folds_expanded_sens <- results_XGB_Folds_expanded_sens %>%
  filter(str_detect(Fold, "Sensitivity"))
# convert Variables column to numeric
results_XGB_Folds_expanded_sens$Variables <- as.numeric(results_XGB_Folds_expanded_sens$Variables)
# calculate the mean Sensitivity for each Variables value
results_XGB_Folds_mean_sens <- aggregate(Sensitivity ~ Variables, data = results_XGB_Folds_expanded_sens, mean)
# plot the Sensitivity curve
ggplot(results_XGB_Folds_expanded_sens, aes(x = Variables, y = Sensitivity)) +
  geom_point(aes(color = Fold), alpha = 0.5) +
  geom_line(data = results_XGB_Folds_mean_sens, aes(x = Variables, y = Sensitivity), size = 1.2) +
  labs(x = "Number of Variables", y = "Sensitivity") +
  theme_classic() +
  theme(plot.title = element_blank(),
        axis.text.x = element_text(size = 24, color = "black", face = "bold"),
        axis.text.y = element_text(size = 24, color = "black", face = "bold"),
        axis.title.x = element_text(size = 14, color = "black", face = "bold"),
        axis.title.y = element_text(size = 14, color = "black", face = "bold"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(size = 1)) +
  guides(color = FALSE)+
  scale_x_continuous(expand = c(0, 0))





# Specificity
# create a long format of results_RF_Folds
results_XGB_Folds_expanded_Specificity <- reshape2::melt(results_XGB_Folds, id.vars = "Variables", 
                                                         variable.name = "Fold", value.name = "Specificity")
# filter only the Specificity columns
results_XGB_Folds_expanded_Specificity <- results_XGB_Folds_expanded_Specificity %>%
  filter(str_detect(Fold, "Specificity"))
# convert Variables column to numeric
results_XGB_Folds_expanded_Specificity$Variables <- as.numeric(results_XGB_Folds_expanded_Specificity$Variables)
# calculate the mean Specificity for each Variables value
results_XGB_Folds_mean_scaled <- aggregate(Specificity ~ Variables, data = results_XGB_Folds_expanded_Specificity, mean)
# plot the Specificity curve
ggplot(results_XGB_Folds_expanded_Specificity, aes(x = Variables, y = Specificity)) +
  geom_point(aes(color = Fold), alpha = 0.5) +
  geom_line(data = results_XGB_Folds_mean_scaled, aes(x = Variables, y = Specificity), size = 1.2) +
  labs(x = "Number of Variables", y = "Specificity",
       title = "Specificity for Clinical Features (XGB)") +
  theme_classic() +
  theme(plot.title = element_blank(),
        axis.text.x = element_text(size = 24, color = "black", face = "bold"),
        axis.text.y = element_text(size = 24, color = "black", face = "bold"),
        axis.title.x = element_text(size = 14, color = "black", face = "bold"),
        axis.title.y = element_text(size = 14, color = "black", face = "bold"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(size = 1)) +
  guides(color = FALSE) +
  scale_x_continuous(expand = c(0, 0))


#DECIDE BEST SUBSET
results_df_XGB$Average <- (results_df_XGB$AUC_train + results_df_XGB$Sensitivity_train + results_df_XGB$Specificity_train)/3






#Create Master  TABLE
Master_table <- left_join(results_df, results_SVM_df, by = "Variables") %>%
  left_join(results_df_XGB, by = "Variables")
Master_table <- Master_table[,c(-5,-9)]

Master_table <- Master_table %>%
  select(Variables,
         contains("AUC"),
         contains("Specificity"),
         contains("Sensitivity"))

master_colnames <- c("Clinical Features", "RF AUC","SVM AUC","XGB AUC", "RF Spec", "SVM Spec", "XGB Spec", "RF Sens", "SVM Sens", "XGB Sens")
colnames(Master_table) <- master_colnames


Master_table[is.na(Master_table)] <- 0

Master_table$average <- rowSums(Master_table[,2:10])/9

write_xlsx(Master_table,"..\\Master Table Clinical.xlsx")





###########################                  Predicting in UNSEEN TEST SET          #######################################################3

#WE Select Best Number of Features from the Master Table
# Best are 4 (imbalanced towards Spec so we also do lowest balanced number which is 11)

# Try 4 and 11

#4

#RANDOM FOREST
best_4 <- combined_df$Variable[1:4]
Top_4_df <- Train_m1[, c("lymph_class", best_4)]

set.seed(42)
rf_model_4 <- randomForest(lymph_class ~ ., data = Top_4_df, ntree = 500, importance = TRUE)

# Make predictions on test set for AUC, sensitivity, and specificity
predict_4 <- predict(rf_model_4, Test_m1, type = "prob")
predict_4_df <- data.frame(prob = predict_4[, 2], class = Test_m1$lymph_class)
auc_4 <- roc(predict_4_df$class, predict_4_df$prob)$auc
sens_4 <- sum(predict_4[Test_m1$lymph_class == 1, 2] >= 0.5)/sum(Test_m1$lymph_class == 1)
spec_4 <- sum(predict_4[Test_m1$lymph_class == 0, 2] < 0.5)/sum(Test_m1$lymph_class == 0)

roc_4 <- roc(predict_4_df$class, predict_4_df$prob)
par(pty = "s")
plot(roc_4, xlim = c(1, 0), col = "blue", main = "ROC using Top 4 Clinical Predictors in Unseen Test", auc.polygon = TRUE, legacy.axes = TRUE, asp = NA)
text(0.3, 0.3, paste("AUC =", round(auc_4, 2)), col = "blue", cex = 1)
text(0.3, 0.24, paste("Sensitivity =", round(sens_4, 2)), col = "blue", cex = 1)
text(0.3, 0.18, paste("Specificity =", round(spec_4, 2)), col = "blue", cex = 1)


#SVM

# Train SVM model
set.seed(42)
svm_model_4 <- svm(lymph_class ~ ., data = Top_4_df, kernel = "sigmoid", 
                   type = "C-classification", probability = TRUE,
                   maxiter = 5000)
# Make predictions on test set for AUC, sensitivity, and specificity
predict_4 <- predict(svm_model_4, Test_m1, probability = TRUE)
predict_4_df <- data.frame(prob = attr(predict_4, "probabilities")[,2], class = Test_m1$lymph_class)
auc_svm_4 <- roc(predict_4_df$class, predict_4_df$prob)$auc
sens_svm_4 <- sum(predict_4[Test_m1$lymph_class == 1] == 1)/sum(Test_m1$lymph_class == 1)
spec_svm_4 <- sum(predict_4[Test_m1$lymph_class == 0] == 0)/sum(Test_m1$lymph_class == 0)

roc_svm_4 <- roc(predict_4_df$class, predict_4_df$prob)
par(pty = "s")
plot(roc_4, xlim = c(1, 0), col = "blue", main = "ROC using Top 4 Clinical Predictors in Unseen Test (SVM)", auc.polygon = TRUE, legacy.axes = TRUE, asp = NA)
text(0.3, 0.3, paste("AUC =", round(auc_svm_4, 2)), col = "blue", cex = 1)
text(0.3, 0.24, paste("Sensitivity =", round(sens_svm_4, 2)), col = "blue", cex = 1)
text(0.3, 0.18, paste("Specificity =", round(spec_svm_4, 2)), col = "blue", cex = 1)



#XGBOOST
best_4 <- combined_df$Variable[1:4]
Top_4_df_XGB <- Train_m1[, c("lymph_class", best_4)]
valid_set_4 <- Test_m1[c("lymph_class", best_4)]

Top_4_df_XGB[,-1] <- as.numeric(ifelse(Top_4_df_XGB[,-1] == 1, 1, 0))
valid_set_4[,-1] <- as.numeric(ifelse(valid_set_4[,-1] == 1, 1, 0))

Top_4_df_XGB$lymph_class <- ifelse(Top_4_df_XGB$lymph_class == 1, 1, 0)
valid_set_4$lymph_class <- ifelse(valid_set_4$lymph_class == 1, 1, 0)

# Train XGBoost model on training set
xgb_model <- xgboost(data = as.matrix(Top_4_df_XGB[,-1]), label = Top_4_df_XGB$lymph_class, nrounds = 100, objective = "binary:logistic")

# Make predictions on validation set for AUC, sensitivity, and specificity

train_prediction_probs <- predict(xgb_model, as.matrix(valid_set_4[,-1]))
train_prediction_df <- data.frame(prob = train_prediction_probs, class = valid_set_4$lymph_class)

roc_xgb_4 <- roc(train_prediction_df$class, train_prediction_df$prob)

auc_xgb_4 <- roc(train_prediction_df$class, train_prediction_df$prob)$auc
sens_xgb_4 <- sum(train_prediction_df$prob[Test_m1$lymph_class == 1] >= 0.5)/sum(Test_m1$lymph_class == 1)
spec_xgb_4 <- sum(train_prediction_df$prob[Test_m1$lymph_class == 0] < 0.5)/sum(Test_m1$lymph_class == 0)



#MUltiPlot
par(pty = "s")
plot(roc_4, col = "#1f77b4", main = "4 Feature Prediction ROC", auc.polygon = FALSE, legacy.axes = TRUE, asp = NA, legacy.show.adj = FALSE, xlim = c(1, 0), ylim = c(0, 1))
plot(roc_svm_4, col = "#ff7f0e", add = TRUE)
plot(roc_xgb_4, col = "#2ca02c", add = TRUE )

# Add AUC values
legend("bottomright", legend = c(paste0("RF AUC = ", round(auc_4, 2)),
                                 paste0("SVM AUC = ", round(auc_svm_4, 2)),
                                 paste0("XGB AUC = ", round(auc_xgb_4, 2))),
       col = c("#1f77b4", "#ff7f0e", "#2ca02c"), lty = 1)



#11
best_11 <- combined_df$Variable[1:11]
Top_11_df <- Train_m1[, c("lymph_class", best_11)]

set.seed(42)
rf_model_11 <- randomForest(lymph_class ~ ., data = Top_11_df, ntree = 500, importance = TRUE)

# Make predictions on test set for AUC, sensitivity, and specificity
predict_11 <- predict(rf_model_11, Test_m1, type = "prob")
predict_11_df <- data.frame(prob = predict_11[, 2], class = Test_m1$lymph_class)
auc_11 <- roc(predict_11_df$class, predict_11_df$prob)$auc
sens_11 <- sum(predict_11[Test_m1$lymph_class == 1, 2] >= 0.5)/sum(Test_m1$lymph_class == 1)
spec_11 <- sum(predict_11[Test_m1$lymph_class == 0, 2] < 0.5)/sum(Test_m1$lymph_class == 0)

roc_11 <- roc(predict_11_df$class, predict_11_df$prob)
par(pty = "s")
plot(roc_11, xlim = c(1, 0), col = "blue", main = "11 Feature Prediction ROC", auc.polygon = TRUE, legacy.axes = TRUE, asp = NA)
text(0.3, 0.3, paste("AUC =", round(auc_11, 2)), col = "blue", cex = 1)
text(0.3, 0.24, paste("Sensitivity =", round(sens_11, 2)), col = "blue", cex = 1)
text(0.3, 0.18, paste("Specificity =", round(spec_11, 2)), col = "blue", cex = 1)


#SVM
# Train SVM model
set.seed(42)
svm_model_11 <- svm(lymph_class ~ ., data = Top_11_df, kernel = "sigmoid", 
                    type = "C-classification", probability = TRUE,
                    maxiter = 5000)
# Make predictions on test set for AUC, sensitivity, and specificity
predict_11 <- predict(svm_model_11, Test_m1, probability = TRUE)
predict_11_df <- data.frame(prob = attr(predict_11, "probabilities")[,2], class = Test_m1$lymph_class)
auc_svm_11 <- roc(predict_11_df$class, predict_11_df$prob)$auc
sens_svm_11 <- sum(predict_11[Test_m1$lymph_class == 1] == 1)/sum(Test_m1$lymph_class == 1)
spec_svm_11 <- sum(predict_11[Test_m1$lymph_class == 0] == 0)/sum(Test_m1$lymph_class == 0)

roc_svm_11 <- roc(predict_11_df$class, predict_11_df$prob)
par(pty = "s")
plot(roc_svm_11, xlim = c(1, 0), col = "red", main = "11 Feature Prediction ROC", auc.polygon = TRUE, legacy.axes = TRUE, asp = NA)
text(0.3, 0.3, paste("AUC =", round(auc_svm_11, 2)), col = "blue", cex = 1)
text(0.3, 0.24, paste("Sensitivity =", round(sens_svm_11, 2)), col = "blue", cex = 1)
text(0.3, 0.18, paste("Specificity =", round(spec_svm_11, 2)), col = "blue", cex = 1)

#XGBOOST
best_11 <- combined_df$Variable[1:11]
Top_11_df_XGB <- Train_m1[, c("lymph_class", best_11)]
valid_set_11 <- Test_m1[c("lymph_class", best_11)]

Top_11_df_XGB[,-1] <- as.numeric(ifelse(Top_11_df_XGB[,-1] == 1, 1, 0))
valid_set_11[,-1] <- as.numeric(ifelse(valid_set_11[,-1] == 1, 1, 0))

Top_11_df_XGB$lymph_class <- ifelse(Top_11_df_XGB$lymph_class == 1, 1, 0)
valid_set_11$lymph_class <- ifelse(valid_set_11$lymph_class == 1, 1, 0)

# Train XGBoost model on training set
xgb_model <- xgboost(data = as.matrix(Top_11_df_XGB[,-1]), label = Top_11_df_XGB$lymph_class, nrounds = 100, objective = "binary:logistic")

# Make predictions on validation set for AUC, sensitivity, and specificity

train_prediction_probs <- predict(xgb_model, as.matrix(valid_set_11[,-1]))
train_prediction_df <- data.frame(prob = train_prediction_probs, class = valid_set_11$lymph_class)

roc_xgb_11 <- roc(train_prediction_df$class, train_prediction_df$prob)

auc_xgb_11 <- roc(train_prediction_df$class, train_prediction_df$prob)$auc
sens_xgb_11 <- sum(train_prediction_df$prob[Test_m1$lymph_class == 1] >= 0.5)/sum(Test_m1$lymph_class == 1)
spec_xgb_11 <- sum(train_prediction_df$prob[Test_m1$lymph_class == 0] < 0.5)/sum(Test_m1$lymph_class == 0)



#MUltiPlot
par(pty = "s")
plot(roc_11, col = "#1f77b4", main = "11 Feature Prediction ROC", auc.polygon = FALSE, legacy.axes = TRUE, asp = NA, legacy.show.adj = FALSE, xlim = c(1, 0), ylim = c(0, 1))
plot(roc_svm_11, col = "#ff7f0e", add = TRUE)
plot(roc_xgb_11, col = "#2ca02c", add = TRUE )

# Add AUC values
legend("bottomright", legend = c(paste0("RF AUC = ", round(auc_11, 2)),
                                 paste0("SVM AUC = ", round(auc_svm_11, 2)),
                                 paste0("XGB AUC = ", round(auc_xgb_11, 2))),
       col = c("#1f77b4", "#ff7f0e", "#2ca02c"), lty = 1)






###################               TRYING ALL POSSIBLE COMBINATIONS OF FEATURES UP TO 4 ############3



library(randomForest)
library(pROC)
install.packages("combinat")
library(combinat)

# Define variables
num_vars <- 3
target_auc <- 0.75

# Initialize empty data frames to store results
results_df <- data.frame(matrix(ncol=4, nrow= 0))
colnames(results_df) <- c("Variables", "AUC_train", "Sensitivity_train", "Specificity_train")

results_RF_Folds <- data.frame(matrix(ncol=16, nrow= 0))
colnames(results_RF_Folds) <- c("Variables", 
                                paste(rep(c("AUC", "Sensitivity", "Specificity"), each=5), 1:5, sep="_"))


# Step 1: Extract row names from imp_df_all
rownames_imp <- rownames(imp_df_all)
# Step 2: Subset columns of Train_m1
Train_m1 <- Train_m1[,c("lymph_class",rownames_imp)]
# Step 3: Update column names
colnames(Train_m1)[2:78] <- rownames_imp

set.seed(42)
# Loop over number of variables to include
for (i in 4:4) {
  # Get all combinations of i variables
  comb <- combinat::combn(names(Train_m1)[-1], 4)
  
  # Loop over all combinations
  for (j in 1:ncol(comb)) {
    selected_var <- comb[, j]
    
    # Create Train data frame with response variable and selected genes
    Training <- Train_m1[, c("lymph_class", selected_var)]
    
    # Set up K-fold cross-validation
    set.seed(42)
    K <- 5
    folds <- sample(rep(1:K, length.out=nrow(Training)))
    
    # Initialize empty vectors to store predictions and labels for each fold
    train_predictions <- numeric(nrow(Training))
    train_labels <- numeric(nrow(Training))
    auc_vec <- numeric(K)
    sens_vec <- numeric(K)
    spec_vec <- numeric(K)
    
    # Train and predict using K-fold cross-validation
    for (k in 1:K) {
      # Split data into training and validation sets
      train_set <- Training[folds != k, ]
      valid_set <- Training[folds == k, ]
      
      # Train random forest model
      rf_model <- randomForest(lymph_class ~ ., data=train_set, ntree=500, importance=TRUE)
      
      # Make predictions on training set for AUC, sensitivity, and specificity
      train_prediction_probs <- predict(rf_model, valid_set, type = "prob")
      train_prediction_df <- data.frame(prob = train_prediction_probs[,2], class = valid_set$lymph_class)
      
      # Calculate AUC, sensitivity, and specificity of ROC curve for validation and training sets
      auc_vec[k] <- roc(train_prediction_df$class, train_prediction_df$prob)$auc
      sens_vec[k] <- sum(train_prediction_probs[valid_set$lymph_class==1,2] >= 0.5)/sum(valid_set$lymph_class==1)
      spec_vec[k] <- sum(train_prediction_probs[valid_set$lymph_class==0,2] < 0.5)/sum(valid_set$lymph_class==0)
    }
    
    # Calculate average AUC, sensitivity, and specificity across all folds
    avg_auc_train <- mean(auc_vec)
    avg_sens_train <- mean(sens_vec)
    avg_spec_train <- mean(spec_vec)
    
    # Store results in results_df
    results_df[j, "Variables"] <- j
    results_df[j, "AUC_train"] <- avg_auc_train
    results_df[j, "Sensitivity_train"] <- avg_sens_train
    results_df[j, "Specificity_train"] <- avg_spec_train
    
    # Store AUC, sensitivity, and specificity for each fold in results_RF_Folds
    results_RF_Folds[j, 1] <- paste(selected_var, collapse = ", ")
    results_RF_Folds[j, c(paste0("AUC_", 1:5))] <- auc_vec
    results_RF_Folds[j, c(paste0("Sensitivity_", 1:5))] <- sens_vec
    results_RF_Folds[j, c(paste0("Specificity_", 1:5))] <- spec_vec
    
    if (avg_auc_train > 0.75) {
      break
    }
  }
  if (avg_auc_train > 0.75) {
    break
  }  
}




### FULL EXPRESSION SCRIPT ####
BRCA_PT_data = readRDS(file = "SE_PT_Samples.RDS")
Count_Matrix_TPM <- assay(BRCA_PT_data, "tpm_unstrand")
Count_Matrix_TPM <- as.data.frame(Count_Matrix_TPM)
write.csv(Count_Matrix_TPM,"Honours Project sample files\\Counts_TPM.csv", row.names = TRUE)

#Import the counts and clinical data (see Survival Analysis for changes made to this data before importing)
setwd("Lymph nodes predictor")
library(readxl)
setwd("Lymph nodes predictor/")
Clinical_lymph <- read_excel("Clinical_lymph_1.xlsx")
Clinical_lymph <- as.data.frame(Clinical_lymph)
View(Clinical_lymph) # n = 1044 (patients with <30 day follow up or time to death have been removed)


setwd("..")
library(readr)
Revised_Counts_Matrix <- read_csv("Counts_TPM.csv")
Revised_Counts_Matrix <- as.data.frame(Revised_Counts_Matrix)
rownames(Revised_Counts_Matrix) <- Revised_Counts_Matrix[,1]
Revised_Counts_Matrix <- Revised_Counts_Matrix[,-1]
View(Revised_Counts_Matrix)
library(stringr) 
colnames(Revised_Counts_Matrix) <- substr(colnames(Revised_Counts_Matrix),1,12)

#Remove <30 day follow up or <30 days to death from Counts
counts <- t(Revised_Counts_Matrix)
row_names_df_to_remove<-c( "TCGA-PL-A8LV", "TCGA-A8-A094", "TCGA-A8-A09Z", "TCGA-C8-A133", "TCGA-AC-A3W6", "TCGA-A8-A096", "TCGA-A8-A09G", "TCGA-D8-A1JK", "TCGA-A8-A083", "TCGA-C8-A12T", "TCGA-E2-A9RU", "TCGA-C8-A12K", "TCGA-C8-A26Y", "TCGA-A8-A081", "TCGA-A8-A08H", "TCGA-LL-A6FP", "TCGA-A8-A06N", "TCGA-A8-A090", "TCGA-A2-A25D", "TCGA-B6-A0IC", "TCGA-AC-A23H", "TCGA-E9-A245", "TCGA-C8-A275", "TCGA-AC-A7VC", "TCGA-BH-A1F8", "TCGA-C8-A1HJ", "TCGA-PL-A8LX", "TCGA-AN-A0AM", "TCGA-AN-A041", "TCGA-PL-A8LY", "TCGA-AN-A0XP", "TCGA-AN-A0XS", "TCGA-AN-A03X", "TCGA-AN-A0FV", "TCGA-AN-A0FZ", "TCGA-AN-A0AT", "TCGA-AN-A0FY", "TCGA-AN-A046", "TCGA-AN-A0FX", "TCGA-AN-A0XR", "TCGA-AN-A0AS", "TCGA-AN-A0AR", "TCGA-AN-A0XN", "TCGA-AN-A0XT", "TCGA-AN-A0XU", "TCGA-AN-A03Y", "TCGA-AN-A0FW", "TCGA-AN-A0G0","TCGA-AN-A049","TCGA-E9-A244","TCGA-E9-A5FL","TCGA-E9-A5FL")
counts <- counts[!(row.names(counts) %in% row_names_df_to_remove),]
counts <- as.data.frame(counts)
View(counts)

#Need to remove patients with no lymph node data (NX), use the code from other laptop
#Removing rows with no data or NX for lymph nodes from clinical and counts
Rows_to_remove <- c(which(Clinical_lymph == "NX", arr.ind=TRUE))
Rows_to_remove <- head(Rows_to_remove, -17)
Clinical_lymph_edit <- Clinical_lymph[!(row.names(Clinical_lymph) %in% Rows_to_remove),]

table(is.na(Clinical_lymph_edit$ajcc_pathologic_n))
which(is.na(Clinical_lymph_edit$ajcc_pathologic_n))
Clinical_lymph_edit <- Clinical_lymph_edit[-1021,]  #remove patient with NA value

Removing <- c(which(Clinical_lymph$ajcc_pathologic_n == "NX", arr.ind = TRUE))
counts <- counts[-Removing,]
counts <- counts[-1021,]
all(row.names(counts) == Clinical_lymph_edit$patient)  #True, so all the counts match the patients and are in the same order

################## NEED TO REMOVE MALES FROM CLINICAL AND EXPRESION DATA
removing2 <- c(which(Clinical_lymph_edit$gender == "male", arr.ind = TRUE))
counts <- counts[-removing2,]

row.names(Clinical_lymph_edit) <- NULL
remove_male <- c(which(Clinical_lymph_edit$gender == "male", arr.ind = TRUE))
Clinical_lymph_edit <- Clinical_lymph_edit[!(row.names(Clinical_lymph_edit) %in% remove_male),]
table(Clinical_lymph_edit$gender)
all(row.names(counts) == Clinical_lymph_edit$patient)  #True, so all the counts match the patients and are in the same order

#Removing genes from the count matrix that have a TPM <1 in 70% or more of the samples
threshold <- 1 # set the threshold TPM count
percent_cutoff <- 0.7 # set the percentage of samples with TPM below threshold
low_count_genes <- colSums(counts < threshold) / nrow(counts) >= percent_cutoff
low_count_genes <- as.data.frame(low_count_genes)
colum_names_to_remove <- rownames(low_count_genes)[which(low_count_genes == "TRUE" , arr.ind = TRUE)[, 1]]
counts_filtered <- counts[,!(colnames(counts) %in% colum_names_to_remove)] 
#This left us with 18976 genes based on if 70% or more of the samples has a TPM value of less than 1, the gene was removed
write.csv(counts_filtered,"..\\Expression data frame for combined analyis.csv", row.names = TRUE)


###################################################################################3

#Perform stratified split of the data into 70% 30% Training and test repectively.
library(caret)
set.seed(32) #random seed for reproducibility
trainIndex <- createDataPartition(Clinical_lymph_edit$lymph_class, p = .7,
                                  list = FALSE,
                                  times = 1)
Train <- Clinical_lymph_edit[ trainIndex,]
Test <- Clinical_lymph_edit[-trainIndex,]

Train_matching_counts <- counts_filtered[ trainIndex,]
Test_matching_counts <- counts_filtered[-trainIndex,]


##################################  Removing near 0 variance genes

nzv <- nearZeroVar(Train_matching_counts, saveMetrics= TRUE)
table(nzv$nzv) # All false, no need to remove

#Removing correlated predictors (using cutoff 0.9 or 0.95, depends how much is left)
descrCor <- cor(Train_matching_counts)
highlyCorDescr <- findCorrelation(descrCor, cutoff = .9) #450 highly correlated genes
Train_matching_counts <- Train_matching_counts[,-highlyCorDescr] #removing 1 of the highly correlated pairs
summary(Train_matching_counts[upper.tri(Train_matching_counts)])



#Create the data frame with the genes and the response variable
Train_data_frame <- cbind(Train$lymph_class,Train_matching_counts)
colnames(Train_data_frame)[colnames(Train_data_frame) == "Train$lymph_class"] <- "lymph_class"
#Need to convert N0 and N1 to 0 and 1
Train_data_frame$lymph_class <- ifelse(Train_data_frame$lymph_class == "N1", 1, 0) #1 = metastasis

Train_data_frame$lymph_class <- as.factor(Train_data_frame$lymph_class)




###########             Robust RFE USINF STRATIFIED SAMPLE SUBSETS   
library(caret)

#Creating a table for subsets of samples which ranks feature importance scores. Then we can get the total feature importance score for each feature across all the subsets and generate ROC curves.
set.seed(42)

# create empty list to store results
results_importance <- list()

# create vector of percentages
percentages <- seq(0.8, 0.1, by = -0.1)

# stratify based on lymph_class column
lymph_0 <- Train_data_frame_original[Train_data_frame_original$lymph_class == 0, ]
lymph_1 <- Train_data_frame_original[Train_data_frame_original$lymph_class == 1, ]

library(progress)

# Initialize progress bar
pb <- progress_bar$new(total = length(percentages))

# loop over each percentage
for (p in percentages) {
  # randomly sample p% of patients for each lymph class
  print(p)
  n_0 <- round(nrow(lymph_0)*p)
  n_1 <- round(nrow(lymph_1)*p)
  
  sample_rows_0 <- sample(nrow(lymph_0), n_0)
  sampled_data_0 <- lymph_0[sample_rows_0, ]
  
  sample_rows_1 <- sample(nrow(lymph_1), n_1)
  sampled_data_1 <- lymph_1[sample_rows_1, ]
  
  # combine samples
  sampled_data <- rbind(sampled_data_0, sampled_data_1)
  
  # set up control function for rfe
  ctrl <- rfeControl(
    method = "repeatedcv",
    repeats = 5,
    verbose = FALSE,
    returnResamp = "final",
    functions = rfFuncs,
  )
  
  # perform rfe
  rfe_result <- rfe(
    x = sampled_data[, -1],
    y = sampled_data[, 1],
    sizes = c(1),
    rfeControl = ctrl
  )
  
  # get feature importance scores
  importance <- varImp(rfe_result, scale = FALSE)  #Scale set to false, worth checking is TRUE gives new ranking (may be better due to mixed data)
  
  # create data frame with feature importance scores
  importance_df <- data.frame(importance[, 1], row.names = rownames(importance))
  colnames(importance_df) <- c("Feature Importance")
  
  # store importance data frame in list of data frames
  results_importance[[as.character(0.1*100)]] <- importance_df
  
  pb$tick()
  
}

# print results
print(results_importance)

#GENERATING THE TABLE WITH IMPORTANCE SCORES FOR EACH GENE ACROSS ALL SUBSETS
# create empty data frame with one column for each subset
# combine all gene names into one vector and remove duplicates
library(tidyverse)

# Set rownames as a regular column in each data frame
results_importance <- lapply(results_importance, function(x) {
  rownames_to_column(x, "Gene")  
})


install.packages("purrr")
library(purrr)
library(dplyr)

# Add a new column to each data frame in the list with the same name as the data frame
names(results_importance) <- paste0("df_", seq_along(results_importance))

# Combine all data frames into one
combined_df <- reduce(results_importance, full_join, by = "Gene")

new_colnames <- c("Gene", "80%", "70%", "60%", "50%", "40%", "30%", "20%", "10%")

colnames(combined_df) <- new_colnames

library(dplyr)

combined_df$total <- rowMeans(combined_df[, 2:9], na.rm = TRUE)

combined_df <- combined_df %>% arrange(desc(total))



### CREATE FEATURE IMPORTANCE PLOT

library(ggplot2)

# create a new data frame with the desired ordering
combined_ordered <- combined_df[order(combined_df$total, decreasing = TRUE), ]


ggplot(combined_ordered, aes(x = total, y = reorder(Gene, total))) +
  geom_bar(stat = "identity", fill = "steelblue") +
  xlab("Average Score") +
  ylab("") +
  scale_y_discrete(labels = NULL) +
  theme(axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 16)) +  # center title
  ggtitle("Variable Importance Scores")



########################################## Training model

library(pROC)
library(randomForest)
library(caret)
#Define number of genes
num_genes <- 200

results_df_RF <- data.frame(matrix(ncol = 4, nrow =num_genes))
colnames(results_df_RF) <- c("Genes", "AUC_train","Sensitivity_train","Specificity_train")


results_RF_Folds <- data.frame(matrix(ncol=16, nrow=num_genes))
colnames(results_RF_Folds) <- c("Genes", 
                                paste(rep(c("AUC", "Sensitivity", "Specificity"), each=5), 1:5, sep="_"))


#Loop over number of genes to include
for (i in 1:num_genes) {
  #selected top i genes
  selected_genes <- combined_df$Gene[1:i]
  
  #Create Train data frame with response variable and selected genes
  Training <- Train_data_frame_original[, c("lymph_class", selected_genes)]
  
  #Set up K-fold cross-Validation
  set.seed(42) #was 42 and got up to 0.75 average    seed 45 also good result, less dispersion
  K <- 5
  folds <- createFolds(Training$lymph_class, K, list = TRUE)
  
  #Initialize empty vectors to store predictions and labels for each fold
  train_predictions <- numeric(nrow(Training))
  train_labels <- numeric(nrow(Training))
  auc_vec <- numeric(K)
  sens_vec <- numeric(K)
  spec_vec <- numeric(K)
  
  #Train and predict using K-fold cross validation
  for (k in 1:K) {
    # Split the data into training and training test
    train_indices <- unlist(folds[-k])
    valid_indices <- unlist(folds[k])
    
    train_set <- Training[train_indices, ]
    valid_set <- Training[valid_indices, ]
    
    #Train random forest model on training set
    rf_model <- randomForest(lymph_class ~ ., data = train_set, ntree = 500, importance = TRUE)
    
    #Make predictions on training set for AUC, Spec and Sens
    train_predictions_probs <- predict(rf_model, valid_set, type = "prob")
    train_prediction_df <- data.frame(prob = train_predictions_probs[,2], class = valid_set$lymph_class)
    
    # Calculate AUC, sensitivity and specificity across all folds
    auc_vec[k] <- roc(train_prediction_df$class, train_prediction_df$prob)$auc
    sens_vec[k] <- sum(train_predictions_probs[valid_set$lymph_class==1,2] >= 0.5)/sum(valid_set$lymph_class==1)
    spec_vec[k] <- sum(train_predictions_probs[valid_set$lymph_class==0,2] < 0.5)/sum(valid_set$lymph_class==0)
  }
  
  # Calculate average AUC, sensitivity, and specificity across all folds
  avg_auc_train <- mean(auc_vec)
  avg_sens_train <- mean(sens_vec)
  avg_spec_train <- mean(spec_vec)
  
  # Store results in results_df
  results_df_RF[i, "Genes"] <- i
  results_df_RF[i, "AUC_train"] <- avg_auc_train
  results_df_RF[i, "Sensitivity_train"] <- avg_sens_train
  results_df_RF[i, "Specificity_train"] <- avg_spec_train
  
  # Store AUC, sensitivity, and specificity for each fold in results_RF_Folds
  results_RF_Folds[i, c(paste0("AUC_", 1:5))] <- auc_vec
  results_RF_Folds[i, c(paste0("Sensitivity_", 1:5))] <- sens_vec
  results_RF_Folds[i, c(paste0("Specificity_", 1:5))] <- spec_vec
}
results_RF_Folds$Genes <- seq(1, num_genes)


library(ggplot2)
library(tidyverse)
#Everything plot

# Reshape the results_RF_Folds dataframe to long format for easier plotting
results_RF_Folds_long <- reshape2::melt(results_RF_Folds, id.vars=c("Genes"), variable.name="Fold", value.name="Value")

# Plot the results using ggplot
ggplot() +
  # Add the AUC values for each fold as points
  geom_point(data=results_RF_Folds_long, aes(x= Genes, y=Value, color=Fold)) +
  # Add the average AUC for each number of variables as a line
  geom_line(data=results_df, aes(x=Genes, y=AUC_train)) +
  # Add axis labels and title
  labs(x="Number of genes", y="AUC", title="Random Forest AUC vs Number of Genes")
#Individual plots

#AUC
library(ggplot2)

# create a long format of results_RF_Folds
results_RF_Folds_long_AUC <- reshape2::melt(results_RF_Folds, id.vars = "Genes", 
                                            variable.name = "Fold", value.name = "AUC")
# filter only the AUC columns
results_RF_Folds_long_AUC <- results_RF_Folds_long_AUC %>%
  filter(str_detect(Fold, "AUC"))
# convert Variables column to numeric
results_RF_Folds_long_AUC$Genes <- as.numeric(results_RF_Folds_long_AUC$Genes)
# calculate the mean AUC for each Variables value
results_RF_Folds_mean <- aggregate(AUC ~ Genes, data = results_RF_Folds_long_AUC, mean)
# plot the AUC curve
ggplot(results_RF_Folds_long_AUC, aes(x = Genes, y = AUC)) +
  geom_point(aes(color = Fold), alpha = 0.5) +
  geom_line(data = results_RF_Folds_mean, aes(x = Genes, y = AUC), size = 1.2) +
  labs(x = "Number of Genes", y = "AUC") +
  theme_classic() +
  theme(plot.title = element_blank(),
        axis.text.x = element_text(size = 24, color = "black", face = "bold"),
        axis.text.y = element_text(size = 24, color = "black", face = "bold"),
        axis.title.x = element_text(size = 14, color = "black", face = "bold"),
        axis.title.y = element_text(size = 14, color = "black", face = "bold"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(size = 1),
        plot.margin = unit(c(1, 1, 0, 1), "cm")) +  # increase left and right margins
  guides(color = FALSE) +
  scale_x_continuous(expand = c(0, 0))+
  scale_y_continuous(limits = c(0.4, 0.9))



geom_vline(xintercept = 26, linetype = "dotted", color = "red") +
  geom_segment(aes(x = 26, y = 0.5, xend = 26, yend = 0.85), linetype = "dotted", color = "red")+
  annotate("text", x = 42, y = 0.5, label = "X = 26", color = "red", vjust = -1, hjust = 1) 


#Sensitivity
# create a long format of results_RF_Folds
results_RF_Folds_long_Sens <- reshape2::melt(results_RF_Folds, id.vars = "Genes", 
                                             variable.name = "Fold", value.name = "Sensitivity")
# filter only the Sensitivity columns
results_RF_Folds_long_Sens <- results_RF_Folds_long_Sens %>%
  filter(str_detect(Fold, "Sensitivity"))
# convert Variables column to numeric
results_RF_Folds_long_Sens$Genes <- as.numeric(results_RF_Folds_long_Sens$Genes)
# calculate the mean Sensitivity for each Variables value
results_RF_Folds_mean_Sens <- aggregate(Sensitivity ~ Genes, data = results_RF_Folds_long_Sens, mean)
# plot the Sensitivity curve
ggplot(results_RF_Folds_long_Sens, aes(x = Genes, y = Sensitivity)) +
  geom_point(aes(color = Fold), alpha = 0.5) +
  geom_line(data = results_RF_Folds_mean_Sens, aes(x = Genes, y = Sensitivity), size = 1.2) +
  labs(x = "Number of Genes", y = "Sensitivity") +
  theme_classic() +
  theme(plot.title = element_blank(),
        axis.text.x = element_text(size = 24, color = "black", face = "bold"),
        axis.text.y = element_text(size = 24, color = "black", face = "bold"),
        axis.title.x = element_text(size = 14, color = "black", face = "bold"),
        axis.title.y = element_text(size = 14, color = "black", face = "bold"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(size = 1),
        plot.margin = unit(c(0, 1, 0, 1), "cm"))+ # increase left and right margins
  guides(color = FALSE)+
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(limits = c(0.4, 1))
geom_vline(xintercept = 26, linetype = "dotted", color = "red") +
  annotate("text", x = 42, y = 0.4, label = "X = 26", color = "red", vjust = -1, hjust = 1) +
  geom_segment(aes(x = 26, y = 0.4, xend = 26, yend = min(results_RF_Folds_mean_Sens$Sensitivity[results_RF_Folds_mean_Sens$Genes == 26], 0.74)), linetype = "dotted", color = "red")


#Specificity
# create a long format of results_RF_Folds
results_RF_Folds_long_Spec <- reshape2::melt(results_RF_Folds, id.vars = "Genes", 
                                             variable.name = "Fold", value.name = "Specificity")
# filter only the Specificity columns
results_RF_Folds_long_Spec <- results_RF_Folds_long_Spec %>%
  filter(str_detect(Fold, "Specificity"))
# convert Variables column to numeric
results_RF_Folds_long_Spec$genes <- as.numeric(results_RF_Folds_long_Spec$Genes)
# calculate the mean Specificity for each Variables value
results_RF_Folds_mean_Spec <- aggregate(Specificity ~ Genes, data = results_RF_Folds_long_Spec, mean)
# plot the Specificity curve
ggplot(results_RF_Folds_long_Spec, aes(x = Genes, y = Specificity)) +
  geom_point(aes(color = Fold), alpha = 0.5) +
  geom_line(data = results_RF_Folds_mean_Spec, aes(x = Genes, y = Specificity), size = 1.2) +
  labs(x = "Number of Genes", y = "Specificity") +
  theme_classic() +
  theme(plot.title = element_blank(),
        axis.text.x = element_text(size = 24, color = "black", face = "bold"),
        axis.text.y = element_text(size = 24, color = "black", face = "bold"),
        axis.title.x = element_text(size = 14, color = "black", face = "bold"),
        axis.title.y = element_text(size = 14, color = "black", face = "bold"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(size = 1),
        plot.margin = unit(c(0, 1, 0, 1), "cm")) + # increase left and right margins
  guides(color = FALSE)+
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(limits = c(0.2, 0.8))
geom_vline(xintercept = 26, linetype = "dotted", color = "red") +
  annotate("text", x = 42, y = 0.4, label = "X = 26", color = "red", vjust = -1, hjust = 1) +
  geom_segment(aes(x = 26, y = 0.4, xend = 26, yend = min(results_RF_Folds_mean_Spec$Specificity[results_RF_Folds_mean_Spec$Genes == 26], 0.74)), linetype = "dotted", color = "red")






########################                 SVM TRAINING MODEL                 #####################################33


library(pROC)
library(caret)

# Define number of variables to use
num_var <- 200


# Initialize empty data frames to store results
results_SVM_df <- data.frame(matrix(ncol=4, nrow=num_var))
colnames(results_SVM_df) <- c("Variables", "AUC_train", "Sensitivity_train", "Specificity_train")

results_SVM_Folds_scaled <- data.frame(matrix(ncol=16, nrow=num_var))
colnames(results_SVM_Folds_scaled) <- c("Variables", 
                                        paste(rep(c("AUC", "Sensitivity", "Specificity"), each=5), 1:5, sep="_"))

# Loop over number of variables to include
set.seed(42)
for (i in 1:num_var) {
  # Select top i genes
  selected_var_scaled <- combined_df$Gene[1:i]
  
  # Create Train data frame with response variable and selected genes
  Training <- Train_data_frame_original[, c("lymph_class", selected_var_scaled)]
  
  # Set up K-fold cross-validation
  set.seed(42)
  K <- 5
  folds <- createFolds(Training$lymph_class, K, list = TRUE)
  
  # Initialize empty vectors to store predictions and labels for each fold
  train_predictions <- numeric(nrow(Training))
  train_labels <- numeric(nrow(Training))
  auc_vec <- numeric(K)
  sens_vec <- numeric(K)
  spec_vec <- numeric(K)
  
  # Train and predict using K-fold cross-validation
  for (k in 1:K) {
    # Split data into training and validation sets
    train_indices <- unlist(folds[-k])
    valid_indices <- unlist(folds[k])
    
    train_set <- Training[train_indices, ]
    valid_set <- Training[valid_indices, ]
    
    
    # Train SVM model on training set
    svm_model <- svm(lymph_class ~ ., data = train_set, kernel = "sigmoid", 
                     type = "C-classification", probability = TRUE,
                     maxiter = 10000)
    
    # Make predictions on training set for AUC, sensitivity, and specificity
    train_prediction_probs <- predict(svm_model, valid_set, type = "prob")
    train_prediction_probs <- as.data.frame(train_prediction_probs)
    train_prediction_df <- data.frame(prob = train_prediction_probs[,1], class = valid_set$lymph_class)
    
    # Calculate AUC, sensitivity, and specificity of ROC curve for validation and training sets
    auc_vec[k] <- roc(as.numeric(train_prediction_df$class),as.numeric(train_prediction_df$prob))$auc
    sens_vec[k] <- sum(train_prediction_df$prob == "1" & train_prediction_df$class == "1")/sum(train_prediction_df$class==1)
    spec_vec[k] <- sum(train_prediction_df$prob == "0" & train_prediction_df$class == "0")/sum(train_prediction_df$class==0)
  }
  
  # Calculate average AUC, sensitivity, and specificity across all folds
  avg_auc_train <- mean(auc_vec)
  avg_sens_train <- mean(sens_vec)
  avg_spec_train <- mean(spec_vec)
  
  # Store results in results_df
  results_SVM_df[i, "Variables"] <- i
  results_SVM_df[i, "AUC_train"] <- avg_auc_train
  results_SVM_df[i, "Sensitivity_train"] <- avg_sens_train
  results_SVM_df[i, "Specificity_train"] <- avg_spec_train
  
  # Store AUC, sensitivity, and specificity for each fold in results_RF_Folds
  results_SVM_Folds_scaled[i, c(paste0("AUC_", 1:5))] <- auc_vec
  results_SVM_Folds_scaled[i, c(paste0("Sensitivity_", 1:5))] <- sens_vec
  results_SVM_Folds_scaled[i, c(paste0("Specificity_", 1:5))] <- spec_vec
}


#plot
results_SVM_Folds_scaled$Variables <- seq(1, num_var)


# Reshape the results_RF_Folds dataframe to long format for easier plotting
results_SVM_Folds_expanded <- reshape2::melt(results_SVM_Folds_scaled, id.vars=c("Variables"), variable.name="Fold", value.name="Value")

# Plot the results using ggplot
ggplot() +
  # Add the AUC values for each fold as points
  geom_point(data=results_SVM_Folds_expanded, aes(x=Variables, y=Value, color=Fold), alpha = 0.5) +
  # Add the average AUC for each number of variables as a line
  geom_line(data=results_SVM_df, aes(x=Variables, y=AUC_train), cex = 1.5) +
  # Add axis labels and title
  labs(x="Number of variables", y="AUC", title="SVM AUC vs Number of Variables")


#AUC
# create a long format of results_RF_Folds
results_SVM_Folds_expanded_AUC <- reshape2::melt(results_SVM_Folds_scaled, id.vars = "Variables", 
                                                 variable.name = "Fold", value.name = "AUC")
# filter only the AUC columns
results_SVM_Folds_expanded_AUC <- results_SVM_Folds_expanded_AUC %>%
  filter(str_detect(Fold, "AUC"))
# convert Variables column to numeric
results_SVM_Folds_expanded_AUC$Variables <- as.numeric(results_SVM_Folds_expanded_AUC$Variables)
# calculate the mean AUC for each Variables value
results_SVM_Folds_mean_scaled <- aggregate(AUC ~ Variables, data = results_SVM_Folds_expanded_AUC, mean)
# plot the AUC curve
ggplot(results_SVM_Folds_expanded_AUC, aes(x = Variables, y = AUC)) +
  geom_point(aes(color = Fold), alpha = 0.5) +
  geom_line(data = results_SVM_Folds_mean_scaled, aes(x = Variables, y = AUC), size = 1.2) +
  labs(x = "Number of Genes", y = "AUC") +
  theme_classic() +
  theme(plot.title = element_blank(),
        axis.text.x = element_text(size = 24, color = "black", face = "bold"),
        axis.text.y = element_text(size = 24, color = "black", face = "bold"),
        axis.title.x = element_text(size = 14, color = "black", face = "bold"),
        axis.title.y = element_text(size = 14, color = "black", face = "bold"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(size = 1),
        plot.margin = unit(c(1, 1, 0, 1), "cm")) +  # increase left and right margins
  guides(color = FALSE) +
  scale_x_continuous(expand = c(0, 0))+
  scale_y_continuous(limits = c(0.4, 0.9))


#leave out for now
geom_vline(xintercept = 54, linetype = "dotted", color = "red") +
  geom_segment(aes(x = 54, y = 0.5, xend = 54, yend = 0.75), linetype = "dotted", color = "red")+
  annotate("text", x = 69, y = 0.5, label = "X = 54", color = "red", vjust = -1, hjust = 1)


# Specificity
# create a long format of results_SVM_Folds
results_SVM_Folds_expanded_spec <- reshape2::melt(results_SVM_Folds_scaled, id.vars = "Variables", 
                                                  variable.name = "Fold", value.name = "Specificity")
# filter only the Specificity columns
results_SVM_Folds_expanded_spec <- results_SVM_Folds_expanded_spec %>%
  filter(str_detect(Fold, "Specificity"))
# convert Variables column to numeric
results_SVM_Folds_expanded_spec$Variables <- as.numeric(results_SVM_Folds_expanded_spec$Variables)
# calculate the mean Specificity for each Variables value
results_SVM_Folds_mean_spec <- aggregate(Specificity ~ Variables, data = results_SVM_Folds_expanded_spec, mean)
# plot the Specificity curve
ggplot(results_SVM_Folds_expanded_spec, aes(x = Variables, y = Specificity)) +
  geom_point(aes(color = Fold), alpha = 0.5) +
  geom_line(data = results_SVM_Folds_mean_spec, aes(x = Variables, y = Specificity), size = 1.2) +
  labs(x = "Number of Genes", y = "Specificity") +
  theme_classic() +
  theme(plot.title = element_blank(),
        axis.text.x = element_text(size = 24, color = "black", face = "bold"),
        axis.text.y = element_text(size = 24, color = "black", face = "bold"),
        axis.title.x = element_text(size = 14, color = "black", face = "bold"),
        axis.title.y = element_text(size = 14, color = "black", face = "bold"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(size = 1),
        plot.margin = unit(c(0, 1, 0, 1), "cm")) + # increase left and right margins
  guides(color = FALSE) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(limits = c(0.2, 0.8))
#leave out for now
geom_vline(xintercept = 54, linetype = "dotted", color = "red") +
  geom_segment(aes(x = 54, y = 0.5, xend = 54, yend = 0.7), linetype = "dotted", color = "red") +
  annotate("text", x = 69, y = 0.5, label = "X = 54", color = "red", vjust = -1, hjust = 1)



# Sensitivity
# create a long format of results_SVM_Folds
results_SVM_Folds_expanded_sens <- reshape2::melt(results_SVM_Folds_scaled, id.vars = "Variables", 
                                                  variable.name = "Fold", value.name = "Sensitivity")
# filter only the Sensitivity columns
results_SVM_Folds_expanded_sens <- results_SVM_Folds_expanded_sens %>%
  filter(str_detect(Fold, "Sensitivity"))
# convert Variables column to numeric
results_SVM_Folds_expanded_sens$Variables <- as.numeric(results_SVM_Folds_expanded_sens$Variables)
# calculate the mean Sensitivity for each Variables value
results_SVM_Folds_mean_sens <- aggregate(Sensitivity ~ Variables, data = results_SVM_Folds_expanded_sens, mean)
# plot the Sensitivity curve
ggplot(results_SVM_Folds_expanded_sens, aes(x = Variables, y = Sensitivity)) +
  geom_point(aes(color = Fold), alpha = 0.5) +
  geom_line(data = results_SVM_Folds_mean_sens, aes(x = Variables, y = Sensitivity), size = 1.2) +
  labs(x = "Number of Genes", y = "Sensitivity") +
  theme_classic() +
  theme(plot.title = element_blank(),
        axis.text.x = element_text(size = 24, color = "black", face = "bold"),
        axis.text.y = element_text(size = 24, color = "black", face = "bold"),
        axis.title.x = element_text(size = 14, color = "black", face = "bold"),
        axis.title.y = element_text(size = 14, color = "black", face = "bold"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(size = 1),
        plot.margin = unit(c(0, 1, 0, 1), "cm")) + # increase left and right margins
  guides(color = FALSE) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(limits = c(0.4, 1))


#Leave out for now
geom_vline(xintercept = 54, linetype = "dotted", color = "red") +
  geom_segment(aes(x = 54, y = 0.5, xend = 54, yend = 0.7), linetype = "dotted", color = "red") +
  annotate("text", x = 69, y = 0.5, label = "X = 54", color = "red", vjust = -1, hjust = 1)





##########################                XGBOOST   TRAINING MODEL             ##############################33333

# Load required packages
library(xgboost)
library(pROC)

# Define number of variables to use
num_var <- 200

# Initialize empty data frames to store results
results_df_XGB <- data.frame(matrix(ncol=4, nrow=num_var))
colnames(results_df_XGB) <- c("Variables", "AUC_train", "Sensitivity_train", "Specificity_train")

results_XGB_Folds <- data.frame(matrix(ncol=16, nrow=num_var))
colnames(results_XGB_Folds) <- c("Variables", 
                                 paste(rep(c("AUC", "Sensitivity", "Specificity"), each=5), 1:5, sep="_"))

# Loop over number of variables to include
set.seed(42)
for (i in 1:num_var) {
  # Select top i genes
  selected_var_scaled <- combined_df$Gene[1:i]
  
  # Create Train data frame with response variable and selected genes
  Training <- Train_data_frame_original[, c("lymph_class", selected_var_scaled)]
  
  # Set up K-fold cross-validation
  set.seed(42)
  K <- 5
  folds <- createFolds(Training$lymph_class, K, list = TRUE)
  
  # Initialize empty vectors to store predictions and labels for each fold
  train_predictions <- numeric(nrow(Training))
  train_labels <- numeric(nrow(Training))
  auc_vec <- numeric(K)
  sens_vec <- numeric(K)
  spec_vec <- numeric(K)
  
  # Train and predict using K-fold cross-validation
  for (k in 1:K) {
    # Split data into training and validation sets
    train_indices <- unlist(folds[-k])
    valid_indices <- unlist(folds[k])
    
    train_set <- Training[train_indices, ]
    valid_set <- Training[valid_indices, ]
    
    train_set$lymph_class <- ifelse(train_set$lymph_class == 1, 1, 0)
    valid_set$lymph_class <- ifelse(valid_set$lymph_class == 1, 1, 0)
    
    # Train XGBoost model on training set
    xgb_model <- xgboost(data = as.matrix(train_set[,-1]), label = train_set$lymph_class, nrounds = 100, objective = "binary:logistic")
    
    # Make predictions on validation set for AUC, sensitivity, and specificity
    
    train_prediction_probs <- predict(xgb_model, as.matrix(valid_set[,-1]))
    train_prediction_df <- data.frame(prob = train_prediction_probs, class = valid_set$lymph_class)
    
    # Calculate AUC, sensitivity, and specificity of ROC curve for validation and training sets
    auc_vec[k] <- roc(train_prediction_df$class,train_prediction_df$prob)$auc
    sens_vec[k] <- sum(train_prediction_probs[valid_set$lymph_class==1] >= 0.5)/sum(valid_set$lymph_class==1)
    spec_vec[k] <- sum(train_prediction_probs[valid_set$lymph_class==0] < 0.5)/sum(valid_set$lymph_class==0)
  }
  
  # Calculate average AUC, sensitivity, and specificity across all folds
  avg_auc_train <- mean(auc_vec)
  avg_sens_train <- mean(sens_vec)
  avg_spec_train <- mean(spec_vec)
  
  # Store results in results_df
  results_df_XGB[i, "Variables"] <- i
  results_df_XGB[i, "AUC_train"] <- avg_auc_train
  results_df_XGB[i, "Sensitivity_train"] <- avg_sens_train
  results_df_XGB[i, "Specificity_train"] <- avg_spec_train
  
  # Store AUC, sensitivity, and specificity for each fold in results_RF_Folds
  results_XGB_Folds[i, c(paste0("AUC_", 1:5))] <- auc_vec
  results_XGB_Folds[i, c(paste0("Sensitivity_", 1:5))] <- sens_vec
  results_XGB_Folds[i, c(paste0("Specificity_", 1:5))] <- spec_vec
}



#plot
results_XGB_Folds$Variables <- seq(1, num_var)


# Reshape the results_RF_Folds dataframe to long format for easier plotting
results_XGB_Folds_expanded <- reshape2::melt(results_XGB_Folds, id.vars=c("Variables"), variable.name="Fold", value.name="Value")

# Plot the results using ggplot
ggplot() +
  # Add the AUC values for each fold as points
  geom_point(data=results_XGB_Folds_expanded, aes(x=Variables, y=Value, color=Fold), alpha = 0.5) +
  # Add the average AUC for each number of variables as a line
  geom_line(data=results_df_XGB, aes(x=Variables, y=AUC_train), cex = 1.5) +
  # Add axis labels and title
  labs(x="Number of variables", y="AUC", title="XGB AUC vs Number of Variables")






#AUC
# create a long format of results_RF_Folds
results_XGB_Folds_expanded_AUC <- reshape2::melt(results_XGB_Folds, id.vars = "Variables", 
                                                 variable.name = "Fold", value.name = "AUC")
# filter only the AUC columns
results_XGB_Folds_expanded_AUC <- results_XGB_Folds_expanded_AUC %>%
  filter(str_detect(Fold, "AUC"))
# convert Variables column to numeric
results_XGB_Folds_expanded_AUC$Variables <- as.numeric(results_XGB_Folds_expanded_AUC$Variables)
# calculate the mean AUC for each Variables value
results_XGB_Folds_mean_scaled <- aggregate(AUC ~ Variables, data = results_XGB_Folds_expanded_AUC, mean)
# plot the AUC curve
ggplot(results_XGB_Folds_expanded_AUC, aes(x = Variables, y = AUC)) +
  geom_point(aes(color = Fold), alpha = 0.5) +
  geom_line(data = results_XGB_Folds_mean_scaled, aes(x = Variables, y = AUC), size = 1.2) +
  labs(x = "Number of Genes", y = "AUC") +
  theme_classic() +
  theme(plot.title = element_blank(),
        axis.text.x = element_text(size = 24, color = "black", face = "bold"),
        axis.text.y = element_text(size = 24, color = "black", face = "bold"),
        axis.title.x = element_text(size = 14, color = "black", face = "bold"),
        axis.title.y = element_text(size = 14, color = "black", face = "bold"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(size = 1),
        plot.margin = unit(c(1, 1, 0, 1), "cm")) + # increase left and right margins
  guides(color = FALSE) +
  scale_x_continuous(expand = c(0, 0))+
  scale_y_continuous(limits = c(0.4, 0.9))
#leave out for now
geom_vline(xintercept = 54, linetype = "dotted", color = "red") +
  geom_segment(aes(x = 54, y = 0.5, xend = 54, yend = 0.75), linetype = "dotted", color = "red")+
  annotate("text", x = 69, y = 0.5, label = "X = 54", color = "red", vjust = -1, hjust = 1)




#Sensitivity
# create a long format of results_RF_Folds
results_XGB_Folds_expanded_sens <- reshape2::melt(results_XGB_Folds, id.vars = "Variables", 
                                                  variable.name = "Fold", value.name = "Sensitivity")
# filter only the Sensitivity columns
results_XGB_Folds_expanded_sens <- results_XGB_Folds_expanded_sens %>%
  filter(str_detect(Fold, "Sensitivity"))
# convert Variables column to numeric
results_XGB_Folds_expanded_sens$Variables <- as.numeric(results_XGB_Folds_expanded_sens$Variables)
# calculate the mean Sensitivity for each Variables value
results_XGB_Folds_mean_sens <- aggregate(Sensitivity ~ Variables, data = results_XGB_Folds_expanded_sens, mean)
# plot the Sensitivity curve
ggplot(results_XGB_Folds_expanded_sens, aes(x = Variables, y = Sensitivity)) +
  geom_point(aes(color = Fold), alpha = 0.5) +
  geom_line(data = results_XGB_Folds_mean_sens, aes(x = Variables, y = Sensitivity), size = 1.2) +
  labs(x = "Number of Genes", y = "Sensitivity") +
  theme_classic() +
  theme(plot.title = element_blank(),
        axis.text.x = element_text(size = 24, color = "black", face = "bold"),
        axis.text.y = element_text(size = 24, color = "black", face = "bold"),
        axis.title.x = element_text(size = 14, color = "black", face = "bold"),
        axis.title.y = element_text(size = 14, color = "black", face = "bold"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(size = 1),
        plot.margin = unit(c(0, 1, 0, 1), "cm")) + # increase left and right margins
  guides(color = FALSE) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(limits = c(0.4, 1))





# Specificity
# create a long format of results_RF_Folds
results_XGB_Folds_expanded_Specificity <- reshape2::melt(results_XGB_Folds, id.vars = "Variables", 
                                                         variable.name = "Fold", value.name = "Specificity")
# filter only the Specificity columns
results_XGB_Folds_expanded_Specificity <- results_XGB_Folds_expanded_Specificity %>%
  filter(str_detect(Fold, "Specificity"))
# convert Variables column to numeric
results_XGB_Folds_expanded_Specificity$Variables <- as.numeric(results_XGB_Folds_expanded_Specificity$Variables)
# calculate the mean Specificity for each Variables value
results_XGB_Folds_mean_scaled <- aggregate(Specificity ~ Variables, data = results_XGB_Folds_expanded_Specificity, mean)
# plot the Specificity curve
ggplot(results_XGB_Folds_expanded_Specificity, aes(x = Variables, y = Specificity)) +
  geom_point(aes(color = Fold), alpha = 0.5) +
  geom_line(data = results_XGB_Folds_mean_scaled, aes(x = Variables, y = Specificity), size = 1.2) +
  labs(x = "Number of Genes", y = "Specificity") +
  theme_classic() +
  theme(plot.title = element_blank(),
        axis.text.x = element_text(size = 24, color = "black", face = "bold"),
        axis.text.y = element_text(size = 24, color = "black", face = "bold"),
        axis.title.x = element_text(size = 14, color = "black", face = "bold"),
        axis.title.y = element_text(size = 14, color = "black", face = "bold"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(size = 1),
        plot.margin = unit(c(0, 1, 0, 1), "cm")) + # increase left and right margins
  guides(color = FALSE) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(limits = c(0.2, 0.8))



#DECIDE BEST SUBSET
results_df_XGB$Average <- (results_df_XGB$AUC_train + results_df_XGB$Sensitivity_train + results_df_XGB$Specificity_train)/3





##############                      GENERATE MASTER TABLE             #############################################3

#Create Master  TABLE
SVM_Master <-results_SVM_df
colnames(SVM_Master)[1] <- "Genes"

XGB_Master <- results_df_XGB
colnames(XGB_Master)[1] <- "Genes"

Master_table <- left_join(results_df_RF, SVM_Master, by = "Genes") %>%
  left_join(XGB_Master, by = "Genes")

Master_table <- Master_table %>%
  select(Genes,
         contains("AUC"),
         contains("Specificity"),
         contains("Sensitivity"))

master_colnames <- c("Transcripts", "RF AUC","SVM AUC","XGB AUC", "RF Spec", "SVM Spec", "XGB Spec", "RF Sens", "SVM Sens", "XGB Sens")
colnames(Master_table) <- master_colnames


Master_table$average <- rowSums(Master_table[,2:10])/9

write_xlsx(Master_table,"..\\Master Table Expression.xlsx")





############################                  Predicting in UNSEEN TEST SET          #######################################################3

#WE Select Best Number of Features from the Master Table

#Best is 149 but there is a slight reduction in accuracy for 54 features
#RANDOM FOREST
best_149 <- combined_df$Gene[1:149]
Top_149_df <- Train_data_frame_original[, c("lymph_class", best_149)]

set.seed(42)
rf_model_149 <- randomForest(lymph_class ~ ., data = Top_149_df, ntree = 500, importance = TRUE)

# Make predictions on test set for AUC, sensitivity, and specificity
predict_149 <- predict(rf_model_149, Test_data_frame, type = "prob")
predict_149_df <- data.frame(prob = predict_149[, 2], class = Test_data_frame$lymph_class)
auc_149 <- roc(predict_149_df$class, predict_149_df$prob)$auc
sens_149 <- sum(predict_149[Test_data_frame$lymph_class == 1, 2] >= 0.5)/sum(Test_data_frame$lymph_class == 1)
spec_149 <- sum(predict_149[Test_data_frame$lymph_class == 0, 2] < 0.5)/sum(Test_data_frame$lymph_class == 0)

roc_149 <- roc(predict_149_df$class, predict_149_df$prob)
par(pty = "s")
plot(roc_149, xlim = c(1, 0), col = "blue", main = "149 Transcript Prediction ROC", auc.polygon = TRUE, legacy.axes = TRUE, asp = NA)
text(0.3, 0.3, paste("AUC =", round(auc_149, 2)), col = "blue", cex = 1)
text(0.3, 0.24, paste("Sensitivity =", round(sens_149, 2)), col = "blue", cex = 1)
text(0.3, 0.18, paste("Specificity =", round(spec_149, 2)), col = "blue", cex = 1)




#SVM

# Train SVM model
set.seed(42)
svm_model_149 <- svm(lymph_class ~ ., data = Top_149_df, kernel = "sigmoid", 
                     type = "C-classification", probability = TRUE,
                     maxiter = 1000)
# Make predictions on test set for AUC, sensitivity, and specificity
predict_149 <- predict(svm_model_149, Test_data_frame, probability = TRUE)
predict_149_df <- data.frame(prob = attr(predict_149, "probabilities")[,2], class = Test_data_frame$lymph_class)
auc_svm_149 <- roc(predict_149_df$class, predict_149_df$prob)$auc
sens_svm_149 <- sum(predict_149[Test_data_frame$lymph_class == 1] == 1)/sum(Test_data_frame$lymph_class == 1)
spec_svm_149 <- sum(predict_149[Test_data_frame$lymph_class == 0] == 0)/sum(Test_data_frame$lymph_class == 0)

roc_svm_149 <- roc(predict_149_df$class, predict_149_df$prob)
par(pty = "s")
plot(roc_149, xlim = c(1, 0), col = "blue", main = "149 Transcript Prediction ROC (SVM)", auc.polygon = TRUE, legacy.axes = TRUE, asp = NA)
text(0.3, 0.3, paste("AUC =", round(auc_svm_149, 2)), col = "blue", cex = 1)
text(0.3, 0.24, paste("Sensitivity =", round(sens_svm_149, 2)), col = "blue", cex = 1)
text(0.3, 0.18, paste("Specificity =", round(spec_svm_149, 2)), col = "blue", cex = 1)


#XGBOOST

Top_149_XGB <- Train_data_frame_original[,c("lymph_class", best_149) ]
valid_set_149 <- Test_data_frame[,c("lymph_class", best_149) ]

Top_149_XGB$lymph_class <- ifelse(Top_149_XGB$lymph_class == 1, 1, 0)
valid_set_149$lymph_class <- ifelse(valid_set_149$lymph_class == 1, 1, 0)

# Train XGBoost model on training set
xgb_model <- xgboost(data = as.matrix(Top_149_XGB[,-1]), label = Top_149_XGB$lymph_class, nrounds = 100, objective = "binary:logistic")

# Make predictions on validation set for AUC, sensitivity, and specificity

train_prediction_probs <- predict(xgb_model, as.matrix(valid_set_149[,-1]))
train_prediction_df <- data.frame(prob = train_prediction_probs, class = valid_set_149$lymph_class)

roc_xgb_149 <- roc(train_prediction_df$class, train_prediction_df$prob)

auc_xgb_149 <- roc(train_prediction_df$class, train_prediction_df$prob)$auc
sens_xgb_149 <- sum(train_prediction_df$prob[Test_data_frame$lymph_class == 1] >= 0.5)/sum(Test_data_frame$lymph_class == 1)
spec_xgb_149 <- sum(train_prediction_df$prob[Test_data_frame$lymph_class == 0] < 0.5)/sum(Test_data_frame$lymph_class == 0)


#MUltiPlot
par(pty = "s")
plot(roc_149, col = "#1f77b4", main = "149 Transcripts Prediction ROC", auc.polygon = FALSE, legacy.axes = TRUE, asp = NA, xlim = c(1, 0), ylim = c(0, 1))
plot(roc_svm_149, col = "#ff7f0e", add = TRUE)
plot(roc_xgb_149, col = "#2ca02c", add = TRUE )

# Add AUC values
legend("bottomright", legend = c(paste0("RF AUC = ", round(auc_149, 2)),
                                 paste0("SVM AUC = ", round(auc_svm_149, 2)),
                                 paste0("XGB AUC = ", round(auc_xgb_149, 2))),
       col = c("#1f77b4", "#ff7f0e", "#2ca02c"), lty = 1)



#54
#RANDOM FOREST
best_54 <- combined_df$Gene[1:54]
Top_54_df <- Train_data_frame_original[, c("lymph_class", best_54)]

set.seed(42)
rf_model_54 <- randomForest(lymph_class ~ ., data = Top_54_df, ntree = 500, importance = TRUE)

#Make predictions on test set for AUC, sensitivity, and specificity
predict_54 <- predict(rf_model_54, Test_data_frame, type = "prob")
predict_54_df <- data.frame(prob = predict_54[, 2], class = Test_data_frame$lymph_class)
auc_54 <- roc(predict_54_df$class, predict_54_df$prob)$auc
sens_54 <- sum(predict_54[Test_data_frame$lymph_class == 1, 2] >= 0.5)/sum(Test_data_frame$lymph_class == 1)
spec_54 <- sum(predict_54[Test_data_frame$lymph_class == 0, 2] < 0.5)/sum(Test_data_frame$lymph_class == 0)

roc_54 <- roc(predict_54_df$class, predict_54_df$prob)
par(pty = "s")
plot(roc_54, xlim = c(1, 0), col = "blue", main = "54 Transcript Prediction ROC", auc.polygon = TRUE, legacy.axes = TRUE, asp = NA)
text(0.3, 0.3, paste("AUC =", round(auc_54, 2)), col = "blue", cex = 1)
text(0.3, 0.24, paste("Sensitivity =", round(sens_54, 2)), col = "blue", cex = 1)
text(0.3, 0.18, paste("Specificity =", round(spec_54, 2)), col = "blue", cex = 1)

# Train SVM model
set.seed(42)
svm_model_54 <- svm(lymph_class ~ ., data = Top_54_df, kernel = "sigmoid", 
                    type = "C-classification", probability = TRUE,
                    maxiter = 1000)
# Make predictions on test set for AUC, sensitivity, and specificity
predict_54 <- predict(svm_model_54, Test_data_frame, probability = TRUE)
predict_54_df <- data.frame(prob = attr(predict_54, "probabilities")[,2], class = Test_data_frame$lymph_class)
auc_svm_54 <- roc(predict_54_df$class, predict_54_df$prob)$auc
sens_svm_54 <- sum(predict_54[Test_data_frame$lymph_class == 1] == 1)/sum(Test_data_frame$lymph_class == 1)
spec_svm_54 <- sum(predict_54[Test_data_frame$lymph_class == 0] == 0)/sum(Test_data_frame$lymph_class == 0)

roc_svm_54 <- roc(predict_54_df$class, predict_54_df$prob)
par(pty = "s")
plot(roc_54, xlim = c(1, 0), col = "blue", main = "54 Transcript Prediction ROC (SVM)", auc.polygon = TRUE, legacy.axes = TRUE, asp = NA)
text(0.3, 0.3, paste("AUC =", round(auc_svm_54, 2)), col = "blue", cex = 1)
text(0.3, 0.24, paste("Sensitivity =", round(sens_svm_54, 2)), col = "blue", cex = 1)
text(0.3, 0.18, paste("Specificity =", round(spec_svm_54, 2)), col = "blue", cex = 1)

# XGBoost
best_54 <- combined_df$Gene[1:54]
Top_54_XGB <- Train_data_frame_original[,c("lymph_class", best_54) ]
valid_set_54 <- Test_data_frame[,c("lymph_class", best_54) ]

Top_54_XGB$lymph_class <- ifelse(Top_54_XGB$lymph_class == 1, 1, 0)
valid_set_54$lymph_class <- ifelse(valid_set_54$lymph_class == 1, 1, 0)

# Train XGBoost model on training set
xgb_model <- xgboost(data = as.matrix(Top_54_XGB[,-1]), label = Top_54_XGB$lymph_class, nrounds = 100, objective = "binary:logistic")

# Make predictions on validation set for AUC, sensitivity, and specificity
train_prediction_probs <- predict(xgb_model, as.matrix(valid_set_54[,-1]))
train_prediction_df <- data.frame(prob = train_prediction_probs, class = valid_set_54$lymph_class)

roc_xgb_54 <- roc(train_prediction_df$class, train_prediction_df$prob)

auc_xgb_54 <- roc(train_prediction_df$class, train_prediction_df$prob)$auc
sens_xgb_54 <- sum(train_prediction_df$prob[Test_data_frame$lymph_class == 1] >= 0.5)/sum(Test_data_frame$lymph_class == 1)
spec_xgb_54 <- sum(train_prediction_df$prob[Test_data_frame$lymph_class == 0] < 0.5)/sum(Test_data_frame$lymph_class == 0)

# MultiPlot
par(pty = "s")
plot(roc_54, col = "#1f77b4", main = "54 Transcripts Prediction ROC", auc.polygon = FALSE, legacy.axes = TRUE, asp = NA, xlim = c(1, 0), ylim = c(0, 1))
plot(roc_svm_54, col = "#ff7f0e", add = TRUE)
plot(roc_xgb_54, col = "#2ca02c", add = TRUE )

# Add AUC values
legend("bottomright", legend = c(paste0("RF AUC = ", round(auc_54, 2)),
                                 paste0("SVM AUC = ", round(auc_svm_54, 2)),
                                 paste0("XGB AUC = ", round(auc_xgb_54, 2))),
       col = c("#1f77b4", "#ff7f0e", "#2ca02c"), lty = 1)





#####################                               GENE ANALYSIS                          ############################3

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("clusterProfiler", dependencies = TRUE, force = TRUE)
BiocManager::install("BiomaRt", dependancies = TRUE)
library(biomaRt)
library(clusterProfiler)

#54 
mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

Top149 <- combined_df$Gene[1:149]
Top149 <- sapply(strsplit(Top149, split = "\\."), "[", 1)
Top149 <- as.data.frame(Top149)
write_xlsx(Top149,"..\\Top 149 Expression.xlsx")


##### FULL COMBINED SCRIPT ######

#Expression
library(readr)
Expression_data_frame_for_combined_analyis <- read_csv("Expression data frame for combined analyis.csv")
View(Expression_data_frame_for_combined_analyis)
Expression_data_frame <- as.data.frame(Expression_data_frame_for_combined_analyis)
rownames(Expression_data_frame) <- Expression_data_frame[,1]
Expression_data_frame <- Expression_data_frame[,-1]

#Clinical   (from Method 1 clinical code)
library(readxl)
Clinical_Info_after_pre_processing_method_1 <- read_excel("Clinical Info after pre-processing method 1.xlsx")
View(Clinical_Info_after_pre_processing_method_1)
Clinical_data_frame <- as.data.frame(Clinical_Info_after_pre_processing_method_1)
rownames(Clinical_data_frame) <- Clinical_data_frame$patients_IDs
Clinical_data_frame <- Clinical_data_frame[,-79]



#CREATE NEW MERGED DATA FRAME

common_rows <- intersect(row.names(Clinical_data_frame), row.names(Expression_data_frame))

# Extract the common rows from each data frame using indexing
dfExpression_common <- Expression_data_frame[common_rows, ]
dfClinical_common <- Clinical_data_frame[common_rows, ]

# Merge the two data frames using cbind() to combine the columns
merged_data_frame <- cbind(dfClinical_common, dfExpression_common)
merged_data_frame <- merged_data_frame %>% mutate_if(is.character, as.factor)


#########################         SPLIT INTO TRAIN AND TEST SET USING STRATIFIED SPLIT              ######################3

#Perform stratified split of the data into 70% 30% Training and test repectively.
library(caret)
set.seed(42) #random seed for reproducibility
trainIndex <- createDataPartition(merged_data_frame$lymph_class, p = .7,
                                  list = FALSE,
                                  times = 1)
Train_combined <- merged_data_frame[ trainIndex,]
Test_combined <- merged_data_frame[-trainIndex,]  ##DONT TOUCH THIS UNTIL THE END







###ROBUST RFE 

set.seed(17)

# create empty list to store results
results_importance_scaled <- list()

# create vector of percentages
percentages <- seq(0.8, 0.1, by = -0.1)

# stratify based on lymph_class column
lymph_0_scaled <- Train_combined[Train_combined$lymph_class == 0, ]
lymph_1_scaled <- Train_combined[Train_combined$lymph_class == 1, ]

library(progress)

# Initialize progress bar
pb <- progress_bar$new(total = length(percentages))

# loop over each percentage
for (p in percentages) {
  # randomly sample p% of patients for each lymph class
  print(p)
  n_0 <- round(nrow(lymph_0_scaled)*p)
  n_1 <- round(nrow(lymph_1_scaled)*p)
  
  sample_rows_0 <- sample(nrow(lymph_0_scaled), n_0)
  sampled_data_0 <- lymph_0_scaled[sample_rows_0, ]
  
  sample_rows_1 <- sample(nrow(lymph_1_scaled), n_1)
  sampled_data_1 <- lymph_1_scaled[sample_rows_1, ]
  
  # combine samples
  sampled_data_scaled <- rbind(sampled_data_0, sampled_data_1)
  
  # set up control function for rfe
  ctrl <- rfeControl(
    method = "repeatedcv",
    repeats = 5,
    verbose = FALSE,
    returnResamp = "final",
    functions = rfFuncs,
  )
  
  # perform rfe
  rfe_result <- rfe(
    x = sampled_data_scaled[, -4],
    y = sampled_data_scaled[, 4],
    sizes = c(1),
    rfeControl = ctrl
  )
  
  # get feature importance scores
  importance_scaled <- varImp(rfe_result, scale = TRUE)  
  # get optimal subset of variables
  optimal_subset_scaled <- rfe_result$optVariables
  
  # create data frame with feature importance scores
  importance_df_scaled <- data.frame(importance_scaled[, 1], row.names = rownames(importance_scaled))
  colnames(importance_df_scaled) <- c("Feature Importance")
  
  # store importance data frame in list of data frames
  results_importance_scaled[[as.character(p*100)]] <- importance_df_scaled
  
  pb$tick()
  
}

###################      GENERATING THE TABLE WITH IMPORTANCE SCORES FOR EACH Variable ACROSS ALL SUBSETS   ############3

# create empty data frame with one column for each subset
install.packages("tidyverse")
library(tidyverse)

# Set rownames as a regular column in each data frame
results_importance_scaled <- lapply(results_importance_scaled, function(x) {
  rownames_to_column(x, "Variable")  
})


install.packages("purrr")
library(purrr)
library(dplyr)

# Add a new column to each data frame in the list with the same name as the data frame
names(results_importance_scaled) <- paste0("df_", seq_along(results_importance_scaled))

# Combine all data frames into one
combined_df_scaled <- reduce(results_importance_scaled, full_join, by = "Variable")

new_colnames <- c("Variable", "80%", "70%", "60%", "50%", "40%", "30%", "20%", "10%")

colnames(combined_df_scaled) <- new_colnames

library(dplyr)

combined_df_scaled$average <- rowMeans(combined_df_scaled[, 2:9], na.rm = TRUE)

combined_df_scaled <- combined_df_scaled %>% arrange(desc(average))



##########    PLOT FEATURES IMPORTANCE
clinical_variables <- colnames(Clinical_data_frame)[-4]

combined_df_scaled_ordered <- combined_df_scaled[order(combined_df_scaled$average, decreasing = TRUE),]

combined_df_scaled_ordered$color <- ifelse(combined_df_scaled_ordered$Variable %in% clinical_variables, "red", "steelblue")


# Plot the data with colored lines
ggplot(combined_df_scaled_ordered, aes(x = average, y = reorder(Variable, average), fill = color)) +
  geom_bar(stat = "identity", width = 0.7, size = 0.5, show.legend = FALSE,
           aes(color = ifelse(color == "red", "Clinical", "Non-clinical"))) +
  scale_color_manual(values = c("Clinical" = "red", "Non-clinical" = "steelblue")) +
  xlab("Average Score") +
  ylab("") +
  scale_y_discrete(labels = NULL) +
  scale_fill_identity() +
  theme(axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 16)) +
  ggtitle("Variable Importance Scores") +
  guides(color = FALSE) +  # remove legend for line widths
  theme(panel.grid = element_blank(),  # remove grid lines
        panel.background = element_blank(),  # remove background color
        axis.line.x = element_line(size = 1)) +  # adjust x-axis line width
  theme(legend.position = "bottom") +  # move legend to bottom
  # Add dotted lines for clinical variables
  geom_segment(data = combined_df_scaled_ordered %>% filter(color == "red"),
               aes(x = average, xend = average + 0.1, y = Variable, yend = Variable), linetype = "dotted") +
  # Add labels for clinical variables
  geom_label(data = combined_df_scaled_ordered %>% filter(color == "red"),
             aes(x = average + 0.1, y = Variable, label = Variable), size = 4, fill = "white", color = "red", show.legend = FALSE)


############################   TRAIN MODEL ON TRAIN SET AND MAKE PREDICTIONS USING CROSS VALIDATION MDOEL TRAINING  #######3



#                               RF CROSS VALIDATION (5 FOLDS) TRAINING MODEL
library(randomForest)
library(pROC)
library(caret)

# Define number of variables to use
num_var <- 200

# Initialize empty data frames to store results
results_df_scaled <- data.frame(matrix(ncol=4, nrow=num_var))
colnames(results_df_scaled) <- c("Variables", "AUC_train", "Sensitivity_train", "Specificity_train")

results_RF_Folds_scaled <- data.frame(matrix(ncol=16, nrow=num_var))
colnames(results_RF_Folds_scaled) <- c("Variables", 
                                       paste(rep(c("AUC", "Sensitivity", "Specificity"), each=5), 1:5, sep="_"))

# Loop over number of variables to include
set.seed(42)
for (i in 1:num_var) {
  # Select top i genes
  selected_var_scaled <- combined_df_scaled$Variable[1:i]
  
  # Create Train data frame with response variable and selected genes
  Training <- merged_data_frame[, c("lymph_class", selected_var_scaled)]
  
  # Set up K-fold cross-validation
  set.seed(42)
  K <- 5
  folds <- createFolds(Training$lymph_class, K, list = TRUE)
  
  # Initialize empty vectors to store predictions and labels for each fold
  train_predictions <- numeric(nrow(Training))
  train_labels <- numeric(nrow(Training))
  auc_vec <- numeric(K)
  sens_vec <- numeric(K)
  spec_vec <- numeric(K)
  
  # Train and predict using K-fold cross-validation
  for (k in 1:K) {
    # Split data into training and validation sets
    train_indices <- unlist(folds[-k])
    valid_indices <- unlist(folds[k])
    
    train_set <- Training[train_indices, ]
    valid_set <- Training[valid_indices, ]
    
    # Train random forest model on training set
    rf_model <- randomForest(lymph_class ~ ., data=train_set, ntree=500, importance=TRUE)
    
    
    # Make predictions on training set for AUC, sensitivity, and specificity
    train_prediction_probs <- predict(rf_model, valid_set, type = "prob")
    train_prediction_df <- data.frame(prob = train_prediction_probs[,2], class = valid_set$lymph_class)
    
    
    # Calculate AUC, sensitivity, and specificity of ROC curve for validation and training sets
    auc_vec[k] <- roc(train_prediction_df$class,train_prediction_df$prob)$auc
    sens_vec[k] <- sum(train_prediction_probs[valid_set$lymph_class==1,2] >= 0.5)/sum(valid_set$lymph_class==1)
    spec_vec[k] <- sum(train_prediction_probs[valid_set$lymph_class==0,2] < 0.5)/sum(valid_set$lymph_class==0)
  }
  
  # Calculate average AUC, sensitivity, and specificity across all folds
  avg_auc_train <- mean(auc_vec)
  avg_sens_train <- mean(sens_vec)
  avg_spec_train <- mean(spec_vec)
  
  # Store results in results_df
  results_df_scaled[i, "Variables"] <- i
  results_df_scaled[i, "AUC_train"] <- avg_auc_train
  results_df_scaled[i, "Sensitivity_train"] <- avg_sens_train
  results_df_scaled[i, "Specificity_train"] <- avg_spec_train
  
  # Store AUC, sensitivity, and specificity for each fold in results_RF_Folds
  results_RF_Folds_scaled[i, c(paste0("AUC_", 1:5))] <- auc_vec
  results_RF_Folds_scaled[i, c(paste0("Sensitivity_", 1:5))] <- sens_vec
  results_RF_Folds_scaled[i, c(paste0("Specificity_", 1:5))] <- spec_vec
}


#plot
results_RF_Folds_scaled$Variables <- seq(1, num_var)

library(ggplot2)

# Reshape the results_RF_Folds dataframe to long format for easier plotting
results_RF_Folds_expanded <- reshape2::melt(results_RF_Folds_scaled, id.vars=c("Variables"), variable.name="Fold", value.name="Value")

# Plot the results using ggplot
ggplot() +
  # Add the AUC values for each fold as points
  geom_point(data=results_RF_Folds_expanded, aes(x=Variables, y=Value, color=Fold)) +
  # Add the average AUC for each number of variables as a line
  geom_line(data=results_df_scaled, aes(x=Variables, y=AUC_train)) +
  # Add axis labels and title
  labs(x="Number of variables", y="AUC", title="Random Forest AUC vs Number of Variables")


#Individual plots

#AUC
library(ggplot2)

# create a long format of results_RF_Folds
results_RF_Folds_expanded_AUC <- reshape2::melt(results_RF_Folds_scaled, id.vars = "Variables", 
                                                variable.name = "Fold", value.name = "AUC")
# filter only the AUC columns
results_RF_Folds_expanded_AUC <- results_RF_Folds_expanded_AUC %>%
  filter(str_detect(Fold, "AUC"))
# convert Variables column to numeric
results_RF_Folds_expanded_AUC$Variables <- as.numeric(results_RF_Folds_expanded_AUC$Variables)
# calculate the mean AUC for each Variables value
results_RF_Folds_mean_scaled <- aggregate(AUC ~ Variables, data = results_RF_Folds_expanded_AUC, mean)
# plot the AUC curve
ggplot(results_RF_Folds_expanded_AUC, aes(x = Variables, y = AUC)) +
  geom_point(aes(color = Fold), alpha = 0.5) +
  geom_line(data = results_RF_Folds_mean_scaled, aes(x = Variables, y = AUC), size = 1.2) +
  labs(x = "Number of Variables", y = "AUC") +
  theme_classic() +
  theme(plot.title = element_blank(),
        axis.text.x = element_text(size = 24, color = "black", face = "bold"),
        axis.text.y = element_text(size = 24, color = "black", face = "bold"),
        axis.title.x = element_text(size = 14, color = "black", face = "bold"),
        axis.title.y = element_text(size = 14, color = "black", face = "bold"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(size = 1),
        plot.margin = unit(c(1, 1, 0, 1), "cm")) + # increase left and right margins
  guides(color = FALSE) +
  scale_x_continuous(expand = c(0, 0))+
  scale_y_continuous(limits = c(0.4, 0.9), expand = c(0, 0))
geom_vline(xintercept = 54, linetype = "dotted", color = "red") +
  geom_segment(aes(x = 54, y = 0.5, xend = 54, yend = 0.85), linetype = "dotted", color = "red")+
  annotate("text", x = 69, y = 0.5, label = "X = 54", color = "red", vjust = -1, hjust = 1)




#Sensitivity
# create a long format of results_RF_Folds
results_RF_Folds_expanded_sens <- reshape2::melt(results_RF_Folds_scaled, id.vars = "Variables", 
                                                 variable.name = "Fold", value.name = "Sensitivity")
# filter only the Sensitivity columns
results_RF_Folds_expanded_sens <- results_RF_Folds_expanded_sens %>%
  filter(str_detect(Fold, "Sensitivity"))
# convert Variables column to numeric
results_RF_Folds_expanded_sens$Variables <- as.numeric(results_RF_Folds_expanded_sens$Variables)
# calculate the mean Sensitivity for each Variables value
results_RF_Folds_mean_Sens_scaled <- aggregate(Sensitivity ~ Variables, data = results_RF_Folds_expanded_sens, mean)
# plot the Sensitivity curve
ggplot(results_RF_Folds_expanded_sens, aes(x = Variables, y = Sensitivity)) +
  geom_point(aes(color = Fold), alpha = 0.5) +
  geom_line(data = results_RF_Folds_mean_Sens_scaled, aes(x = Variables, y = Sensitivity), size = 1.2) +
  labs(x = "Number of Variables", y = "Sensitivity") +
  theme_classic() +
  theme(plot.title = element_blank(),
        axis.text.x = element_text(size = 24, color = "black", face = "bold"),
        axis.text.y = element_text(size = 24, color = "black", face = "bold"),
        axis.title.x = element_text(size = 14, color = "black", face = "bold"),
        axis.title.y = element_text(size = 14, color = "black", face = "bold"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(size = 1),
        plot.margin = unit(c(0, 1, 0, 1), "cm")) + # increase left and right margins
  guides(color = FALSE) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(limits = c(0.4, 1))
geom_vline(xintercept = 54, linetype = "dotted", color = "red") +
  geom_segment(aes(x = 54, y = 0.5, xend = 54, yend = 0.85), linetype = "dotted", color = "red")+
  annotate("text", x = 69, y = 0.5, label = "X = 54", color = "red", vjust = -1, hjust = 1)



#Specificity
# create a long format of results_RF_Folds
results_RF_Folds_expanded_spec <- reshape2::melt(results_RF_Folds_scaled, id.vars = "Variables", 
                                                 variable.name = "Fold", value.name = "Specificity")
# filter only the Specificity columns
results_RF_Folds_expanded_spec <- results_RF_Folds_expanded_spec %>%
  filter(str_detect(Fold, "Specificity"))
# convert Variables column to numeric
results_RF_Folds_expanded_spec$Variables <- as.numeric(results_RF_Folds_expanded_spec$Variables)
# calculate the mean Specificity for each Variables value
results_RF_Folds_mean_Spec_scaled <- aggregate(Specificity ~ Variables, data = results_RF_Folds_expanded_spec, mean)
# plot the Specificity curve
ggplot(results_RF_Folds_expanded_spec, aes(x = Variables, y = Specificity)) +
  geom_point(aes(color = Fold), alpha = 0.5) +
  geom_line(data = results_RF_Folds_mean_Spec_scaled, aes(x = Variables, y = Specificity), size = 1.2) +
  labs(x = "Number of Variables", y = "Specificity") +
  theme_classic() +
  theme(plot.title = element_blank(),
        axis.text.x = element_text(size = 24, color = "black", face = "bold"),
        axis.text.y = element_text(size = 24, color = "black", face = "bold"),
        axis.title.x = element_text(size = 14, color = "black", face = "bold"),
        axis.title.y = element_text(size = 14, color = "black", face = "bold"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(size = 1),
        plot.margin = unit(c(0, 1, 0, 1), "cm")) + # increase left and right margins
  guides(color = FALSE) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(limits = c(0.3, 0.8))
geom_vline(xintercept = 54, linetype = "dotted", color = "red") +
  geom_segment(aes(x = 54, y = 0.5, xend = 54, yend = 0.75), linetype = "dotted", color = "red")+
  annotate("text", x = 69, y = 0.4, label = "X = 54", color = "red", vjust = -1, hjust = 1)







#############################             NEURAL NETWORK MODEL   (Doesnt Work well)              ########################################3

# Load required packages
library(neuralnet)
library(pROC)


# Define number of variables to use
num_var <-50

# Initialize empty data frames to store results
results_df_NN <- data.frame(matrix(ncol=4, nrow=num_var))
colnames(results_df_NN) <- c("Variables", "AUC_train", "Sensitivity_train", "Specificity_train")

results_NN_Folds <- data.frame(matrix(ncol=16, nrow=num_var))
colnames(results_NN_Folds) <- c("Variables", 
                                paste(rep(c("AUC", "Sensitivity", "Specificity"), each=5), 1:5, sep="_"))

# Loop over number of variables to include
set.seed(42)
for (i in 1:num_var) {
  # Select top i genes
  selected_var_scaled <- combined_df_scaled$Variable[1:(i+1)]
  
  # Create Train data frame with response variable and selected genes
  Training <- merged_data_frame[, c("lymph_class", selected_var_scaled)]
  
  # Set up K-fold cross-validation
  set.seed(42)
  K <- 5
  folds <- createFolds(Training$lymph_class, K, list = TRUE)
  
  # Initialize empty vectors to store predictions and labels for each fold
  train_predictions <- numeric(nrow(Training))
  train_labels <- numeric(nrow(Training))
  auc_vec <- numeric(K)
  sens_vec <- numeric(K)
  spec_vec <- numeric(K)
  
  # Train and predict using K-fold cross-validation
  for (k in 1:K) {
    # Split data into training and validation sets
    train_indices <- unlist(folds[-k])
    valid_indices <- unlist(folds[k])
    
    train_set <- Training[train_indices, ]
    valid_set <- Training[valid_indices, ]
    
    
    
    # Train neural network model on training set
    nn_model <- neuralnet(lymph_class ~ ., data = train_set, hidden = 10, linear.output = TRUE, threshold = 0.1, stepmax = 1e+08, act.fct = "logistic")
    
    # Make predictions on validation set for AUC, sensitivity, and specificity
    train_prediction_probs <- predict(nn_model, valid_set[,-1])
    train_prediction_df <- data.frame(prob = train_prediction_probs[,2], class = valid_set$lymph_class)
    
    # Calculate AUC, sensitivity, and specificity of ROC curve for validation and training sets
    auc_vec[k] <- roc(train_prediction_df$class,train_prediction_df$prob)$auc
    sens_vec[k] <- sum(train_prediction_probs[valid_set$lymph_class==1,2] >= 0.5)/sum(valid_set$lymph_class==1)
    spec_vec[k] <- sum(train_prediction_probs[valid_set$lymph_class==0,2] < 0.5)/sum(valid_set$lymph_class==0)
  }
  
  # Calculate average AUC, sensitivity, and specificity across all folds
  avg_auc_train <- mean(auc_vec)
  avg_sens_train <- mean(sens_vec)
  avg_spec_train <- mean(spec_vec)
  
  # Store results in results_df
  results_df_NN[i, "Variables"] <- i
  results_df_NN[i, "AUC_train"] <- avg_auc_train
  results_df_NN[i, "Sensitivity_train"] <- avg_sens_train
  results_df_NN[i, "Specificity_train"] <- avg_spec_train
  
  # Store AUC, sensitivity, and specificity for each fold in results_RF_Folds
  results_NN_Folds[i, c(paste0("AUC_", 1:5))] <- auc_vec
  results_NN_Folds[i, c(paste0("Sensitivity_", 1:5))] <- sens_vec
  results_NN_Folds[i, c(paste0("Specificity_", 1:5))] <- spec_vec
}



#plot
results_NN_Folds$Variables <- seq(1, num_var)

library(ggplot2)

# Reshape the results_RF_Folds dataframe to long format for easier plotting
results_NN_Folds_expanded <- reshape2::melt(results_NN_Folds, id.vars=c("Variables"), variable.name="Fold", value.name="Value")

# Plot the results using ggplot
ggplot() +
  # Add the AUC values for each fold as points
  geom_point(data=results_NN_Folds_expanded, aes(x=Variables, y=Value, color=Fold)) +
  # Add the average AUC for each number of variables as a line
  geom_line(data=results_df_NN, aes(x=Variables, y=AUC_train)) +
  # Add axis labels and title
  labs(x="Number of variables", y="AUC", title="Neural Network AUC vs Number of Variables")












################################3                          SVM MODEL             ############################################################################3

library(pROC)
library(caret)

# Define number of variables to use
num_var <- 200

# Initialize empty data frames to store results
results_SVM_df <- data.frame(matrix(ncol=4, nrow=num_var))
colnames(results_SVM_df) <- c("Variables", "AUC_train", "Sensitivity_train", "Specificity_train")

results_SVM_Folds_scaled <- data.frame(matrix(ncol=16, nrow=num_var))
colnames(results_SVM_Folds_scaled) <- c("Variables", 
                                        paste(rep(c("AUC", "Sensitivity", "Specificity"), each=5), 1:5, sep="_"))

# Loop over number of variables to include
set.seed(42)
for (i in 1:num_var) {
  # Select top i genes
  selected_var_scaled <- combined_df_scaled$Variable[1:i]
  
  # Create Train data frame with response variable and selected genes
  Training <- merged_data_frame[, c("lymph_class", selected_var_scaled)]
  
  # Set up K-fold cross-validation
  set.seed(42)
  K <- 5
  folds <- createFolds(Training$lymph_class, K, list = TRUE)
  
  # Initialize empty vectors to store predictions and labels for each fold
  train_predictions <- numeric(nrow(Training))
  train_labels <- numeric(nrow(Training))
  auc_vec <- numeric(K)
  sens_vec <- numeric(K)
  spec_vec <- numeric(K)
  
  # Train and predict using K-fold cross-validation
  for (k in 1:K) {
    # Split data into training and validation sets
    train_indices <- unlist(folds[-k])
    valid_indices <- unlist(folds[k])
    
    train_set <- Training[train_indices, ]
    valid_set <- Training[valid_indices, ]
    
    
    # Train SVM model on training set
    svm_model <- svm(lymph_class ~ ., data = train_set, kernel = "sigmoid", 
                     type = "C-classification", probability = TRUE,
                     maxiter = 10000)
    
    # Make predictions on training set for AUC, sensitivity, and specificity
    train_prediction_probs <- predict(svm_model, valid_set, type = "prob")
    train_prediction_probs <- as.data.frame(train_prediction_probs)
    train_prediction_df <- data.frame(prob = train_prediction_probs[,1], class = valid_set$lymph_class)
    
    # Calculate AUC, sensitivity, and specificity of ROC curve for validation and training sets
    auc_vec[k] <- roc(as.numeric(train_prediction_df$class),as.numeric(train_prediction_df$prob))$auc
    sens_vec[k] <- sum(train_prediction_df$prob == "1" & train_prediction_df$class == "1")/sum(train_prediction_df$class==1)
    spec_vec[k] <- sum(train_prediction_df$prob == "0" & train_prediction_df$class == "0")/sum(train_prediction_df$class==0)
  }
  
  # Calculate average AUC, sensitivity, and specificity across all folds
  avg_auc_train <- mean(auc_vec)
  avg_sens_train <- mean(sens_vec)
  avg_spec_train <- mean(spec_vec)
  
  # Store results in results_df
  results_SVM_df[i, "Variables"] <- i
  results_SVM_df[i, "AUC_train"] <- avg_auc_train
  results_SVM_df[i, "Sensitivity_train"] <- avg_sens_train
  results_SVM_df[i, "Specificity_train"] <- avg_spec_train
  
  # Store AUC, sensitivity, and specificity for each fold in results_RF_Folds
  results_SVM_Folds_scaled[i, c(paste0("AUC_", 1:5))] <- auc_vec
  results_SVM_Folds_scaled[i, c(paste0("Sensitivity_", 1:5))] <- sens_vec
  results_SVM_Folds_scaled[i, c(paste0("Specificity_", 1:5))] <- spec_vec
}


#plot
results_SVM_Folds_scaled$Variables <- seq(1, num_var)


# Reshape the results_RF_Folds dataframe to long format for easier plotting
results_SVM_Folds_expanded <- reshape2::melt(results_SVM_Folds_scaled, id.vars=c("Variables"), variable.name="Fold", value.name="Value")

# Plot the results using ggplot
ggplot() +
  # Add the AUC values for each fold as points
  geom_point(data=results_SVM_Folds_expanded, aes(x=Variables, y=Value, color=Fold), alpha = 0.5) +
  # Add the average AUC for each number of variables as a line
  geom_line(data=results_SVM_df, aes(x=Variables, y=AUC_train), cex = 1.5) +
  # Add axis labels and title
  labs(x="Number of variables", y="AUC", title="SVM AUC vs Number of Variables")


#AUC
# create a long format of results_RF_Folds
results_SVM_Folds_expanded_AUC <- reshape2::melt(results_SVM_Folds_scaled, id.vars = "Variables", 
                                                 variable.name = "Fold", value.name = "AUC")
# filter only the AUC columns
results_SVM_Folds_expanded_AUC <- results_SVM_Folds_expanded_AUC %>%
  filter(str_detect(Fold, "AUC"))
# convert Variables column to numeric
results_SVM_Folds_expanded_AUC$Variables <- as.numeric(results_SVM_Folds_expanded_AUC$Variables)
# calculate the mean AUC for each Variables value
results_SVM_Folds_mean_scaled <- aggregate(AUC ~ Variables, data = results_SVM_Folds_expanded_AUC, mean)
# plot the AUC curve
ggplot(results_SVM_Folds_expanded_AUC, aes(x = Variables, y = AUC)) +
  geom_point(aes(color = Fold), alpha = 0.5) +
  geom_line(data = results_SVM_Folds_mean_scaled, aes(x = Variables, y = AUC), size = 1.2) +
  labs(x = "Number of Variables", y = "AUC") +
  theme_classic() +
  theme(plot.title = element_blank(),
        axis.text.x = element_text(size = 24, color = "black", face = "bold"),
        axis.text.y = element_text(size = 24, color = "black", face = "bold"),
        axis.title.x = element_text(size = 14, color = "black", face = "bold"),
        axis.title.y = element_text(size = 14, color = "black", face = "bold"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(size = 1),
        plot.margin = unit(c(1, 1, 0, 1), "cm")) + # increase left and right margins
  guides(color = FALSE) +
  scale_x_continuous(expand = c(0, 0))+
  scale_y_continuous(limits = c(0.4, 0.9), expand = c(0, 0))

#leave out for now
geom_vline(xintercept = 54, linetype = "dotted", color = "red") +
  geom_segment(aes(x = 54, y = 0.5, xend = 54, yend = 0.75), linetype = "dotted", color = "red")+
  annotate("text", x = 69, y = 0.5, label = "X = 54", color = "red", vjust = -1, hjust = 1)


# Specificity
# create a long format of results_SVM_Folds
results_SVM_Folds_expanded_spec <- reshape2::melt(results_SVM_Folds_scaled, id.vars = "Variables", 
                                                  variable.name = "Fold", value.name = "Specificity")
# filter only the Specificity columns
results_SVM_Folds_expanded_spec <- results_SVM_Folds_expanded_spec %>%
  filter(str_detect(Fold, "Specificity"))
# convert Variables column to numeric
results_SVM_Folds_expanded_spec$Variables <- as.numeric(results_SVM_Folds_expanded_spec$Variables)
# calculate the mean Specificity for each Variables value
results_SVM_Folds_mean_spec <- aggregate(Specificity ~ Variables, data = results_SVM_Folds_expanded_spec, mean)
# plot the Specificity curve
ggplot(results_SVM_Folds_expanded_spec, aes(x = Variables, y = Specificity)) +
  geom_point(aes(color = Fold), alpha = 0.5) +
  geom_line(data = results_SVM_Folds_mean_spec, aes(x = Variables, y = Specificity), size = 1.2) +
  labs(x = "Number of Variables", y = "Specificity") +
  theme_classic() +
  theme(plot.title = element_blank(),
        axis.text.x = element_text(size = 24, color = "black", face = "bold"),
        axis.text.y = element_text(size = 24, color = "black", face = "bold"),
        axis.title.x = element_text(size = 14, color = "black", face = "bold"),
        axis.title.y = element_text(size = 14, color = "black", face = "bold"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(size = 1),
        plot.margin = unit(c(0, 1, 0, 1), "cm")) + # increase left and right margins
  guides(color = FALSE) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(limits = c(0.3, 0.8))
#leave out for now
geom_vline(xintercept = 54, linetype = "dotted", color = "red") +
  geom_segment(aes(x = 54, y = 0.5, xend = 54, yend = 0.7), linetype = "dotted", color = "red") +
  annotate("text", x = 69, y = 0.5, label = "X = 54", color = "red", vjust = -1, hjust = 1)



# Sensitivity
# create a long format of results_SVM_Folds
results_SVM_Folds_expanded_sens <- reshape2::melt(results_SVM_Folds_scaled, id.vars = "Variables", 
                                                  variable.name = "Fold", value.name = "Sensitivity")
# filter only the Sensitivity columns
results_SVM_Folds_expanded_sens <- results_SVM_Folds_expanded_sens %>%
  filter(str_detect(Fold, "Sensitivity"))
# convert Variables column to numeric
results_SVM_Folds_expanded_sens$Variables <- as.numeric(results_SVM_Folds_expanded_sens$Variables)
# calculate the mean Sensitivity for each Variables value
results_SVM_Folds_mean_sens <- aggregate(Sensitivity ~ Variables, data = results_SVM_Folds_expanded_sens, mean)
# plot the Sensitivity curve
ggplot(results_SVM_Folds_expanded_sens, aes(x = Variables, y = Sensitivity)) +
  geom_point(aes(color = Fold), alpha = 0.5) +
  geom_line(data = results_SVM_Folds_mean_sens, aes(x = Variables, y = Sensitivity), size = 1.2) +
  labs(x = "Number of Variables", y = "Sensitivity",
       title = "Sensitivity across Top 200 Variables (SVM)") +
  theme_classic() +
  theme(plot.title = element_blank(),
        axis.text.x = element_text(size = 24, color = "black", face = "bold"),
        axis.text.y = element_text(size = 24, color = "black", face = "bold"),
        axis.title.x = element_text(size = 14, color = "black", face = "bold"),
        axis.title.y = element_text(size = 14, color = "black", face = "bold"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(size = 1),
        plot.margin = unit(c(0, 1, 0, 1), "cm")) + # increase left and right margins
  guides(color = FALSE) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(limits = c(0.4, 1))
#Leave out for now
geom_vline(xintercept = 54, linetype = "dotted", color = "red") +
  geom_segment(aes(x = 54, y = 0.5, xend = 54, yend = 0.7), linetype = "dotted", color = "red") +
  annotate("text", x = 69, y = 0.5, label = "X = 54", color = "red", vjust = -1, hjust = 1)






#DECIDE BEST SUBSET
results_SVM_df$Average <- (results_SVM_df$AUC_train + results_SVM_df$Sensitivity_train + results_SVM_df$Specificity_train)/3










################################                   XGBOOST MODEL              #####################################333

# Load required packages
library(xgboost)
library(pROC)

# Define number of variables to use
num_var <- 200

# Initialize empty data frames to store results
results_df_XGB <- data.frame(matrix(ncol=4, nrow=num_var))
colnames(results_df_XGB) <- c("Variables", "AUC_train", "Sensitivity_train", "Specificity_train")

results_XGB_Folds <- data.frame(matrix(ncol=16, nrow=num_var))
colnames(results_XGB_Folds) <- c("Variables", 
                                 paste(rep(c("AUC", "Sensitivity", "Specificity"), each=5), 1:5, sep="_"))

# Loop over number of variables to include
set.seed(42)
for (i in 1:num_var) {
  # Select top i genes
  selected_var_scaled <- combined_df_scaled$Variable[1:i]
  
  # Create Train data frame with response variable and selected genes
  Training <- merged_data_frame[, c("lymph_class", selected_var_scaled)]
  
  # Set up K-fold cross-validation
  set.seed(42)
  K <- 5
  folds <- createFolds(Training$lymph_class, K, list = TRUE)
  
  # Initialize empty vectors to store predictions and labels for each fold
  train_predictions <- numeric(nrow(Training))
  train_labels <- numeric(nrow(Training))
  auc_vec <- numeric(K)
  sens_vec <- numeric(K)
  spec_vec <- numeric(K)
  
  # Train and predict using K-fold cross-validation
  for (k in 1:K) {
    # Split data into training and validation sets
    train_indices <- unlist(folds[-k])
    valid_indices <- unlist(folds[k])
    
    train_set <- Training[train_indices, ]
    valid_set <- Training[valid_indices, ]
    
    train_set$lymph_class <- ifelse(train_set$lymph_class == 1, 1, 0)
    valid_set$lymph_class <- ifelse(valid_set$lymph_class == 1, 1, 0)
    
    # Train XGBoost model on training set
    xgb_model <- xgboost(data = as.matrix(train_set[,-1]), label = train_set$lymph_class, nrounds = 100, objective = "binary:logistic")
    
    # Make predictions on validation set for AUC, sensitivity, and specificity
    
    train_prediction_probs <- predict(xgb_model, as.matrix(valid_set[,-1]))
    train_prediction_df <- data.frame(prob = train_prediction_probs, class = valid_set$lymph_class)
    
    # Calculate AUC, sensitivity, and specificity of ROC curve for validation and training sets
    auc_vec[k] <- roc(train_prediction_df$class,train_prediction_df$prob)$auc
    sens_vec[k] <- sum(train_prediction_probs[valid_set$lymph_class==1] >= 0.5)/sum(valid_set$lymph_class==1)
    spec_vec[k] <- sum(train_prediction_probs[valid_set$lymph_class==0] < 0.5)/sum(valid_set$lymph_class==0)
  }
  
  # Calculate average AUC, sensitivity, and specificity across all folds
  avg_auc_train <- mean(auc_vec)
  avg_sens_train <- mean(sens_vec)
  avg_spec_train <- mean(spec_vec)
  
  # Store results in results_df
  results_df_XGB[i, "Variables"] <- i
  results_df_XGB[i, "AUC_train"] <- avg_auc_train
  results_df_XGB[i, "Sensitivity_train"] <- avg_sens_train
  results_df_XGB[i, "Specificity_train"] <- avg_spec_train
  
  # Store AUC, sensitivity, and specificity for each fold in results_RF_Folds
  results_XGB_Folds[i, c(paste0("AUC_", 1:5))] <- auc_vec
  results_XGB_Folds[i, c(paste0("Sensitivity_", 1:5))] <- sens_vec
  results_XGB_Folds[i, c(paste0("Specificity_", 1:5))] <- spec_vec
}



#plot
results_XGB_Folds$Variables <- seq(1, num_var)


# Reshape the results_RF_Folds dataframe to long format for easier plotting
results_XGB_Folds_expanded <- reshape2::melt(results_XGB_Folds, id.vars=c("Variables"), variable.name="Fold", value.name="Value")

# Plot the results using ggplot
ggplot() +
  # Add the AUC values for each fold as points
  geom_point(data=results_XGB_Folds_expanded, aes(x=Variables, y=Value, color=Fold), alpha = 0.5) +
  # Add the average AUC for each number of variables as a line
  geom_line(data=results_df_XGB, aes(x=Variables, y=AUC_train), cex = 1.5) +
  # Add axis labels and title
  labs(x="Number of variables", y="AUC", title="XGB AUC vs Number of Variables")






#AUC
# create a long format of results_RF_Folds
results_XGB_Folds_expanded_AUC <- reshape2::melt(results_XGB_Folds, id.vars = "Variables", 
                                                 variable.name = "Fold", value.name = "AUC")
# filter only the AUC columns
results_XGB_Folds_expanded_AUC <- results_XGB_Folds_expanded_AUC %>%
  filter(str_detect(Fold, "AUC"))
# convert Variables column to numeric
results_XGB_Folds_expanded_AUC$Variables <- as.numeric(results_XGB_Folds_expanded_AUC$Variables)
# calculate the mean AUC for each Variables value
results_XGB_Folds_mean_scaled <- aggregate(AUC ~ Variables, data = results_XGB_Folds_expanded_AUC, mean)
# plot the AUC curve
ggplot(results_XGB_Folds_expanded_AUC, aes(x = Variables, y = AUC)) +
  geom_point(aes(color = Fold), alpha = 0.5) +
  geom_line(data = results_XGB_Folds_mean_scaled, aes(x = Variables, y = AUC), size = 1.2) +
  labs(x = "Number of Variables", y = "AUC",
       title = "AUC curve across Top 200 Variables (XGB)") +
  theme_classic() +
  theme(plot.title = element_blank(),
        axis.text.x = element_text(size = 24, color = "black", face = "bold"),
        axis.text.y = element_text(size = 24, color = "black", face = "bold"),
        axis.title.x = element_text(size = 14, color = "black", face = "bold"),
        axis.title.y = element_text(size = 14, color = "black", face = "bold"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(size = 1),
        plot.margin = unit(c(1, 1, 0, 1), "cm")) + # increase left and right margins
  guides(color = FALSE) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(limits = c(0.4, 0.9))
#leave out for now
geom_vline(xintercept = 54, linetype = "dotted", color = "red") +
  geom_segment(aes(x = 54, y = 0.5, xend = 54, yend = 0.75), linetype = "dotted", color = "red")+
  annotate("text", x = 69, y = 0.5, label = "X = 54", color = "red", vjust = -1, hjust = 1)




#Sensitivity
# create a long format of results_RF_Folds
results_XGB_Folds_expanded_sens <- reshape2::melt(results_XGB_Folds, id.vars = "Variables", 
                                                  variable.name = "Fold", value.name = "Sensitivity")
# filter only the Sensitivity columns
results_XGB_Folds_expanded_sens <- results_XGB_Folds_expanded_sens %>%
  filter(str_detect(Fold, "Sensitivity"))
# convert Variables column to numeric
results_XGB_Folds_expanded_sens$Variables <- as.numeric(results_XGB_Folds_expanded_sens$Variables)
# calculate the mean Sensitivity for each Variables value
results_XGB_Folds_mean_sens <- aggregate(Sensitivity ~ Variables, data = results_XGB_Folds_expanded_sens, mean)
# plot the Sensitivity curve
ggplot(results_XGB_Folds_expanded_sens, aes(x = Variables, y = Sensitivity)) +
  geom_point(aes(color = Fold), alpha = 0.5) +
  geom_line(data = results_XGB_Folds_mean_sens, aes(x = Variables, y = Sensitivity), size = 1.2) +
  labs(x = "Number of Variables", y = "Sensitivity",
       title = "Sensitivity curve across Top 200 Variables(XGB)") +
  theme_classic() +
  theme(plot.title = element_blank(),
        axis.text.x = element_text(size = 24, color = "black", face = "bold"),
        axis.text.y = element_text(size = 24, color = "black", face = "bold"),
        axis.title.x = element_text(size = 14, color = "black", face = "bold"),
        axis.title.y = element_text(size = 14, color = "black", face = "bold"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(size = 1),
        plot.margin = unit(c(0, 1, 0, 1), "cm")) + # increase left and right margins
  guides(color = FALSE) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(limits = c(0.4, 1))





# Specificity
# create a long format of results_RF_Folds
results_XGB_Folds_expanded_Specificity <- reshape2::melt(results_XGB_Folds, id.vars = "Variables", 
                                                         variable.name = "Fold", value.name = "Specificity")
# filter only the Specificity columns
results_XGB_Folds_expanded_Specificity <- results_XGB_Folds_expanded_Specificity %>%
  filter(str_detect(Fold, "Specificity"))
# convert Variables column to numeric
results_XGB_Folds_expanded_Specificity$Variables <- as.numeric(results_XGB_Folds_expanded_Specificity$Variables)
# calculate the mean Specificity for each Variables value
results_XGB_Folds_mean_scaled <- aggregate(Specificity ~ Variables, data = results_XGB_Folds_expanded_Specificity, mean)
# plot the Specificity curve
ggplot(results_XGB_Folds_expanded_Specificity, aes(x = Variables, y = Specificity)) +
  geom_point(aes(color = Fold), alpha = 0.5) +
  geom_line(data = results_XGB_Folds_mean_scaled, aes(x = Variables, y = Specificity), size = 1.2) +
  labs(x = "Number of Variables", y = "Specificity") +
  theme_classic() +
  theme(plot.title = element_blank(),
        axis.text.x = element_text(size = 24, color = "black", face = "bold"),
        axis.text.y = element_text(size = 24, color = "black", face = "bold"),
        axis.title.x = element_text(size = 14, color = "black", face = "bold"),
        axis.title.y = element_text(size = 14, color = "black", face = "bold"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(size = 1),
        plot.margin = unit(c(0, 1, 0, 1), "cm")) + # increase left and right margins
  guides(color = FALSE) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(limits = c(0.3, 0.8))




#DECIDE BEST SUBSET
results_df_XGB$Average <- (results_df_XGB$AUC_train + results_df_XGB$Sensitivity_train + results_df_XGB$Specificity_train)/3





###############                MASTER TABLE               #################################################3
#Create Master  TABLE
Master_table <- left_join(results_df, results_SVM_df, by = "Variables") %>%
  left_join(results_df_XGB, by = "Variables")
Master_table <- Master_table[,c(-8,-12)]

Master_table <- Master_table %>%
  select(Variables,
         contains("AUC"),
         contains("Specificity"),
         contains("Sensitivity"))

master_colnames <- c("Features", "RF AUC","SVM AUC","XGB AUC", "RF Sens", "SVM Sens", "XGB Sens", "RF Spec", "SVM Spec", "XGB Spec")
colnames(Master_table) <- master_colnames


Master_table[is.na(Master_table)] <- 0

Master_table$average <- rowSums(Master_table[,2:10])/9

write_xlsx(Master_table,"..\\Master Table Combined.xlsx")



######################        MAKE PREDICTIONS USING TEST SET         ################################################3

#DECIDE BEST SUBSET FROM MASTER TABLE

#Best is 120 features

best_120 <- combined_df_scaled$Variable[1:120]
Top_120_df <- Train_combined[, c("lymph_class", best_120)]

set.seed(42)
rf_model_120 <- randomForest(lymph_class ~ ., data = Top_120_df, ntree = 500, importance = TRUE)

# Make predictions on test set for AUC, sensitivity, and specificity
predict_120 <- predict(rf_model_120, Test_combined, type = "prob")
predict_120_df <- data.frame(prob = predict_120[, 2], class = Test_combined$lymph_class)
auc_120 <- roc(predict_120_df$class, predict_120_df$prob)$auc
sens_120 <- sum(predict_120[Test_combined$lymph_class == 1, 2] >= 0.5)/sum(Test_combined$lymph_class == 1)
spec_120 <- sum(predict_120[Test_combined$lymph_class == 0, 2] < 0.5)/sum(Test_combined$lymph_class == 0)

roc_120 <- roc(predict_120_df$class, predict_120_df$prob)
par(pty = "s")
plot(roc_120, xlim = c(1, 0), col = "blue", main = "120 Transcript Prediction ROC", auc.polygon = TRUE, legacy.axes = TRUE, asp = NA)
text(0.3, 0.3, paste("AUC =", round(auc_120, 2)), col = "blue", cex = 1)
text(0.3, 0.24, paste("Sensitivity =", round(sens_120, 2)), col = "blue", cex = 1)
text(0.3, 0.18, paste("Specificity =", round(spec_120, 2)), col = "blue", cex = 1)



# Train SVM model
set.seed(42)
svm_model_120 <- svm(lymph_class ~ ., data = Top_120_df, kernel = "sigmoid", 
                     type = "C-classification", probability = TRUE,
                     maxiter = 1000)
# Make predictions on test set for AUC, sensitivity, and specificity
predict_120 <- predict(svm_model_120, Test_combined, probability = TRUE)
predict_120_df <- data.frame(prob = attr(predict_120, "probabilities")[,2], class = Test_combined$lymph_class)
auc_svm_120 <- roc(predict_120_df$class, predict_120_df$prob)$auc
sens_svm_120 <- sum(predict_120[Test_combined$lymph_class == 1] == 1)/sum(Test_combined$lymph_class == 1)
spec_svm_120 <- sum(predict_120[Test_combined$lymph_class == 0] == 0)/sum(Test_combined$lymph_class == 0)

roc_svm_120 <- roc(predict_120_df$class, predict_120_df$prob)
par(pty = "s")
plot(roc_120, xlim = c(1, 0), col = "blue", main = "120 Transcript Prediction ROC (SVM)", auc.polygon = TRUE, legacy.axes = TRUE, asp = NA)
text(0.3, 0.3, paste("AUC =", round(auc_svm_120, 2)), col = "blue", cex = 1)
text(0.3, 0.24, paste("Sensitivity =", round(sens_svm_120, 2)), col = "blue", cex = 1)
text(0.3, 0.18, paste("Specificity =", round(spec_svm_120, 2)), col = "blue", cex = 1)


#XGBOOST

Top_120_XGB <- Train_combined[,c("lymph_class", best_120) ]
valid_set_120 <- Test_combined[,c("lymph_class", best_120) ]

Top_120_XGB$lymph_class <- ifelse(Top_120_XGB$lymph_class == 1, 1, 0)
valid_set_120$lymph_class <- ifelse(valid_set_120$lymph_class == 1, 1, 0)

# Train XGBoost model on training set
xgb_model <- xgboost(data = as.matrix(Top_120_XGB[,-1]), label = Top_120_XGB$lymph_class, nrounds = 100, objective = "binary:logistic")

# Make predictions on validation set for AUC, sensitivity, and specificity

train_prediction_probs <- predict(xgb_model, as.matrix(valid_set_120[,-1]))
train_prediction_df <- data.frame(prob = train_prediction_probs, class = valid_set_120$lymph_class)

roc_xgb_120 <- roc(train_prediction_df$class, train_prediction_df$prob)

auc_xgb_120 <- roc(train_prediction_df$class, train_prediction_df$prob)$auc
sens_xgb_120 <- sum(train_prediction_df$prob[Test_combined$lymph_class == 1] >= 0.5)/sum(Test_combined$lymph_class == 1)
spec_xgb_120 <- sum(train_prediction_df$prob[Test_combined$lymph_class == 0] < 0.5)/sum(Test_combined$lymph_class == 0)


#MUltiPlot
par(pty = "s")
plot(roc_120, col = "#1f77b4", main = "120 Transcripts Prediction ROC", auc.polygon = FALSE, legacy.axes = TRUE, asp = NA, xlim = c(1, 0), ylim = c(0, 1))
plot(roc_svm_120, col = "#ff7f0e", add = TRUE)
plot(roc_xgb_120, col = "#2ca02c", add = TRUE )

# Add AUC values
legend("bottomright", legend = c(paste0("RF AUC = ", round(auc_120, 2)),
                                 paste0("SVM AUC = ", round(auc_svm_120, 2)),
                                 paste0("XGB AUC = ", round(auc_xgb_120, 2))),
       col = c("#1f77b4", "#ff7f0e", "#2ca02c"), lty = 1)



#94
#RF
best_94 <- combined_df_scaled$Variable[1:94]
Top_94_df <- Train_combined[, c("lymph_class", best_94)]

set.seed(42)
rf_model_94 <- randomForest(lymph_class ~ ., data = Top_94_df, ntree = 500, importance = TRUE)

# Make predictions on test set for AUC, sensitivity, and specificity
predict_94 <- predict(rf_model_94, Test_combined, type = "prob")
predict_94_df <- data.frame(prob = predict_94[, 2], class = Test_combined$lymph_class)
auc_94 <- roc(predict_94_df$class, predict_94_df$prob)$auc
sens_94 <- sum(predict_94[Test_combined$lymph_class == 1, 2] >= 0.5)/sum(Test_combined$lymph_class == 1)
spec_94 <- sum(predict_94[Test_combined$lymph_class == 0, 2] < 0.5)/sum(Test_combined$lymph_class == 0)

roc_94 <- roc(predict_94_df$class, predict_94_df$prob)
par(pty = "s")
plot(roc_94, xlim = c(1, 0), col = "blue", main = "94 Transcript Prediction ROC", auc.polygon = TRUE, legacy.axes = TRUE, asp = NA)
text(0.3, 0.3, paste("AUC =", round(auc_94, 2)), col = "blue", cex = 1)
text(0.3, 0.24, paste("Sensitivity =", round(sens_94, 2)), col = "blue", cex = 1)
text(0.3, 0.18, paste("Specificity =", round(spec_94, 2)), col = "blue", cex = 1)


# Train SVM model
set.seed(42)
svm_model_94 <- svm(lymph_class ~ ., data = Top_94_df, kernel = "sigmoid", 
                    type = "C-classification", probability = TRUE,
                    maxiter = 1000)
# Make predictions on test set for AUC, sensitivity, and specificity
predict_94 <- predict(svm_model_94, Test_combined, probability = TRUE)
predict_94_df <- data.frame(prob = attr(predict_94, "probabilities")[,2], class = Test_combined$lymph_class)
auc_svm_94 <- roc(predict_94_df$class, predict_94_df$prob)$auc
sens_svm_94 <- sum(predict_94[Test_combined$lymph_class == 1] == 1)/sum(Test_combined$lymph_class == 1)
spec_svm_94 <- sum(predict_94[Test_combined$lymph_class == 0] == 0)/sum(Test_combined$lymph_class == 0)

roc_svm_94 <- roc(predict_94_df$class, predict_94_df$prob)
par(pty = "s")
plot(roc_94, xlim = c(1, 0), col = "blue", main = "94 Transcript Prediction ROC (SVM)", auc.polygon = TRUE, legacy.axes = TRUE, asp = NA)
text(0.3, 0.3, paste("AUC =", round(auc_svm_94, 2)), col = "blue", cex = 1)
text(0.3, 0.24, paste("Sensitivity =", round(sens_svm_94, 2)), col = "blue", cex = 1)
text(0.3, 0.18, paste("Specificity =", round(spec_svm_94, 2)), col = "blue", cex = 1)



#XGBOOST

Top_94_XGB <- Train_combined[,c("lymph_class", best_94) ]
valid_set_94 <- Test_combined[,c("lymph_class", best_94) ]

Top_94_XGB$lymph_class <- ifelse(Top_94_XGB$lymph_class == 1, 1, 0)
valid_set_94$lymph_class <- ifelse(valid_set_94$lymph_class == 1, 1, 0)

#Train XGBoost model on training set
xgb_model <- xgboost(data = as.matrix(Top_94_XGB[,-1]), label = Top_94_XGB$lymph_class, nrounds = 100, objective = "binary:logistic")

#Make predictions on validation set for AUC, sensitivity, and specificity
train_prediction_probs <- predict(xgb_model, as.matrix(valid_set_94[,-1]))
train_prediction_df <- data.frame(prob = train_prediction_probs, class = valid_set_94$lymph_class)

roc_xgb_94 <- roc(train_prediction_df$class, train_prediction_df$prob)

auc_xgb_94 <- roc(train_prediction_df$class, train_prediction_df$prob)$auc
sens_xgb_94 <- sum(train_prediction_df$prob[Test_combined$lymph_class == 1] >= 0.5)/sum(Test_combined$lymph_class == 1)
spec_xgb_94 <- sum(train_prediction_df$prob[Test_combined$lymph_class == 0] < 0.5)/sum(Test_combined$lymph_class == 0)

#MUltiPlot
par(pty = "s")
plot(roc_94, col = "#1f77b4", main = "94 Transcripts Prediction ROC", auc.polygon = FALSE, legacy.axes = TRUE, asp = NA, xlim = c(1, 0), ylim = c(0, 1))
plot(roc_svm_94, col = "#ff7f0e", add = TRUE)
plot(roc_xgb_94, col = "#2ca02c", add = TRUE )

#Add AUC values
legend("bottomright", legend = c(paste0("RF AUC = ", round(auc_94, 2)),
                                 paste0("SVM AUC = ", round(auc_svm_94, 2)),
                                 paste0("XGB AUC = ", round(auc_xgb_94, 2))),
       col = c("#1f77b4", "#ff7f0e", "#2ca02c"), lty = 1)




#############################    GENE PATHWAY AD FURTHER DOWNSTREAM ANALYSIS    #################################3
Top_54_genes <- combined_df_scaled$Variable[1:54]
Top_54_genes <- sapply(strsplit(Top_54_genes, split = "\\."), "[",1)
Top_54_genes <- as.data.frame(Top_54_genes)
write_xlsx(Top_54_genes,"..\\Top_54_genes.xlsx")


Top_200_genes <- combined_df_scaled$Variable[1:200]
Top_200_genes <- sapply(strsplit(Top_200_genes, split = "\\."), "[",1)
Top_200_genes <- as.data.frame(Top_200_genes)
write_xlsx(Top_200_genes,"..\\Top_200_genes.xlsx")
