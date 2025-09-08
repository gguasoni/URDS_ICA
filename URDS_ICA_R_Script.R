#Installing and loading any necessary packages

#Make sure the relevant files are in the working directory
path <- setwd(getwd())

library(ggplot2)
library(tidyverse)
library(pheatmap)
library(tools)

#The code was adapted from the following site(s):
#https://stackoverflow.com/questions/17499013/how-do-i-make-a-list-of-data-frames

# Methods: Defensive Coding
## Loading the Data Files
# A function was created to ensure that all the relevant files were converted to data frames that could later be manipulated in RStudio. If no .csv or .txt files were present in the directory, the function would stop and return a message. 
#Loading the data sets and creating variables
#Creating a list vector of all the .csv and .txt files
relevant_data  <- list.files(pattern = "\\.(csv|txt)$")

#Creating a function that will read the files in 'relevant_data' and turn them into data frames
create_dataframes <- function(relevant_data){
  if (length(relevant_data) == 0){
    stop("No data to analyse.")
  } else if (length(relevant_data) > 0) {
    for (i in seq_along(relevant_data)){
      if (grepl("\\.csv", relevant_data[i])) {
        message(paste0(relevant_data[i], " is a .csv") )
        assign(paste0(file_path_sans_ext(relevant_data[i]), "_df"), read.csv(relevant_data[i]), envir = .GlobalEnv)
      } else if (grepl("\\.txt", relevant_data[i])){
        message(paste0(relevant_data[i], "is a .txt") )
        assign(paste0(file_path_sans_ext(relevant_data[i]), "_df"), read.delim(relevant_data[i]), envir = .GlobalEnv)
      }}}
dataframes_list <- ls(pattern = "_df", envir = .GlobalEnv)
dataframes_list <- assign("dataframes_list", dataframes_list, envir = .GlobalEnv)
message("All the relevant data has been loaded. The data frames for this project are: ", paste(dataframes_list, collapse = ", "), ".")
  }

#Running the function
create_dataframes(relevant_data) 
```

##Preparing for Data Manipulation
# To ensure that the recently converted data frames are ready for data manipulation, two 'for' loops check for the presence of missing values and duplicate rows. If any are found, the data frames are cleaned and an appropriate message is sent to the console.

#Ensuring the data frames contain no empty cells
for (dataframe in dataframes_list){
  if (sum(is.na(dataframe)) > 0){
    message(sum(is.na(dataframe)), " missing value(s) have been found and the relevant columns/rows have been removed. The data frame is now ready for manipulation. Proceed if you are satisfied.")
    #Remove columns containing missing values
    dataframe <- dataframe[ , colSums(is.na(dataframe))==0]
  } else if (sum(is.na(dataframe)) == 0){
    message(dataframe, " is ready for manipulation.")
  }
}

#Ensuring the dataframes don't contain duplicates
for (name_df in dataframes_list){
  df = get(name_df)
  duplicates = sum(duplicated(df))
  if (duplicates > 0){
    message(sum(duplicated(df)), " duplicate(s) have been found in ", name_df, " and the relevant rows have been removed.", "\n", "The data frame is now ready for manipulation. Proceed if you are satisfied.")
    df <- unique(df)
    assign(name_df, df, envir = .GlobalEnv)
  } else if (duplicates == 0){
    message(name_df, " is ready for manipulation.")
  }
}

#Creating a new data frame with all the information
#The following code was adapted from the following websites: 
#https://dplyr.tidyverse.org/reference/across.html
#https://www.rdocumentation.org/packages/base/versions/3.6.2/topics/numeric

#Binding together the data_all and gene_annotations data frames
total <- merge(gene_annotation_df, data_all_df, by="X") 
#Sorting the assigned genes in numerical order for ease of reading and manipulation
gene_list_sorted <- genelist_9_df[order(genelist_9_df$x), ]
#Matching the gene names in gene_list to the names in the 'total' data frame and then creating a new data frame with all the relevant information
matched_genes <- match(gene_list_sorted, total$X)
sample <- total[matched_genes, ]
#Remove the X column from the 'sample' data frame since it contains the same values as the Gene column
sample = subset(sample, select = -c(X))
#Set the values in the Gene column as characters (relevant for next step)
sample$Gene = as.character(sample$Gene)
#log-transform the values in the treatment columns
sample_logged <- mutate(sample, across(where(is.numeric), log10))

#Changing the data formatting
#The following code was adapted from the following website:
#https://tidyr.tidyverse.org/reference/pivot_longer.html

#create a vector containing the column names for the treatments
treatment <- c("A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L")
#converting the data type
sample_logged_long <- pivot_longer(sample_logged, cols = all_of(treatment), names_to = "treatment", values_to = "expression")
#Adding the TreatmentGroup column to the new data frame. We're repeating the column 40 times because we have 480 rows in the sample_logged_long data frame.
treatment_group <- rep(sample_annotation_df$TreatmentGroup, times = 40)
sample_logged_long$TreatmentGroup <- as.character(treatment_group)

#Plotting the data
#The code was adapted from the following websites
#https://r-charts.com/correlation/scatter-plot-group-ggplot2/
sample_logged_long$comb <- paste0(sample_logged_long$Type, sample_logged_long$LongName)

ggplot(sample_logged_long, aes(x = comb, y = expression, colour = Type, shape = TreatmentGroup)) +
  geom_point() +
  scale_color_manual(values = c("XA" = "#39FF39", 
                                "XB" = "#66CCFF", 
                                "XC" = "#FF79C6")) +
  theme(legend.position = "right") +
          theme_bw() +
          theme(panel.grid.major = element_blank()) +
         theme(panel.grid.minor = element_blank()) +
          theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
          labs(
            title = "log(expression) (vertical) for each gene (horizontal)."
            ) +
          xlab("Treatment Group") +
          ylab("log(Expression)") +
          theme(plot.title = element_text(size=12))

#Create two heatmaps:
#a) one where both the genes and the samples are clustered 

#Creating a matrix file only containing the log(expression) values and the LongName column as row names
no_Gene <- subset(sample_logged, select = -c(Gene, Type))
rownames(no_Gene) <- no_Gene$LongName
no_Gene <- subset(no_Gene, select = -c(LongName))
test <- as.matrix(no_Gene)

#Creating metadata data frames to annotate the heatmaps with Treatment Type and Gene Type.
metadata_c <- as.data.frame(c(sample_logged_long$TreatmentGroup), row.names = sample_logged_long$treatment)
names(metadata_c)[1] <- 'Treatment'
metadata_r <- as.data.frame(c(sample_logged_long$Type), row.names = sample_logged_long$LongName)
names(metadata_r)[1] <- 'Gene Type'

#Colour-coding the Treatment and Gene Type values
annotationscolours <- list( "Treatment" = c("1" = "#B572C7",
                                            "2" = "#9B4D96", 
                                            "3" = "#6A1E9C",
                                            "4" = "#5A137C"),
                            "Gene Type" = c("XA" = "#39FF39", 
                                            "XB" = "#66CCFF", 
                                            "XC" = "#FF79C6"))
#Creating the heatmap
pheatmap(test, 
         cluster_rows = TRUE, 
         cluster_cols = TRUE,
         main = "Heatmap Where Both Genes and Treatment are Clustered",
         annotation_legend = TRUE,
         annotation_col = metadata_c,
         annotation_row = metadata_r,
         fontsize = 7,
         annotation_names_row = FALSE,
         annotation_names_col = FALSE,
         annotation_colors = annotationscolours)

## Clustering Data Only Based on Genes
#b) one where only the genes are clustered.
pheatmap(test,
         cluster_rows = TRUE,
         cluster_cols = FALSE,
         main = "Heatmap Where Only Genes are Clustered",
         annotation_legend = TRUE,
         annotation_col = metadata_c,
         annotation_row = metadata_r,
         fontsize = 7,
         annotation_names_row = FALSE,
         annotation_names_col = FALSE,
         annotation_colors = annotationscolours)


## Representing Distributions via Boxplot
#Create a boxplot with an overlaying jitter-plot to represent the data
colours <- c("1" = "#B572C7", "2" = "#9B4D96", "3" = "#6A1E9C", "4" = "#5A137C")

boxplot(sample_logged_long$expression ~ sample_logged_long$TreatmentGroup, data = sample_logged_long, col = colours, main = "Figure 4: log(expression) of genes in response to Treatments 1, 2, 3, and 4. Jitterplots of gene types XA (green), XB (blue), and XC (pink) are overlayed on top of the boxplots", par(cex.main = 0.85), xlab = "Treatment Group", ylab ="log(Expression)")

substrip <- sample_logged_long[sample_logged_long$Type == "XA",]
stripchart(substrip$expression ~ substrip$TreatmentGroup, 
           data = substrip,
           method = "jitter",
           pch = 20,
           col = "#39FF39",
           vertical = TRUE,
           add = TRUE,
           at = c(0.7, 1.7, 2.7, 3.7))

substrip <- sample_logged_long[sample_logged_long$Type == "XB",]
stripchart(substrip$expression ~ substrip$TreatmentGroup, 
           data = substrip,
           method = "jitter",
           pch = 20,
           col = "#66CCFF",
           vertical = TRUE,
           add = TRUE,
           at = c(1, 2, 3, 4))

substrip <- sample_logged_long[sample_logged_long$Type == "XC",]
stripchart(substrip$expression ~ substrip$TreatmentGroup, 
           data = substrip,
           method = "jitter",
           pch = 20,
           col = "#FF79C6",
           vertical = TRUE,
           add = TRUE,
           at = c(1.3, 2.3, 3.3, 4.3))
