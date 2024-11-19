## ----global_options, include=TRUE, echo=FALSE---------------------------------
knitr::opts_chunk$set(
    echo = TRUE,
    warning = TRUE,
    message = TRUE,
    error = FALSE)

## ----'install1', eval=FALSE, message=FALSE, warning=FALSE---------------------
# # Installing the BiocManager package if needed
# if (!require("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")
# # Installing the fgsea package with the BiocManager
# BiocManager::install("fgsea")

## ----'install2', eval=FALSE, message=FALSE, warning=FALSE---------------------
# install.packages("mulea")

## ----'install3', eval=FALSE, message=FALSE, warning=FALSE---------------------
# # Installing the devtools package if needed
# if (!require("devtools", quietly = TRUE))
#     install.packages("devtools")
# 
# # Installing the mulea package from GitHub
# devtools::install_github("https://github.com/ELTEbioinformatics/mulea")

## ----'calling1', eval=FALSE---------------------------------------------------
# library(mulea)
# library(tidyverse)

## ----'calling2', echo=FALSE---------------------------------------------------
suppressMessages(library(mulea))
suppressMessages(library(tidyverse))

## ----'read_GMT1', eval=FALSE--------------------------------------------------
# # Reading the mulea GMT file locally
# tf_ontology <- read_gmt("Transcription_factor_RegulonDB_Escherichia_coli_GeneSymbol.gmt")

## ----'read_GMT2'--------------------------------------------------------------
# Reading the GMT file from the GitHub repository
tf_ontology <- read_gmt("https://raw.githubusercontent.com/ELTEbioinformatics/GMT_files_for_mulea/main/GMT_files/Escherichia_coli_83333/Transcription_factor_RegulonDB_Escherichia_coli_GeneSymbol.gmt")

## ----'read_GMT3', eval=FALSE--------------------------------------------------
# # Reading the Enrichr GMT file locally
# tf_enrichr_ontology <- read_gmt("TRRUST_Transcription_Factors_2019.txt")
# 
# # The ontology_name is empty, therefore we need to fill it with the ontology_id
# tf_enrichr_ontology$ontology_name <- tf_enrichr_ontology$ontology_id
# 

## ----'read_GMT4', eval=FALSE--------------------------------------------------
# # Reading the MsigDB GMT file locally
# tf_msigdb_ontology <- read_gmt("c3.tft.v2023.2.Hs.symbols.gmt")

## ----'read_muleaData', eval=FALSE---------------------------------------------
# # Installing the ExperimentHub package from Bioconductor
# BiocManager::install("ExperimentHub")
# 
# # Calling the ExperimentHub library.
# library(ExperimentHub)
# 
# # Downloading the metadata from ExperimentHub.
# eh <- ExperimentHub()
# 
# # Creating the muleaData variable.
# muleaData <- query(eh, "muleaData")
# 
# # Looking for the ExperimentalHub ID of the ontology.
# EHID <- mcols(muleaData) %>%
#   as.data.frame() %>%
#   dplyr::filter(title == "Transcription_factor_RegulonDB_Escherichia_coli_GeneSymbol.rds") %>%
#   rownames()
# 
# # Reading the ontology from the muleaData package.
# tf_ontology <- muleaData[[EHID]]
# 
# # Change the header
# tf_ontology <- tf_ontology %>%
#   rename(ontology_id = "ontologyId",
#          ontology_name = "ontologyName",
#          list_of_values = "listOfValues")

## ----'exclude_ontology'-------------------------------------------------------
# Filtering the ontology
tf_ontology_filtered <- filter_ontology(gmt = tf_ontology,
                                        min_nr_of_elements = 3,
                                        max_nr_of_elements = 400)

## ----'save_gmt', eval=FALSE---------------------------------------------------
# # Saving the ontology to GMT file
# write_gmt(gmt = tf_ontology_filtered,
#           file = "Filtered.gmt")

## ----'list_to_gmt_example', eval=FALSE----------------------------------------
# # Creating a list of gene sets
# ontology_list <- list(gene_set1 = c("gene1", "gene2", "gene3"),
#                       gene_set2 = c("gene4", "gene5", "gene6"))
# 
# # Converting the list to a ontology (GMT) object
# new_ontology_df <- list_to_gmt(ontology_list)

## ----'reading_target_bg'------------------------------------------------------
# Taget set
target_set <- readLines("https://raw.githubusercontent.com/ELTEbioinformatics/mulea/master/inst/extdata/target_set.txt")

# Background set
background_set  <- readLines("https://raw.githubusercontent.com/ELTEbioinformatics/mulea/master/inst/extdata/background_set.txt")

## ----'ora'--------------------------------------------------------------------
# Creating the ORA model using the GMT variable
ora_model <- ora(gmt = tf_ontology_filtered, 
                 # Test set variable
                 element_names = target_set, 
                 # Background set variable
                 background_element_names = background_set, 
                 # p-value adjustment method
                 p_value_adjustment_method = "eFDR", 
                 # Number of permutations
                 number_of_permutations = 10000,
                 # Number of processor threads to use
                 nthreads = 2, 
                 # Setting a random seed for reproducibility
                 random_seed = 1) 

# Running the ORA
ora_results <- run_test(ora_model)

## ----'ora_size'---------------------------------------------------------------
ora_results %>%
  # Rows where the eFDR < 0.05
  filter(eFDR < 0.05) %>% 
  # Number of such rows
  nrow()

## ----'print_ora', eval=FALSE--------------------------------------------------
# ora_results %>%
#   # Arrange the rows by the eFDR values
#   arrange(eFDR) %>%
#   # Rows where the eFDR < 0.05
#   filter(eFDR < 0.05)

## ----'print_ora2', echo=FALSE-------------------------------------------------
ora_results %>%
  # Arrange the rows by the eFDR values
  arrange(eFDR) %>% 
  # Rows where the eFDR < 0.05
  filter(eFDR < 0.05) %>% 
  knitr::kable()

## ----'init_plot_ora'----------------------------------------------------------
# Reshapeing the ORA results for visualisation
ora_reshaped_results <- reshape_results(model = ora_model, 
                                        model_results = ora_results, 
                                        # Choosing which column to use for the
                                        #     indication of significance
                                        p_value_type_colname = "eFDR")

## ----'lollipop_plot_ora'------------------------------------------------------
plot_lollipop(reshaped_results = ora_reshaped_results,
              # Column containing the names we wish to plot
              ontology_id_colname = "ontology_id",
              # Upper threshold for the value indicating the significance
              p_value_max_threshold = 0.05,
              # Column that indicates the significance values
              p_value_type_colname = "eFDR")

## ----'bar_plot_ora'-----------------------------------------------------------
plot_barplot(reshaped_results = ora_reshaped_results,
              # Column containing the names we wish to plot
              ontology_id_colname = "ontology_id",
              # Upper threshold for the value indicating the significance
              p_value_max_threshold = 0.05,
              # Column that indicates the significance values
              p_value_type_colname = "eFDR")

## ----'network_plot_ora'-------------------------------------------------------
plot_graph(reshaped_results = ora_reshaped_results,
           # Column containing the names we wish to plot
           ontology_id_colname = "ontology_id",
           # Upper threshold for the value indicating the significance
           p_value_max_threshold = 0.05,
           # Column that indicates the significance values
           p_value_type_colname = "eFDR")

## ----'heatmap_ora'------------------------------------------------------------
plot_heatmap(reshaped_results = ora_reshaped_results,
             # Column containing the names we wish to plot
             ontology_id_colname = "ontology_id",
             # Column that indicates the significance values
             p_value_type_colname = "eFDR")

## ----'ora_bh_bonferroni'------------------------------------------------------
# Creating the ORA model using the Benjamini-Hochberg p-value correction method
BH_ora_model <- ora(gmt = tf_ontology_filtered, 
                 # Test set variable
                 element_names = target_set, 
                 # Background set variable
                 background_element_names = background_set, 
                 # p-value adjustment method
                 p_value_adjustment_method = "BH",
                 # Number of processor threads to use
                 nthreads = 2) 

# Running the ORA
BH_results <- run_test(BH_ora_model)

# Creating the ORA model using the Bonferroni p-value correction method
Bonferroni_ora_model <- ora(gmt = tf_ontology_filtered, 
                            # Test set variable
                            element_names = target_set, 
                            # Background set variable
                            background_element_names = background_set, 
                            # p-value adjustment method
                            p_value_adjustment_method = "bonferroni",
                            # Number of processor threads to use
                            nthreads = 2) 

# Running the ORA
Bonferroni_results <- run_test(Bonferroni_ora_model)

## ----'compare_p.adj'----------------------------------------------------------
# Merging the Benjamini-Hochberg and eFDR results
merged_results <- BH_results %>% 
  # Renaming the column
  rename(BH_adjusted_p_value = adjusted_p_value) %>% 
  # Selecting the necessary columns
  select(ontology_id, BH_adjusted_p_value) %>%
  # Joining with the eFDR results
  left_join(ora_results, ., by = "ontology_id") %>% 
  # Converting the data.frame to a tibble
  tibble()

# Merging the Bonferroni results with the merged results
merged_results <- Bonferroni_results %>% 
  # Renaming the column
  rename(Bonferroni_adjusted_p_value = adjusted_p_value) %>% 
  # Selecting the necessary columns
  select(ontology_id, Bonferroni_adjusted_p_value) %>%
  # Joining with the eFDR results
  left_join(merged_results, ., by = "ontology_id") %>% 
  # Arranging by the p-value
  arrange(p_value)

# filter the p-value < 0.05 results
merged_results_filtered <- merged_results %>% 
  filter(p_value < 0.05) %>% 
  # remove the unnecessary columns
  select(-ontology_id, -nr_common_with_tested_elements, 
         -nr_common_with_background_elements)

## ----'print_compare_p.adj', echo=FALSE----------------------------------------
merged_results_filtered %>% 
  knitr::kable()

## ----'reading_ordered'--------------------------------------------------------
# Reading the tsv containing the ordered set
ordered_set <- read_tsv("https://raw.githubusercontent.com/ELTEbioinformatics/mulea/master/inst/extdata/ordered_set.tsv")

## ----'gsea', warning=FALSE, message=FALSE-------------------------------------
# Creating the GSEA model using the GMT variable
gsea_model <- gsea(gmt = tf_ontology_filtered,
                   # Names of elements to test
                   element_names = ordered_set$Gene.symbol,
                   # LogFC-s of elements to test
                   element_scores = ordered_set$logFC,
                   # Consider elements having positive logFC values only
                   element_score_type = "pos",
                   # Number of permutations
                   number_of_permutations = 10000)

# Running the GSEA
gsea_results <- run_test(gsea_model)

## ----'gsea_size'--------------------------------------------------------------
gsea_results %>%
  # rows where the adjusted_p_value < 0.05
  filter(adjusted_p_value < 0.05) %>% 
  # the number of such rows
  nrow()

## ----'print_gsea', eval=FALSE-------------------------------------------------
# gsea_results %>%
#   # arrange the rows by the adjusted_p_value values
#   arrange(adjusted_p_value) %>%
#   # rows where the adjusted_p_value < 0.05
#   filter(adjusted_p_value < 0.05)

## ----'print_gsea2', echo=FALSE------------------------------------------------
gsea_results %>%
  # arrange the rows by the adjusted_p_value values
  arrange(adjusted_p_value) %>% 
  # rows where the adjusted_p_value < 0.05
  filter(adjusted_p_value < 0.05) %>% 
  knitr::kable()

## ----'init_plot_gsea'---------------------------------------------------------
# Reshaping the GSEA results for visualisation
gsea_reshaped_results <- reshape_results(model = gsea_model, 
                                         model_results = gsea_results, 
                                         # choosing which column to use for the
                                         # indication of significance
                                         p_value_type_colname = "adjusted_p_value")

## ----'network_plot_gsea'------------------------------------------------------
plot_graph(reshaped_results = gsea_reshaped_results,
           # the column containing the names we wish to plot
           ontology_id_colname = "ontology_id",
           # upper threshold for the value indicating the significance
           p_value_max_threshold = 0.05,
           # column that indicates the significance values
           p_value_type_colname = "adjusted_p_value")

## ----'DE1', eval=TRUE, message=FALSE, warning=FALSE---------------------------
# Importing necessary libraries and reading the DE results table
geo2r_result_tab <- read_tsv("https://raw.githubusercontent.com/ELTEbioinformatics/mulea/master/inst/extdata/GSE55662.table_wt_non_vs_cipro.tsv")

## ----'print_geo1', eval=FALSE-------------------------------------------------
# # Printing the first few rows of the data frame
# geo2r_result_tab %>%
#   head(3)

## ----'print_geo2', echo=FALSE-------------------------------------------------
# Printing the first few rows of the data frame
geo2r_result_tab %>%  
  head(3) %>% 
  knitr::kable()

## ----'format_geo'-------------------------------------------------------------
# Formatting the data frame
geo2r_result_tab <- geo2r_result_tab %>% 
  # Extracting the primary gene symbol and removing extraneous information
  mutate(Gene.symbol = str_remove(string = Gene.symbol,
                                  pattern = "\\/.*")) %>% 
  # Filtering out rows with NA gene symbols
  filter(!is.na(Gene.symbol)) %>% 
  # Sorting by logFC
  arrange(desc(logFC))

## ----'print_geo_formatted1', eval=FALSE---------------------------------------
# # Printing the first few rows of the formatted data frame
# geo2r_result_tab %>%
#   head(3)

## ----'print_geo_formatted2', echo=FALSE---------------------------------------
# Printing the first few rows of the formatted data frame
geo2r_result_tab %>%  
  head(3) %>% 
  knitr::kable()

## ----'target_set'-------------------------------------------------------------
target_set <- geo2r_result_tab %>% 
  # Filtering for adjusted p-value < 0.05 and logFC > 1
  filter(adj.P.Val < 0.05 & logFC > 1) %>% 
  # Selecting the Gene.symbol column
  select(Gene.symbol) %>% 
  # Converting the tibble to a vector
  pull() %>% 
  # Removing duplicates
  unique()

## ----'target_head'------------------------------------------------------------
target_set %>% 
  head(10)

## ----'target_gene_nr'---------------------------------------------------------
target_set %>% 
  length()

## ----'background_set'---------------------------------------------------------
background_set <- geo2r_result_tab %>% 
  # Selecting the Gene.symbol column
  select(Gene.symbol) %>% 
  # Converting the tibble to a vector
  pull() %>% 
  # Removing duplicates
  unique()

## ----'background_gene_nr'-----------------------------------------------------
background_set %>% 
  length()

## ----'save_target_bg', eval=FALSE---------------------------------------------
# # Save taget set to text file
# target_set %>%
#   writeLines("target_set.txt")
# 
# # Save background set to text file
# background_set %>%
#   writeLines("inst/extdata/background_set.txt")

## ----'gsea_input'-------------------------------------------------------------
# If there are duplicated Gene.symbols keep the first one only
ordered_set <- geo2r_result_tab %>% 
  # Grouping by Gene.symbol to be able to filter
  group_by(Gene.symbol) %>%
  # Keeping the first row for each Gene.symbol from rows with the same 
  #     Gene.symbol
  filter(row_number()==1) %>% 
  # Ungrouping
  ungroup() %>% 
  # Arranging by logFC in descending order
  arrange(desc(logFC)) %>%
  select(Gene.symbol, logFC)

## ----'ordered_genes_length'---------------------------------------------------
ordered_set %>% 
  nrow()

## ----'save_ordered', eval=FALSE-----------------------------------------------
# # Save ordered set to text file
# ordered_set %>%
#   write_tsv("ordered_set.tsv")

## ----'session_info'-----------------------------------------------------------
sessionInfo()

