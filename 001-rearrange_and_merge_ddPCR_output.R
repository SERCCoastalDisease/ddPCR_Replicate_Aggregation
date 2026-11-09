#--------------------------------------------------------------------------#
#-------------Script to combine Digital Droplet PCR replicates-------------#
#--------------------------------------------------------------------------#
#
# Authors: 
# Emma Palmer - SERC
# Formulas taken from the Droplet Digital PCR Applications Guide (ch 3)
#   -link: https://www.bio-rad.com/webroot/web/pdf/lsr/literature/Bulletin_6407.pdf
#
# Details:
# - Input: csv files of the ddPCR output results 
#
# - Script:
#     - Ultimate goal: calculate final statistics (concentration, copy number) from ddPCR results across multiple replicates
#     - This script can additionally add the metadata to the combined results, making the data ready for write-up or further analysis

# - Output: xlsx file of the ddPCR output results with metadata
# 
# Before Running Script: 
#   - Make sure to threshold the output BEFORE running this script; Failure to do so will result in inaccurate results
#   - Ensure that the ddPCR results are named using the following format
#       - If using another format, make sure to change line 52 (plate_names) to accurately reflect the naming convention
# 
# TODO:
#     - Ultimate goal: calculate final statistics (concentration, copy number) from ddPCR results across multiple replicates
#--------------------------------------------------------------------------#

# Libraries
library(tidyverse)  # for pipelines
library(openxlsx)   # for creating final result

source("000-custom_functions.R")

options(scipen = 999) # turns off scientific notation usage

# -------- Gathering the Data

# ----- Project Parameters

# TODO: Update the values in this section to match the project
project_name <- "ProjectName"                                 # Name of the project
folder_OG_data <- paste0("ddPCR_output/", project_name)       # folder path to ddPCR output files
folder_results <- "adjusted_results/"                         # folder in which to place final results
vSample <- 5                                                  # Volume of sample added to ddPCR reaction (µL)
vTotal <- 25                                                  # Total reaction volume (µL)

# ----- Metadata

# TODO: Update path to metadata file
# IMPORTANT: make sure the names across the metadata and ddPCR results match
# metadata <- read.xlsx("metadata/Metadata_ProjectName_Date.xlsx", sheet = 2)
metadata <- read.csv("metadata/Metadata_ProjectName_Date.csv")

# ----- Pulling all data files

# Pulls names of all csv files in project folder
all_OG_plates <- sort(list.files(folder_OG_data, pattern = "*.csv", full.names = T))

# Change if plate names are not in the last section of the file, creates unique names for the plates
plate_names <- str_remove(str_split_i(all_OG_plates, "_", -1), "\\.csv")  

# Creates a single dataframe of all csv files and adds the associated plate names
ddPCR_output <- lapply(all_OG_plates, read.csv, fileEncoding='latin1') %>%
  set_names(., str_remove(str_remove(all_OG_plates, paste0("\\.csv")), paste0(folder_OG_data, "/"))) %>%
  map2_dfr(., plate_names, ~.x %>% mutate(Plate = .y, 
                            across(everything(), as.character))) %>% 
  bind_rows(.)

# ----- Trim and summarize data

# Aim: Combine replicate ddPCR runs into a single result per sample
#     - filter(str_detect(Status.Reason, "<", negate = TRUE) | is.na(Status.Reason))): Removes failed runs (< 10000 droplets)
#     - clean_counts: Transforms runs with <3 positive droplets to negative runs
#     - summarize_counts: Calculates the summary metrics (positive droplets, concentration, copies per µL)
# TODO: change the last argument in join_by to the proper joining column if adding metadata  
#         Note: the only select certain columns or change the order of the final output 

adjusted_output <- ddPCR_output %>% 
  filter(str_detect(Status.Reason, "<", negate = TRUE) | is.na(Status.Reason)) %>%              
  select(Plate, Sample.description.1:Sample.description.3, Conc.copies.µL., Positives, Negatives) %>% 
  dplyr::rename(SampleID = Sample.description.1, Replicate = Sample.description.2, Marker = Sample.description.3, Conc_copies_uL = Conc.copies.µL.) %>% 
  clean_counts() %>% 
  summarize_counts(., vSample, vTotal) %>% 
  select(FullID = SampleID, Marker, Plate, PositiveDroplets:PresAbs, Conc_copies_uL, contains("Poisson"), Copies_uL) %>% 
  mutate(SampleID = ifelse(str_detect(FullID, "PosControl"), "PosControl", FullID)) %>% 
  left_join(metadata, by = join_by(FullID == JoinName)) %>%           # TODO: change "JoinName" to share column in metadata and ddPCR_output
  select(FullID, Metadata1:Metadata2, Plate:Copies_uL)                # TODO: change Metadata1:Metadata2 to column names that you wish to keep from metadata

# ----- Creating Final Workbook

# File containing formatted ReadME file added to final result
read_me_template <- read.xlsx("TEMPLATE_README_file.xlsx", colNames = FALSE, skipEmptyRows = FALSE)

# Replaces placeholders in template with project specific information
readme <- createREADme(read_me_template)

# Adjusts xlsx formatting
header_style <- createStyle(
  textDecoration = "bold",
  border = "bottom")

title_style <- createStyle(textDecoration = "bold")

# Creates workbook of the final results
ddPCR_adjusted_wb <- createWorkbook()
addWorksheet(ddPCR_adjusted_wb, "READme")
addWorksheet(ddPCR_adjusted_wb, "ddPCR_data")

writeData(ddPCR_adjusted_wb, "READme", readme, startCol = 1, startRow = 1, colNames = FALSE)
addStyle(ddPCR_adjusted_wb, "READme", title_style, rows = c(3, 6, 9), cols = 1)
addStyle(ddPCR_adjusted_wb, "READme", header_style, rows = 12, cols = 1:2)

writeData(ddPCR_adjusted_wb, "ddPCR_data", adjusted_output, startCol = 1, startRow = 1)
addStyle(ddPCR_adjusted_wb, "ddPCR_data", header_style, rows = 1, cols = 1:41)

# TODO: Uncomment line below to save final workbook
# saveWorkbook(ddPCR_adjusted_wb, paste0(folder_results, project_name, "-adjusted_ddPCR_output.xlsx"), overwrite = TRUE)
