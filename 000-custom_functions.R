#------------------------------------------------------------------------------#
#------Functions used in data cleaning or analysis of metabarcoding data-------#
#------------------------------------------------------------------------------#
#
# AUTHOR:
# Emma Palmer - Smithsonian Environmental Research Center
# email: palmerem@si.edu
#
# CITATION: there is not a DOI yet to cite this work, please do not distribute
#
# Details:
# 
# How to use: see 001-rearrange_and_merge_ddPCR_output.R file
#------------------------------------------------------------------------------#

clean_counts <- function(data){
  # Purpose: remove failed (< 10000 droplets) and identify false positive runs 
  # Args
  #   data: complete data frame of the runs associated with a ddPCR project 
  # 
  # Output: data with two extra columns and modified Droplet columns:
  #   pseudoConc: recodes "No Call" to make the column numeric for future calculations
  #   totalDroplets: total number of droplets read from the run
  #   Negatives and Positives: if Positives is <= 3, identify the run as a false positive (all reads are negative)
  # set up columns we need
  data %>%
    mutate(across(Conc_copies_uL:Negatives, as.numeric),
           pseudoConc = if_else(Conc_copies_uL == "No Call", 0, as.numeric(Conc_copies_uL), missing = 0),
           totalDroplets = Positives + Negatives,
           Negatives = ifelse(Positives <= 3, totalDroplets, Negatives),
           Positives = ifelse(Positives <= 3, 0, Positives)) %>%
    filter(totalDroplets >= 10000)
}

summarize_counts <- function(data, vSample, vTotal){
  # Purpose: calculate summary metrics of the ddPCR runs 
  #   All formulas within the function are derived from the ddPCR Applications Guide
  # Args
  #   data: complete data frame of the runs associated with a ddPCR project 
  #   vSample: Volume of sample added to ddPCR reaction (µL)
  #   vTotal: Total reaction volume (µL)
  # 
  # Output: summarized dataframe including the following new columns
  #   pseudoConc: recodes "No Call" to make the column numeric for future calculations
  #   totalDroplets: total number of droplets read from the run
  #   NegativeDroplets and PositiveDroplets: if Positives is <= 3, identify the run as a false positive (all reads are negative)

  data %>%
    group_by(SampleID, Marker) %>%
    mutate(totalDroplets = Positives + Negatives) %>%
    summarize(flagged_Conc = paste(unique(Conc_copies_uL), collapse = "-"),
              n_Pos = sum(Positives > 0),
              Positives = round(mean(Positives), 1),
              Negatives = round(mean(Negatives), 1),
              AcceptedDroplets = Positives + Negatives,
              n_Rep = length(Replicate),
              Plate = paste(unique(Plate), collapse = "-"),) %>%
    relocate(n_Pos, .after = n_Rep) %>%
    mutate(PresAbs = case_when(n_Pos == 1 & ceiling(n_Rep/2) == 1 ~ "Low Concentration",
                               n_Pos >= ceiling(n_Rep/2) | n_Pos == n_Rep ~ "Present",
                               n_Pos < ceiling(n_Rep/2) & n_Pos != 0 ~ "Low Concentration",
                               flagged_Conc == "No Call" & n_Pos == 0 ~ "No Call",
                               n_Pos == 0 ~ "Absent",
                               T ~ "TO_CHECK"),
           Conc_copies_uL = -log(Negatives/AcceptedDroplets)/0.000795,
           p = Positives/AcceptedDroplets,
           sdNeg = sqrt(AcceptedDroplets * p * (1-p)),
           NnegMin = Negatives - (1.96 * sdNeg),
           NnegMax = Negatives + (1.96 * sdNeg),
           PoissonConfMin = -log(NnegMax/AcceptedDroplets)/0.000795,
           PoissonConfMin = ifelse(PoissonConfMin < 0, 0, PoissonConfMin),
           PoissonConfMax = -log(NnegMin/AcceptedDroplets)/0.000795,
           Copies_uL = (Conc_copies_uL * vTotal)/vSample) %>%
    ungroup() %>% 
    dplyr::rename(PositiveDroplets = Positives, NegativeDroplets = Negatives)
}

createREADme <- function(read_me) {
  # Purpose: replace ReadMe placeholders with project specific information
  # Args
  #   readme: template ReadME file associated with the Project, can be modified for projects/users
  # 
  # Output: readme with the system data and project name 
  read_me %>%
    mutate(X1 = str_replace(X1, "XXXX-XX-XX", as.character(Sys.Date())),
           X2 = str_replace(X2, "ProjectNameHere", paste0(project_name, "-adjusted_ddPCR_output")))
}