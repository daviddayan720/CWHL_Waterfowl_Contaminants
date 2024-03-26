# --------------------------------------------
# WFCS PFAS processing
#
# Code Author: David Dayan
# Code QAQC: Brenda Hanley
# Project Supervisor: Krysten Schuler
#
# Project: Waterfowl Contaminants in the NE US. 
# 
# Location: Cornell Wildlife Health Lab
# Date: February 2024
# License: MIT
# 
# Description: This code processes PFAS data and makes 
# a new csv for processed PFAS data.
# 
# Funding: This work was supported under the project 
# name 'Contaminant Loads in Waterfowl of the Northeast Atlantic
# Flyway: New Threats and Outdated Advisories' of the 
# Multistate Conservation Grant Program.
#
# --------------------------------------------

# Load required pacakges 
library(tidyverse)

# Read in study metadata
metadata <- read_csv("Waterfowl Contaminant Study Sample Collection Metadata 2021-2022.csv")

# Read in PFAS data from SGS AXYS
PFAS <- read_csv("Waterfowl Muscle Tissue Per- and Polyfluoroalkyl Substance Data from SGS AXYS 2021-2022.csv")

# Split the PFAS data based on SDG (sample batch identifier)
PFAS_list <- split(PFAS, PFAS$SDG, drop = TRUE)

# Generate names for the list elements
names(PFAS_list) <- paste0("PFAS", seq_along(PFAS_list))

PFAS1 <- PFAS_list[[1]]
PFAS2 <- PFAS_list[[2]]
PFAS3 <- PFAS_list[[3]]
PFAS4 <- PFAS_list[[4]]
PFAS5 <- PFAS_list[[5]]
PFAS6 <- PFAS_list[[6]]

# Create a vector of data frame names
PFAS_dataframes <- c("PFAS1", "PFAS2", "PFAS3",
                     "PFAS4", "PFAS5", "PFAS6")

# Loop through each PFAS dataframe file.
for (i in seq_along(PFAS_dataframes)) {
  
  # Access data frame by name
  raw_data <- get(PFAS_dataframes[i])
  
  # Subset columns, remove % Lipid and Moisture results, rename Sample_ID to ID.
  raw_data <- raw_data %>% select(Sample_Type, Analyte, Result, Result_Qualifier, Detect, LOQ, Reporting_Limit, Sample_ID, Result_Unit ) %>% 
    filter(!Analyte %in% c("% Lipid", "% Moisture")) %>% rename(ID = Sample_ID)
  
  # Subset to get only blanks data with detectable signals.
  blanks_data <- raw_data %>% 
    filter(Sample_Type == "BLANK") %>% 
    filter(Detect == "Y") %>% 
    select(Analyte, Result) %>% 
    rename(blank_result = Result)
  
  # Remove non-samples from raw data.
  raw_data2 <- raw_data %>% filter(Sample_Type == "Sample")
  
  # Merge blank and sample data.
  sample_blank_merged <- merge(raw_data2, blanks_data,by = "Analyte", all.x = TRUE)
  
  # Subtract detectable blank values from data, correct detectable/undetectable.
  blank_corrected <- sample_blank_merged %>% 
    mutate(Result = if_else(is.na(blank_result), Result, Result - blank_result)) %>% 
    mutate(Detect = if_else(Result < Reporting_Limit | Detect == "N","N", "Y")) %>% 
    mutate(Result = if_else(Result < Reporting_Limit | is.na(Result), 0, Result)) 
  
  # Read in reference doses.
  rfds <- read_csv("reference_doses.csv")
  rfds2 <- rfds %>% select(Analyte, rfd)
  
  # Read in cancer slope factors.
  CSFs <- read_csv("cancer_slope_factors.csv")
  CSFs2 <- CSFs %>% select(Analyte, CSF)
  
  # Merge blank corrected data with reference doses.
  PFAS_merged <- merge(merge(blank_corrected, rfds2, by = "Analyte", all.x = TRUE),CSFs2, by = "Analyte", all.x = TRUE)
  
  # Subset columns.
  PFAS_merged <- PFAS_merged %>% select(ID, Analyte, Result, rfd, CSF) 
  
  # Sort columns.
  PFAS_merged <- PFAS_merged %>% arrange(ID, Analyte)
  
  # Create a new dataframe to store the results.
  PFAS_results <- data.frame(ID = unique(PFAS_merged$ID))
  
  # Loop through each analyte and perform the required operations.
  analytes <- unique(PFAS_merged$Analyte)
  for (analyte in analytes) {
    # Select rows corresponding to the current analyte
    analyte_data <- PFAS_merged %>% filter(Analyte == analyte)
    # Extract the result and reference dose columns
    result_col <- analyte_data$Result
    ref_dose <- unique(analyte_data$rfd)
    slope_factor <- unique(analyte_data$CSF)
    # Create a new column in the result_df for the current analyte
    PFAS_results[[analyte]] <- result_col
    # Create a new column for HQ_Analyte 2 for meals/month only if ref_dose is available
    #(concentration ng/g * 2 meals/month * 227 g/meal)/(rfd ng/kgBWday * 80 kgBW * 30.4 days/month)
    if (!is.na(ref_dose)) {
      PFAS_results[[paste0("HQ_", analyte)]] <- 
        signif((result_col * 2 * 227) / (ref_dose * 80 * 30.4), digits = 3)
    }
    # Create a new column for CR_Analyte for 2 meals/month only if slope_factor is available
    # (concentration ng/g * 2 meals/month * 227 g/meal * CSF kgBWday/ng)/(80 kgBW * 30.4 days/month)
    if (!is.na(slope_factor)) {
      PFAS_results[[paste0("CR_", analyte)]] <- 
        signif((result_col * 2 * 227 * slope_factor) / (80 * 30.4), digits = 3)
    }
    # Create a new column for EPA_CR_Analyte only if slope_factor is available
    # 24 meals per yr over 26 yrs with 70 yr averaging time
    # (concentration ng/g * 227 g/meal * 24 meals/year * 26 years * CSF kgBWday/ng)/(80 kgBW * 25550 days)
    if (!is.na(slope_factor)) {
      PFAS_results[[paste0("EPA_CR_", analyte)]] <- 
        signif((result_col*24*26*227 * slope_factor) / (80 * 25550), digits = 3)
    }
  } # End loop through each analyte and perform the required operations.
  
  # Create a new dataframe with the desired name in the global environment
  assign(paste0('PFAS_full', i), PFAS_results)
  
} # End loop through each CSV file.

# Combine processed wide PFAS dataframes into final dataframe.
PFAS_processed <- bind_rows(PFAS_full1,PFAS_full2,PFAS_full3,PFAS_full4,PFAS_full5,PFAS_full6,.id = "original dataset")

# Make combined PFAS HQ, CR, and EPA CR.
PFAS_processed2 <- PFAS_processed %>% 
  mutate(Sum_PFAS_HQ = rowSums(select(., starts_with("HQ_")), na.rm = TRUE)) %>% 
  mutate(Sum_PFAS_CR = rowSums(select(., starts_with("CR_")), na.rm = TRUE)) %>% 
  mutate(Sum_PFAS_EPACR = rowSums(select(., starts_with("EPA_CR_")), na.rm = TRUE)) %>% 
  select(-`original dataset`)

# Check IDs in processed data for mismatches with metadata
# Check that SGS IDs match metadata IDs
metadata_PFAStest <- metadata %>% 
  filter(PFAS_test == "Y")
ID_mismatch <- PFAS_processed2 %>% 
  filter(!ID %in% metadata_PFAStest$ID) 
  #NJ-AGWT-02-NJ instead of NJ_AGWT_02_NJ

#rename NJ-AGWT-02-NJ to NJ_AGWT_02_NJ
PFAS_processed3 <- PFAS_processed2 %>% 
  mutate(ID = if_else(ID == "NJ-AGWT-02-NJ", "NJ_AGWT_02_NJ", ID))

# Check IDs in processed data for mismatches with metadata
# Check that SGS IDs match metadata IDs
ID_mismatch <- PFAS_processed3 %>% 
  filter(!ID %in% metadata_PFAStest$ID) # no mismatches

# Write the data.
write_csv(PFAS_processed3, "WFCS_PFAS_processed.csv")
