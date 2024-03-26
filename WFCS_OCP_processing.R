# --------------------------------------------
# OCP Processing
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
# Description: This code processes OCP data and makes 
# a new csv for processed OCP data.
# 
# Funding: This work was supported under the project 
# name 'Contaminant Loads in Waterfowl of the Northeast Atlantic
# Flyway: New Threats and Outdated Advisories' of the 
# Multistate Conservation Grant Program.
#
# --------------------------------------------

#Load required packages.
library(tidyverse)
library(knitr)

# Read in study metadata
metadata <- read_csv("Waterfowl Contaminant Study Sample Collection Metadata 2021-2022.csv")

# Read in E1 OCP data from SGS AXYS
E1OCPs <- read_csv("Waterfowl Muscle Tissue E1 Organochlorine Pesticide Data from SGS AXYS 2021-2022.csv")

# Split the E1 OCP data based on SDG (sample batch identifier)
E1OCPs_list <- split(E1OCPs, E1OCPs$SDG, drop = TRUE)

# Generate names for the list elements
names(E1OCPs_list) <- paste0("E1OCPs", seq_along(E1OCPs_list))

E1OCPs1 <- E1OCPs_list[[1]]
E1OCPs2 <- E1OCPs_list[[2]]
E1OCPs3 <- E1OCPs_list[[3]]
E1OCPs4 <- E1OCPs_list[[4]]
E1OCPs5 <- E1OCPs_list[[5]]

# Create a vector of data frame names
E1OCPs_dataframes <- c("E1OCPs1","E1OCPs2", "E1OCPs3", "E1OCPs4", "E1OCPs5" )

# Loop through each E1OCP dataframe.
for (i in seq_along(E1OCPs_dataframes)) {
    # Access data frame by name
    raw_data <- get(E1OCPs_dataframes[i])
  
  # Subset columns, remove % Lipid and Moisture results, rename Sample_ID to ID.
  raw_data <- raw_data %>% 
    select(Sample_Type, Analyte, Result, Result_Qualifier, Detect, LOQ, Reporting_Limit, Sample_ID, Result_Unit ) %>% 
    filter(!Analyte %in% c("% Lipid", "% Moisture")) %>% 
    rename(ID = Sample_ID)
  
  # BLANK CORRECTION--------------------------------------------------------------
  # Subset to get only blanks data with detectable signals.
  blanks_data <- raw_data %>% 
    filter(Sample_Type == "BLANK") %>% 
    filter(Detect == "Y") %>% 
    select(Analyte, Result) %>% 
    rename(blank_result = Result)
  
  # Remove non-samples from raw data.
  raw_data2 <- raw_data %>% filter(Sample_Type == "Sample")
  
  # Merge blank and sample data.
  sample_blank_merged <- merge(raw_data2, blanks_data, by = "Analyte", all.x = TRUE)
  
  # Subtract detectable blank values from data, correct detectable/undetectable.
  blank_corrected <- sample_blank_merged %>% 
    mutate(Result = if_else(is.na(blank_result), Result, Result - blank_result)) %>% 
    mutate(Detect = if_else(Result < Reporting_Limit | Detect == "N","N", "Y")) %>% 
    mutate(Result = if_else(Result < Reporting_Limit | is.na(Result), 0, Result)) 
  
  # HAZARD QUOTIENT AND CANCER RISK CALCULATIONS---------------------------------
  # Read in reference doses.
  rfds <- read_csv("reference_doses.csv")
  rfds2 <- rfds %>% select(Analyte, rfd)
  
  # Read in cancer slope factors.
  CSFs <- read_csv("cancer_slope_factors.csv")
  CSFs2 <- CSFs %>% select(Analyte, CSF)
  
  # Merge blank corrected data with reference doses and CSFs.
  E1_OCPs_merged <- merge(merge(blank_corrected, rfds2, by = "Analyte", all.x = TRUE),CSFs2, by = "Analyte", all.x = TRUE)
  
  # Subset columns.
  E1_OCPs_merged <- E1_OCPs_merged %>% select(ID, Analyte, Result, rfd, CSF)
  
  # Sort columns.
  E1_OCPs_merged <- E1_OCPs_merged %>% arrange(ID, Analyte)
  
  # Create a new dataframe to store the results.
  E1_OCP_results <- data.frame(ID = unique(E1_OCPs_merged$ID))
  
  # Loop through each analyte and perform the required operations.
  analytes <- unique(E1_OCPs_merged$Analyte)  
  for (analyte in analytes) {
    # Select rows corresponding to the current analyte.
    analyte_data <- E1_OCPs_merged %>% filter(Analyte == analyte)
    # Extract the result and reference dose columns.
    result_col <- analyte_data$Result
    ref_dose <- unique(analyte_data$rfd)
    slope_factor <- unique(analyte_data$CSF)
    # Create a new column in the results data frame for the current analyte.
    E1_OCP_results[[analyte]] <- result_col
    # Create a new column for HQ_Analyte for 2 meals/month.
    #(concentration ng/g * 2 meals/month * 227 g/meal)/(rfd ng/kgBWday * 80 kgBW * 30.4 days/month)
    E1_OCP_results[[paste0("HQ_", analyte)]] <- signif((result_col*2*227)/(ref_dose*80*30.4), digits = 3)
    # Create a new column for CR_Analyte for 2 meals/month.
    # (concentration ng/g * 2 meals/month * 227 g/meal * CSF kgBWday/ng)/(80 kgBW * 30.4 days/month)
    E1_OCP_results[[paste0("CR_", analyte)]] <- signif((result_col*2*227*slope_factor)/(80*30.4), digits = 3)
    # Create a new column for EPA_CR_Analyte for 24 meals per yr over 26 yrs with 70 yr averaging time:
    # (concentration ng/g * 227 g/meal * 24 meals/year * 26 years * CSF kgBWday/ng)/(80 kgBW * 25550 days)
    E1_OCP_results[[paste0("EPA_CR_", analyte)]] <- signif((result_col*24*227*26*slope_factor)/(80*25550), digits = 3)
  } # End analyte in analyte
  
  # Create a new dataframe with the desired name in the global environment.
  assign(paste0('E1_OCPs_full', i), E1_OCP_results)
  
} # End loop through each dataframe.


# Combine processed wide E1_OCP data frames into final dataframe.
E1_OCPs_processed <- bind_rows(E1_OCPs_full1, E1_OCPs_full2,E1_OCPs_full3, E1_OCPs_full4, E1_OCPs_full5, .id = "original dataset")
E1_OCPs_processed <- E1_OCPs_processed %>% select(-`original dataset`)

# Check IDs in processed data for mismatches with metadata
# Check that SGS IDs match metadata IDs
metadata_PCBtest <- metadata %>% 
  filter(PCBOCP_test == "Y")
ID_mismatch <- E1_OCPs_processed %>% 
  filter(!ID %in% metadata_PCBtest$ID) # no mismatched IDs


# Read in E2 OCP data from SGS AXYSdata:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAABwAAAAcCAYAAAByDd+UAAAC10lEQVR42u2W20vTYRjHva/+jQkSEzrQRdB9VwXdVEvYxAJnUSDzwlDotBZqTrYfLdq0aefVDMNp22JbOxVGB6OMrY3pclP7CSUTB7/s7f1OFmv+jtpV9MLD7/x8fu/zfp/neWtq/g86ahsd22p11gMqndVI7b5Kx/ioeVQ6ixP38AzvbBpUp2e2UodNO1rsYb0tkLjqmWIHItmVWy/nV/ujuR+93vRyh+s9q+n1J+qb7WG8u2Ew/nr7CVuobTCaevR6gXO/Zcndia9kID5ProfzhAnliBUWzJG+QI5cHs9wmr5Auq7JFsK3CkPIaPYZhuK3Y9nCk8lFcmdigdgopOfZLDE9zZKLY1lyYWymdDTSa5M3S654v5Su21yfCrtPO+PwIQumarTs39vqjD18lS8+frdI7NE50u1fc3beMyNo+AEYztuHk8Vdp/pj8CUKUzd0bUEYneHMEmAIHWZUdgQLJr6JgsvWeu/DEsILHQjPTsscN9yMphBGR2yuFKpqRxhyoUfMwRR8Ckp/Z4s94qYCwZr1+Gf/mFklUC60cyTDqZtvRHiVq9JaDuptwQTUeO15XnDNKocc6KEubwK+ecRivdTt+ciWZkfVKOSgekhBTzrfsCgO64G0gtCkLiLPIBS5QClouztZhG8eIONDBYEyxVJAaAhBz41Or8I33wxHHbRcoYLwiUUKGEp+Vwpkhsy+9DLKldIZCsHEQ0pF0+GaZFEbjQrWUAwmLhoqXY3Zn0Ah5kt4PqAUTDQtkJz1NElN49McCrHQOiqBiSZ+ubQd7Qumyt1ACCgHJlnayg0XBdfwYKog5EQuDG1KsniXG++eM4Mvzg5/LspxzKtM+i18yG7EKq21AU0Uf6kU9rsB6yzHFG8xEJLD5kC6cyTNSQskzeHdDW0xKpWLjZGabpAgb+QUTeQVWj1+wnCOe3im3uwmah2Y5lLVNtG3dk7v0Wd/BfRPjF/sOXqT33GGYwAAAABJRU5ErkJggg==
E2OCPs <- read_csv("Waterfowl Muscle Tissue E2 Organochlorine Pesticide Data from SGS AXYS 2021-2022.csv")

# Split the E2 OCP data based on SDG (sample batch identifier)
E2OCPs_list <- split(E2OCPs, E2OCPs$SDG, drop = TRUE)

# Generate names for the list elements
names(E2OCPs_list) <- paste0("E2OCPs", seq_along(E2OCPs_list))

E2OCPs1 <- E2OCPs_list[[1]]
E2OCPs2 <- E2OCPs_list[[2]]
E2OCPs3 <- E2OCPs_list[[3]]
E2OCPs4 <- E2OCPs_list[[4]]
E2OCPs5 <- E2OCPs_list[[5]]

# Create a vector of data frame names
E2OCPs_dataframes <- c("E2OCPs1","E2OCPs2", "E2OCPs3", "E2OCPs4", "E2OCPs5")

# Loop through each E2OCP dataframe.
for (i in seq_along(E2OCPs_dataframes)) {
  # Access data frame by name
  raw_data <- get(E2OCPs_dataframes[i])
  
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
  
  # Merge blank and sample data
  sample_blank_merged <- merge(raw_data2, blanks_data,by = "Analyte", all.x = TRUE)
  
  # Subtract detectable blank values from data, correct detectable/undetectable
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
  E2_OCPs_merged <- merge(merge(blank_corrected, rfds2, by = "Analyte", all.x = TRUE),CSFs2, by = "Analyte", all.x = TRUE)
  
  # Subset columns.
  E2_OCPs_merged <- E2_OCPs_merged %>% select(ID, Analyte, Result, rfd, CSF)
  
  # Sort columns.
  E2_OCPs_merged <- E2_OCPs_merged %>% arrange(ID, Analyte)
  
  # Create a new dataframe to store the results.
  E2_OCP_results <- data.frame(ID = unique(E2_OCPs_merged$ID))
  
  # Loop through each analyte and perform the required operations.
  analytes <- unique(E2_OCPs_merged$Analyte)
  for (analyte in analytes) {
    # Select rows corresponding to the current analyte.
    analyte_data <- E2_OCPs_merged %>% filter(Analyte == analyte)
    # Extract the result and reference dose columns.
    result_col <- analyte_data$Result
    ref_dose <- unique(analyte_data$rfd)
    slope_factor <- unique(analyte_data$CSF)
    # Create a new column in the result data frame for the current analyte.
    E2_OCP_results[[analyte]] <- result_col
    # If reference dose available, create a new column for HQ_Analyte for 2 meals/month:
    #(concentration ng/g * 2 meals/month * 227 g/meal)/(rfd ng/kgBWday * 80 kgBW * 30.4 days/month)
    if (!is.na(ref_dose)) {E2_OCP_results[[paste0("HQ_", analyte)]] <- signif((result_col * 2 * 227) / (ref_dose * 80 * 30.4), digits = 3)}
    # If reference dose available, create a new column for CR_Analyte for 2 meals/month:
    # (concentration ng/g * 2 meals/month * 227 g/meal * CSF kgBWday/ng)/(80 kgBW * 30.4 days/month)
    if (!is.na(slope_factor)) {E2_OCP_results[[paste0("CR_", analyte)]] <- signif((result_col * 2 * 227 * slope_factor) / (80 * 30.4), digits = 3)}
    # If reference dose available, create a new column for EPA_CR_Analyte for 24 meals per yr over 26 yrs with 70 yr averaging time:
    # (concentration ng/g * 227 g/meal * 24 meals/year * 26 years * CSF kgBWday/ng)/(80 kgBW * 25550 days) 
    if (!is.na(slope_factor)) {E2_OCP_results[[paste0("CR_", analyte)]] <- signif((result_col*24*26*227 * slope_factor) / (80 * 25550), digits = 3)}
  } # End loop through each analyte and perform the required operations.
  
  # Create a new dataframe with the desired name in the global environment
  assign(paste0('E2_OCPs_full', i), E2_OCP_results)
  
} # End loop through each CSV file.

# Combine processed wide E2_OCP dataframes into final dataframe.
E2_OCPs_processed <- bind_rows(E2_OCPs_full1, E2_OCPs_full2, E2_OCPs_full3, E2_OCPs_full4, E2_OCPs_full5, .id = "original dataset")
E2_OCPs_processed <- E2_OCPs_processed %>% select(-`original dataset`)

# Check IDs in processed data for mismatches with metadata
# Check that SGS IDs match metadata IDs
metadata_PCBtest <- metadata %>% 
  filter(PCBOCP_test == "Y")
ID_mismatch2 <- E2_OCPs_processed %>% 
  filter(!ID %in% metadata_PCBtest$ID) # no mismatched IDs

# Merge the two data frames.
OCPs_processed <- merge(E1_OCPs_processed, E2_OCPs_processed,by = "ID", all.x = TRUE)

# Process the data frame. 
OCPs_processed2 <- OCPs_processed %>% 
  mutate(Sum_OCPs_HQ = rowSums(select(., starts_with("HQ_")), na.rm = TRUE)) %>% 
  mutate(Sum_OCPs_CR = rowSums(select(., starts_with("CR_")), na.rm = TRUE)) %>% 
  mutate(Sum_OCPs_EPACR = rowSums(select(., starts_with("EPA_CR_")), na.rm = TRUE))

# Write the csv file.
write_csv(OCPs_processed2, "WFCS_OCPs_processed.csv")

