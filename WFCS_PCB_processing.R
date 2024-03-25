# --------------------------------------------
# WFCS PCB Processing
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
# Description: This code processes PCB data from SGS AXYS and makes 
# a new csv for processed PCB data.
# 
# Funding: This work was supported under the project 
# name 'Contaminant Loads in Waterfowl of the Northeast Atlantic
# Flyway: New Threats and Outdated Advisories' of the 
# Multistate Conservation Grant Program.
#
# --------------------------------------------

# Load required package.
library(tidyverse)

# Read in study metadata
metadata <- read_csv("Waterfowl Contaminant Study Sample Collection Metadata 2021-2022.csv")

# Read in PCB data from SGS AXYS
PCBs <- read_csv("Waterfowl Muscle Tissue Polychlorinated Biphenyl Data from SGS AXYS 2021-2022.csv")

# Split the PCB data based on SDG (sample batch identifier)
PCB_list <- split(PCBs, PCBs$SDG, drop = TRUE)

# Generate names for the list elements
names(PCB_list) <- paste0("PCB", seq_along(PCB_list))

# Now, PCB_list has names like "PCB1", "PCB2", ..., "PCB5"

PCB1 <- PCB_list[[1]]
PCB2 <- PCB_list[[2]]
PCB3 <- PCB_list[[3]]
PCB4 <- PCB_list[[4]]
PCB5 <- PCB_list[[5]]

# Create a vector of data frame names
PCB_dataframes <- c("PCB1", "PCB2", "PCB3", "PCB4", "PCB5")

# Loop through each PCB Dataframe.
for (i in seq_along(PCB_dataframes)) {
  # Access data frame by name
  raw_data <- get(PCB_dataframes[i])
  
  # Subset columns, remove % Lipid and Moisture results, rename Sample_ID to ID
  raw_data <- raw_data %>% select(Sample_Type, Analyte, Result, Result_Qualifier, Detect, LOQ, Reporting_Limit, Sample_ID, Result_Unit ) %>% filter(!Analyte %in% c("% Lipid", "% Moisture")) %>% rename(ID = Sample_ID)
  
  # Subset to get only blanks data with detectable signals
  blanks_data <- raw_data %>% filter(Sample_Type == "BLANK") %>% filter(Detect == "Y") %>% select(Analyte, Result) %>% rename(blank_result = Result)
  
  # Remove non-samples (including duplicates) from raw data.
  raw_data2 <- raw_data %>% filter(Sample_Type == "Sample")
  
  # Merge blank and sample data.
  sample_blank_merged <- merge(raw_data2, blanks_data,by = "Analyte", all.x = TRUE)
  
  # Subtract detectable blank values from data, correct detectable/undetectable.
  blank_corrected <- sample_blank_merged %>% 
    mutate(Result = if_else(is.na(blank_result), Result, Result - blank_result)) %>% 
    mutate(Detect = if_else(Result < Reporting_Limit | Detect == "N","N", "Y"))
  
  # Separate congeners.
  PCB_congeners <- blank_corrected %>% filter(Analyte != "TOTAL PCBs")
  
  # Separate total PCBs.
  TPCBs <- blank_corrected %>% filter(Analyte == "TOTAL PCBs")
  
  # Read in reference doses.
  rfds <- read_csv("reference_doses.csv")
  rfds2 <- rfds %>% select(Analyte, rfd)
  
  # Read in cancer slope factors.
  CSFs <- read_csv("cancer_slope_factors.csv")
  CSFs2 <- CSFs %>% select(Analyte, CSF)
  
  # Merge TOTAL PCBs with rfd and csf.
  TPCBs_merged <- merge(merge(TPCBs, rfds2, by = "Analyte", all.x = TRUE),CSFs2, by = "Analyte", all.x=TRUE)
  
  # Subset columns
  # make HQ_TPCB for 2 meals/month:
  # (concentration ng/g * 2 meals/month * 227 g/meal)/(rfd ng/kgBWday * 80 kgBW * 30.4 days/month)
  # make CR_TPCB for 2 meals/month:
  # (concentration ng/g * 2 meals/month * 227 g/meal * CSF kgBWday/ng)/(80 kgBW * 30.4 days/month)
  # make EPA_CR_TPCB for 24 meals per yr over 26 yrs with 70 yr averaging time:
  # (concentration ng/g * 227 g/meal * 24 meals/year * 26 years * CSF kgBWday/ng)/(80 kgBW * 25550 days)
  TPCBs_index <- TPCBs_merged %>% 
    select(ID, Analyte, Result, rfd, CSF) %>% 
    mutate(HQ_TPCB = (Result*2*227)/(rfd * 80 * 30.4)) %>% 
    mutate(HQ_TPCB = signif(HQ_TPCB, digits = 3)) %>%
    mutate(CR_TPCB = (Result*2*227*CSF)/(80 * 30.4)) %>% 
    mutate(CR_TPCB = signif(CR_TPCB, digits = 3)) %>%
    mutate(EPA_CR_TPCB = (Result*227*24*26*CSF)/(80*25550)) %>% 
    mutate(EPA_CR_TPCB = signif(EPA_CR_TPCB, digits = 3)) %>% 
    select(-rfd, -CSF)
  
  # Pivot TPCB_HQ data to wide format.
  TPCBs_index_wide <- TPCBs_index %>% pivot_wider(names_from = Analyte, values_from = Result)
  
  # Change congener NA values and values < DL to 0.
  # data protocol involved setting values below detection limits to 0 ng/g
  # NA values existed for a few samples that received testing but were not quantified. These were also set to 0
  PCB_congeners_noNA <- PCB_congeners %>% mutate(Result = if_else(Result < Reporting_Limit | is.na(Result), 0, Result)) 
  
  # Subset columns from congener data and pivot to wide format.
  PCB_congeners_wide <- PCB_congeners_noNA %>% select(ID, Analyte, Result) %>% pivot_wider(names_from = Analyte, values_from = Result)
  
  # Join wide format TPCB/HQ and congener data.
  PCBs_full <- full_join(TPCBs_index_wide, PCB_congeners_wide,by = "ID")
  
  # Create a new dataframe with the desired name in the global environment.
  assign(paste0('PCBs_full', i), PCBs_full)
  
} # End loop through each CSV file.

# Combine processed wide PCB dataframes into final dataframe.
PCBs_processed <- bind_rows(PCBs_full1,PCBs_full2,PCBs_full3,PCBs_full4, PCBs_full5, .id = "original dataset") %>% rename("TPCBs" = "TOTAL PCBs")

# Process the data.
PCBs_processed <- PCBs_processed %>% select(-`original dataset`)

# Write the csv file.
write_csv(PCBs_processed, "WFCS_PCBs_processed.csv")

