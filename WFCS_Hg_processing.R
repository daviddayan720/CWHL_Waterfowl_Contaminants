# --------------------------------------------
# WFCS Hg processing
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
# Description: This code processes Hg data from Hale Creek Field Station
# and makes a new csv for processed Hg data.
# 
# Funding: This work was supported under the project 
# name 'Contaminant Loads in Waterfowl of the Northeast Atlantic
# Flyway: New Threats and Outdated Advisories' of the 
# Multistate Conservation Grant Program.

# --------------------------------------------

# Load required packages.
library(tidyverse) # Version 2.0.0

# Read in study metadata.
metadata <- read_csv("Waterfowl Contaminant Study Sample Collection Metadata 2021-2022.csv")

# Read in raw Hg data from Hale Creek Field Station (HFCS).
Hg_raw <- read_csv("Waterfowl Muscle Tissue Total Mercury Data from Hale Creek Field Station 2021-22.csv")

# Remove blank rows from Hg_raw (i.e., rows without a TAGNO).
mercury <- Hg_raw %>% 
  filter(!is.na(TAGNO))

# Check TAGNO values from HFCS data.
unique_TAGNO <- c(unique(mercury$TAGNO))
# 103 samples were analyzed, but 2 samples mistakenly had identical TAGNO
# based on SDATE and LOCATION:
# LABNO 21-1548-H corresponds to metadata ID NJ_MALL_01_AD
# LABNO 21-1556-H corresponds to metadata ID NJ_MALL_40_AD
# Rename TAGNO for LABNO 21-1556-H to NJ_MALL_40_AD
mercury2 <- mercury %>% 
  mutate(TAGNO = if_else(LABNO == "21-1556-H","NJ_MALL_40_AD", TAGNO))
# Check to make sure there are 103 unique TAGNO.
unique_TAGNO2 <- c(unique(mercury2$TAGNO)) #TRUE

# Determine which TAGNO do not match their corresponding metadata ID.
TAGNO_mismatch <- mercury2 %>% 
  filter(!TAGNO %in% metadata$ID)

# TAGNO from NJ are missing a '_' compared to metadata ID.
# example: TAGNO NJ_WODU_11AD corresponds to ID NJ_WODU_11_AD

# TAGNO starting with B correspond to metadata ID by last digit.
# examples: 
# B921078_1 corresponds to ID 1
# B921078_2 corresponds to ID 2
wrong_TAGNO <- c(TAGNO_mismatch$TAGNO)
right_ID <- c("NJ_MALL_60_AD", "NJ_MALL_10_AD", "NJ_ABDU_02_LC", 
              "NJ_AGWT_01_LC", "NJ_ABDU_30_AD", "NJ_CAGO_30_AD", 
              "NJ_WODU_03_AD", "NJ_WODU_02_AD", "NJ_MALL_11_AD", 
              "NJ_WODU_01_JK", "NJ_CAGO_01_AD", "NJ_MALL_01_AD", 
              "NJ_AGWT_02_AD", "NJ_WODU_01_AD", "NJ_MALL_30_AD", 
              "NJ_WODU_11_AD", "NJ_WODU_03_LC", "NJ_ABDU_01_AD", 
              "NJ_AGWT_01_AD", "NJ_WODU_01_TN", "NJ_AGWT_02_NJ", 
              "1", "2", "6", "7", "8")
replacement_rules <- tibble(
  wrong_TAGNOs = wrong_TAGNO,  # Insert the wrong values here
  right_IDs = right_ID)    # Insert the corresponding correct values here

# Replace incorrect TAGNO with corresponding ID.
mercury3 <- mercury2 %>%
  mutate(TAGNO = case_when(
    TAGNO %in% replacement_rules$wrong_TAGNOs ~ 
      replacement_rules$right_IDs[match(TAGNO, replacement_rules$wrong_TAGNOs)],
    TRUE ~ TAGNO  # Keep original value if it doesn't match any replacement rule
  ))

# Check that there are no more mismatched TAGNO and ID.
TAGNO_mismatch2 <- mercury3 %>% 
  filter(!TAGNO %in% metadata$ID) #TRUE

# Select only TAGNO and Hg columns.
mercury4 <- mercury3 %>% 
  select(TAGNO, Hg)

# Pivot data to long format.
mercury5 <- mercury4 %>%
  pivot_longer(cols = -TAGNO, names_to = "Analyte", values_to = "Result")

# Make a new column with the MDL (method detection limit).
# The MDL from HCFS was 0.004 µg/g.
mercury6 <- mercury5 %>% 
  mutate(MDL = 0.004)

# Change Results < MDL to 0, change Result from µg/g to ng/g (*1000), and remove MDL.
mercury7 <- mercury6 %>% 
  mutate(Result = if_else(Result < MDL, 0, Result)) %>% 
  mutate(Result = Result * 1000) %>% 
  select(TAGNO, Analyte, Result)

# Read in reference doses.
rfds <- read_csv("reference_doses.csv")
rfds2 <- rfds %>% select(Analyte, rfd)

# Merge Hg data with reference doses.
mercury_merged <- merge(mercury7, rfds2, by = "Analyte", all.x = TRUE)

# Make HQ_Hg to 3 sigfigs for 2 meals/month: 
# (concentration ng/g * 2 meals/month * 227 g/meal)/(rfd ng/kgBWday * 80 kgBW * 30.4 days/month)
Hg_HQ <- mercury_merged %>% 
  mutate(HQ_Hg = (Result*2*227)/(rfd * 80 * 30.4)) %>%
  mutate(HQ_Hg = signif(HQ_Hg, digits = 3)) %>% 
  select(-rfd)

# Rename TAGNO to ID and pivot to wide format.
Hg_processed <- Hg_HQ %>%
  rename(ID = TAGNO) %>% 
  pivot_wider(names_from = Analyte, values_from = Result)

# Write the csv file.
write_csv(Hg_processed, "WFCS_Hg_processed.csv")
