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
# Description: This code processes Hg data and makes 
# a new csv for processed Hg data.
# 
# Funding: This work was supported under the project 
# name 'Contaminant Loads in Waterfowl of the Northeast Atlantic
# Flyway: New Threats and Outdated Advisories' of the 
# Multistate Conservation Grant Program.

# --------------------------------------------

# Load required packages.
library(tidyverse)
library(stringr)

# Read in study metadata.
metadata <- read_csv("Waterfowl Contaminant Study Sample Collection Metadata 2021-2022.csv")

# Read in raw Hg data from Hale Creek Field Station (HFCS).
Hg_raw <- read_csv("Waterfowl Muscle Tissue Total Mercury Data from Hale Creek Field Station 2021-22.csv")

# Remove blank rows from Hg_raw (i.e., rows without a TAGNO).
mercury <- Hg_raw %>% 
  filter(!is.na(TAGNO))

# Check TAGNO values from HFCS data
unique_TAGNO <- c(unique(mercury$TAGNO))
# 103 samples were analyzed, but 2 samples mistakenly had identical TAGNO
  # based on SDATE and LOCATION:
    # LABNO 21-1548-H corresponds to metadata ID NJ_MALL_01_AD
    # LABNO 21-1556-H corresponds to metadata ID NJ_MALL_40_AD
# rename TAGNO for LABNO 21-1556-H to NJ_MALL_40_AD
mercury2 <- mercury %>% 
  mutate(TAGNO = if_else(LABNO == "21-1556-H","NJ_MALL_40_AD", TAGNO))

# determine which TAGNO do not match their corresponding metadata ID
TAGNO_mismatch <- mercury2 %>% 
  filter(!TAGNO %in% metadata$ID)
# TAGNO from NJ are missing a '_' compared to metadata ID
# TAGNO starting with B correspond to metadata ID by last digit
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
              "NJ_AGWT_01_AD", "NJ_WODU_01_TN", "NJ_AGWT_02N_J", 
              "1", "2", "6", "7", "8")

replacement_rules <- tibble(
  wrong_TAGNOs = wrong_TAGNO,  # Insert the wrong values here
  right_IDs = right_ID)    # Insert the corresponding correct values here


# Assuming mercury2 is your data frame
mercury3 <- mercury2 %>%
  mutate(TAGNO = case_when(
    TAGNO %in% replacement_rules$wrong_TAGNOs ~ 
      replacement_rules$right_IDs[match(TAGNO, replacement_rules$wrong_TAGNOs)],
    TRUE ~ TAGNO  # Keep original value if it doesn't match any replacement rule
  ))
