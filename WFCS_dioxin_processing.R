# --------------------------------------------
# WFCS Dioxin processing
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
# Description: This code processes dioxin data from Pace Analytical
# and makes a new csv for processed dioxin data.
# 
# Funding: This work was supported under the project 
# name 'Contaminant Loads in Waterfowl of the Northeast Atlantic
# Flyway: New Threats and Outdated Advisories' of the 
# Multistate Conservation Grant Program.
#
# --------------------------------------------

# Load required pacakges.
library(tidyverse)

# Read in study metadata.
metadata <- read_csv("Waterfowl Contaminant Study Sample Collection Metadata 2021-2022.csv")

# Read in raw dioxin data from Pace Analytical.
dioxins <- read_csv("Waterfowl Muscle Tissue Dioxin and Furan Data from Pace Analytical 2021-22.csv")

# Remove duplicate analyses and Percent Moisture and Percent Lipid
dioxins2 <- dioxins %>% 
  filter(!grepl("DUP", ID)) %>% 
  filter(!Analyte %in% c("Moisture_Percent", "Lipid_Percent"))

# Check ID values from Pace Analytical
dioxin_IDs <- unique(dioxins2$ID) #104 samples were analyzed, 104 unique IDs

# Check that Pace IDs match metadata IDs
metadata_dioxintest <- metadata %>% 
  filter(dioxin_test == "Y")
ID_mismatch <- dioxins2 %>% 
  filter(!ID %in% metadata_dioxintest$ID) # no mismatched IDs

# Change ND Results to 0
dioxins3 <- dioxins2 %>%
  mutate(Result = if_else(Result == "ND", "0", Result))
dioxins3$Result = as.numeric(dioxins3$Result)

# check to make sure no NAs were introduced
sum(is.na(dioxins3$Result)) #no NAs

# Change ng/kg to ng/g (/1000).
dioxins4 <- dioxins3 %>% 
  mutate(Result = Result / 1000)

# Separate raw data into TEQ and congener data frames.
congeners <- dioxins4 %>% filter(Analyte != "TEQ")
TEQs <- dioxins4 %>% filter(Analyte == "TEQ")

# Read in reference doses.
rfds <- read_csv("reference_doses.csv")
rfds2 <- rfds %>% select(Analyte, rfd, unit)

# Read in cancer slope factors.
CSFs <- read_csv("cancer_slope_factors.csv")
CSFs2 <- CSFs %>% select(Analyte, CSF)

# Merge TEQ data with reference doses and CSFs.
TEQ_merged <- merge(merge(TEQs, rfds2, by = "Analyte", all.x = TRUE),CSFs2, by = "Analyte", all.x = TRUE)

# Remove result qualifiers, rfd units, and EDL
  #Make HQ_TEQ for 2 meals/month:
    # (concentration ng/g * 2 meals/month * 227 g/meal)/(rfd ng/kgBWday * 80 kgBW * 30.4 days/month)
  # Make CR_TEQ for 2 meals/month:
    # (concentration ng/g * 2 meals/month * 227 g/meal * CSF kgBWday/ng)/(80 kgBW * 30.4 days/month)
  # Make EPA_CR_TEQ for 24 meals per yr over 26 yrs with 70 yr averaging time:
    # (concentration ng/g * 227 g/meal * 24 meals/year * 26 years * CSF kgBWday/ng)/(80 kgBW * 25550 days)
TEQ_index <- TEQ_merged %>% 
  select(-Result_Qualifier, -EDL, -unit) %>% 
  mutate(HQ_TEQ = (Result*2*227)/(rfd * 80 * 30.4)) %>%
  mutate(HQ_TEQ = signif(HQ_TEQ, digits = 3)) %>% 
  mutate(CR_TEQ = (Result*2*227*CSF)/ (80 * 30.4)) %>%
  mutate(CR_TEQ = signif(CR_TEQ, digits = 3)) %>% 
  mutate(EPA_CR_TEQ = (Result*227*24*26*CSF)/(80*25550)) %>% 
  mutate(EPA_CR_TEQ = signif(EPA_CR_TEQ, digits = 3)) %>% 
  select(-rfd, - CSF)

# Pivot TEQ and TEQ HQ data to wide format.
TEQ_wide <- TEQ_index %>% pivot_wider(names_from = Analyte, values_from = Result)

# remove result qualifier and EDL from congener data and pivot to wide
congeners_wide <- congeners %>% 
  select(-EDL, -Result_Qualifier) %>% 
  pivot_wider(names_from = Analyte, values_from = Result)

# Join wide format TEQ/HQ and wide format congener data.
dioxins_processed <- full_join(TEQ_wide, congeners_wide, by = "ID")

# Write the dioxins data.
write_csv(dioxins_processed, "WFCS_dioxins_processed.csv")
