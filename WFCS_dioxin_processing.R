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

# Load required package.
library(tidyverse) # Version 2.0.0

# Load study metadata.
metadata <- read_csv("Waterfowl Contaminant Study Sample Collection Metadata 2021-2022.csv")

# Load raw dioxin data from Pace Analytical.
dioxins <- read_csv("Waterfowl Muscle Tissue Dioxin and Furan Data from Pace Analytical 2021-22.csv")

# Remove duplicate analyses and Percent Moisture and Percent Lipid.
dioxins2 <- dioxins %>% 
  filter(!grepl("DUP", ID)) %>% 
  filter(!Analyte %in% c("Moisture_Percent", "Lipid_Percent"))

# Check ID values from Pace Analytical.
dioxin_IDs <- unique(dioxins2$ID) #104 samples were analyzed, 104 unique IDs

# Check that Pace IDs match metadata IDs.
metadata_dioxintest <- metadata %>% 
  filter(dioxin_test == "Y")
ID_mismatch <- dioxins2 %>% 
  filter(!ID %in% metadata_dioxintest$ID) # no mismatched IDs

# Change ND Results to 0.
dioxins3 <- dioxins2 %>%
  mutate(Result = if_else(Result == "ND", "0", Result))
dioxins3$Result = as.numeric(dioxins3$Result)

# Check to make sure no NAs were introduced.
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

# Remove result qualifiers, rfd units, and EDL.
# Make HQ_TEQ for 2 meals/month:
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

# Remove result qualifier and EDL from congener data and pivot to wide.
congeners_wide <- congeners %>% 
  select(-EDL, -Result_Qualifier) %>% 
  pivot_wider(names_from = Analyte, values_from = Result)

# Join wide format TEQ/HQ and wide format congener data.
dioxins_processed <- full_join(TEQ_wide, congeners_wide, by = "ID")

# Write the dioxins data.
write_csv(dioxins_processed, "WFCS_dioxins_processed.csv")


# EMPC and 2024 TEF Sensitivity Testing -----------------------------------

toxic_congeners <- unique(congeners$Analyte)[!grepl("Total", unique(congeners$Analyte))]
TEFs <- data.frame(Analyte = toxic_congeners,
                   TEF_2005 = c(.1, 1,.03,
                                .3, 1,.1, 
                                .1,.1,.1,
                                .1,.1,.1,
                                .01,.01,.01,
                                .0003, .0003),
                   TEF_2024 = c(.07, 1, .01,
                                .1,.4,.3,
                                .09,.1,.2,
                                .09,.07,.05,
                                .02,.1,.05,
                                .002, .001))
congeners_TEF <- full_join(congeners, TEFs, by = "Analyte") %>% filter(Analyte %in% toxic_congeners)

congeners_equiv <- congeners_TEF %>% 
  mutate(equiv_2024 = Result * TEF_2024,
         equiv_2005 = Result * TEF_2005,
         equiv_2005_no_empc = if_else(Result_Qualifier %in% c("IJ", "PJ", "I", "P"), 0, Result * TEF_2005)
         )
manual_TEQs <- congeners_equiv %>% 
  group_by(ID) %>% 
  summarize(TEQ_2024 = sum(equiv_2024),
            TEQ_2005 = sum(equiv_2005),
            TEQ_2005_empc = sum(equiv_2005_no_empc))
manual_TEQs_check <- full_join(manual_TEQs, TEQ_wide %>% select(ID, TEQ), by = "ID") %>% 
  mutate(manual_check = TEQ_2005/TEQ,
         EMPC_check = 100- 100*(TEQ_2005_empc / TEQ_2005),
         TEF_check = 100* ((TEQ_2024 - TEQ_2005)/TEQ_2005))
manual_TEQs_check %>% 
  summarize(
    min_EMPC_contribution = min(EMPC_check, na.rm = TRUE),
    max_EMPC_contribution = max(EMPC_check, na.rm = TRUE),
    mean_EMPC_contribution = mean(EMPC_check, na.rm = TRUE),
    median_EMPC_contribution = median(EMPC_check, na.rm = TRUE),
    min_TEF_change = min(TEF_check, na.rm = TRUE),
    max_TEF_change = max(TEF_check, na.rm = TRUE),
    mean_TEF_change = mean(TEF_check, na.rm = TRUE),
    median_TEF_change = median(TEF_check, na.rm = TRUE),
  ) %>% 
  pivot_longer(cols = 1:8,
               names_to = "Summary", 
               values_to = "Value")
