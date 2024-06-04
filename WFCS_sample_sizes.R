# --------------------------------------------
# WFCS Sample Sizes
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
# Description: This code computes sample sizes.
# 
# Funding: This work was supported under the project 
# name 'Contaminant Loads in Waterfowl of the Northeast Atlantic
# Flyway: New Threats and Outdated Advisories' of the 
# Multistate Conservation Grant Program.
#
# --------------------------------------------

# Load required package.
library(tidyverse) # version 2.0.0.

# Load metadata.
metadata <- read_csv("Waterfowl Contaminant Study Sample Collection Metadata 2021-2022.csv")

# Group ABDU/MALL hybrid as ABDU (mallards are over represented and hybrids are similar to both).
metadata <- metadata %>% mutate(species = if_else(species == "ABDU/MALL", "ABDU", species))

# Make a table of sample sizes for collection (IDs = N), analysis of each contaminant, and analysis of all contaminants.
result_table <- metadata %>%
  group_by(species) %>%
  summarise(IDs = n(),
            Hg_test_Y = sum(Hg_test == "Y"),
            dioxin_test_Y = sum(dioxin_test == "Y"),
            PFAS_test_Y = sum(PFAS_test == "Y"),
            PCBOCP_test_Y = sum(PCBOCP_test == "Y"),
            all_tests_Y = sum(Hg_test == "Y" & dioxin_test == "Y" & PFAS_test == "Y" & PCBOCP_test == "Y")
  )

# View the table.
view(result_table)

# Make a dataframe with samples that had at least 1 contaminant test.
metadata2 <- metadata %>%
  filter_at(vars(2:5), any_vars(. == "Y"))
metadata2 %>%
  group_by(ecoregion, sample_group, species) %>% 
  summarize(num_samples = n()) %>% view()
