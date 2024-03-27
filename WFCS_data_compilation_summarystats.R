# --------------------------------------------
# WFCS Data Compilation and Summary Stats
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
# Description: This code compiles the processed data and creates
# the summary statistics.
# 
# Funding: This work was supported under the project 
# name 'Contaminant Loads in Waterfowl of the Northeast Atlantic
# Flyway: New Threats and Outdated Advisories' of the 
# Multistate Conservation Grant Program.
#
# --------------------------------------------

# Load required packages.
library(tidyverse)
library(knitr)

#load metadata
metadata <- read_csv("Waterfowl Contaminant Study Sample Collection Metadata 2021-2022.csv")

# Load processed data.
Hg_df <- read_csv("WFCS_Hg_processed.csv")
dioxin_df <- read_csv("WFCS_dioxins_processed.csv")
PCB_df <- read_csv("WFCS_PCBs_processed.csv")
OCP_df <- read_csv("WFCS_OCPs_processed.csv")
PFAS_df <- read_csv("WFCS_PFAS_processed.csv")

# Merge processed data.
list_to_merge <- list(Hg_df, dioxin_df, PCB_df, OCP_df, PFAS_df)
contaminant_df <- Reduce(function(x, y) merge(x, y, by = "ID", all = TRUE), list_to_merge)

# Make a data frame with only concentrations.
processed_concentrations <- contaminant_df %>% select(-matches("HQ|EPA|CR"))

#merge concentration data with metadata
concentrations_metadata <- merge(metadata, processed_concentrations, by = "ID", all = TRUE)

# Print the csv file of all concentrations.
write_csv(concentrations_metadata, "Waterfowl Contaminant Study Final Processed Data 2021-2022.csv")

# Make a data frame for summary data for all concentrations.
summary_stats_conc <- processed_concentrations %>%
  pivot_longer(cols = -ID, names_to = "Contaminant", values_to = "Concentration") %>%
  group_by(Contaminant) %>%
  summarize(
    Min = signif(min(Concentration, na.rm = TRUE), digits = 3),
    Max = signif(max(Concentration, na.rm = TRUE), digits = 3),
    Median = signif(median(Concentration, na.rm = TRUE), digits = 3),
    Mean = signif(mean(Concentration, na.rm = TRUE), digits = 3),
    DF = signif(mean(Concentration > 0, na.rm = TRUE) * 100, digits = 3),
    N_tested = sum(!is.na(Concentration))
  )

# Print the csv file for summary of concentrations.
write_csv(summary_stats_conc, "WFCS_concentration_summary_stats.csv")

#make a dataframe for summary stats of dioxins
summary_dioxins <-  processed_concentrations %>% 
  select(1, 3:28) %>% 
  pivot_longer(cols = -ID, names_to = "Contaminant", values_to = "Concentration") %>%
  group_by(Contaminant) %>%
  summarize(
    Min = signif(min(Concentration, na.rm = TRUE), digits = 3),
    Max = signif(max(Concentration, na.rm = TRUE), digits = 3),
    Median = signif(median(Concentration, na.rm = TRUE), digits = 3),
    Mean = signif(mean(Concentration, na.rm = TRUE), digits = 3),
    DF = signif(mean(Concentration > 0, na.rm = TRUE) * 100, digits = 3),
    N_tested = sum(!is.na(Concentration))
  )

write.csv(summary_dioxins, "Dioxin summary stats.csv")

#make a dataframe for summary stats of PCBs
summary_PCBs <- processed_concentrations %>% 
  select(1, 29:206) %>% 
  pivot_longer(cols = -ID, names_to = "Contaminant", values_to = "Concentration") %>%
  group_by(Contaminant) %>%
  summarize(
    Min = signif(min(Concentration, na.rm = TRUE), digits = 3),
    Max = signif(max(Concentration, na.rm = TRUE), digits = 3),
    Median = signif(median(Concentration, na.rm = TRUE), digits = 3),
    Mean = signif(mean(Concentration, na.rm = TRUE), digits = 3),
    DF = signif(mean(Concentration > 0, na.rm = TRUE) * 100, digits = 3),
    N_tested = sum(!is.na(Concentration))
  )

write.csv(summary_PCBs, "PCB summary stats.csv")  

#make a dataframe for summary stats of OCPs
summary_OCPs <- processed_concentrations %>% 
  select(1, 207:234) %>% 
  pivot_longer(cols = -ID, names_to = "Contaminant", values_to = "Concentration") %>%
  group_by(Contaminant) %>%
  summarize(
    Min = signif(min(Concentration, na.rm = TRUE), digits = 3),
    Max = signif(max(Concentration, na.rm = TRUE), digits = 3),
    Median = signif(median(Concentration, na.rm = TRUE), digits = 3),
    Mean = signif(mean(Concentration, na.rm = TRUE), digits = 3),
    DF = signif(mean(Concentration > 0, na.rm = TRUE) * 100, digits = 3),
    N_tested = sum(!is.na(Concentration))
  )

write.csv(summary_OCPs, "OCP summary stats.csv")

# make a dataframe for summary stats of PFAS
summary_PFAS <- processed_concentrations %>% 
  select(1, 235:274) %>% 
  pivot_longer(cols = -ID, names_to = "Contaminant", values_to = "Concentration") %>%
  group_by(Contaminant) %>%
  summarize(
    Min = signif(min(Concentration, na.rm = TRUE), digits = 3),
    Max = signif(max(Concentration, na.rm = TRUE), digits = 3),
    Median = signif(median(Concentration, na.rm = TRUE), digits = 3),
    Mean = signif(mean(Concentration, na.rm = TRUE), digits = 3),
    DF = signif(mean(Concentration > 0, na.rm = TRUE) * 100, digits = 3),
    N_tested = sum(!is.na(Concentration))
  )

write.csv(summary_PFAS, "PFAS summary stats.csv")


# Make a vector with contaminants for the table in the main text.
contaminants_to_select <- c("Hg", "TEQ", "TPCBs", "HCH, gamma", "Heptachlor Epoxide", "Aldrin", "PFOS", "PFOA", "PFNA")

# Make a data frame for summary data for contaminants for main text.
select_summary_stats <- summary_stats_conc %>% filter(Contaminant %in% contaminants_to_select)

# Make a data frame with group HQs, CRs, and EPACRs.
WFCS_risk_levels <- contaminant_df %>% 
  select(ID,
         HQ_Hg, 
         HQ_TEQ, 
         CR_TEQ, 
         EPA_CR_TEQ,
         HQ_TPCB, 
         CR_TPCB, 
         EPA_CR_TPCB,
         Sum_OCPs_HQ, 
         Sum_OCPs_CR, 
         Sum_OCPs_EPACR,
         Sum_PFAS_HQ, 
         Sum_PFAS_CR, 
         Sum_PFAS_EPACR)

# Make total HQ, CR, and EPACR values in WFCS_risk_levels data frame.
WFCS_risk_levels <- WFCS_risk_levels %>% 
  mutate(
    Total_HI = signif(rowSums(select(., contains("HQ")), na.rm = FALSE), digits = 3),
    Total_EPACR = signif(rowSums(select(., contains("EPA")), na.rm = FALSE), digits = 3),
    Total_CR = signif(rowSums(select(., contains("CR") & !contains("EPA")), na.rm = FALSE), digits = 3)
  )

# Write a csv file for overall risk level data.
write_csv(WFCS_risk_levels, "WFCS_group_risk_levels.csv")

