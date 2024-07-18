# --------------------------------------------
# WFCS Determinisitc risk assessment stat analysis
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
# Description: This code includes the deterministic 
# risk assessment, statistical analyses, and visualizations.
# 
# Funding: This work was supported under the project 
# name 'Contaminant Loads in Waterfowl of the Northeast Atlantic
# Flyway: New Threats and Outdated Advisories' of the 
# Multistate Conservation Grant Program.
#
# --------------------------------------------

# Load required packages.
library(tidyverse) # Version 2.0.0
library(ggbreak) # Version 0.1.2
library(ggpubr) # Version 0.6.0
library(knitr) # Version 1.43
library(cowplot) # Version 1.1.3
library(dunn.test) # Version 1.3.5

# Load risk data.
risk_data <- read_csv("WFCS_group_risk_levels.csv")

# Load metadata.
metadata <- read_csv("Waterfowl Contaminant Study Sample Collection Metadata 2021-2022.csv")

# Transform ABDU/MALL hybrid into ABDU.
metadata <- metadata %>% mutate(species = if_else(species == "ABDU/MALL", "ABDU", species))

# Merge metadata and risk data.
risk_data_merged <- merge(risk_data, metadata, by = "ID", all.x = TRUE)

# Factor the species in merged risk data.
risk_data_merged$species <- factor(risk_data_merged$species,levels = c("MALL", "ABDU", "AGWT", "CAGO", "WODU"))

# Non-Cancer Risk by Species Stat Comparison ----------------------------------------

# Create a list of variables.
variables_to_test_HQ <- c("HQ_Hg", "HQ_TEQ", "HQ_TPCB", "Sum_OCPs_HQ", "Sum_PFAS_HQ", "Total_HI")

# Compare the variables
for (variable in variables_to_test_HQ) {
  result <- dunn.test(x = risk_data_merged[[variable]], g = risk_data_merged$species, method = "holm")
}

# ---------------------------------------------------------------------------
# Total HI by Species Histogram ---------------------------------------------

max_value_HI <- max(risk_data_merged$Total_HI, na.rm = TRUE) + 0.5

# Calculate the number of non-NA values for each species.
species_counts_HI <- risk_data_merged %>% group_by(species) %>% summarize(count = sum(!is.na(Total_HI)))

# Create custom labeller function for facet labels.
label_species_HI <- function(variable) {
  count <- species_counts_HI$count[species_counts_HI$species == variable]
  return(paste(variable, " (n=", count, ")", sep = ""))
} # End label species HI function.

# Create helper function to generate breaks and labels.
generate_breaks_labels <- function(max_value, bin_width, label_interval) {
  breaks <- seq(0, max_value, by = bin_width)
  labels <- ifelse(breaks %% label_interval == 0, as.character(breaks), "")
  return(list(breaks = breaks, labels = labels))
} # End generate breaks function.

# Use the generate breaks labels function on the data.
breaks_labels_HI <- generate_breaks_labels(max_value_HI, 0.25, 0.5)

# Make a data frame of factors with significance labels from Dunn's Test.
sig_labels_HI <- data.frame(species = factor(c("MALL", "ABDU", "AGWT", "CAGO", "WODU")), label = c("a", "a", "a", "b", "b"))

# Produce the plot.
Total_HI_histo <- risk_data_merged %>%
  filter(!is.na(Total_HI)) %>%
  ggplot(aes(x = Total_HI)) +
  geom_histogram(breaks = breaks_labels_HI$breaks) +
  geom_vline(xintercept = 1, linetype = "solid", color = "red") +
  scale_x_continuous(breaks = breaks_labels_HI$breaks, labels = breaks_labels_HI$labels) +
  facet_wrap(~ species, ncol = 1, strip.position = "left", labeller = as_labeller(label_species_HI)) +
  labs(x = "Total Hazard Index",y = "Sample Count") +
  scale_x_break(c(9.9, 36.5)) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text = element_text(hjust = 0.5, size = 8),
        axis.title.x = element_text(margin = margin(t = -8, unit = "pt")),
        axis.title.y = element_text(margin = margin(r = -8, unit = "pt")))+
  geom_text(aes(x = 37.25 , y = 12, label = label), data = sig_labels_HI, size = 5)

# Save the histogram.
ggsave("Total_HI_by_species.png", plot = Total_HI_histo, dpi = 1000, width = 7.8, height = 5.83)

# ---------------------------------------------------------------------------
# Hg HQ by Species Histogram ------------------------------------------------

# Calculate the number of non-NA values for each species.
species_counts_Hg <- risk_data_merged %>% group_by(species) %>% summarize(count = sum(!is.na(HQ_Hg)))

# Make a data frame of factors with significance labels from Dunn's Test.
sig_labels_Hg <- data.frame(species = factor(c("MALL", "ABDU", "AGWT", "CAGO", "WODU")), label = c("bc", "a", "ab", "d", "c"))

# Create custom labeller function for facet labels.
label_species_Hg <- function(variable) {
  count <- species_counts_Hg$count[species_counts_Hg$species == variable]
  return(paste(variable, " (n=", count, ")", sep = ""))
} # End custom labeller function for facet labels.

# Produce the plot.
Hg_HQ_histo <- risk_data_merged %>% 
  filter(!is.na(HQ_Hg)) %>% 
  ggplot(aes(x = HQ_Hg)) +
  geom_histogram(breaks = seq(0,1.1, by = .05))+
  geom_vline(xintercept = 1, linetype = "solid", color = "red") +
  scale_x_continuous(sec.axis = dup_axis(name = NULL), breaks = seq(0, 1.1, by = .05), labels = seq(0,1.1, by = .05)) +  facet_wrap(~ species, ncol = 1, strip.position = "left", labeller = as_labeller(label_species_Hg)) +
  labs(x = "Mercury Hazard Quotient",y = "Sample Count") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text = element_text(hjust = 0.5),
        axis.title.x = element_text(margin = margin(t = 2, unit = "pt")),
        axis.title.y = element_text(margin = margin(r = 2, unit = "pt"))) +
  theme(strip.text = element_text(size = 8))+
  geom_text(aes(x = 1.1 , y = 12, label = label), data = sig_labels_Hg, size = 5)

# Save the plot.
ggsave("Hg_HQ_by_species.png", plot = Hg_HQ_histo, dpi = 1000, width = 7.8, height = 5.83)

# ---------------------------------------------------------------------------
# TEQ HQ by Species Histogram -----------------------------------------------

max_value_TEQ <- max(risk_data_merged$HQ_TEQ, na.rm = TRUE) + 0.5

# Calculate the number of non-NA values for each species.
species_counts_TEQ <- risk_data_merged %>% group_by(species) %>% summarize(count = sum(!is.na(HQ_TEQ)))

# Create custom labeler function for facet labels.
label_species_TEQ <- function(variable) {
  count <- species_counts_TEQ$count[species_counts_TEQ$species == variable]
  return(paste(variable, " (n=", count, ")", sep = ""))
} # End create custom labeler.

# Create helper function to generate breaks and labels.
generate_breaks_labels_TEQ <- function(max_value, bin_width, label_interval) {
  breaks <- seq(0, max_value, by = bin_width)
  labels <- ifelse(breaks %% label_interval == 0, as.character(breaks), "")
  return(list(breaks = breaks, labels = labels))
} # End create helper function.

# Apply helper function to data. 
breaks_labels_TEQ <- generate_breaks_labels_TEQ(max_value_TEQ, 0.25, 0.5)

# Make a data frame of factors with significance labels from Dunn's Test.
sig_labels_TEQ <- data.frame(species = factor(c("MALL", "ABDU", "AGWT", "CAGO", "WODU")), label = c("a", "a", "ab", "bc", "c"))

# Produce the plot. 
TEQ_HQ_histo <- risk_data_merged %>%
  filter(!is.na(HQ_TEQ)) %>%
  ggplot(aes(x = HQ_TEQ)) +
  geom_histogram(breaks = breaks_labels_TEQ$breaks) +
  geom_vline(xintercept = 1, linetype = "solid", color = "red") +
  scale_x_continuous(breaks = breaks_labels_TEQ$breaks, labels = breaks_labels_TEQ$labels) +
  facet_wrap(~ species, ncol = 1, strip.position = "left", labeller = as_labeller(label_species_TEQ)) +
  labs(x = "TCDD Toxic Equivalency Hazard Quotient",y = "Sample Count") +
  scale_x_break(c(4.9, 13.5)) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text = element_text(hjust = 0.5),
        axis.title.x = element_text(margin = margin(t = -8, unit = "pt")),
        axis.title.y = element_text(margin = margin(r = -8, unit = "pt"))) +
  theme(strip.text = element_text(size = 8))+
  geom_text(aes(x = 14 , y = 15, label = label), data = sig_labels_TEQ, size = 5)

# Save the plot. 
ggsave("TEQ_HQ_by_species.png", plot = TEQ_HQ_histo, dpi = 1000, width = 7.8, height = 5.83)

# ---------------------------------------------------------------------------
# TPCB HQ by Species Histogram ----------------------------------------------

max_value_TPCB <- max(risk_data_merged$HQ_TPCB, na.rm = TRUE) + 0.5

# Calculate the number of non-NA values for each species.
species_counts_TPCB <- risk_data_merged %>% group_by(species) %>% summarize(count = sum(!is.na(HQ_TPCB)))

# Create custom labeller function for facet labels.
label_species_TPCB <- function(variable) {
  count <- species_counts_TPCB$count[species_counts_TPCB$species == variable]
  return(paste(variable, " (n=", count, ")", sep = ""))
} # End custom labeler.

# Create helper function to generate breaks and labels.
generate_breaks_labels_TPCB <- function(max_value, bin_width, label_interval) {
  breaks <- seq(0, max_value, by = bin_width)
  labels <- ifelse(breaks %% label_interval == 0, as.character(breaks), "")
  return(list(breaks = breaks, labels = labels))
} # End helper function.

# Make a data frame of factors with significance labels from Dunn's Test
sig_labels_TPCB <- data.frame(species = factor(c("MALL", "ABDU", "AGWT", "CAGO", "WODU")), label = c("a", "a", "ab", "b", "b"))

# Use the helper function on data.
breaks_labels_TPCB <- generate_breaks_labels_TPCB(max_value_TPCB, 0.25, 0.5)

# Produce the plot. 
PCB_HQ_histo <- risk_data_merged %>%
  filter(!is.na(HQ_TPCB)) %>%
  ggplot(aes(x = HQ_TPCB)) +
  geom_histogram(breaks = breaks_labels_TPCB$breaks) +
  geom_vline(xintercept = 1, linetype = "solid", color = "red") +
  scale_x_continuous(breaks = breaks_labels_TPCB$breaks, labels = breaks_labels_TPCB$labels) +
  facet_wrap(~ species, ncol = 1, strip.position = "left", labeller = as_labeller(label_species_TPCB)) +
  labs(x = expression("\u03A3PCB Hazard Quotient"),y = "Sample Count") +
  scale_x_break(c(4.5, 20.5)) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text = element_text(hjust = 0.5),
        axis.title.x = element_text(margin = margin(t = -8, unit = "pt")),
        axis.title.y = element_text(margin = margin(r = -8, unit = "pt"))) +
  theme(strip.text = element_text(size = 8))+
  geom_text(aes(x = 21.25 , y = 15, label = label), data = sig_labels_TPCB, size = 5)

# Save the plot. 
ggsave("PCB_HQ_by_species.png", plot = PCB_HQ_histo, dpi = 1000, width = 7.8, height = 5.83)

# ---------------------------------------------------------------------------
# OCPs HI by Species Histogram ----------------------------------------------

# Calculate the number of non-NA values for each species.
species_counts_OCP <- risk_data_merged %>% group_by(species) %>% summarize(count = sum(!is.na(Sum_OCPs_HQ)))

# Create custom labeler function for facet labels.
label_species_OCP <- function(variable) {
  count <- species_counts_OCP$count[species_counts_OCP$species == variable]
  return(paste(variable, " (n=", count, ")", sep = ""))
} # End labeler 

# Make a data frame of factors with significance labels from Dunn's Test.
sig_labels_OCPs <- data.frame(species = factor(c("MALL", "ABDU", "AGWT", "CAGO", "WODU")), label = c("ab", "a", "ab", "ab", "b"))

# Produce the plot.
OCP_HI_histo <- risk_data_merged %>% 
  filter(!is.na(Sum_OCPs_HQ)) %>% 
  ggplot(aes(x = Sum_OCPs_HQ)) +
  geom_histogram(breaks = seq(0, 1.1, by = .05))+
  geom_vline(xintercept = 1, linetype = "solid", color = "red") +
  scale_x_continuous(sec.axis = dup_axis(name = NULL), breaks = seq(0, 1.1, by = .05), labels = seq(0,1.1, by = .05)) +
  facet_wrap(~ species, ncol = 1, strip.position = "left", labeller = as_labeller(label_species_OCP)) +
  labs(x = "OCPs Hazard Index",y = "Sample Count") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text = element_text(hjust = 0.5),
        axis.title.x = element_text(margin = margin(t = 2, unit = "pt")),
        axis.title.y = element_text(margin = margin(r = 2, unit = "pt"))) +
  theme(strip.text = element_text(size = 8))+
  geom_text(aes(x = 1.1 , y = 15, label = label), data = sig_labels_OCPs, size = 5)

# Save the plot.
ggsave("OCPs_HI_by_species.png", plot = OCP_HI_histo, dpi = 1000, width = 7.8, height = 5.83)

# ---------------------------------------------------------------------------
# PFAS HI by Species Histogram ----------------------------------------------

max_value_PFAS <- max(risk_data_merged$Sum_PFAS_HQ, na.rm = TRUE) + 0.5

# Calculate the number of non-NA values for each species.
species_counts_PFAS <- risk_data_merged %>% group_by(species) %>% summarize(count = sum(!is.na(Sum_PFAS_HQ)))

# Create custom labeller function for facet labels.
label_species_PFAS <- function(variable) {
  count <- species_counts_PFAS$count[species_counts_PFAS$species == variable]
  return(paste(variable, " (n=", count, ")", sep = ""))
} # End create custom labeler.

# Create helper function to generate breaks and labels.
generate_breaks_labels_PFAS <- function(max_value, bin_width, label_interval) {
  breaks <- seq(0, max_value, by = bin_width)
  labels <- ifelse(breaks %% label_interval == 0, as.character(breaks), "")
  return(list(breaks = breaks, labels = labels))
} # End create helper function. 

# Make a data frame of factors with significance labels from Dunn's Test.
sig_labels_PFAS <- data.frame(species = factor(c("MALL", "ABDU", "AGWT", "CAGO", "WODU")), label = c("a", "ab", "a", "c", "bc"))

# Apply the helper function to data.
breaks_labels_PFAS <- generate_breaks_labels_PFAS(max_value_PFAS, 0.25, 0.5)

# Produce the plot. 
PFAS_HI_histo <- risk_data_merged %>%
  filter(!is.na(Sum_PFAS_HQ)) %>%
  ggplot(aes(x = Sum_PFAS_HQ)) +
  geom_histogram(breaks = breaks_labels_PFAS$breaks) +
  geom_vline(xintercept = 1, linetype = "solid", color = "red") +
  scale_x_continuous(breaks = breaks_labels_PFAS$breaks, labels = breaks_labels_PFAS$labels) +
  facet_wrap(~ species, ncol = 1, strip.position = "left", labeller = as_labeller(label_species_PFAS)) +
  labs(x = "PFAS Hazard Index",y = "Sample Count") +
  theme_bw() +
  scale_x_break(c(2.3, 4.25))+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text = element_text(hjust = 0.5),
        axis.title.x = element_text(margin = margin(t = -8, unit = "pt")),
        axis.title.y = element_text(margin = margin(r = -8, unit = "pt"))) +
  theme(strip.text = element_text(size = 8))+
  geom_text(aes(x = 4.9 , y = 15, label = label), data = sig_labels_PFAS, size = 5)

# Save the plot. 
ggsave("PFAS_HI_by_species.png", plot = PFAS_HI_histo, dpi = 1000, width = 7.8, height = 5.83)

# ---------------------------------------------------------------------------
# Non-Cancer Deterministic Summary Stats ------------------------------------

# List of variables.
variables_of_interest <- c("HQ_Hg", "HQ_TEQ", "HQ_TPCB", "Sum_OCPs_HQ", "Sum_PFAS_HQ", "Total_HI") 

# Create a new dataframe to store the statistics.
non_cancer_deterministic_summary <- data.frame(
  min = numeric(length(variables_of_interest)),
  max = numeric(length(variables_of_interest)),
  mean = numeric(length(variables_of_interest)),
  median = numeric(length(variables_of_interest)),
  confidence_interval = character(length(variables_of_interest)),
  percentile_95 = numeric(length(variables_of_interest))
)

# Calculate statistics for each variable.
for (i in seq_along(variables_of_interest)) {
  variable <- variables_of_interest[i]
  # Filter out missing values or NaNs.
  non_missing_values <- na.omit(risk_data[[variable]])
  # Calculate basic statistics.
  stats <- summary(non_missing_values)
  # Create a linear model with an intercept term.
  l_model <- lm(non_missing_values ~ 1)
  # Calculate the confidence interval.
  confidence_interval <- confint(l_model, level = 0.95)
  # Extract the lower and upper bounds.
  lower_bound <- confidence_interval[1]
  upper_bound <- confidence_interval[2]
  # Ensure that lower_bound and upper_bound are numeric.
  lower_bound <- as.numeric(lower_bound)
  upper_bound <- as.numeric(upper_bound)
  # Calculate 95th percentile.
  percentile_95 <- quantile(non_missing_values, 0.95)
  # Update the statistics dataframe.
  non_cancer_deterministic_summary[i, ] <- c(stats[1], stats[6], mean(non_missing_values), median(non_missing_values), paste(round(lower_bound, 3), "-", round(upper_bound, 3)), percentile_95)
} # End loop for variables of interest in sequence.

# Set row names to variable names.
row.names(non_cancer_deterministic_summary) <- variables_of_interest

# Display the resulting statistics dataframe.
view(non_cancer_deterministic_summary)

# Save the file. 
write_csv(non_cancer_deterministic_summary, "non_cancer_deterministic_summary.csv")

# -----------------------------------------------------------------------------------
# Total CR by Species Histogram ---------------------------------------------

# Calculate the number of non-NA values for each species.
species_counts_CR <- risk_data_merged %>% group_by(species) %>% summarize(count = sum(!is.na(Total_CR)))

# Create custom labeler function for facet labels.
label_species_CR <- function(variable) {
  count <- species_counts_CR$count[species_counts_CR$species == variable]
  return(paste(variable, " (n=", count, ")", sep = ""))
} # End create custom labeler. 

# Create the faceted histogram plot.
Total_CR_histo <- risk_data_merged %>% 
  filter(!is.na(Total_CR)) %>% 
  ggplot() +
  geom_histogram(mapping = aes(x = Total_CR)) +
  facet_wrap(~ species, ncol = 1, strip.position = "left", labeller = as_labeller(label_species_CR)) +
  scale_x_log10(sec.axis = dup_axis(name = NULL)) +
  labs(x = "Total Cancer Risk", y = "Sample Count") +
  theme_bw()+
  geom_vline(xintercept = 1e-6, linetype = "dashed", color = "orange") +
  geom_vline(xintercept = 1e-4, linetype = "solid", color = "red") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text = element_text(hjust = 0.5),
        axis.title.x = element_text(margin = margin(t = 1, unit = "pt")),
        axis.title.y = element_text(margin = margin(r = 2, unit = "pt"))) +
  theme(strip.text = element_text(size = 8))

# Save the plot. 
ggsave("Total_CR_by_species.png", plot = Total_CR_histo, dpi = 1000, width = 7.8, height = 5.83)

# ----------------------------------------------------------------------------
# TCDD TEQ CR by Species Histogram -------------------------------------------

# Calculate the number of non-NA values for each species.
species_counts_CRTEQ <- risk_data_merged %>% group_by(species) %>% summarize(count = sum(!is.na(CR_TEQ)))

# Create custom labeller function for facet labels.
label_species_CRTEQ <- function(variable) {
  count <- species_counts_CRTEQ$count[species_counts_CRTEQ$species == variable]
  return(paste(variable, " (n=", count, ")", sep = ""))
} # End create labeller. 

# Create the faceted histogram plot.
TEQ_CR_histo <- risk_data_merged %>% 
  filter(!is.na(CR_TEQ)) %>% 
  ggplot() +
  geom_histogram(mapping = aes(x = CR_TEQ)) +
  facet_wrap(~ species, ncol = 1, strip.position = "left", labeller = as_labeller(label_species_CRTEQ)) +
  scale_x_log10(sec.axis = dup_axis(name = NULL), breaks = c(1e-8, 1e-7, 1e-6, 1e-5, 1e-4, 1e-3))+
  labs(x = "TCDD Toxic Equivalency Cancer Risk", y = "Sample Count") +
  theme_bw()+
  geom_vline(xintercept = 1e-6, linetype = "dashed", color = "orange") +
  geom_vline(xintercept = 1e-4, linetype = "solid", color = "red") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text = element_text(hjust = 0.5),
        axis.title.x = element_text(margin = margin(t = 1, unit = "pt")),
        axis.title.y = element_text(margin = margin(r = 2, unit = "pt"))) +
  theme(strip.text = element_text(size = 8))

# Make barplots of zeros and non-zeros as logscale histogram excludes zeros.
# Get number of samples per species with CR_TEQ = 0 and not = 0.
risk_data_merged %>% group_by(species) %>% summarise(TEQ_CR_0 = sum(CR_TEQ == 0, na.rm = TRUE),TEQ_CR_not_0 = sum(CR_TEQ != 0, na.rm = TRUE))

# Make new variable "CR_TEQ_above_zero" (Y or N) to make barplots.
risk_data_merged_TEQCRzeros <- risk_data_merged %>% mutate(CR_TEQ_above_zero = ifelse(is.na(CR_TEQ), NA,ifelse(CR_TEQ > 0, "Y", "N")))

# Make faceted barplots.
# Calculate counts for each category and species.
count_data <- risk_data_merged_TEQCRzeros %>% count(species, CR_TEQ_above_zero, name = "count")

# Create a data frame with all possible combinations of species and CR_TEQ_above_zero.
all_combinations <- expand.grid(species = unique(risk_data_merged_TEQCRzeros$species),CR_TEQ_above_zero = c("N", "Y"))

# Merge with count_data to include counts for all combinations.
count_data_complete <- merge(all_combinations, count_data, by = c("species", "CR_TEQ_above_zero"), all.x = TRUE)

# Replace NA counts with 0.
count_data_complete[is.na(count_data_complete$count), "count"] <- 0

# Filter out rows with count > 0.
count_data_positive <- count_data_complete %>% filter(count > 0)

# Plot bars.
TEQ_CR_bar <- count_data_positive %>%
  ggplot(aes(x = CR_TEQ_above_zero, y = count)) +
  geom_bar(stat = "identity", width = 0.7, position = "dodge", color = "gray") +
  ylim(0, 40) +
  facet_wrap(~species, ncol = 1, strip.position = "left", labeller = as_labeller(label_species_CRTEQ)) +  
  scale_x_discrete(labels = c("N" = "TCDD-TEQ CR = 0", "Y" = "TCDD-TEQ CR > 0")) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  labs(x = NULL, y = "Sample Count")

# Add text labels for all counts.
TEQ_CR_bar <- TEQ_CR_bar +
  geom_text(data = count_data_complete, aes(label = ifelse(count == 0, "0", count)), vjust = -0.5, position = position_dodge(width = 0.7))

# Combine plots together.
TEQ_CR_combined <- plot_grid(TEQ_CR_bar, TEQ_CR_histo,labels = c("a)", "b)"),
                             ncol = 2,rel_widths = c(1, 2)  # Adjust widths of the plots
)
plot(TEQ_CR_combined)

# Save the plot. 
ggsave("TCDD_CR_by_species.png", plot = TEQ_CR_combined, dpi = 1000, width = 10.5, height = 5.83)

# ---------------------------------------------------------------------------
# PCB CR by Species Histogram -----------------------------------------------

# Calculate the number of non-NA values for each species.
species_counts_CRTPCB <- risk_data_merged %>% group_by(species) %>% summarize(count = sum(!is.na(CR_TPCB)))

# Create custom labeller function for facet labels.
label_species_CRTPCB <- function(variable) {
  count <- species_counts_CRTPCB$count[species_counts_CRTPCB$species == variable]
  return(paste(variable, " (n=", count, ")", sep = ""))
} # End create custom labeler. 

# Create the faceted histogram plot.
PCB_CR_histo <- risk_data_merged %>% 
  filter(!is.na(CR_TPCB)) %>% 
  ggplot() +
  geom_histogram(bins = 40,mapping = aes(x = CR_TPCB)) +
  facet_wrap(~ species, ncol = 1, strip.position = "left", labeller = as_labeller(label_species_CRTPCB)) +
  scale_x_log10(sec.axis = dup_axis(name = NULL), breaks = c(1e-8, 1e-7, 1e-6, 1e-5, 1e-4, 1e-3))+
  labs(x = expression("\u03A3PCB Cancer Risk"), y = "Sample Count") +
  theme_bw()+
  geom_vline(xintercept = 1e-6, linetype = "dashed", color = "orange") +
  geom_vline(xintercept = 1e-4, linetype = "solid", color = "red") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text = element_text(hjust = 0.5),
        axis.title.x = element_text(margin = margin(t = 1, unit = "pt")),
        axis.title.y = element_text(margin = margin(r = 2, unit = "pt"))) +
  theme(strip.text = element_text(size = 8))

# Save the plot.
ggsave("PCB_CR_by_species.png", plot = PCB_CR_histo, dpi = 1000, width = 7.8, height = 5.83)

# ---------------------------------------------------------------------------
# OCPs CR by Species Histogram ----------------------------------------------

# Calculate the number of non-NA values for each species.
species_counts_CROCP <- risk_data_merged %>% group_by(species) %>% summarize(count = sum(!is.na(Sum_OCPs_CR)))

# Create custom labeller function for facet labels.
label_species_CROCP <- function(variable) {
  count <- species_counts_CROCP$count[species_counts_CROCP$species == variable]
  return(paste(variable, " (n=", count, ")", sep = ""))
} # End create custom labeler. 

# Create the plot.
OCP_CR_histo <- risk_data_merged %>% 
  filter(!is.na(Sum_OCPs_HQ)) %>% 
  ggplot() +
  geom_histogram(mapping = aes(x = Sum_OCPs_CR)) +
  facet_wrap(~ species, ncol = 1, strip.position = "left", labeller = as_labeller(label_species_CROCP)) +
  scale_x_log10(sec.axis = dup_axis(name = NULL), breaks = c(1e-8, 1e-7, 1e-6, 1e-5, 1e-4, 1e-3))+
  labs(x = "OCPs Cancer Risk", y = "Sample Count") +
  theme_bw()+
  geom_vline(xintercept = 1e-6, linetype = "dashed", color = "orange") +
  geom_vline(xintercept = 1e-4, linetype = "solid", color = "red") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text = element_text(hjust = 0.5),
        axis.title.x = element_text(margin = margin(t = 1, unit = "pt")),
        axis.title.y = element_text(margin = margin(r = 2, unit = "pt"))) +
  theme(strip.text = element_text(size = 8))

# Save the plot. 
ggsave("OCP_CR_by_species.png", plot = OCP_CR_histo, dpi = 1000, width = 7.8, height = 5.83)

# ---------------------------------------------------------------------------
# PFAS CR by Species Histogram ----------------------------------------------

# Calculate the number of non-NA values for each species.
species_counts_CRPFAS <- risk_data_merged %>% group_by(species) %>% summarize(count = sum(!is.na(Sum_PFAS_CR)))

# Create custom labeler function for facet labels.
label_species_CRPFAS <- function(variable) {
  count <- species_counts_CRPFAS$count[species_counts_CRPFAS$species == variable]
  return(paste(variable, " (n=", count, ")", sep = ""))
} # End custom labeler. 

# Create the faceted histogram plot.
PFAS_CR_histo <- risk_data_merged %>% 
  filter(!is.na(Sum_PFAS_CR)) %>% 
  ggplot() +
  geom_histogram(mapping = aes(x = Sum_PFAS_CR)) +
  facet_wrap(~ species, ncol = 1, strip.position = "left", labeller = as_labeller(label_species_CRPFAS)) +
  scale_x_log10(sec.axis = dup_axis(name = NULL),breaks = c(1e-8, 1e-7, 1e-6, 1e-5, 1e-4, 1e-3))+
  labs(x = "PFAS Cancer Risk", y = "Sample Count") +
  theme_bw()+
  geom_vline(xintercept = 1e-6, linetype = "dashed", color = "orange") +
  geom_vline(xintercept = 1e-4, linetype = "solid", color = "red") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text = element_text(hjust = 0.5),
        axis.title.x = element_text(margin = margin(t = 1, unit = "pt")),
        axis.title.y = element_text(margin = margin(r = 2, unit = "pt"))) +
  theme(strip.text = element_text(size = 8))

# Make barplots of zeros and non-zeros as logscale histogram excludes zeros.
# Get number of samples per species with  Sum_PFAS_Cr = 0 and not = 0.
risk_data_merged %>% group_by(species) %>% summarise(PFAS_CR_0 = sum(Sum_PFAS_CR == 0, na.rm = TRUE), PFAS_CR_not_0 = sum(Sum_PFAS_CR != 0, na.rm = TRUE))

# Make new variable "CR_PFAS_above_zero" (Y or N) to make barplots.
risk_data_merged_PFASCRzeros <- risk_data_merged %>% mutate(CR_PFAS_above_zero = ifelse(is.na(Sum_PFAS_CR), NA, ifelse(Sum_PFAS_CR > 0, "Y", "N")))

# Make faceted barplots.
# Calculate counts for each category and species.
count_data_PFAS <- risk_data_merged_PFASCRzeros %>% count(species, CR_PFAS_above_zero, name = "count")

# Create a data frame with all possible combinations of species and CR_TEQ_above_zero.
all_combinations_PFAS <- expand.grid(species = unique(risk_data_merged_PFASCRzeros$species),CR_PFAS_above_zero = c("N", "Y"))

# Merge with count_data to include counts for all combinations.
count_data_complete_PFAS <- merge(all_combinations_PFAS, count_data_PFAS, by = c("species", "CR_PFAS_above_zero"), all.x = TRUE)

# Replace NA counts with 0.
count_data_complete_PFAS[is.na(count_data_complete_PFAS$count), "count"] <- 0

# Filter out rows with count > 0.
count_data_positive_PFAS <- count_data_complete_PFAS %>%
  filter(count > 0)

# Plot bars.
PFAS_CR_bar <- count_data_positive_PFAS %>%
  ggplot(aes(x = CR_PFAS_above_zero, y = count)) +
  geom_bar(stat = "identity", width = 0.7, position = "dodge", color = "gray") +
  ylim(0, 40) +
  facet_wrap(~species, ncol = 1, strip.position = "left", labeller = as_labeller(label_species_CRPFAS)) +  
  scale_x_discrete(labels = c("N" = "PFAS CR = 0", "Y" = "PFAS CR > 0")) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  labs(x = NULL, y = "Sample Count")

# Add text labels for all counts.
PFAS_CR_bar <- PFAS_CR_bar + geom_text(data = count_data_complete_PFAS, aes(label = ifelse(count == 0, "0", count)), vjust = -0.5, position = position_dodge(width = 0.7))


# Combine plots.
PFAS_CR_combined <- plot_grid(
  PFAS_CR_bar, PFAS_CR_histo,
  labels = c("a)", "b)"),ncol = 2,rel_widths = c(1, 2)  # Adjust widths of the plots
)
plot(PFAS_CR_combined)

# Save plot. 
ggsave("PFAS_CR_by_species.png", plot = PFAS_CR_combined, dpi = 1000, width = 10.5, height = 5.83)

# ---------------------------------------------------------------------------
# Cancer Deterministic Summary Stats ----------------------------------------

# List of variables.
variables_of_interest_CR <- c("CR_TEQ", "CR_TPCB", "Sum_OCPs_CR", "Sum_PFAS_CR", "Total_CR") 

# Create a new dataframe to store the statistics.
cancer_deterministic_summary <- data.frame(
  min = numeric(length(variables_of_interest_CR)),
  max = numeric(length(variables_of_interest_CR)),
  mean = numeric(length(variables_of_interest_CR)),
  median = numeric(length(variables_of_interest_CR)),
  confidence_interval = character(length(variables_of_interest_CR)),
  percentile_95 = numeric(length(variables_of_interest_CR))
)

# Calculate statistics for each variable.
for (i in seq_along(variables_of_interest_CR)) {
  variable <- variables_of_interest_CR[i]
  # Filter out missing values or NaNs.
  non_missing_values <- na.omit(risk_data[[variable]])
  # Calculate basic statistics.
  stats <- summary(non_missing_values)
  # Create a linear model with an intercept term.
  l_model <- lm(non_missing_values ~ 1)
  # Calculate the confidence interval.
  confidence_interval <- confint(l_model, level = 0.95)
  # Extract the lower and upper bounds and format them.
  lower_bound <- format(as.numeric(confidence_interval[1]), scientific = TRUE)
  upper_bound <- format(as.numeric(confidence_interval[2]), scientific = TRUE)
  # Calculate 95th percentile.
  percentile_95 <- format(quantile(non_missing_values, 0.95), scientific = TRUE)
  # Update the statistics dataframe.
  cancer_deterministic_summary[i, ] <- c(
    format(stats[1], scientific = TRUE), 
    format(stats[6], scientific = TRUE), 
    format(mean(non_missing_values), scientific = TRUE), 
    format(median(non_missing_values), scientific = TRUE), 
    paste(lower_bound, "-", upper_bound), 
    percentile_95)
} # End statistics computation loop. 

# Set row names to variable names.
row.names(cancer_deterministic_summary) <- variables_of_interest_CR

# Display the resulting statistics dataframe.
view(cancer_deterministic_summary)

# Save the file. 
write_csv(cancer_deterministic_summary, "cancer_deterministic_summary.csv")




