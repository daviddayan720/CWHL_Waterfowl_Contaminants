# --------------------------------------------
# WFCS Multivariate Analysis of Contaminant Concentrations
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
# Description: This code includes the multivariate 
# statistical analyses and and visualizations.
# 
# Funding: This work was supported under the project 
# name 'Contaminant Loads in Waterfowl of the Northeast Atlantic
# Flyway: New Threats and Outdated Advisories' of the 
# Multistate Conservation Grant Program.
#
# --------------------------------------------
# Load required packages.
library(tidyverse) #version 2.0.0
library(vegan) #version 2.6.8
library(ConfidenceEllipse) #version 1.0.0
library(viridis) #version 0.6.5

#Load contaminant data
waterfowl_raw <- read_csv("Waterfowl Contaminant Study Final Processed Data 2021-2022.csv")

#set seed for reproducibility
set.seed(10706)

# Process the raw data
# Only include birds tested for all contaminants
# Select columns representing unique contaminant data
# Group the ABDU/MALL hybrid as ABDU as in previous analyses
# Remove columns representing contaminants that were never detected
waterfowl1 <- waterfowl_raw %>% 
  filter(if_all(all_of(c("Hg_test", "dioxin_test", "PCBOCP_test", "PFAS_test")), ~ . == "Y")) %>% 
  select(ID, species, 
         Hg, #mercury
         24, 26, 28, 29, 31, 33:36, 38:40, 42, 43, 45, 47, 48, #dioxin congeners
         50:209, #PCB congeners
         227:254, #OCPs
         255:294) %>% #PFAS
  mutate(species = if_else(species == "ABDU/MALL", "ABDU", species)) %>% 
  select(ID, species, where(~ any(. != 0)))

#relativize pollutant values by dividing by pollutant maximum
waterfowl_relativized <- waterfowl1 %>% 
  mutate(across(-c(ID, species), ~ . / max(.))) 

permanova1 <- adonis2(waterfowl_relativized[3:180] ~ species, 
                      data = waterfowl_relativized,
                      method = "bray",
                      permutations = 9999)


#Pairwise tests following PERMANOVA
#make the function as per "https://github.com/jon-bakker/appliedmultivariatestatistics/blob/main/pairwise.adonis2.R"
pairwise.adonis2 <- function(resp, fact, p.method = "none", nperm = 999) {
  require(vegan)
  resp <- as.matrix(resp)
  fact <- factor(fact)
  fun.p <- function(i, j) {
    fact2 <- droplevels(fact[as.numeric(fact) %in% c(i, j)])
    index <- which(fact %in% levels(fact2))
    resp2 <- as.dist(resp[index, index])
    result <- adonis2(resp2 ~ fact2, permutations = nperm)
    result$`Pr(>F)`[1]
  }
  multcomp <- pairwise.table(fun.p, levels(fact), p.adjust.method = p.method)
  return(list(fact = levels(fact), p.value = multcomp, p.adjust.method = p.method))
}

#make a matrix from the concentrations
relativized_matrix <- vegdist(waterfowl_relativized[3:180], method = "bray") #distance matrix on pollutant concentrations

#Run pairwise tests
pairwise.adonis2(resp = relativized_matrix, 
                 fact = waterfowl_relativized$species,
                 p.method = "none",
                 nperm = 9999)

permdisp1 <- betadisper(
  relativized_matrix,
  waterfowl_relativized$species,
  type = "median",
  bias.adjust = FALSE,
  sqrt.dist = FALSE,
  add = FALSE
)
#Visualize dispersions by species
boxplot(permdisp1)
#test differences in dispersions by species using PERMANOVA
adonis2(formula = dist(permdisp1$distances)~permdisp1$group)

#test differences in dispersions by species using ANOVA
aov(permdisp1$distances ~ permdisp1$group) %>% summary()

#principal coordinates analysis (PCoA) to identify # of dimensions for NMDS
PCoA1 <- wcmdscale(relativized_matrix, eig = TRUE)

#Non-metric MultiDimensional Scaling (NMDS)
NMDS1 <- metaMDS(comm = waterfowl_relativized[3:180],
                 autotransform = FALSE,
                 distance = "bray",
                 k = 9,
                 weakties = TRUE,
                 model = "global",
                 maxit = 1000,
                 try = 40,
                 trymax = 100)

#stress was 0.0635 

#join with species
NMDS1_species <- as.data.frame(NMDS1$points) %>% 
  mutate(species = waterfowl_relativized$species,
         species = as.factor(species))

ellipse_grp <- confidence_ellipse(NMDS1_species, x = MDS1, y = MDS2, .group_by = species, conf_level = .95)

#Visualize NMDS by species
NMDS_figure <- NMDS1_species %>% ggplot() +
  geom_point(data = NMDS1_species, mapping = aes(x = MDS1, y = MDS2, color = species), size = 3, alpha = 0.7) +  # Add points
  geom_polygon(data = ellipse_grp, aes(x = x, y = y, color = species),fill = NA, alpha = .45, linetype = "solid") +
  labs(
    x = "NMDS Axis 1 Score",
    y = "NMDS Axis 2 Score",
    color = "Species"
  ) +
  scale_color_viridis_d() + # Discrete color scale from viridis
  annotate("text", x = -1.5, y = -1.5, label = "Stress = 6.35", hjust = 0, size = 4) + # Add annotation
  theme_light() +
  theme_light() +
  theme(legend.position = "top",
        panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank())

#save the plot
ggsave("NMDS.png", plot = NMDS_figure, dpi = 1000, width = 7.8, height = 5.83)
