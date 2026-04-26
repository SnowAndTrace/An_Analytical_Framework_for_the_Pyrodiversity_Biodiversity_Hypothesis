# ==============================================================================
# Generate Figure 3.4 - Parametric GLM abundance functions
# ==============================================================================

library(tidyverse)
library(ggplot2)

# 1. Setup the Figures directory
if(!dir.exists("Figures")) {
  dir.create("Figures", recursive = TRUE)
}

# 2. Load the calculated GLM parameters
data_power <- read_csv("data/species_power_dist.csv", show_col_types = FALSE)

species_list <- data_power$species
alpha_list <- data_power$alpha
beta_list <- data_power$beta

# 3. Generate theoretical predictions for the time sequence
t_list <- c(1:1000)

data_abun <- data.frame(tsf = t_list, stringsAsFactors = FALSE)

# Generate abundance values using your defined power function alpha * t^beta
for (i in 1:nrow(data_power)){
  predictions <- alpha_list[i] * (t_list ^ beta_list[i])
  data_abun[[species_list[i]]] <- predictions
}

# 4. Reshape data to long format for ggplot2
data_long <- data_abun %>%
  pivot_longer(
    cols = -tsf,
    names_to = "species",
    values_to = "Abundance"
  )

# Convert species to a factor to maintain consistent legend ordering
data_long$species <- factor(data_long$species, levels = species_list)

# 5. Generate and save Figure 3.4
p_glm <- ggplot(data_long, aes(x = tsf, y = Abundance, color = species)) +
  geom_line(linewidth = 1) +
  theme_minimal() +
  labs(
    title = "Abundance over Time",
    x = "Time Since Fire (tsf)",
    y = "Abundance"
  ) +
  theme(
    legend.position = "right",
    legend.title = element_text(face = "bold")
  )

ggsave("Figures/Figure_3_4_GLM_Abundance.png", plot = p_glm, width = 9, height = 6, dpi = 300)
cat("Saved: Figures/Figure_3_4_GLM_Abundance.png\n")