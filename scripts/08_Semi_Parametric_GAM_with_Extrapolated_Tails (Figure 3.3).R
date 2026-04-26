# ==============================================================================
# Generate Figure 3.3 - Semi-Parametric GAM with extrapolated tails
# ==============================================================================

library(tidyverse)
library(mgcv)
library(readxl)
library(ggplot2)

# 1. Setup the Figures directory
if(!dir.exists("Figures")) {
  dir.create("Figures", recursive = TRUE)
}

# 2. Load the empirical reptile data and power distribution data
data <- read_excel("data/data.xlsx")
species_names <- names(data)[3:ncol(data)]
data_power <- read_csv("data/species_power_dist.csv", show_col_types = FALSE)

# 3. Define and fit the GAM models
gam_fun <- function(species){
  gam(data[[species]] ~ s(tsf), family = binomial(link = "cloglog"), method = "REML", data = data)
}

gam_models <- lapply(species_names, gam_fun)
names(gam_models) <- species_names

# 4. Calculate breaking points and dilation factors for the tail
breaking_point_gam <- data.frame(species = species_names, value = NA_real_, stringsAsFactors = FALSE)

for (species in species_names){
  prediction <- -log(1 - predict(gam_models[[species]], newdata = data.frame(tsf = 154), type = "response"))
  breaking_point_gam$value[breaking_point_gam$species == species] <- prediction
}

power_abun_tlast <- data.frame(species = species_names, stringsAsFactors = FALSE)
power_abun_tlast$value <- data_power$alpha * 154^data_power$beta

dilation_factors <- breaking_point_gam$value / power_abun_tlast$value
data_power$modified_alpha <- data_power$alpha * dilation_factors

# 5. Build prediction grids (GAM range and Tail range)
t_gam <- seq(min(data$tsf, na.rm = TRUE), max(data$tsf, na.rm = TRUE))   # 0 to 154
t_tail <- seq(max(data$tsf, na.rm = TRUE) + 1, 200, by = 1)              # 155 to 200

predict_gam_trans <- function(sp, tgrid) {
  preds <- predict(gam_models[[sp]], newdata = data.frame(tsf = tgrid), type = "response")
  preds[preds >= 1] <- 1 - 1e-12
  preds[preds <= 0] <- 1e-12
  -log(1 - preds)
}

# Build GAM dataframe
gam_list <- lapply(species_names, function(sp) {
  data.frame(tsf = t_gam, species = sp, value = predict_gam_trans(sp, t_gam), 
             segment = "GAM", stringsAsFactors = FALSE)
})
gam_df <- bind_rows(gam_list)

# Build Extrapolated Tail dataframe
tail_list <- lapply(species_names, function(sp) {
  alpha <- data_power$modified_alpha[data_power$species == sp]
  beta  <- data_power$beta[data_power$species == sp]
  vals <- alpha * (t_tail ^ beta)
  data.frame(tsf = t_tail, species = sp, value = vals, 
             segment = "tail", stringsAsFactors = FALSE)
})
tail_df <- bind_rows(tail_list)

# Combine and clean data
full_df <- bind_rows(gam_df, tail_df)
full_df$value[full_df$value <= 0] <- 1e-7
full_df$species <- factor(full_df$species, levels = species_names)

# 6. Generate and save Figure 3.3
p_facet <- ggplot(full_df, aes(x = tsf, y = value)) +
  geom_line(aes(linetype = segment), linewidth = 0.8) +
  facet_wrap(~species, scales = "fixed", ncol = 3) +
  labs(title = "Per-species abundance with extrapolation",
       x = "tsf", y = "Abundance (transformed)") +
  theme_minimal() +
  scale_linetype_manual(values = c(GAM = "solid", tail = "dashed"))

ggsave("Figures/Figure_3_3_GAM_Abundance.png", plot = p_facet, width = 10, height = 8, dpi = 300)
cat("Saved: Figures/Figure_3_3_GAM_Abundance.png\n")