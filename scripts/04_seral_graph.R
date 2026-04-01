library(tidyverse)
library(ggplot2)
library(viridis)

# ==========================================
# 1. SETUP DIRECTORIES & LOAD DATA
# ==========================================

# Create the nested directory structure if it doesn't exist yet
if(!dir.exists("Figures/Seral_Figures")) {
  dir.create("Figures/Seral_Figures", recursive = TRUE)
}

# Load the fitted parameters from the GLM step
fits <- read_csv("data/seral_model_fits.csv")

# ==========================================
# 2. DEFINE THE PARAMETER GRID
# ==========================================

# The specific rate parameters (b) you want to generate individual graphs for
b_vals <- c(0.005, 0.01, 0.5, 0.9, 1.0, 3.0)

# ==========================================
# 3. CALCULATE & PLOT IN A LOOP
# ==========================================

for (b_v in b_vals) {
  
  grid <- fits %>% mutate(b = b_v)
  
  # Step B: Calculate the indices
  results <- grid %>%
    mutate(
      # Expected abundance in the log domain (Equation 3.8)
      log_d_i = log(alpha) + (1/k)*log(b) - (beta/k)*log(b + gamma) + 
        lgamma(beta/k) - lgamma(1/k)
    ) %>%
    # Group strictly by k to aggregate the 11 species into a single community metric
    group_by(k) %>% 
    summarise(
      log_G = mean(log_d_i), 
      .groups = "drop"
    ) %>%
    mutate(
      # Re-attach the current b value and calculate Pyrodiversity (P)
      b = b_v, 
      pyro_index = (1/k) * (1 - log(b_v)) + lgamma((1/k) + 1)
    )
  
  # Step C: Generate the ggplot for this specific b value
  p <- ggplot(results, aes(x = pyro_index, y = log_G, color = k)) +
    geom_path(linewidth = 1.2, arrow = arrow(length = unit(0.3, "cm"), type = "closed")) +
    scale_color_viridis_c(name = "Shape (k)", option = "viridis") +
    labs(
      x = "Pyrodiversity Index (P)",
      y = "log(Geometric Species Index) [ln(G)]",
      title = paste("Trajectory Plot for Empirical Seral Community (b =", b_v, ")")
    ) +
    theme_minimal()
  
  # Step D: Save the plot with a LaTeX-safe file name
  b_str <- gsub("\\.", "_", as.character(b_v))
  file_name <- paste0("Figures/Seral_Figures/Pyro_Bio_Seral_b_", b_str, ".png")
  
  ggsave(filename = file_name, plot = p, width = 7, height = 5, dpi = 300, bg = "white")
  cat("Successfully saved:", file_name, "\n")
}