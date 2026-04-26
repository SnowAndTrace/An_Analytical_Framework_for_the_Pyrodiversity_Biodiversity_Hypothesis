library(tidyverse)
library(ggplot2)
library(viridis)

#Figure 3.5 code

# ==========================================
# 1. SETUP DIRECTORIES & LOAD DATA
# ==========================================

if(!dir.exists("Figures/Seral_Figures")) {
  dir.create("Figures/Seral_Figures", recursive = TRUE)
}

fits <- read_csv("data/seral_model_fits.csv")

# ==========================================
# 2. DEFINE THE PARAMETER GRID
# ==========================================

k_vals_to_plot <- c(0.5, 1.0, 1.5, 2.0) 
b_seq <- seq(0.01, 3.0, length.out = 300)

# ==========================================
# 3. CALCULATE & PLOT IN A LOOP
# ==========================================

for (k_v in k_vals_to_plot) {
  
  k_fits <- fits %>% filter(abs(k - k_v) < 1e-5)
  
  if(nrow(k_fits) == 0) {
    cat("No data found for k =", k_v, "- skipping.\n")
    next
  }
  
  grid <- crossing(k_fits, b = b_seq)
  
  results <- grid %>%
    mutate(
      log_d_i = log(alpha) + (1/k)*log(b) - (beta/k)*log(b + gamma) + 
        lgamma(beta/k) - lgamma(1/k)
    ) %>%
    group_by(b, k) %>% 
    summarise(
      log_G = mean(log_d_i), 
      .groups = "drop"
    ) %>%
    mutate(
      pyro_index = (1/k) * (1 - log(b)) + lgamma((1/k) + 1)
    )
  
  # Step C: Calculate the Pearson correlation (r-value) AND Linear Gradient
  r_val <- cor(results$pyro_index, results$log_G, use = "complete.obs")
  
  # Fit a linear model to extract the slope (gradient) of the best-fit line
  lm_fit <- lm(log_G ~ pyro_index, data = results)
  grad_val <- coef(lm_fit)["pyro_index"] # Extracts the gradient/slope
  
  # Step D: Generate the ggplot
  p <- ggplot(results, aes(x = pyro_index, y = log_G, color = b)) +
    geom_path(linewidth = 1.2, arrow = arrow(length = unit(0.3, "cm"), type = "closed")) +
    # Optional: uncomment the line below if you want to visually draw the linear trendline on the graph
    # geom_smooth(method = "lm", color = "black", linetype = "dashed", se = FALSE, linewidth = 0.5) +
    scale_color_viridis_c(name = "Rate (b)", option = "plasma", direction = -1) +
    labs(
      x = "Pyrodiversity Index (P)",
      y = "log(Geometric Species Index) [L]",
      title = paste("Trajectory Plot for Empirical Seral Community (Fixed k =", k_v, ")"),
      # Subtitle now dynamically prints both values
      subtitle = sprintf("Pearson r = %.4f", r_val) 
    ) +
    theme_minimal()
  
  # Step E: Save the plot
  k_str <- gsub("\\.", "_", as.character(k_v))
  file_name <- paste0("Figures/Seral_Figures/Pyro_Bio_Seral_k_", k_str, ".png")
  
  ggsave(filename = file_name, plot = p, width = 7, height = 5, dpi = 300, bg = "white")
  cat(sprintf("Saved: %s (r = %.4f)\n", file_name, r_val))
}