# ==============================================================================
# Compare the flexible GAM approach against the Seral parametric model.
# ==============================================================================

library(tidyverse)
library(mgcv)
library(readxl)
library(ggplot2)
library(viridis)

# --- 1. Load Data & Fits ---
data <- read_excel("data/data.xlsx")
species_names <- names(data)[3:ncol(data)]
data_power <- read_csv("data/species_power_dist.csv", show_col_types = FALSE)
seral_fits <- read_csv("data/seral_model_fits.csv", show_col_types = FALSE)

# --- 2. Fit GAM Models & Calculate Tail Dilation (Approach 2) ---
# This ensures the GAM has the correct extrapolated tails for comparison
gam_models <- lapply(species_names, function(sp) {
  gam(data[[sp]] ~ s(tsf), family = binomial(link = "cloglog"), method = "REML", data = data)
})
names(gam_models) <- species_names

breaking_point_gam <- sapply(species_names, function(sp) {
  -log(1 - predict(gam_models[[sp]], newdata = data.frame(tsf = 154), type = "response"))
})

power_abun_tlast <- data_power$alpha * 154^data_power$beta
data_power$modified_alpha <- data_power$alpha * (breaking_point_gam / power_abun_tlast)

# --- 3. Integration Functions ---
trapz <- function(x, y) sum((y[-1] + y[-length(y)])/2 * diff(x))

# Function to calculate ln(G) for the GAM across b and k
calc_logG_gam <- function(b, k) {
  # Range [0, 154]
  t1 <- seq(0, 154, length.out = 500)
  g1 <- b^(1/k) * exp(-b * t1^k) / gamma(1/k + 1)
  areas1 <- sapply(species_names, function(sp) {
    preds <- predict(gam_models[[sp]], newdata = data.frame(tsf = t1), type = "response")
    trapz(t1, g1 * (-log(1 - pmax(pmin(preds, 1-1e-12), 1e-12))))
  })
  
  # Tail [154, 1000]
  t2 <- seq(154, 1000, length.out = 500)
  g2 <- b^(1/k) * exp(-b * t2^k) / gamma(1/k + 1)
  areas2 <- sapply(species_names, function(sp) {
    alpha <- data_power$modified_alpha[data_power$species == sp]
    beta  <- data_power$beta[data_power$species == sp]
    trapz(t2, g2 * (alpha * t2^beta))
  })
  
  mean(log(areas1 + areas2))
}

# --- 4. Evaluate Trajectories ---
k_vals_to_plot <- c(0.5, 1.0, 1.5, 2.0)
b_seq <- seq(0.01, 3.0, length.out = 100)

cat("Calculating trajectories...\n")
plot_results <- list()

for (k_v in k_vals_to_plot) {
  # Get Seral coefficients for this k
  k_fits <- seral_fits %>% filter(abs(k - k_v) < 1e-5)
  
  # Calculate per-point indices 
  temp_df <- data.frame(b = b_seq, k = k_v)
  # GLM value using the formula directly
  temp_df$logG_seral <- sapply(b_seq, function(b_v) {
    # Equation 3.8: Analytic solution for the Seral Model
    log_d_i <- log(k_fits$alpha) + (1/k_v)*log(b_v) - (k_fits$beta/k_v)*log(b_v + k_fits$gamma) + 
      lgamma(k_fits$beta/k_v) - lgamma(1/k_v)
    mean(log_d_i)
  })
  
  # GAM value using the trapezoidal rule approximation
  temp_df$logG_gam <- sapply(b_seq, function(b_v) calc_logG_gam(b_v, k_v))
  
  # Calculate R2 for this specific k trajectory
  r2 <- cor(temp_df$logG_seral, temp_df$logG_gam)^2
  temp_df$r2_label <- sprintf("k = %.1f", k_v)
  # temp_df$r2_label <- sprintf("k = %.1f (R² = %.4f)", k_v, r2)
  
  plot_results[[as.character(k_v)]] <- temp_df
}

final_df <- bind_rows(plot_results)

# --- 5. Generate and Save Plot ---
if(!dir.exists("Figures")) dir.create("Figures")

p <- ggplot(final_df, aes(x = logG_seral, y = logG_gam, color = b)) +
  geom_path(linewidth = 1.2, arrow = arrow(length = unit(0.3, "cm"), type = "closed")) +
  geom_abline(slope = 1, intercept = 0, color = "black", linetype = "dashed", alpha = 0.5) +
  facet_wrap(~r2_label, scales = "free") +
  scale_color_viridis_c(name = "Rate (b)", option = "plasma", direction = -1) +
  labs(
    x = "Seral Parametric (L)",
    y = "GAM Semi-Parametric (L)",
    title = "Comparison of Biodiversity Estimates across Fire Regimes",
    subtitle = "Trajectories show increasing fire frequency (increasing b)"
  ) +
  theme_minimal()

ggsave("Figures/Compare_Seral_GAM_Trajectories.png", p, width = 10, height = 8, dpi = 300, bg = "white")
cat("Saved: Figures/Compare_Seral_GAM_Trajectories.png\n")