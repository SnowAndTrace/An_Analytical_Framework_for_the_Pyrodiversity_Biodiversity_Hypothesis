# ==============================================================================
# Generate Figure 3.2 - Comparison between GLM and GAM
# ==============================================================================

rm(list = ls())
library(tidyverse)
library(mgcv)
library(readxl)
library(ggplot2)

# --- 1. Set up and Data Loading ---
data <- read_excel("data/data.xlsx")
species_names = names(data[3:length(names(data))])

# fit the GAM model
gam_fun = function(species){
  gam_model = gam(data[[species]] ~ s(tsf), family = binomial(link = "cloglog"), method = "REML", data = data)
  return(gam_model)
}

# Pre-fit GAM models once
gam_models <- lapply(species_names, gam_fun)
names(gam_models) <- species_names

# --- 2. Calculate Dilation Factors ---
breaking_point_gam <- data.frame(species = species_names, value = NA_real_, stringsAsFactors = FALSE)
for (species in species_names){
  prediction = -log(1 - predict(gam_fun(species), newdata = data.frame(tsf = 154), type = "response"))
  breaking_point_gam$value[breaking_point_gam$species == species] = prediction
}

data_power = read_csv("data/species_power_dist.csv", show_col_types = FALSE)
species_list = data_power$species
alpha_list = data_power$alpha
beta_list = data_power$beta
n = length(data_power$species)

power_abun_tlast <- data.frame(species = species_names, stringsAsFactors = FALSE)
power_abun_tlast$value = alpha_list * 154^beta_list

dilation_factors = breaking_point_gam$value/power_abun_tlast$value
data_power$modified_alpha = data_power$alpha * dilation_factors

# --- 3. Integration Functions ---
trapz <- function(x, y) {
  sum((y[-1] + y[-length(y)])/2 * diff(x))
}

gam_area_trap <- function(b, k, species_names, gam_models, n = 1000) {
  t_grid <- seq(0, 154, length.out = n)
  g_vals <- b^(1/k) * exp(-b * t_grid^k) / gamma(1/k + 1)
  areas <- sapply(species_names, function(sp) {
    preds <- predict(gam_models[[sp]], newdata = data.frame(tsf = t_grid), type = "response")
    resp_vals <- -log(1 - preds)   
    trapz(t_grid, g_vals * resp_vals)
  })
  return(areas)
}

glm_area_trap <- function(b, k, species_names, data_power, t_max = 1000, n = 1000) {
  t_grid <- seq(154, t_max, length.out = n)
  g_vals <- b^(1/k) * exp(-b * t_grid^k) / gamma(1/k + 1)
  areas <- sapply(species_names, function(sp) {
    alpha <- data_power$modified_alpha[data_power$species == sp]
    beta  <- data_power$beta[data_power$species == sp]
    power_vals <- alpha * t_grid^beta
    trapz(t_grid, g_vals * power_vals)
  })
  return(areas)
}

bio_index_trap <- function(b, k, species_names, gam_models, data_power, n = 1000, t_max = 1000) {
  gam_areas <- gam_area_trap(b, k, species_names, gam_models, n)
  glm_areas <- glm_area_trap(b, k, species_names, data_power, t_max = 1000, n = 1000)
  total_area <- gam_areas + glm_areas
  G_log <- sum(log(total_area)) / length(total_area)
  return(G_log)  
}

bio_index_thm = function(b, k){
  d = alpha_list * b^(-beta_list/k) * gamma((beta_list + 1)/k) / gamma(1/k)
  power = sum(log(d))/n
  G = exp(power)
  return (log(G))
}
bio_index_glm <- bio_index_thm

# --- 4. Main Grid Evaluation ---
b_list_grid <- seq(0.2, 3, length.out = 50)
k_list_grid <- seq(0.6, 2, length.out = 30)
grid <- expand.grid(b = b_list_grid, k = k_list_grid)

grid$G_glm <- mapply(bio_index_glm, grid$b, grid$k)
grid$G_gam <- mapply(function(b, k) {
  bio_index_trap(b, k, species_names, gam_models, data_power, n = 1000, t_max = 1000)
}, grid$b, grid$k)

grid$abs_diff <- abs(grid$G_gam - grid$G_glm)

# --- PANEL A: Scatterplot ---
fit <- lm(G_gam ~ G_glm, data = grid)
coef_fit <- coef(fit)
slope <- round(coef_fit[2], 3)
intercept <- round(coef_fit[1], 3)
r2 <- round(summary(fit)$r.squared, 3)
corr <- round(cor(grid$G_glm, grid$G_gam), 3)

stat_label <- paste0("R² = ", r2, ", r = ", corr)

pA <- ggplot(grid, aes(x = G_glm, y = G_gam)) +
  geom_point(alpha = 0.5) +
  geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
  annotate("text", x = min(grid$G_glm), y = max(grid$G_gam),
           label = stat_label, hjust = 0, vjust = 1, size = 5, color = "blue") +
  labs(title = "log(G) index, GAM vs GLM with linear fit",
       x = "GLM log(G) index",
       y = "GAM log(G) index")
print(pA)

# --- PANEL B: Heatmap ---
pB <- ggplot(grid, aes(x = b, y = k, fill = abs_diff)) +
  geom_tile() +
  scale_fill_gradient(low = "yellow", high = "red") +
  labs(title = "Absolute Difference log(G) across (b,k)",
       x = "b", y = "k", fill = "|log(G_GAM) - log(G_GLM)|") +
  theme_minimal()
print(pB)

# --- PANEL C: Fixed b ---
b_list_fix <- c(0.2, 0.5, 1, 2, 3)
k_list_fix <- seq(0.6, 2, length.out = 30)

grid_fix_b <- expand.grid(b = b_list_fix, k = k_list_fix)
grid_fix_b$G_glm <- mapply(bio_index_glm, grid_fix_b$b, grid_fix_b$k)
grid_fix_b$G_gam <- mapply(function(b, k) {
  bio_index_trap(b, k, species_names, gam_models, data_power, n = 1000, t_max = 1000)
}, grid_fix_b$b, grid_fix_b$k)

pC <- ggplot(grid_fix_b, aes(x = G_glm, y = G_gam, color = factor(b), group = b)) +
  geom_point() +
  geom_line() +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  labs(
    title = "Comparison of log(G) index: GAM vs GLM (fixed b)",
    x = "GLM log(G)",
    y = "GAM log(G)",
    color = "b"
  ) +
  theme_minimal(base_size = 14)
print(pC)

# --- PANEL D: Fixed k ---
k_list_fix2 <- c(0.6, 0.7, 1, 1.5, 2)
b_list_fix2 <- seq(0.2, 3, length.out = 50)

grid_fix_k <- expand.grid(b = b_list_fix2, k = k_list_fix2)
grid_fix_k$G_glm <- mapply(bio_index_glm, grid_fix_k$b, grid_fix_k$k)
grid_fix_k$G_gam <- mapply(function(b, k) {
  bio_index_trap(b, k, species_names, gam_models, data_power, n = 1000, t_max = 1000)
}, grid_fix_k$b, grid_fix_k$k)

pD <- ggplot(grid_fix_k, aes(x = G_glm, y = G_gam, color = factor(k), group = k)) +
  geom_point() +
  geom_line() +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  labs(
    title = "Comparison of log(G) index: GAM vs GLM (fixed k)",
    x = "GLM log(G)",
    y = "GAM log(G)",
    color = "k"
  ) +
  theme_minimal(base_size = 14)
print(pD)

# 1. Ensure directory exists
if(!dir.exists("Figures")) {
  dir.create("Figures", recursive = TRUE)
}

# 2. Put your plots into a named list
# 1. Ensure directory exists
if(!dir.exists("Figures")) {
  dir.create("Figures", recursive = TRUE)
}

# 2. Put your plots into a named list matching the LaTeX filenames
my_plots <- list(
  "Compare_GLM_GAM_Linear_NoLabel" = pA,
  "Compare_GLM_GAM_HeatMap" = pB,
  "Compare_GLM_GAM_Linear_fixed_b" = pC,
  "Compare_GLM_GAM_Linear_fixed_k" = pD
)

# 3. Loop through and save them
for (plot_name in names(my_plots)) {
  file_path <- paste0("Figures/", plot_name, ".png")
  
  ggsave(file_path, plot = my_plots[[plot_name]], width = 9, height = 6, dpi = 300)
  cat("Saved:", file_path, "\n")
}