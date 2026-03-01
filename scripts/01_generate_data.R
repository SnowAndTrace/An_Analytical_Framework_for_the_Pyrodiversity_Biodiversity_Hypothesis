# scripts/01_generate_data.R
library(tidyverse)
source("src/trajectory_functions.R")

# 1. Load raw data
data <- read_csv("data/species_power_dist.csv") #
alpha_list <- data$alpha #

# 2. Define Parameters for the grid
k_vals <- seq(0.3, 3, length.out = 10000) #
b_fixed_vals <- c(0.005, 0.01, 0.5, 0.9, 1, 3) 

# Create an empty dataframe to store results
all_trajectories <- data.frame()

# 3. Loop through both scenarios
scenarios <- c("positive", "negative")

for (scenario in scenarios) {
  
  # Set the modified beta list based on the scenario
  if (scenario == "positive") {
    beta_list_mod <- data$beta + 1 #
    is_pos <- TRUE
  } else {
    beta_list_mod <- -data$beta 
    is_pos <- FALSE
  }
  
  # Loop through each fixed 'b' value
  for (b_val in b_fixed_vals) {
    
    df_temp <- expand.grid(b = b_val, k = k_vals)
    df_temp$scenario <- scenario
    
    # Calculate indices
    df_temp$pyro <- mapply(pyro_index, df_temp$b, df_temp$k) #
    df_temp$logBio <- mapply(function(b, k) {
      bio_index_thm(b, k, alpha_list, beta_list_mod, is_positive_scenario = is_pos)
    }, df_temp$b, df_temp$k)
    
    all_trajectories <- bind_rows(all_trajectories, df_temp)
  }
}

# 4. Save the computed data
saveRDS(all_trajectories, "output/trajectory_data.rds")
cat("Data generation complete. Saved to output/trajectory_data.rds\n")