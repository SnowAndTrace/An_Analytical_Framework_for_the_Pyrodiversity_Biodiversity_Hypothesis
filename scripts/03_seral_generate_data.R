library(tidyverse)
library(readxl)
# 1. Load the actual dataset
data <- read_excel("data/data.xlsx")

# Create a clean log_tsf variable. 
# (Your data already has a 'log(tsf)' column, but naming it cleanly 
# prevents syntax errors in the glm formula).
data <- data %>% mutate(log_tsf = log(tsf))

# Extract the 11 species names (starting from column 3)
species_names <- names(data)[3:ncol(data)]

# 2. Define the fitting function
fit_seral_species <- function(species, k_val, df) {
  
  # Create the t^k term dynamically based on the current loop's k value
  df$t_k <- df$tsf^k_val
  
  # Try fitting the Seral Model
  seral_model <- glm(df[[species]] ~ log_tsf + t_k, 
                     family = binomial(link = "cloglog"), 
                     data = df)
  
  coefs <- coef(seral_model)
  alpha_est <- exp(coefs[1])     # Intercept corresponds to alpha
  beta_est <- coefs[2] + 1       # log_tsf coefficient corresponds to beta
  gamma_est <- -coefs[3]         # t_k coefficient corresponds to gamma
  
  # 3.Gamma must be positive for a valid seral "hump-shaped" model.
  if (!is.na(gamma_est) && gamma_est > 0) {
    model_type <- "Seral"
  } else {
    # If gamma is negative or zero, fallback to the standard Power Function Model
    power_model <- glm(df[[species]] ~ log_tsf, 
                       family = binomial(link = "cloglog"), 
                       data = df)
    
    coefs_power <- coef(power_model)
    alpha_est <- exp(coefs_power[1])
    beta_est <- coefs_power[2] + 1
    gamma_est <- 0  # Gamma is 0 in the power function
    model_type <- "Power (Fallback)"
  }
  
  # Return the row of results
  return(data.frame(
    species = species,
    k = k_val,
    alpha = alpha_est,
    beta = beta_est,
    gamma = gamma_est,
    model_used = model_type,
    stringsAsFactors = FALSE
  ))
}

# 4. Execute the loop across the grid of k values
# (e.g., 0.5 to 2.0 as discussed for realistic fire regimes)
k_values <- seq(0.5, 2.0, by = 0.1) 
all_results <- list()

for (k in k_values) {
  for (sp in species_names) {
    tryCatch({
      res <- fit_seral_species(sp, k, data)
      all_results[[length(all_results) + 1]] <- res
    }, error = function(e) {
      cat("Error fitting", sp, "at k =", k, "\n")
    })
  }
}

# 5. Combine and save the results
final_model_fits <- bind_rows(all_results)
write_csv(final_model_fits, "data/seral_model_fits.csv")
cat("Successfully saved to: data/seral_model_fits.csv\n")

# View the top few rows
head(final_model_fits)