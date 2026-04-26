library(dplyr)

# Load fitted data
fits <- read.csv("data/seral_model_fits.csv")

# Constants
power_baseline <- 0.0505
n_species <- 11
k_vals_to_test <- c(0.5, 1.0, 1.5, 2.0)

cat("--- Finding the Crossover b Value ---\n\n")

for (k_v in k_vals_to_test) {
  # Filter for the specific k value
  k_fits <- fits %>% filter(abs(k - k_v) < 1e-5)
  
  if(nrow(k_fits) == 0) next
  
  # Equation 3.15
  # Define a function that calculates the difference between the seral 
  # gradient (at a given b) and the target 0.0505 baseline.
  # We want to find where this function equals 0. 
  gradient_difference <- function(b_val) {
    sum_term <- sum(k_fits$beta / (b_val + k_fits$gamma))
    local_gradient <- -1 + (b_val / n_species) * sum_term
    return(local_gradient - power_baseline)
  }
  
  # Use uniroot to computationally solve for b
  # We search within a wide ecological range for b (e.g., 0.0001 to 100)
  result <- tryCatch({
    uniroot(gradient_difference, interval = c(0.0001, 100))
  }, error = function(e) {
    return(NULL) # Failsafe if no crossover exists in that range
  })
  
  if(!is.null(result)) {
    crossover_b <- result$root
    cat(sprintf("For k = %.1f:\n", k_v))
    cat(sprintf("  The gradients perfectly match at b = %.4f\n\n", crossover_b))
  } else {
    cat(sprintf("For k = %.1f: No crossover point exists in the search interval.\n\n", k_v))
  }
}