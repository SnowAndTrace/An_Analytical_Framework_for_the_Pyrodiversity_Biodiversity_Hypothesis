# src/trajectory_functions.R

# --- Pyro index ---
pyro_index <- function(b, k) {
  1/k * (1 - log(b)) + lgamma(1/k + 1) #
}

# --- Biodiversity index (GLM Theoretical) ---
bio_index_thm <- function(b, k, alpha_list, beta_list_mod, is_positive_scenario) {
  n <- length(alpha_list)
  
  if (is_positive_scenario) {
    # Formula for positive beta star
    d <- alpha_list * b^(-(beta_list_mod - 1)/k) * gamma(beta_list_mod / k) / gamma(1/k) #
  } else {
    # Formula for negative beta star
    d <- alpha_list * b^(-beta_list_mod/k) * gamma((beta_list_mod + 1)/k) / gamma(1/k)
  }
  
  power <- sum(log(d)) / n
  G <- exp(power)
  return(log(G))
}