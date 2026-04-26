# ==============================================================================
# Fit the parametric GLM (Approach 1) using a complementary log-log link, 
# extract the alpha and beta parameters for all species, and save the 
# results to generate Table 3.1.
# ==============================================================================

library(tidyverse)
library(lme4)
library(readxl)

# 1. Load the empirical data
data <- read_excel("data/data.xlsx")

# 2. Extract the 11 species names
species_names <- names(data)[3:length(names(data))]

# 3. Initialise an empty dataframe to store the parameter results
results <- data.frame(
  species = character(),  
  alpha = numeric(),      
  beta = numeric(),       
  stringsAsFactors = FALSE
)

# 4. Loop through each species to fit the GLM and extract parameters
for (species in species_names){
  
  # Define the model formula
  formula <- as.formula(paste(species, "~ log(tsf)"))
  
  # Fit the GLM with a cloglog link function
  model1 <- glm(formula, family = binomial(link = "cloglog"), data = data)
  
  # Extract coefficients
  B0 <- summary(model1)$coefficients[1] # Intercept
  B1 <- summary(model1)$coefficients[2] # Slope
  
  results <- rbind(results, data.frame(
    species = species,
    alpha = exp(B0),    # alpha = exp(B0)
    beta = B1 + 1       # beta = B1 + 1 
  ))
}


# 5. Save the output to be used by all subsequent scripts
# write.csv(results, "species_power_dist.csv", row.names = FALSE)

# 6. Display the final table (Table 3.1)
print(results)
# cat("Successfully saved parameters to: species_power_dist.csv\n")