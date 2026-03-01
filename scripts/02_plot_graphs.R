# scripts/02_plot_graphs.R
library(tidyverse)
library(ggplot2)
library(viridis)

# 1. Load the pre-computed data
df_all <- readRDS("output/trajectory_data.rds")

# Ensure the Figures directory exists
if(!dir.exists("Figures")) dir.create("Figures")

# 2. Get unique scenarios and b values
scenarios <- unique(df_all$scenario)
b_vals <- unique(df_all$b)

# 3. Generate and save plots dynamically
for (scen in scenarios) {
  for (b_v in b_vals) {
    
    # Filter data for this specific plot
    plot_data <- df_all %>% filter(scenario == scen, b == b_v)
    
    # Set the mathematical title based on the scenario
    beta_title <- ifelse(scen == "positive", "beta* > 0", "beta* < 0")
    full_title <- paste("Trajectory Plot for", beta_title, "and b =", b_v)
    
    # Generate the ggplot
    p <- ggplot(plot_data, aes(x = pyro, y = logBio)) +
      geom_path(aes(color = k), linewidth = 1.2, arrow = arrow(length = unit(0.3, "cm"))) + #
      scale_color_viridis_c() + #
      labs(
        x = "Pyro Index",
        y = "log(Geometric Species Index)",
        title = full_title,
        color = "k" 
      ) +
      theme_minimal() #
    
    # Create the exact filename for Overleaf (e.g., Pyro_Bio_beta_positive_b_0_5.png)
    b_str <- gsub("\\.", "_", as.character(b_v))
    file_name <- paste0("Figures/Pyro_Bio_beta_", scen, "_b_", b_str, ".png")
    
    # Save the plot
    ggsave(filename = file_name, plot = p, width = 6, height = 4, dpi = 300, bg = "white")
    cat("Saved:", file_name, "\n")
  }
}