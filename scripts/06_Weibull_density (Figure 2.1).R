# weibull_density.R
library(ggplot2)
library(dplyr)

# Create the Figures directory if it doesn't exist
if(!dir.exists("Figures")) {
  dir.create("Figures", recursive = TRUE)
}

# Define the Weibull hazard function in terms of b and k
weibull_hazard <- function(t, k, b) {
  b * k * (t^(k - 1))
}

# ===================================================================
# Figure 2.1 (A): Effect of the rate parameter (b)
# ===================================================================

# Generate data for b = 10 and b = 5 (with fixed k = 1.5)
data_b <- data.frame(t = rep(seq(0.1, 60, by = 0.1), 2)) %>%
  mutate(group = factor(rep(c("10", "5"), each = 600), levels = c("10", "5"))) %>%
  mutate(h_t = case_when(
    group == "10" ~ weibull_hazard(t, k = 1.5, b = 10),
    group == "5"  ~ weibull_hazard(t, k = 1.5, b = 5)
  ))

# Create the plot
p_rate <- ggplot(data_b, aes(x = t, y = h_t, color = group)) +
  geom_line(linewidth = 0.8) +
  theme_classic() +
  scale_color_manual(name = "Rate (b)\n(fixed k = 1.5)", values = c("red", "blue")) +
  labs(
    x = expression(paste("time-since-fire (", italic("t"), ")")),
    y = expression(paste("h(", italic("t"), ")"))
  )

# Save to Figures folder
ggsave("Figures/Weibull_Hazard_Rate_b.png", plot = p_rate, width = 6, height = 4, dpi = 300)
cat("Saved: Figures/Weibull_Hazard_Rate_b.png\n")

# ===================================================================
# Figure 2.1 (B): Effect of the shape parameter (k)
# ===================================================================

# Generate data for k = 0.7, 1.0, and 1.5 (with fixed b = 1)
data_k <- data.frame(t = rep(seq(0.1, 60, by = 0.1), 3)) %>%
  mutate(group = factor(rep(c("0.7", "1.0", "1.5"), each = 600), levels = c("0.7", "1.0", "1.5"))) %>%
  mutate(h_t = case_when(
    group == "0.7" ~ weibull_hazard(t, k = 0.7, b = 1),
    group == "1.0" ~ weibull_hazard(t, k = 1.0, b = 1),
    group == "1.5" ~ weibull_hazard(t, k = 1.5, b = 1)
  ))

# Create the plot
p_shape <- ggplot(data_k, aes(x = t, y = h_t, color = group)) +
  geom_line(linewidth = 0.8) +
  theme_classic() +
  scale_color_manual(name = "Shape (k)\n(fixed b = 1)", values = c("black", "blue", "red")) +
  coord_cartesian(ylim = c(0, 12)) + # Constrain y-axis since k=0.7 spikes near t=0
  labs(
    x = expression(paste("time-since-fire (", italic("t"), ")")),
    y = expression(paste("h(", italic("t"), ")"))
  )

# Save to Figures folder
ggsave("Figures/Weibull_Hazard_Shape_k.png", plot = p_shape, width = 6, height = 4, dpi = 300)
cat("Saved: Figures/Weibull_Hazard_Shape_k.png\n")