# species_density.R

# Create the Figures directory if it doesn't exist
if(!dir.exists("Figures")) {
  dir.create("Figures", recursive = TRUE)
}

# Define the species abundance function A_i(t)
A_i <- function(t, alpha, beta) {
  alpha * t^(beta - 1)
}

# Time sequence (starting slightly above 0 to prevent infinity issues when beta < 1)
t <- seq(0.01, 5, length.out = 500)

# Colors for curves
cols <- c("black", "blue")

# ===================================================================
# Figure 2.2 (A): Early-successional response (beta < 1)
# ===================================================================

params_early <- list(
  list(alpha = 1, beta = 0.5),
  list(alpha = 2, beta = 0.5)
)

# Open the PNG device
png("Figures/Species_Density_Early_Successional.png", width = 1800, height = 1200, res = 300)

# Plot empty frame with larger axis elements
plot(t, A_i(t, 1, 1), type = "n", ylim = c(0, 5),
     xlab = "time since the most recent fire (t)", ylab = expression(A[i](t)),
     cex.lab = 1.3)

# Add each curve
for (i in seq_along(params_early)) {
  lines(t, A_i(t, params_early[[i]]$alpha, params_early[[i]]$beta), col = cols[i], lwd = 2)
}

# Add legend
legend("topright", legend = sapply(params_early, function(p) {
  paste0("alpha = ", p$alpha, ", beta = ", p$beta)
}), col = cols, lwd = 2, cex = 1.2)

# Close the device and save
dev.off()
cat("Saved: Figures/Species_Density_Early_Successional.png\n")

# ===================================================================
# Figure 2.2 (B): Late-successional response (beta > 1)
# ===================================================================

params_late <- list(
  list(alpha = 1, beta = 1.5),
  list(alpha = 2, beta = 1.5)
)

# Open the PNG device
png("Figures/Species_Density_Late_Successional.png", width = 1800, height = 1200, res = 300)

# Plot empty frame
plot(t, A_i(t, 1, 1), type = "n", ylim = c(0, 5),
     xlab = "time since the most recent fire (t)", ylab = expression(A[i](t)),
     cex.lab = 1.3)

# Add each curve
for (i in seq_along(params_late)) {
  lines(t, A_i(t, params_late[[i]]$alpha, params_late[[i]]$beta), col = cols[i], lwd = 2)
}

# Add legend (Moved to topleft so it doesn't overlap the rising lines)
legend("topleft", legend = sapply(params_late, function(p) {
  paste0("alpha = ", p$alpha, ", beta = ", p$beta)
}), col = cols, lwd = 2, cex = 1.2)

# Close the device and save
dev.off()
cat("Saved: Figures/Species_Density_Late_Successional.png\n")