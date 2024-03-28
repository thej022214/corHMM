corHMM:::compute_neglnlikelihood

set.seed(123) # For reproducibility
data <- rnorm(100, mean = 50, sd = 10) # 100 random normal data points

# Generate a grid of parameter values
mu_vals <- seq(40, 60, by = 0.5)
sigma_vals <- seq(5, 25, by = 0.5)

# Initialize a matrix to store likelihood values
lnLikelihood_matrix <- matrix(nrow = length(mu_vals), ncol = length(sigma_vals))

# Calculate likelihood for each combination of mu and sigma
for (i in 1:length(mu_vals)) {
  for (j in 1:length(sigma_vals)) {
    lnLikelihood_matrix[i, j] <- sum(dnorm(data, mean = mu_vals[i], sd = sigma_vals[j], log = TRUE))
  }
}

# Convert log-likelihoods to likelihoods for plotting
likelihood_matrix <- exp(lnLikelihood_matrix - max(lnLikelihood_matrix))

# Now plot the contour plot with corrected dimensions
par(mfrow=c(1,2))
contour(mu_vals, sigma_vals, lnLikelihood_matrix, xlab = "Mu", ylab = "Sigma", main = "Contour Likelihood Plot", nlevels = 80, lty = 2)
points(x = mu_vals[23], y = sigma_vals[9], pch = 21, bg="red")
image(mu_vals, sigma_vals, lnLikelihood_matrix)

# Plot the contour plot
contour(sigma_vals, mu_vals, likelihood_matrix, xlab = "Sigma", ylab = "Mu", main = "Contour Likelihood Plot")
# Assuming likelihood_matrix is filled correctly with:
# Rows corresponding to mu_vals (y-axis)
# Columns corresponding to sigma_vals (x-axis)

# Directly use the original matrix with contour plotting
?contour(sigma_vals, mu_vals, likelihood_matrix, xlab = "Sigma", ylab = "Mu", main = "Contour Likelihood Plot")



require(grDevices) # for colours
require(grDevices) # for colours
require(grDevices) # for colours
x <- -6:16
op <- par(mfrow = c(2, 2))
contour(outer(x, x), method = "edge", vfont = c("sans serif", "plain"))
z <- outer(x, sqrt(abs(x)), FUN = "/")
image(x, x, z)
contour(x, x, z, col = "pink", add = TRUE, method = "edge",
        vfont = c("sans serif", "plain"))
contour(x, x, z, ylim = c(1, 6), method = "simple", labcex = 1,
        xlab = quote(x[1]), ylab = quote(x[2]))
contour(x, x, z, ylim = c(-6, 6), nlevels = 20, lty = 2, method = "simple",
        main = "20 levels; \"simple\" labelling method")
par(op)


