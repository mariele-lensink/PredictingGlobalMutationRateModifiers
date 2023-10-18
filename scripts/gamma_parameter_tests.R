library(ggplot2)

# Generate random numbers from gamma distributions with different shape parameters
shape_values <- c(1, 2, 5, 10)
n <- 1000

# Create a list to store the generated random numbers
random_numbers <- list()

# Generate random numbers for each shape parameter value
for (shape in shape_values) {
  random_numbers[[as.character(shape)]] <- rgamma(n, shape = shape, rate = 1)
}

# Combine the random numbers into a single data frame
df <- data.frame(
  shape = rep(shape_values, each = n),
  value = unlist(random_numbers)
)

# Plot the distributions
ggplot(df, aes(x = value, fill = factor(shape))) +
  geom_density(alpha = 0.5) +
  labs(title = "Gamma Distribution with Different Shape Parameters", x = "Value", y = "Density") +
  scale_fill_discrete(name = "Shape Parameter")





# Generate random numbers from gamma distributions with different rate parameters
rate_values <- c(1, 2, 5, 10)
n <- 1000

# Create a list to store the generated random numbers
random_numbers <- list()

# Generate random numbers for each rate parameter value
for (rate in rate_values) {
  random_numbers[[as.character(rate)]] <-rgamma(n, shape = 4, rate = rate)
}

# Combine the random numbers into a single data frame
df <- data.frame(
  rate = rep(rate_values, each = n),
  value = unlist(random_numbers)
)


# Plot the distributions
ggplot(df, aes(x = value, fill = factor(rate))) +
  geom_density(alpha = 0.5) +
  labs(title = "Gamma Distribution with Different Rate Parameters", x = "Value", y = "Density") +
  scale_fill_discrete(name = "Rate Parameter")

