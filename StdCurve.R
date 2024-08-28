library(ggplot2)
library(ggpmisc)
library(dplyr)

data <- data.frame(
  Conc = c(0,0.0125,0.025,0.05,0.1,0.2,0.3),
  Abs = c(0.02133333,0.08933333,0.17966667,0.35166667,0.755,1.57966667,2.1855),
  SD = c(0.008082904,0.033857545,0.016441817,0.013650397,0.011269428,0.106978191,0.064346717)
)

model <- lm(Abs ~ Conc, data = data)

# Extract model coefficients and R^2 value
intercept <- coef(model)[1]
slope <- coef(model)[2]
r_squared <- summary(model)$r.squared

# Create the scatterplot with error bars and best-fit line
plot <- ggplot(data, aes(x = Conc, y = Abs)) +
  geom_point() +  # Plot the points
  geom_errorbar(aes(ymin = Abs - SD, ymax = Abs + SD)) +  # Add error bars
  geom_smooth(method = "lm", color = "cornflowerblue", se = FALSE) +  # Add best-fit line
  theme_bw() +  # Minimal theme
  labs(title = "Proline Standard Curve, Absorbance at 520nm against Proline Standards in mg", 
       x = "Proline Concentration (mg)", 
       y = "Absorbance at 520nm (Abs)" ) + 
  annotate("text", x = 0.1, y = 2.2, 
           label = sprintf("y = %.3f + %.2fx\nR² = %.3f", intercept, slope, r_squared), 
           hjust = 1.1, vjust = 1.5, size = 5, color = "cornflowerblue", fontface = "italic") +
  theme(plot.title = element_text(size=12.5))

# Print the plot
print(plot)

