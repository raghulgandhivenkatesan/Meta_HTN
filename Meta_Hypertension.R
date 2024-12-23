####################Load necessary libraries#######################
library(metafor)
library(meta)

####################Data for meta-analysis########################
hypertension_data <- data.frame(
  Author = c("Mahapatra et al., (2021)", "Santosh et al., (2019)", 
             "Basavanagowdappa et al., (2016)", "Boppana et al., (2023)", 
             "Annie Caroline et al., (2021)", "Bharatia et al., (2016)", 
             "Gupta et al., (2019)", "Mandal et al., (2014)", 
             "Umamaheswari et al., (2019)","Sagarad et al.,(2013)","Garg et al.,(2018)"),
  Year = c(2021, 2019, 2016, 2023, 2021, 2016, 2019, 2014, 2019, 2013, 2018), # Pub. years
  Events = c(31, 67, 248, 100, 130, 922, 525, 70, 18, 30, 4), # Hypertension cases
  N = c(275, 169, 1537, 200, 200, 4725, 3073, 300, 100, 990, 42) # Total sample size
)
# Sorting the data by Year in ascending order
hypertension_data <- hypertension_data[order(hypertension_data$Year), ]
# Display the sorted data
print(hypertension_data)


####################Step 1: Calculate observed proportions######################
# Calculate prevalence (proportions) and their variances
ies <- escalc(measure = "PR", xi = Events, ni = N, data = hypertension_data)
print(ies)

####################Step 2: Perform random-effects meta-analysis##################
# Using DerSimonian-Laird (DL) method for random-effects model
pes <- rma(yi, vi, data = ies, method = "DL")
print(pes, digits = 2)
confint(pes)

####################Step 3: Identifying outliers with residuals######################
stud_res <- rstudent(pes)
abs_z <- abs(stud_res$z)
outliers <- stud_res[order(-abs_z), ]
print(outliers)

####################Step 4: Leave-One-Out Analysis####################################
l1o <- leave1out(pes)
# Extract estimates and confidence intervals
yi <- l1o$estimate             # Leave-one-out effect sizes
sei <- l1o$se                  # Standard errors
ci.lb <- yi - 1.96 * sei       # Lower bounds of 95% CI
ci.ub <- yi + 1.96 * sei       # Upper bounds of 95% CI

# Create Forest Plot for Leave-One-Out Analysis (based on prevalence)
forest(yi, sei = sei,
       slab = hypertension_data$Author, 
       xlab = "Summary Proportions Leaving Out Each Study",
       refline = pes$b,         # Overall summary proportion
       digits = 4,              # Number of digits
       alim = c(min(ci.lb), max(ci.ub)),  # Adjust x-axis limits
       cex = 0.8,               # Text size
       col = "blue",            # Point color
       main = "Leave-One-Out Sensitivity Analysis (Prevalence)")

# Calculate I-squared values for Leave-One-Out Analysis
isq <- sapply(1:nrow(hypertension_data), function(i) {
  temp_pes <- rma(yi = ies$yi[-i], vi = ies$vi[-i], method = "DL") # Subset yi and vi correctly
  return(temp_pes$I2)
})

# Create Forest Plot for Leave-One-Out Analysis (based on I-squared)
forest(isq, sei = rep(0, length(isq)),
       slab = hypertension_data$Author,
       xlab = "I-squared Values Leaving Out Each Study",
       refline = pes$I2,         # Overall I-squared value
       digits = 2,               # Number of digits
       alim = c(min(isq), max(isq)),  # Adjust x-axis limits
       cex = 0.8,                # Text size
       col = "red",             # Point color
       main = "Leave-One-Out Sensitivity Analysis (I-squared)")

####################Step 5: Baujat Plot (Heterogeneity Diagnostics)####################
baujat(pes, main = "Baujat Plot for Heterogeneity Diagnostics")

####################Step 6: Influence Diagnostics######################################
inf <- influence(pes)
print(inf)
# Plot the Influence Diagnostics
par(mar = c(5, 5, 4, 2) + 0.1)
plot(inf, cex = 0.8)  # Reduce text size to 80% of default
# Add a custom title to the plot
mtext("Influence Diagnostics for Hypertension Meta-Analysis", side = 3, line = 1, cex = 1.2, font = 2)

####################Step 7: Remove Outliers (if needed, for example 10,11)############
outlier_removed <- ies[-c(10,11), ]
pes_updated <- rma(yi, vi, data = outlier_removed, method = "DL")
print(pes_updated, digits = 2)

####################Step 8: Forest Plot for Overall Meta-Analysis#####################
# Forest plot for all studies
overall_meta <- metaprop(
  event = hypertension_data$Events,
  n = hypertension_data$N,
  studlab = hypertension_data$Author,
  sm = "PR", method.ci = "NAsm", method.tau = "DL"
)
# Generate forest plot for overall meta-analysis
forest.meta(
  overall_meta,
  leftcols = c("studlab", "event", "n", "w.random"), # Study details and weights
  leftlabs = c("Study", "Cases", "Total", "Weight"),
  rightcols = c("effect", "ci"), # Effect size and confidence intervals
  rightlabs = c("Prevalence", "95% CI"),
  comb.fixed = FALSE,           # Disable common-effect (fixed-effect) model
  comb.random = TRUE,           # Display random-effects model only
  col.square = "blue",          # Color for individual study markers
  col.diamond = "red",          # Color for random-effect model (diamond)
  col.diamond.lines = "black",  # Diamond outline
  type.random = "diamond",      # Random-effect model as a diamond
  prediction = FALSE,            # Show prediction interval
  print.Q = TRUE,               # Display Q statistic for heterogeneity
  print.I2 = TRUE,              # Display I-squared heterogeneity
  print.tau2 = TRUE,            # Display tau-squared heterogeneity
  overall = TRUE,               # Include overall random-effect summary
  fixed = FALSE,                # Completely suppress fixed-effects model
  digits = 3,                   # Precision for numbers
  xlim = c(0, 1),               # Adjust x-axis range
  main = "Forest Plot: Random Effect Model Only"
)
# Generate forest plot for meta-analysis after removing outliers
forest.meta(
  metaprop(
    event = hypertension_data$Events[-c(10,11)],         # Remove outlier events
    n = hypertension_data$N[-c(10,11)],                  # Remove outlier sample sizes
    studlab = hypertension_data$Author[-c(10,11)],       # Remove outlier study labels
    sm = "PR", method.ci = "NAsm", method.tau = "DL"
  ),
  leftcols = c("studlab", "event", "n", "w.random"), # Study details and weights
  leftlabs = c("Study", "Cases", "Total", "Weight"),
  rightcols = c("effect", "ci"), # Effect size and confidence intervals
  rightlabs = c("Prevalence", "95% CI"),
  comb.fixed = FALSE,           # Disable common-effect (fixed-effect) model
  comb.random = TRUE,           # Display random-effects model only
  col.square = "blue",          # Color for individual study markers
  col.diamond = "red",          # Color for random-effect model (diamond)
  col.diamond.lines = "black",  # Diamond outline
  type.random = "diamond",      # Random-effect model as a diamond
  prediction = FALSE,            # Show prediction interval
  print.Q = TRUE,               # Display Q statistic for heterogeneity
  print.I2 = TRUE,              # Display I-squared heterogeneity
  print.tau2 = TRUE,            # Display tau-squared heterogeneity
  overall = TRUE,               # Include overall random-effect summary
  fixed = FALSE,                # Completely suppress fixed-effects model
  digits = 3,                   # Precision for numbers
  xlim = c(0, 1),               # Adjust x-axis range
  main = "Forest Plot: Random Effect Model After Removing Outlier (Study 5)"
)
-----------------------------------------------------------------------------------------
# (Optional)#Display forest plot for overall studies in R console
forest(overall_meta,
       leftcols = c("studlab", "event", "n", "w.random"),
       leftlabs = c("Study", "Cases", "Total", "Weight"),
       rightcols = c("effect", "ci"),
       rightlabs = c("Prevalence", "95% CI"),
       col.square = "blue", col.diamond = "red",
       col.diamond.lines = "black",
       type.random = "diamond",  # Ensure random-effects model is displayed
       common = FALSE,           # Remove common-effect model
       overall = TRUE,           # Display random-effects summary
       print.tau2 = TRUE, print.I2 = TRUE, print.Q = TRUE, # Include heterogeneity stats
       prediction = TRUE,        # Ensure prediction interval is displayed
       digits = 3, 
       xlim = c(0, 1),           # Proper x-axis range
       main = "Forest Plot with Random Effect Model and Prediction Interval")
---------------------------------------------------------------------------------------
#####################Step 9: Funnel Plot for Publication Bias###########################
funnel(pes, main = "Funnel Plot for Hypertension Data")

#####################Step 10: Trim-and-Fill Analysis for Publication Bias###################
pes_trimfill <- trimfill(pes)
print(pes_trimfill)
funnel(pes_trimfill, main = "Trim-and-Fill Funnel Plot")

#####################Step 11: Egger's Regression Test for Small-Study Effects##############
egger_test <- regtest(pes, model = "lm", predictor = "sei")
print(egger_test)

#####################Step 12: Rank Correlation Test######################################
rank_corr <- ranktest(pes)
print(rank_corr)


####################Step 13: Perform meta-regression with sample size as a moderator#################
meta_reg_sample <- rma(yi, vi, mods = ~ N, data = ies, method = "DL")
print(meta_reg_sample, digits = 2)
# Summary of meta-regression
summary(meta_reg_sample)

####################Step 13a: Perform meta-regression with COVID period as a moderator#################
# Add a classification column for Pre- and Post-COVID years
hypertension_data$COVID_Period <- ifelse(hypertension_data$Year >= 2020, "Post-COVID", "Pre-COVID")

# Convert COVID_Period to a factor and set Pre-COVID as the reference level
hypertension_data$COVID_Period <- factor(hypertension_data$COVID_Period, levels = c("Pre-COVID", "Post-COVID"))

# Add COVID_Period to the ies dataframe
ies$COVID_Period <- hypertension_data$COVID_Period

# Perform meta-regression with COVID_Period as a moderator
meta_reg_covid <- rma(yi, vi, mods = ~ COVID_Period, data = ies, method = "DL")
print(meta_reg_covid, digits = 2)
# Summary of meta-regression
summary(meta_reg_covid)

####################Step 13b: Perform bivariate meta-regression with COVID period and sample size#################
meta_reg_bivariate <- rma(yi, vi, mods = ~ COVID_Period + N, data = ies, method = "DL")
print(meta_reg_bivariate, digits = 2)
# Summary of bivariate meta-regression
summary(meta_reg_bivariate)

####################Step 14: Bubble Plot for Meta-Regression with Sample Size########################

# Bubble plot for Sample Size with modified CI and PI lines
bubble_sample <- predict(meta_reg_sample, newmods = ies$N)
y_range_sample <- range(c(bubble_sample$ci.lb, bubble_sample$ci.ub, ies$yi))

plot(ies$N, ies$yi,
     xlab = "Sample Size (N)", 
     ylab = "Observed Effect Size (Log Proportions)", 
     main = "Bubble Plot for Meta-Regression (Sample Size)",
     pch = 21, bg = "lightblue", cex = 1.5, ylim = y_range_sample)

# Add meta-regression line
lines(sort(ies$N), bubble_sample$pred[order(ies$N)], col = "blue", lwd = 2)

# Add 95% confidence intervals (solid black lines)
lines(sort(ies$N), bubble_sample$ci.lb[order(ies$N)], col = "black", lty = 1, lwd = 1.5)
lines(sort(ies$N), bubble_sample$ci.ub[order(ies$N)], col = "black", lty = 1, lwd = 1.5)

# Add legend
legend("topright", 
       legend = c("Meta-Regression Line", "95% CI"), 
       col = c("blue", "black"), 
       lty = c(1, 1), 
       lwd = c(2, 1.5))

####################Step 14a: Bubble Plot for Meta-Regression with COVID Period###########

# Predict values for bubble plot with COVID Period
bubble_covid <- predict(meta_reg_covid, newmods = as.numeric(ies$COVID_Period == "Post-COVID"))

# Create the bubble plot
plot(as.numeric(ies$COVID_Period), ies$yi,
     xlab = "COVID Period (1 = Pre-COVID, 2 = Post-COVID)", 
     ylab = "Observed Effect Size (Log Proportions)", 
     main = "Bubble Plot for Meta-Regression (COVID Period)",
     pch = 21, bg = "lightgreen", cex = 1.5, 
     ylim = range(c(bubble_covid$ci.lb, bubble_covid$ci.ub)))

# Add meta-regression line
lines(as.numeric(ies$COVID_Period), bubble_covid$pred, col = "blue", lwd = 2)

# Add 95% confidence intervals (solid black lines)
lines(as.numeric(ies$COVID_Period), bubble_covid$ci.lb, col = "black", lty = 1, lwd = 1.5)
lines(as.numeric(ies$COVID_Period), bubble_covid$ci.ub, col = "black", lty = 1, lwd = 1.5)

# Add legend
legend("topleft", 
       legend = c("Meta-Regression Line", "95% CI"), 
       col = c("blue", "black"), 
       lty = c(1, 1), 
       lwd = c(2, 1.5))

################### Step 15: Distribution of the True Effect ###################
# Ensure that pes$b and pes$tau2 are numeric and properly defined
if (!is.null(pes$b) && !is.null(pes$tau2) && pes$tau2 > 0) {
  true_effects <- seq(
    from = as.numeric(pes$b - 4 * sqrt(pes$tau2)),
    to = as.numeric(pes$b + 4 * sqrt(pes$tau2)),
    length.out = 1000
  )
} else {
  stop("Error: Invalid 'pes' object or missing 'tau2' value.")
}
# Density of the distribution
density_values <- dnorm(true_effects, mean = as.numeric(pes$b), sd = sqrt(as.numeric(pes$tau2)))

# Plot the distribution
plot(true_effects, density_values, type = "l", lwd = 2, col = "blue",
     xlab = "True Effect Size (Log Proportions)",
     ylab = "Density",
     main = "Distribution of the True Effect")

# Add vertical line for the estimated overall effect size
abline(v = as.numeric(pes$b), col = "red", lwd = 2, lty = 2)

# Add legend
legend("topright", legend = c("True Effect Distribution", "Overall Effect Size"),
       col = c("blue", "red"), lty = c(1, 2), lwd = 2)




