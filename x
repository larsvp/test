############################################
#### Preparation: packages #################
############################################

# List of used packages
packages <- c("ggplot2", "tidyverse", "dplyr", "ggpubr",
              "funModeling", "Hmisc", "knitr")

# Check whether all required packages are installed
pack_installed <- c()
for (i in seq_along(packages)) {
  pack_installed[i] <- packages[i] %in% rownames(installed.packages())
  if (!pack_installed[i]) {
    print(paste(packages[i], " is NOT installed."))
  }
}; if (all(pack_installed)) {
  print("All required packages have already been installed.")
}

# install.packages("ggplot2")
# install.packages("tidyverse")
# install.packages("dplyr")
# install.packages("ggpubr")
# install.packages("funModeling")
# install.packages("Hmisc")
# install.packages("knitr")

library(ggplot2)
library(tidyverse)
library(dplyr)
library(ggpubr)
library(funModeling)
library(Hmisc)
library(knitr)

############################################
#### Preparation: Reading in the data ######
############################################

nasal <- c(15.4, 13.5, 13.3, 12.4, 12.8, 13.5, 14.5, 13.9,
           11.0, 15.0, 17.0, 13.8, 17.4, 16.5, 14.4)
endo <- c(16.5, 13.2, 13.6, 13.6, 14.0, 14.0, 16.0, 14.1,
          11.5, 14.4, 16.0, 13.2, 16.6, 18.5, 14.5)

# Create data frame
data <- data.frame(cbind(nasal, endo))

############################################
#### 1. Exploration ########################
############################################

# Scatterplot Nasal Brushing vs Endobronchial Forceps Biopsy
plot(data$nasal,
     data$endo,
     main="Nasal Brushing vs Endobronchial Forceps Biopsy",
     xlab="Nasal Brushing",
     ylab="Endobronchial Forceps Biopsy",
     pch=19)
# Regression line
abline(lm(endo~nasal), col="red")
# The scatterplot indicates that there’s a relatively strong
# positive linear correlation between Nasal Brushing and
# Endobronchial Forceps Biopsy. Additionally, the graph shows
# that the scales for Nasal Brushing and Endobronchial Forceps Biopsy
# are very similar.

# To be able to interpret the scatterplot better,
# the x and y-axis are put on an equal scale.
ggplot(data, aes(x=nasal, y=endo)) +
  geom_abline(intercept=0, slope=1, color="gray", size=1) +
  geom_point(shape=16, color="blue", size=3) +
  geom_smooth(method=lm,  linetype=1, color="red", fill=NA, size=1.2) +
  ggtitle("Nasal Brushing vs Endobronchial Forceps Biopsy") +
  xlab("Nasal Brushing") +
  ylab("Endobronchial Forceps Biopsy") +
  scale_x_discrete(limits=c(11:18), drop=c(19)) +
  scale_y_discrete(limits=c(11:18), drop=c(19)) +
  coord_cartesian(xlim = c(11, 18.5), ylim=c(11, 18.5))

# Overview of the variables
glimpse(data)

# Examine missing values and outliers
df_status(data)

# Examine distribution
profiling_num(data)

# Mean, sd, percentiles, skewness, kurtosis, ...
plot_num(data)

# Distribution
par(mfrow=c(2,2))
hist(nasal)
qqnorm(nasal, main='Q-Q plot of Nasal'); qqline(nasal, col=2, lwd=2, lty=2)
hist(endo)
qqnorm(endo, main='Q-Q plot of Endo');qqline(endo, col=2, lwd=2, lty=2)
par(mfrow=c(1,1))

# Boxplot
boxplot(nasal, endo, data=data)

############################################
#### 2. Spearman's rank correlation ########
############################################

# Dataframe of ranked variables
data.ranked <- data.frame(cbind(
  rank(data$nasal, ties.method="average"),
  rank(data$endo, ties.method="average")
))
colnames(data.ranked) <- c("nasal", "endo")

# Amount of observations
n <- length(data.ranked$nasal)

# Transposed table of ranks
kable(as.data.frame(t(as.matrix(data.ranked))))

# Amount of ties
n.x_ties <- n - length(unique(data.ranked$nasal))
n.y_ties <- n - length(unique(data.ranked$endo))

# Sum of squares
sum_squares <- sum((data.ranked$nasal - data.ranked$endo)^2)

# Spearman formula
spearman <- function(sum_sq, n) {
  1 - ((6 * sum_sq) / (n^3 - n))
}

# Spearman rank correlation coefficient
spearman_rho_manual <- spearman(sum_squares, n)
spearman_rho_manual

# Spearman with cor.test()
spearmman_rho_cort <- cor.test(x=data$nasal,
                           y=data$endo,
                           method="spearman",
                           exact=TRUE)
spearmman_rho_cort

# Spearman as Pearson correlation coefficient between ranked variables
# -> equal result as cor.test()
cov_sd_spearman <- cov(data.ranked) /
  (sd(data.ranked$nasal) * sd(data.ranked$endo))
cov_sd_spearman[[2]]
as.numeric(spearmman_rho_cort$estimate)

# Pearson of ranked vars: sum of squares
# -> equal to sum of squares from cor.test()
pearson <- cor(x=data.ranked$nasal, y=data.ranked$endo, method='pearson')
s <- (n^3 - n) * (1 - pearson) / 6
s
as.numeric(spearmman_rho_cort$statistic) # S value from cor.test()

# Pearson of ranked vars: p-value
# -> equal to the p-value from cor.test()
t <- pearson * (sqrt((n - 2) / (1 - pearson^2)))
p <- 2 * (1 - pt(t, n - 2))
p
spearmman_rho_cort$p.value # P-value from cor.test()

############################################
#### 3. Permutation test ###################
############################################

# Setting a seed for reproducibility
set.seed(71991)

# Settings
R <- 10000                    # Amount of replicates
simulations_rho_manual <- c() # Result array for manual rho

# Permutations
for (i in 1:R) {
  # Sample and rank both variables
  nasal_r <- sample(data.ranked$nasal, replace=F)
  endo_r <- sample(data.ranked$endo, replace=F)
  sample_r <- data.frame(cbind(nasal_r, endo_r))
  # Calculate spearman
  spearman_r <- cov(sample_r) /
                (sd(sample_r$nasal_r) * sd(sample_r$endo_r))
  simulations_rho_manual[i] <- spearman_r[[2]]
}



# Critical value
critical_value_manual <- quantile(simulations_rho_manual, probs=0.95)
critical_value_manual

# P-value
p.value_manual <- mean(simulations_rho_manual >=
                         as.numeric(spearmman_rho_cort$estimate))
p.value_manual

# Histogram with critical value and the observed spearman
hist(simulations_rho_manual, breaks=20,xlim=c(-1,1),
     main="Exact permutation null distribution manual")
abline(v=critical_value_manual, col="blue") # Critical value
abline(v=spearman_rho_manual, col="red")    # Observed spearman

############################################
#### 4. Asymptotic approximation ###########
############################################

# Standardized rho:
options(digits=5)

# T-value
t <- spearman_rho_manual*(sqrt((15-2)/(1-(spearman_rho_manual)^2)))
t

# Critical values for t-distribution (95%))
criticalValues <- c(qt(0.025, 13), qt(0.975, 13))

# P-value
p <- 2 * (1 - pt(t, n - 2))
p

# Spearman with the Hmisc package
spearman.test(nasal,endo)







