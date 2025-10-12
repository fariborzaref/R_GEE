#########################################################
# Generalized Estimating Equations (GEE) – Dr. Fariborz Aref
# Purpose: Estimate population-average effects in clustered social data
# Packages: geepack, dplyr
#########################################################

# Load libraries
library(geepack)
library(dplyr)

# --- 1. Simulate clustered social data ---
set.seed(2025)

# Suppose we have 60 communities (clusters) each with 10 respondents
n_clusters <- 60
cluster_size <- 10

data <- data.frame(
  community_id = rep(1:n_clusters, each = cluster_size),
  age = rnorm(n_clusters * cluster_size, mean = 40, sd = 10),
  income = rnorm(n_clusters * cluster_size, mean = 55000, sd = 12000),
  education = sample(c("Low", "Medium", "High"), n_clusters * cluster_size, replace = TRUE),
  gender = sample(c("Male", "Female"), n_clusters * cluster_size, replace = TRUE)
)

# Create a binary social-health outcome influenced by predictors
data <- data %>%
  mutate(
    health_score = 0.3 * (age / 10) + 0.00003 * income +
                   ifelse(education == "High", 1.2, ifelse(education == "Medium", 0.6, 0)) +
                   ifelse(gender == "Female", 0.5, 0) +
                   rnorm(n(), 0, 1),
    poor_health = ifelse(health_score > 5.5, 0, 1) # binary outcome (1 = poor health)
  )

# --- 2. Fit GEE model ---
gee_model <- geeglm(
  poor_health ~ age + income + education + gender,
  id = community_id,
  data = data,
  family = binomial(link = "logit"),
  corstr = "exchangeable"
)

# --- 3. Output ---
summary(gee_model)

# --- 4. Extract results ---
coef_df <- data.frame(
  Variable = names(coef(gee_model)),
  Estimate = coef(gee_model),
  Std_Error = summary(gee_model)$coefficients[, "Std.err"],
  P_value = summary(gee_model)$coefficients[, "Pr(>|W|)"]
)

print(coef_df)

# --- 5. Interpret ---
cat("\nInterpretation:\n")
cat("This GEE model estimates population-average effects of age, income, education, and gender\n")
cat("on the likelihood of poor health across clustered communities.\n")
cat("Exchangeable correlation assumes respondents within a community share similar outcomes.\n")
