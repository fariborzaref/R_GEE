############################################################
# GENERALIZED ESTIMATING EQUATIONS (GEE) — DR. FARIBORZ AREF
# Goal: Population-average effects with clustered data
# Features:
#   • Reads real data if present (gee_input.csv); else simulates
#   • Proper factor handling, reference levels
#   • Compares working correlations via QIC (exchangeable / AR-1 / independence)
#   • Odds ratios with robust 95% CI
#   • Sensitivity checks (link, correlation)
#   • Leave-one-cluster-out (LOCO) influence on a target coefficient
#   • Publication-grade plots (coef forest, marginal curve)
############################################################

# ====================== 0) PACKAGES & STYLE ======================
req <- c("geepack","data.table","dplyr","ggplot2","broom","purrr","tidyr")
to_install <- setdiff(req, rownames(installed.packages()))
if(length(to_install)) install.packages(to_install, repos = "https://cloud.r-project.org")
invisible(lapply(req, library, character.only = TRUE))

theme_set(theme_minimal(base_size = 13))
set.seed(2025)

# ====================== 1) DATA: READ OR SIMULATE =================
# If you have real data, create R_GEE/gee_input.csv with columns:
#   community_id, age, income, education, gender, outcome (0/1)
# education in {Low,Medium,High}; gender in {Male,Female}

path <- "R_GEE/gee_input.csv"
if (file.exists(path)) {
  df <- fread(path)
  names(df) <- tolower(names(df))
  stopifnot(all(c("community_id","age","income","education","gender","outcome") %in% names(df)))
  df <- df |> mutate(
    community_id = as.integer(community_id),
    education = factor(education, levels = c("Low","Medium","High")),
    gender    = factor(gender, levels = c("Male","Female")),
    poor_health = as.integer(outcome)
  )
} else {
  # Simulated clustered data (realistic scale)
  n_clusters <- 60; n_i <- 10
  df <- data.frame(
    community_id = rep(seq_len(n_clusters), each = n_i),
    age    = rnorm(n_clusters*n_i, 40, 10),
    income = rnorm(n_clusters*n_i, 55000, 12000),
    education = factor(sample(c("Low","Medium","High"), n_clusters*n_i, TRUE),
                       levels = c("Low","Medium","High")),
    gender = factor(sample(c("Male","Female"), n_clusters*n_i, TRUE),
                    levels = c("Male","Female"))
  )
  # latent logit for poor health (1 = poor)
  lp <- -3.2 + 0.05*(df$age - 40) - 0.00002*(df$income - 55000) -
        0.35*(df$education == "Medium") - 0.8*(df$education == "High") -
        0.25*(df$gender == "Female")
  p  <- 1/(1+exp(-lp))
  df$poor_health <- rbinom(nrow(df), 1, p)
}

# Quick sanity
stopifnot(all(df$poor_health %in% c(0,1)))
df <- df |> mutate(
  education = droplevels(education),
  gender    = droplevels(gender)
)

# ====================== 2) MODEL SPECIFICATION ====================
form <- poor_health ~ age + income + education + gender

# Fit with three working correlations
fit_exch <- geeglm(form, id = community_id, data = df,
                   family = binomial("logit"), corstr = "exchangeable")
fit_ar1  <- geeglm(form, id = community_id, data = df,
                   family = binomial("logit"), corstr = "ar1")
fit_ind  <- geeglm(form, id = community_id, data = df,
                   family = binomial("logit"), corstr = "independence")

# QIC comparison (lower is better)
qic_tbl <- tibble(
  model = c("exchangeable","ar1","independence"),
  QIC   = c(QIC(fit_exch)[1], QIC(fit_ar1)[1], QIC(fit_ind)[1]),
  QICu  = c(QIC(fit_exch)[2], QIC(fit_ar1)[2], QIC(fit_ind)[2])
) |> arrange(QIC)
cat("\n=== QIC (MODEL SELECTION) ===\n"); print(qic_tbl)

# Choose best working correlation by QIC
best_fit <- list(exch = fit_exch, ar1 = fit_ar1, ind = fit_ind)[[c("exch","ar1","ind")[which.min(qic_tbl$QIC)]]]
cat("\nSelected correlation structure:", best_fit$corstr, "\n")

# ====================== 3) RESULTS: OR WITH 95% CI =================
sum_best <- summary(best_fit)
coefs <- as.data.frame(sum_best$coefficients)
coefs$term <- rownames(coefs); rownames(coefs) <- NULL
coefs <- coefs |> select(term, Estimate, Std.err = Std.err, p = `Pr(>|W|)`) |>
  mutate(
    OR   = exp(Estimate),
    lo95 = exp(Estimate - 1.96*Std.err),
    hi95 = exp(Estimate + 1.96*Std.err)
  )

cat("\n=== POPULATION-AVERAGE EFFECTS (LOGIT LINK) ===\n")
print(coefs)

# ====================== 4) SENSITIVITY CHECKS =====================
# (a) Link function sensitivity: probit vs logit
fit_probit <- geeglm(form, id = community_id, data = df,
                     family = binomial("probit"), corstr = best_fit$corstr)
# (b) Correlation sensitivity: second-best QIC
second_idx <- order(qic_tbl$QIC)[2]
second_fit <- list(exch = fit_exch, ar1 = fit_ar1, ind = fit_ind)[[c("exch","ar1","ind")[second_idx]]]

sens_tbl <- tibble(
  spec   = c("Logit + Best corstr", "Probit + Best corstr", paste0("Logit + ", second_fit$corstr)),
  age    = c(coef(best_fit)["age"], coef(fit_probit)["age"], coef(second_fit)["age"]),
  income = c(coef(best_fit)["income"], coef(fit_probit)["income"], coef(second_fit)["income"])
)
cat("\n=== SENSITIVITY (KEY COEFFICIENTS) ===\n"); print(round(sens_tbl, 4))

# ====================== 5) CLUSTER INFLUENCE (LOCO) ===============
# Δ of a target coefficient when dropping each cluster
target_term <- "age"
base_beta   <- unname(coef(best_fit)[target_term])

delta <- df |>
  distinct(community_id) |>
  mutate(delta = map_dbl(community_id, function(cid) {
    sub <- df |> filter(community_id != cid)
    m   <- geeglm(form, id = community_id, data = sub,
                  family = binomial("logit"), corstr = best_fit$corstr)
    unname(coef(m)[target_term]) - base_beta
  }))

cat("\n=== LEAVE-ONE-CLUSTER-OUT Δ ON", target_term, "===\n"); print(arrange(delta, desc(abs(delta)))[1:10,])

# ====================== 6) VISUALIZATIONS =========================
# (a) Coefficient forest (OR with CI)
plot_df <- coefs |> filter(term != "(Intercept)") |>
  mutate(term = gsub("education", "educ:", term),
         term = gsub("gender", "gender:", term))

p_forest <- ggplot(plot_df, aes(x = OR, y = reorder(term, OR))) +
  geom_vline(xintercept = 1, linetype = 3) +
  geom_point(size = 2) +
  geom_errorbarh(aes(xmin = lo95, xmax = hi95), height = 0.15) +
  scale_x_continuous(trans = "log10") +
  labs(title = "Population-Average Odds Ratios (95% CI)",
       x = "Odds Ratio (log scale)", y = NULL)

# (b) Marginal predicted probabilities over age (median covariates)
ref <- df |>
  summarise(
    income = median(income, na.rm = TRUE),
    education = levels(education)[1],
    gender = levels(gender)[1]
  )
age_grid <- tibble(age = seq(20, 70, by = 1))
newd <- crossing(ref, age_grid)
lp   <- predict(best_fit, newdata = newd, type = "link", se.fit = TRUE)
newd$pred <- best_fit$family$linkinv(lp$fit)
newd$lo   <- best_fit$family$linkinv(lp$fit - 1.96*lp$se.fit)
newd$hi   <- best_fit$family$linkinv(lp$fit + 1.96*lp$se.fit)

p_marg <- ggplot(newd, aes(age, pred)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.15) +
  geom_line(size = 1) +
  labs(title = "Predicted Probability of Poor Health vs. Age",
       x = "Age", y = "Predicted probability") +
  coord_cartesian(ylim = c(0, 1))

# (c) Cluster influence plot
p_infl <- ggplot(delta, aes(x = reorder(factor(community_id), delta), y = delta)) +
  geom_hline(yintercept = 0, linetype = 3) +
  geom_point(size = 2) +
  coord_flip() +
  labs(title = paste0("LOCO Δ on ", target_term, " Coefficient"),
       x = "Community ID", y = expression(Delta~beta))

print(p_forest); print(p_marg); print(p_infl)

# ====================== 7) SAVE ARTIFACTS =========================
if (!dir.exists("R_GEE/out"))  dir.create("R_GEE/out", recursive = TRUE)
if (!dir.exists("R_GEE/figs")) dir.create("R_GEE/figs", recursive = TRUE)

fwrite(coefs, "R_GEE/out/gee_coefficients_or.csv")
fwrite(qic_tbl, "R_GEE/out/gee_qic_comparison.csv")
fwrite(delta,  "R_GEE/out/gee_loco_delta_age.csv")

ggsave("R_GEE/figs/gee_or_forest.png", p_forest, width = 7.4, height = 4.8, dpi = 300)
ggsave("R_GEE/figs/gee_marg_age.png",  p_marg,   width = 7.4, height = 4.6, dpi = 300)
ggsave("R_GEE/figs/gee_loco_age.png",  p_infl,   width = 7.4, height = 6.0, dpi = 300)

cat("\n=== DONE ===\nArtifacts written to R_GEE/out and R_GEE/figs.\n")

