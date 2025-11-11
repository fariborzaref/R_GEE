# Generalized Estimating Equations (GEE)

**Author:** Dr. Fariborz Aref  
**Discipline:** Quantitative Sociology & Inequality Research  
**Date:** October 17, 2025  
**License:** MIT  

---

### Purpose
Estimate **population-average effects** in clustered social data to study structural inequality and health disparities.  
Implements **logit GEE models**, **QIC-based model comparison**, and **leave-one-cluster-out (LOCO)** sensitivity checks.

---

### Structure
R_GEE/


---
############################################################
# GENERALIZED ESTIMATING EQUATIONS (GEE)
# Author: Dr. Fariborz Aref | Sociologist & Quantitative Methodologist
# Date: 2025-10-17 | License: MIT
#
# Goal: Estimate population average effects in clustered data for health and inequality work.
# Input: R_GEE/gee_input.csv with variables:
#   community_id, age, income, education{Low,Medium,High}, gender{Male,Female}, outcome{0,1}
#   If missing, synthetic data are generated.
# Output:
#   R_GEE/out/gee_coefficients_or.csv
#   R_GEE/out/gee_qic_comparison.csv
#   R_GEE/out/gee_loco_delta_age.csv
#   R_GEE/figs/gee_or_forest.png
#   R_GEE/figs/gee_marg_age.png
#   R_GEE/figs/gee_loco_age.png
############################################################

# Packages
req <- c("geepack","data.table","dplyr","ggplot2","broom","purrr","tidyr","tibble")
to_install <- setdiff(req, rownames(installed.packages()))
if (length(to_install)) install.packages(to_install, repos = "https://cloud.r-project.org")
suppressPackageStartupMessages(invisible(lapply(req, library, character.only = TRUE)))

# Plot theme (academic, compact)
theme_set(
  theme_minimal(base_size = 11, base_family = "serif") +
    theme(
      plot.title   = element_text(size = 12, face = "bold", hjust = 0.5),
      axis.title   = element_text(size = 10),
      axis.text    = element_text(size = 9),
      strip.text   = element_text(size = 10, face = "bold"),
      legend.title = element_text(size = 10),
      legend.text  = element_text(size = 9),
      panel.grid.minor = element_blank()
    )
)

set.seed(2025)
options(stringsAsFactors = FALSE)

# Paths
in_path   <- "R_GEE/gee_input.csv"
out_dir   <- "R_GEE/out"
fig_dir   <- "R_GEE/figs"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

# Data
if (file.exists(in_path)) {
  df <- data.table::fread(in_path)
  names(df) <- tolower(names(df))
  need <- c("community_id","age","income","education","gender","outcome")
  if (!all(need %in% names(df))) stop("Input file missing required columns.")
  df <- df |>
    dplyr::mutate(
      community_id = as.integer(community_id),
      education    = factor(education, levels = c("Low","Medium","High")),
      gender       = factor(gender,    levels = c("Male","Female")),
      poor_health  = as.integer(outcome)
    )
} else {
  n_clusters <- 60
  n_i <- 10
  N <- n_clusters * n_i
  df <- tibble::tibble(
    community_id = rep(seq_len(n_clusters), each = n_i),
    age          = rnorm(N, 40, 10),
    income       = rnorm(N, 55000, 12000),
    education    = factor(sample(c("Low","Medium","High"), N, TRUE),
                          levels = c("Low","Medium","High")),
    gender       = factor(sample(c("Male","Female"), N, TRUE),
                          levels = c("Male","Female"))
  )
  linpred <- -3.2 +
             0.05*(df$age - 40) -
             0.00002*(df$income - 55000) -
             0.35*(df$education == "Medium") -
             0.80*(df$education == "High") -
             0.25*(df$gender == "Female")
  p <- 1/(1 + exp(-linpred))
  df$poor_health <- rbinom(nrow(df), 1, p)
}

# Checks
if (!all(df$poor_health %in% c(0,1))) stop("Outcome must be 0 or 1.")
df <- df |>
  dplyr::mutate(education = droplevels(education),
                gender    = droplevels(gender))

# Model specification
form <- poor_health ~ age + income + education + gender

# Fit GEE with common working correlation structures
fit_exch <- geepack::geeglm(form, id = community_id, data = df,
                            family = binomial("logit"), corstr = "exchangeable")
fit_ar1  <- geepack::geeglm(form, id = community_id, data = df,
                            family = binomial("logit"), corstr = "ar1")
fit_ind  <- geepack::geeglm(form, id = community_id, data = df,
                            family = binomial("logit"), corstr = "independence")

# QIC comparison
QIC_pair <- function(f) { q <- geepack::QIC(f); c(QIC = unname(q[1]), QICu = unname(q[2])) }
q_exch <- QIC_pair(fit_exch); q_ar1 <- QIC_pair(fit_ar1); q_ind <- QIC_pair(fit_ind)

qic_tbl <- tibble::tibble(
  model = c("exchangeable","ar1","independence"),
  QIC   = c(q_exch["QIC"], q_ar1["QIC"], q_ind["QIC"]),
  QICu  = c(q_exch["QICu"], q_ar1["QICu"], q_ind["QICu"])
) |>
  dplyr::arrange(QIC)

# Select best model
order_map <- c("exchangeable" = "exch", "ar1" = "ar1", "independence" = "ind")
fits <- list(exch = fit_exch, ar1 = fit_ar1, ind = fit_ind)
best_label <- order_map[qic_tbl$model[1]]
best_fit <- fits[[best_label]]

# Odds ratios with confidence
sum_best <- summary(best_fit)
coefs <- as.data.frame(sum_best$coefficients)
coefs$term <- rownames(coefs); rownames(coefs) <- NULL
coefs <- coefs |>
  dplyr::select(term, Estimate, Std.err = Std.err, p = `Pr(>|W|)`) |>
  dplyr::mutate(
    OR   = exp(Estimate),
    lo95 = exp(Estimate - 1.96*Std.err),
    hi95 = exp(Estimate + 1.96*Std.err)
  )

# Within cluster correlation note for exchangeable
icc_note <- NA_real_
if (tolower(best_fit$corstr) == "exchangeable") icc_note <- as.numeric(best_fit$geese$alpha)

# Sensitivity checks
fit_probit <- geepack::geeglm(form, id = community_id, data = df,
                              family = binomial("probit"), corstr = best_fit$corstr)
second_fit <- fits[[ order_map[qic_tbl$model[2]] ]]

sens_tbl <- tibble::tibble(
  spec   = c("Logit + best corstr", "Probit + best corstr", paste0("Logit + ", second_fit$corstr)),
  age    = c(coef(best_fit)["age"],    coef(fit_probit)["age"],    coef(second_fit)["age"]),
  income = c(coef(best_fit)["income"], coef(fit_probit)["income"], coef(second_fit)["income"])
)

# Leave one cluster out influence on target coefficient
target_term <- "age"
base_beta   <- unname(coef(best_fit)[target_term])
delta <- df |>
  dplyr::distinct(community_id) |>
  dplyr::mutate(delta = purrr::map_dbl(community_id, function(cid) {
    sub <- dplyr::filter(df, community_id != cid)
    m   <- geepack::geeglm(form, id = community_id, data = sub,
                           family = binomial("logit"), corstr = best_fit$corstr)
    unname(coef(m)[target_term]) - base_beta
  }))

# Plots
plot_df <- coefs |>
  dplyr::filter(term != "(Intercept)") |>
  dplyr::mutate(term = gsub("education", "educ:", term),
                term = gsub("gender", "gender:", term))

p_forest <- ggplot2::ggplot(plot_df, aes(x = OR, y = reorder(term, OR))) +
  ggplot2::geom_point(size = 1.8) +
  ggplot2::geom_errorbarh(aes(xmin = lo95, xmax = hi95), height = 0.12) +
  ggplot2::scale_x_continuous(trans = "log10") +
  ggplot2::labs(
    title = "Population average odds ratios with 95 percent confidence",
    x = "Odds ratio, log scale", y = NULL
  )

ref <- df |>
  dplyr::summarise(
    income = stats::median(income, na.rm = TRUE),
    education = levels(df$education)[1],
    gender = levels(df$gender)[1]
  )

age_grid <- tibble::tibble(age = seq(20, 70, by = 1))
newd <- tidyr::crossing(ref, age_grid)
lp   <- stats::predict(best_fit, newdata = newd, type = "link", se.fit = TRUE)
newd$pred <- best_fit$family$linkinv(lp$fit)
newd$lo   <- best_fit$family$linkinv(lp$fit - 1.96*lp$se.fit)
newd$hi   <- best_fit$family$linkinv(lp$fit + 1.96*lp$se.fit)

p_marg <- ggplot2::ggplot(newd, aes(age, pred)) +
  ggplot2::geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.12) +
  ggplot2::geom_line(linewidth = 0.6) +
  ggplot2::labs(
    title = "Predicted probability of poor health by age",
    x = "Age", y = "Predicted probability"
  ) +
  ggplot2::coord_cartesian(ylim = c(0, 1))

p_infl <- ggplot2::ggplot(delta, aes(x = reorder(factor(community_id), delta), y = delta)) +
  ggplot2::geom_point(size = 1.8) +
  ggplot2::coord_flip() +
  ggplot2::labs(
    title = paste0("Leave one cluster out delta on ", target_term),
    x = "Community ID", y = expression(Delta~beta)
  )

# Save
data.table::fwrite(coefs,  file.path(out_dir, "gee_coefficients_or.csv"))
data.table::fwrite(qic_tbl, file.path(out_dir, "gee_qic_comparison.csv"))
data.table::fwrite(delta,  file.path(out_dir, "gee_loco_delta_age.csv"))

ggplot2::ggsave(file.path(fig_dir, "gee_or_forest.png"), p_forest, width = 6.5, height = 4.2, dpi = 300)
ggplot2::ggsave(file.path(fig_dir, "gee_marg_age.png"),  p_marg,   width = 6.5, height = 4.0, dpi = 300)
ggplot2::ggsave(file.path(fig_dir, "gee_loco_age.png"),  p_infl,   width = 6.5, height = 5.0, dpi = 300)

# Console summary
cat("\nQIC model comparison (lower is better):\n"); print(qic_tbl)
cat("Selected working correlation: ", best_fit$corstr, "\n", sep = "")
if (!is.na(icc_note)) {
  cat(sprintf("Within cluster correlation alpha (exchangeable): %.3f\n", icc_note))
} else {
  cat("Within cluster correlation not directly reported for this structure.\n")
}
cat("\nDone. Results saved to R_GEE/out and R_GEE/figs.\n")

### Key Methods
- **Modeling:** Generalized Estimating Equations (logit / probit)
- **Comparison:** QIC and QICu metrics  
- **Sensitivity:** Leave-one-cluster-out (LOCO) delta  
- **Visualization:** Predicted probabilities, odds ratios, influence plots  

---

### Quick Example
```r
library(geepack)
fit <- geeglm(poor_health ~ age + income + education + gender,
              id = community_id, data = df,
              family = binomial("logit"), corstr = "exchangeable")
summary(fit)
