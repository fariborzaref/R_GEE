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
