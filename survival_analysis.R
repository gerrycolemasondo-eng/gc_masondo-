# =============================================================================
# SURVIVAL ANALYSIS FOR BREAST CANCER PATIENTS IN ZIMBABWE
# Case: Chitungwiza Central Hospital (2013-2019)
# Author: Gerald Collins Masondo - Midlands State University
# Tools: R 3.5.1 - survival, survminer, prodlim, Publish
# =============================================================================

# 1. SETUP
library(survival)
library(survminer)
library(prodlim)
library(dplyr)

# 2. LOAD DATA - Use relative path for GitHub
BC <- read.csv("breast_cancer_zimbabwe.csv", header = TRUE)

# 3. CLEANING & LABELLING
BC$SEX <- factor(BC$SEX, levels = c(0, 1), labels = c("Male", "Female"))
BC$MARITAL.STATUS <- factor(BC$MARITAL.STATUS, levels = c(1,2,3,4),
                            labels = c("Single","Married","Widowed","Divorced"))
BC$TREATMENT <- factor(BC$TREATMENT, levels = c(0,1), labels = c("No","Yes"))
BC$OCCUP <- factor(BC$OCCUP, levels = c(0,1), labels = c("Unemployed","Employed"))

# 4. KAPLAN-MEIER
km_fit <- survfit(Surv(SURVIVAL_TIME, STATUS) ~ AGE_GROUP, data = BC)
ggsurvplot(km_fit, data = BC, pval = TRUE, risk.table = TRUE,
           title = "KM Survival Curve by Age Group - Zimbabwe")

# 5. LOG-RANK TEST
survdiff(Surv(SURVIVAL_TIME, STATUS) ~ AGE_GROUP, data = BC)

# 6. COX MODEL & PH ASSUMPTION
cox_model <- coxph(Surv(SURVIVAL_TIME, STATUS) ~ AGE_GROUP + SEX + TREATMENT, data = BC)
summary(cox_model)
cox.zph(cox_model)  # PH test - all p > 0.05 = assumption met

# 7. FINAL MODEL - Significant: Age 41-50 (p=0.0038), 61-70 (p=0.0449)
cox_final <- coxph(Surv(SURVIVAL_TIME, STATUS) ~ AGE_GROUP, data = BC)
exp(cbind(HR = coef(cox_final), confint(cox_final)))
