# Document 1

library(tidyverse)
library(MASS)
library(broom)
library(jsonlite)
select <- dplyr::select

data <- read_csv("D:/GVA/Excel Files/Gulf_Final_Merged_MASTER.csv")

data <- data %>%
  mutate(POP_SQMI = E_TOTPOP / AREA_SQMI)

url <- "https://api.census.gov/data/2022/acs/acs5?get=NAME,B01001_001E,B01001_002E&for=county:*"
raw <- fromJSON(url)
sex_df <- as.data.frame(raw[-1,], stringsAsFactors = FALSE)
colnames(sex_df) <- c("NAME", "total_pop", "male_pop", "state", "county")
sex_df <- sex_df %>%
  mutate(
    total_pop = as.numeric(total_pop),
    male_pop  = as.numeric(male_pop),
    EP_MALE   = (male_pop / total_pop) * 100,
    FIPS      = paste0(state, county)
  ) %>%
  dplyr::select(FIPS, EP_MALE)

data <- data %>% left_join(sex_df, by = "FIPS")

#Complete Case Data 
data_clean <- data %>%
  filter(
    !is.na(total_events), !is.na(Homicide),
    !is.na(RUCC_2023),
    !is.na(EP_POV150), !is.na(EP_UNEMP), !is.na(EP_AGE65), !is.na(EP_NOHSDP),
    !is.na(POP_SQMI), !is.na(EP_AFAM), !is.na(EP_HISP), !is.na(EP_MALE),
    !is.na(HeatIndex_2022_Mean), !is.na(Green_2022),
    !is.na(eal_score), !is.na(Total_Disasters_2005_2023),
    E_TOTPOP > 0
  )

cat("N complete cases:", nrow(data_clean), "counties\n")

# ============================================================
# MODEL SET A: DV = total_events (Firearm Total Count)
# ============================================================

# Model A1: Base predictors
modelA1 <- glm(
  total_events ~
    EP_POV150 + EP_UNEMP + EP_AGE65 + EP_NOHSDP + POP_SQMI +
    EP_AFAM + EP_HISP + EP_MALE + RUCC_2023,
  family = "poisson", data = data_clean
)
cat("\n=== MODEL A1 (DV=total_events, Base) ===\n")
print(summary(modelA1))
cat("Dispersion:", deviance(modelA1)/df.residual(modelA1), "\n")
rr_A1 <- exp(cbind(RR = coef(modelA1), confint(modelA1)))
cat("Rate Ratios:\n"); print(rr_A1)

# Model A2: Full model
modelA2 <- glm(
  total_events ~
    EP_POV150 + EP_UNEMP + EP_AGE65 + EP_NOHSDP + POP_SQMI +
    EP_AFAM + EP_HISP + EP_MALE + RUCC_2023 +
    HeatIndex_2022_Mean + Green_2022 + eal_score + Total_Disasters_2005_2023,
  family = "poisson", data = data_clean
)
cat("\n=== MODEL A2 (DV=total_events, Full) ===\n")
print(summary(modelA2))
cat("Dispersion:", deviance(modelA2)/df.residual(modelA2), "\n")
rr_A2 <- exp(cbind(RR = coef(modelA2), confint(modelA2)))
cat("Rate Ratios:\n"); print(rr_A2)


# Model B1: Base predictors + population offset
modelB1 <- glm(
  Homicide ~
    EP_POV150 + EP_UNEMP + EP_AGE65 + EP_NOHSDP + POP_SQMI +
    EP_AFAM + EP_HISP + EP_MALE + RUCC_2023 +
    offset(log(E_TOTPOP)),
  family = "poisson", data = data_clean
)
cat("\n=== MODEL B1 (DV=Homicide Rate, Base) ===\n")
print(summary(modelB1))
cat("Dispersion:", deviance(modelB1)/df.residual(modelB1), "\n")
rr_B1 <- exp(cbind(RR = coef(modelB1), confint(modelB1)))
cat("Rate Ratios:\n"); print(rr_B1)

# Model B2: Full model + population offset
modelB2 <- glm(
  Homicide ~
    EP_POV150 + EP_UNEMP + EP_AGE65 + EP_NOHSDP + POP_SQMI +
    EP_AFAM + EP_HISP + EP_MALE + RUCC_2023 +
    HeatIndex_2022_Mean + Green_2022 + eal_score + Total_Disasters_2005_2023 +
    offset(log(E_TOTPOP)),
  family = "poisson", data = data_clean
)
cat("\n=== MODEL B2 (DV=Homicide Rate, Full) ===\n")
print(summary(modelB2))
cat("Dispersion:", deviance(modelB2)/df.residual(modelB2), "\n")
rr_B2 <- exp(cbind(RR = coef(modelB2), confint(modelB2)))
cat("Rate Ratios:\n"); print(rr_B2)


cat("\n=== MODEL FIT COMPARISON ===\n")
models <- list(
  "A1 (total_events, Base)"  = modelA1,
  "A2 (total_events, Full)"  = modelA2,
  "B1 (Homicide Rate, Base)" = modelB1,
  "B2 (Homicide Rate, Full)" = modelB2
)
for (nm in names(models)) {
  m <- models[[nm]]
  cat(sprintf("\n%s\n", nm))
  cat(sprintf("  Null deviance:     %.2f (df=%d)\n", m$null.deviance, m$df.null))
  cat(sprintf("  Residual deviance: %.2f (df=%d)\n", deviance(m), df.residual(m)))
  cat(sprintf("  AIC:               %.1f\n", AIC(m)))
  cat(sprintf("  Dispersion ratio:  %.4f\n", deviance(m)/df.residual(m)))
  cat(sprintf("  N:                 %d\n", nrow(data_clean)))
}


cat("\n=== CONFIRMING OVERDISPERSION FROM POISSON ===\n")
cat("Model A (total_events) dispersion was: 200.61 and 157.53\n")
cat("Model B (Homicide)     dispersion was:  35.91 and  33.31\n")
cat("All far above 2.0 — Negative Binomial is appropriate.\n")


# ============================================================
# MODEL SET A (NB): DV = total_events (Firearm Total Count)
# ============================================================

# --- Model A1 NB: Base predictors ---
modelA1_nb <- glm.nb(
  total_events ~
    EP_POV150 + EP_UNEMP + EP_AGE65 + EP_NOHSDP + POP_SQMI +
    EP_AFAM + EP_HISP + EP_MALE + RUCC_2023,
  data = data_clean
)
cat("\n=== MODEL A1 NB (DV=total_events, Base) ===\n")
print(summary(modelA1_nb))
cat("Dispersion (theta):", modelA1_nb$theta, "\n")
cat("Dispersion ratio (resid dev / df):", deviance(modelA1_nb)/df.residual(modelA1_nb), "\n")
cat("Rate Ratios:\n")
rr_A1_nb <- exp(cbind(RR = coef(modelA1_nb), confint(modelA1_nb)))
print(rr_A1_nb)

# --- Model A2 NB: Full model ---
modelA2_nb <- glm.nb(
  total_events ~
    EP_POV150 + EP_UNEMP + EP_AGE65 + EP_NOHSDP + POP_SQMI +
    EP_AFAM + EP_HISP + EP_MALE + RUCC_2023 +
    HeatIndex_2022_Mean + Green_2022 + eal_score + Total_Disasters_2005_2023,
  data = data_clean
)
cat("\n=== MODEL A2 NB (DV=total_events, Full) ===\n")
print(summary(modelA2_nb))
cat("Dispersion (theta):", modelA2_nb$theta, "\n")
cat("Dispersion ratio (resid dev / df):", deviance(modelA2_nb)/df.residual(modelA2_nb), "\n")
cat("Rate Ratios:\n")
rr_A2_nb <- exp(cbind(RR = coef(modelA2_nb), confint(modelA2_nb)))
print(rr_A2_nb)

# ============================================================
# MODEL SET B (NB): DV = Homicide (with population offset)
# ============================================================

# --- Model B1 NB: Base predictors + population offset ---
modelB1_nb <- glm.nb(
  Homicide ~
    EP_POV150 + EP_UNEMP + EP_AGE65 + EP_NOHSDP + POP_SQMI +
    EP_AFAM + EP_HISP + EP_MALE + RUCC_2023 +
    offset(log(E_TOTPOP)),
  data = data_clean
)
cat("\n=== MODEL B1 NB (DV=Homicide Rate, Base) ===\n")
print(summary(modelB1_nb))
cat("Dispersion (theta):", modelB1_nb$theta, "\n")
cat("Dispersion ratio (resid dev / df):", deviance(modelB1_nb)/df.residual(modelB1_nb), "\n")
cat("Rate Ratios:\n")
rr_B1_nb <- exp(cbind(RR = coef(modelB1_nb), confint(modelB1_nb)))
print(rr_B1_nb)

# --- Model B2 NB: Full model + population offset ---
modelB2_nb <- glm.nb(
  Homicide ~
    EP_POV150 + EP_UNEMP + EP_AGE65 + EP_NOHSDP + POP_SQMI +
    EP_AFAM + EP_HISP + EP_MALE + RUCC_2023 +
    HeatIndex_2022_Mean + Green_2022 + eal_score + Total_Disasters_2005_2023 +
    offset(log(E_TOTPOP)),
  data = data_clean
)
cat("\n=== MODEL B2 NB (DV=Homicide Rate, Full) ===\n")
print(summary(modelB2_nb))
cat("Dispersion (theta):", modelB2_nb$theta, "\n")
cat("Dispersion ratio (resid dev / df):", deviance(modelB2_nb)/df.residual(modelB2_nb), "\n")
cat("Rate Ratios:\n")
rr_B2_nb <- exp(cbind(RR = coef(modelB2_nb), confint(modelB2_nb)))
print(rr_B2_nb)

# ============================================================
# FULL MODEL FIT COMPARISON: Poisson vs NB side by side
# ============================================================
cat("\n=== FULL MODEL FIT COMPARISON (Poisson vs NB) ===\n")

cat(sprintf("\n%-30s %12s %16s %10s %12s\n",
            "Model", "Resid Dev (df)", "AIC", "Dispersion", "N"))

# Poisson (from previous run — hardcoded for reference)
poisson_ref <- list(
  list(name="A1 Poisson (total_events)",  resdev="81,849 (408)", aic="84,112", disp="200.61"),
  list(name="A2 Poisson (total_events)",  resdev="63,643 (404)", aic="65,913", disp="157.53"),
  list(name="B1 Poisson (Homicide Rate)", resdev="14,651 (408)", aic="16,593", disp=" 35.91"),
  list(name="B2 Poisson (Homicide Rate)", resdev="13,459 (404)", aic="15,409", disp=" 33.31")
)
for (p in poisson_ref) {
  cat(sprintf("%-30s %12s %16s %10s %12s\n", p$name, p$resdev, p$aic, p$disp, "418"))
}

cat("\n--- Negative Binomial ---\n")
nb_models <- list(
  list(name="A1 NB (total_events)",  m=modelA1_nb),
  list(name="A2 NB (total_events)",  m=modelA2_nb),
  list(name="B1 NB (Homicide Rate)", m=modelB1_nb),
  list(name="B2 NB (Homicide Rate)", m=modelB2_nb)
)
for (nm in nb_models) {
  m <- nm$m
  cat(sprintf("%-30s %12s %16s %10s %12s\n",
              nm$name,
              paste0(round(deviance(m), 1), " (", df.residual(m), ")"),
              round(AIC(m), 1),
              round(deviance(m)/df.residual(m), 4),
              nrow(data_clean)
  ))
}

# ============================================================
# MISSING COUNTIES EXPLANATION
# ============================================================
cat("\n=== MISSING COUNTIES BREAKDOWN ===\n")
data_full <- data  # full 534-county dataset with EP_MALE merged

missing_summary <- data_full %>%
  mutate(
    miss_gva     = is.na(total_events) | is.na(Homicide),
    miss_rucc    = is.na(RUCC_2023),
    miss_climate = is.na(HeatIndex_2022_Mean) | is.na(Green_2022),
    miss_nri     = is.na(eal_score),
    miss_any     = miss_gva | miss_rucc | miss_climate | miss_nri
  ) %>%
  summarise(
    Total_counties     = n(),
    Missing_GVA        = sum(miss_gva),
    Missing_RUCC       = sum(miss_rucc),
    Missing_Climate    = sum(miss_climate),
    Missing_NRI        = sum(miss_nri),
    Missing_any        = sum(miss_any),
    Included_in_model  = sum(!miss_any)
  )

print(as.data.frame(missing_summary))

cat("\nCounties missing GVA data by state:\n")
data_full %>%
  filter(is.na(total_events)) %>%
  count(STATE, sort = TRUE) %>%
  print(n = 20)


# Document 2

library(tidyverse)
library(MASS)
library(broom)
library(jsonlite)
library(readxl)
select <- dplyr::select


data <- read_csv("D:/GVA/Excel Files/Gulf_Final_Merged_MASTER.csv")
data <- data %>% mutate(POP_SQMI = E_TOTPOP / AREA_SQMI)

url <- "https://api.census.gov/data/2022/acs/acs5?get=NAME,B01001_001E,B01001_002E&for=county:*"
raw <- fromJSON(url)
sex_df <- as.data.frame(raw[-1,], stringsAsFactors = FALSE)
colnames(sex_df) <- c("NAME", "total_pop", "male_pop", "state", "county")
sex_df <- sex_df %>%
  mutate(
    total_pop = as.numeric(total_pop),
    male_pop  = as.numeric(male_pop),
    EP_MALE   = (male_pop / total_pop) * 100,
    FIPS      = paste0(state, county)
  ) %>%
  dplyr::select(FIPS, EP_MALE)
data <- data %>% left_join(sex_df, by = "FIPS")


rucc_new <- read_excel("C:/Desktop/Research/GVA/Excel Files/Ruralurbancontinuumcodes2023 (1).xlsx")

cat("=== NEW RUCC FILE COLUMNS ===\n")
print(names(rucc_new))
cat("Rows:", nrow(rucc_new), "\n")
print(head(rucc_new))

rucc_new_clean <- rucc_new %>%
  mutate(FIPS = str_pad(as.character(FIPS), width = 5, pad = "0")) %>%
  dplyr::select(FIPS, RUCC_2023)

cat("\nRUCC new file FIPS sample:\n")
print(head(rucc_new_clean))
cat("Urban (RUCC 1-3):", sum(rucc_new_clean$RUCC_2023 <= 3, na.rm = TRUE), "\n")
cat("Rural (RUCC 4-9):", sum(rucc_new_clean$RUCC_2023 >= 4, na.rm = TRUE), "\n")

data <- data %>%
  dplyr::select(-RUCC_2023) %>%
  left_join(rucc_new_clean, by = "FIPS")

cat("\nRUCC NAs after merge:", sum(is.na(data$RUCC_2023)), "\n")
cat("RUCC distribution:\n")
print(table(data$RUCC_2023, useNA = "always"))

cat("\n=== ZERO COUNT CHECK ===\n")
cat("Counties with total_events = 0 (NA replaced with 0):\n")

data <- data %>%
  mutate(
    total_events = ifelse(is.na(total_events), 0, total_events),
    Homicide     = ifelse(is.na(Homicide), 0, Homicide)
  )

cat("total_events == 0:", sum(data$total_events == 0), "\n")
cat("Homicide == 0:", sum(data$Homicide == 0), "\n")
cat("total_events > 0:", sum(data$total_events > 0), "\n")

zero_pct_events   <- mean(data$total_events == 0) * 100
zero_pct_homicide <- mean(data$Homicide == 0) * 100
cat(sprintf("Zero %% for total_events: %.1f%%\n", zero_pct_events))
cat(sprintf("Zero %% for Homicide:     %.1f%%\n", zero_pct_homicide))
cat("If zero % > 30%, zero-inflated models (ZIP/ZINB) should be considered.\n")

base_vars <- c("EP_POV150", "EP_UNEMP", "EP_AGE65", "EP_NOHSDP", "POP_SQMI",
               "EP_AFAM", "EP_HISP", "EP_MALE", "RUCC_2023",
               "HeatIndex_2022_Mean", "Green_2022", "eal_score",
               "Total_Disasters_2005_2023", "E_TOTPOP")

data_with_zeros <- data %>%
  filter(
    !is.na(RUCC_2023),
    !is.na(EP_POV150), !is.na(EP_UNEMP), !is.na(EP_AGE65), !is.na(EP_NOHSDP),
    !is.na(POP_SQMI), !is.na(EP_AFAM), !is.na(EP_HISP), !is.na(EP_MALE),
    !is.na(HeatIndex_2022_Mean), !is.na(Green_2022),
    !is.na(eal_score), !is.na(Total_Disasters_2005_2023),
    E_TOTPOP > 0
  )

data_no_zeros <- data_with_zeros %>%
  filter(total_events > 0, Homicide > 0)

cat("\n=== SAMPLE SIZES ===\n")
cat("With zeros  — N:", nrow(data_with_zeros), "\n")
cat("Without zeros — N:", nrow(data_no_zeros), "\n")


run_nb_pair <- function(dv, offset_var = NULL, dataset, label) {
  formula_base <- as.formula(paste0(
    dv, " ~ EP_POV150 + EP_UNEMP + EP_AGE65 + EP_NOHSDP + POP_SQMI + ",
    "EP_AFAM + EP_HISP + EP_MALE + RUCC_2023",
    if (!is.null(offset_var)) paste0(" + offset(log(", offset_var, "))") else ""
  ))
  formula_full <- as.formula(paste0(
    dv, " ~ EP_POV150 + EP_UNEMP + EP_AGE65 + EP_NOHSDP + POP_SQMI + ",
    "EP_AFAM + EP_HISP + EP_MALE + RUCC_2023 + ",
    "HeatIndex_2022_Mean + Green_2022 + eal_score + Total_Disasters_2005_2023",
    if (!is.null(offset_var)) paste0(" + offset(log(", offset_var, "))") else ""
  ))
  
  cat(sprintf("\n========== %s — Model 1 (Base) ==========\n", label))
  m1 <- glm.nb(formula_base, data = dataset)
  print(summary(m1))
  cat("Theta:", m1$theta, "| Dispersion ratio:", round(deviance(m1)/df.residual(m1), 4), "\n")
  rr1 <- exp(cbind(RR = coef(m1), confint(m1)))
  cat("Rate Ratios:\n"); print(rr1)
  
  cat(sprintf("\n========== %s — Model 2 (Full) ==========\n", label))
  m2 <- glm.nb(formula_full, data = dataset)
  print(summary(m2))
  cat("Theta:", m2$theta, "| Dispersion ratio:", round(deviance(m2)/df.residual(m2), 4), "\n")
  rr2 <- exp(cbind(RR = coef(m2), confint(m2)))
  cat("Rate Ratios:\n"); print(rr2)
  
  cat(sprintf("\nFit — %s:\n", label))
  cat(sprintf("  M1: Null Dev=%.1f(df=%d) | Resid Dev=%.1f(df=%d) | AIC=%.1f | Disp=%.4f | N=%d\n",
              m1$null.deviance, m1$df.null, deviance(m1), df.residual(m1), AIC(m1),
              deviance(m1)/df.residual(m1), nrow(dataset)))
  cat(sprintf("  M2: Null Dev=%.1f(df=%d) | Resid Dev=%.1f(df=%d) | AIC=%.1f | Disp=%.4f | N=%d\n",
              m2$null.deviance, m2$df.null, deviance(m2), df.residual(m2), AIC(m2),
              deviance(m2)/df.residual(m2), nrow(dataset)))
  
  return(list(m1 = m1, m2 = m2))
}

cat("\n\n############################################################\n")
cat("SET A: DV = total_events — WITH ZEROS\n")
cat("############################################################\n")
modA_wz <- run_nb_pair("total_events", offset_var = NULL,
                       dataset = data_with_zeros,
                       label = "NB: total_events (WITH zeros)")

cat("\n\n############################################################\n")
cat("SET B: DV = total_events — WITHOUT ZEROS\n")
cat("############################################################\n")
modA_nz <- run_nb_pair("total_events", offset_var = NULL,
                       dataset = data_no_zeros,
                       label = "NB: total_events (WITHOUT zeros)")

cat("\n\n############################################################\n")
cat("SET C: DV = Homicide Rate — WITH ZEROS\n")
cat("############################################################\n")
modB_wz <- run_nb_pair("Homicide", offset_var = "E_TOTPOP",
                       dataset = data_with_zeros,
                       label = "NB: Homicide Rate (WITH zeros)")

cat("\n\n############################################################\n")
cat("SET D: DV = Homicide Rate — WITHOUT ZEROS\n")
cat("############################################################\n")
modB_nz <- run_nb_pair("Homicide", offset_var = "E_TOTPOP",
                       dataset = data_no_zeros,
                       label = "NB: Homicide Rate (WITHOUT zeros)")

cat("\n=== FINAL MISSING DATA SUMMARY ===\n")
cat("Total Gulf Coast counties in master dataset: 534\n")
cat("Counties with complete predictor data (with zeros):", nrow(data_with_zeros), "\n")
cat("Counties with complete predictor data (without zeros):", nrow(data_no_zeros), "\n")
cat("RUCC NAs remaining:", sum(is.na(data$RUCC_2023)), "\n")

cat("\nZero counts summary:\n")
cat(sprintf("  total_events zeros: %d (%.1f%%)\n",
            sum(data_with_zeros$total_events == 0),
            mean(data_with_zeros$total_events == 0) * 100))
cat(sprintf("  Homicide zeros:     %d (%.1f%%)\n",
            sum(data_with_zeros$Homicide == 0),
            mean(data_with_zeros$Homicide == 0) * 100))


# Document 3


library(tidyverse)
library(MASS)
library(jsonlite)
library(readxl)
select <- dplyr::select

data <- read_csv("D:/GVA/Excel Files/Gulf_Final_Merged_MASTER.csv")
data <- data %>% mutate(POP_SQMI = E_TOTPOP / AREA_SQMI)

url <- "https://api.census.gov/data/2022/acs/acs5?get=NAME,B01001_001E,B01001_002E&for=county:*"
raw <- fromJSON(url)
sex_df <- as.data.frame(raw[-1,], stringsAsFactors = FALSE)
colnames(sex_df) <- c("NAME", "total_pop", "male_pop", "state", "county")
sex_df <- sex_df %>%
  mutate(total_pop = as.numeric(total_pop), male_pop = as.numeric(male_pop),
         EP_MALE = (male_pop / total_pop) * 100, FIPS = paste0(state, county)) %>%
  dplyr::select(FIPS, EP_MALE)
data <- data %>% left_join(sex_df, by = "FIPS")

rucc_new <- read_excel("C:/Desktop/Research/GVA/Excel Files/Ruralurbancontinuumcodes2023 (1).xlsx")
rucc_new_clean <- rucc_new %>%
  mutate(FIPS = str_pad(as.character(FIPS), width = 5, pad = "0")) %>%
  dplyr::select(FIPS, RUCC_2023)
data <- data %>% dplyr::select(-RUCC_2023) %>% left_join(rucc_new_clean, by = "FIPS")

data <- data %>%
  mutate(
    total_events = ifelse(is.na(total_events), 0, total_events),
    Homicide     = ifelse(is.na(Homicide), 0, Homicide),
    Assault      = ifelse(is.na(Assault), 0, Assault)
  )

cat(sprintf("Assault zeros: %d (%.1f%%)\n", sum(data$Assault == 0), mean(data$Assault == 0)*100))

data_with_zeros <- data %>%
  filter(!is.na(RUCC_2023), !is.na(EP_POV150), !is.na(EP_UNEMP), !is.na(EP_AGE65),
         !is.na(EP_NOHSDP), !is.na(POP_SQMI), !is.na(EP_AFAM), !is.na(EP_HISP),
         !is.na(EP_MALE), !is.na(HeatIndex_2022_Mean), !is.na(Green_2022),
         !is.na(eal_score), !is.na(Total_Disasters_2005_2023), E_TOTPOP > 0)

data_no_zeros <- data_with_zeros %>% filter(Assault > 0)

cat("With zeros N:", nrow(data_with_zeros), "\n")
cat("Without zeros (Assault > 0) N:", nrow(data_no_zeros), "\n")

# ============================================================
# MODEL C1-WZ: Base predictors, WITH zeros
# ============================================================
modelC1_wz <- glm.nb(
  Assault ~ EP_POV150 + EP_UNEMP + EP_AGE65 + EP_NOHSDP + POP_SQMI +
    EP_AFAM + EP_HISP + EP_MALE + RUCC_2023 + offset(log(E_TOTPOP)),
  data = data_with_zeros)
cat("\n=== MODEL C1-WZ (Assault Rate, Base, WITH zeros) ===\n")
print(summary(modelC1_wz))
cat("Theta:", modelC1_wz$theta, "| Dispersion:", round(deviance(modelC1_wz)/df.residual(modelC1_wz),4), "\n")
rr_C1_wz <- exp(cbind(RR = coef(modelC1_wz), confint(modelC1_wz)))
cat("Rate Ratios:\n"); print(rr_C1_wz)
cat(sprintf("Fit: NullDev=%.1f(df=%d) | ResidDev=%.1f(df=%d) | AIC=%.1f | Disp=%.4f | N=%d\n",
            modelC1_wz$null.deviance, modelC1_wz$df.null, deviance(modelC1_wz),
            df.residual(modelC1_wz), AIC(modelC1_wz),
            deviance(modelC1_wz)/df.residual(modelC1_wz), nrow(data_with_zeros)))

# ============================================================
# MODEL C2-WZ: Full model, WITH zeros
# ============================================================
modelC2_wz <- glm.nb(
  Assault ~ EP_POV150 + EP_UNEMP + EP_AGE65 + EP_NOHSDP + POP_SQMI +
    EP_AFAM + EP_HISP + EP_MALE + RUCC_2023 +
    HeatIndex_2022_Mean + Green_2022 + eal_score + Total_Disasters_2005_2023 +
    offset(log(E_TOTPOP)),
  data = data_with_zeros)
cat("\n=== MODEL C2-WZ (Assault Rate, Full, WITH zeros) ===\n")
print(summary(modelC2_wz))
cat("Theta:", modelC2_wz$theta, "| Dispersion:", round(deviance(modelC2_wz)/df.residual(modelC2_wz),4), "\n")
rr_C2_wz <- exp(cbind(RR = coef(modelC2_wz), confint(modelC2_wz)))
cat("Rate Ratios:\n"); print(rr_C2_wz)
cat(sprintf("Fit: NullDev=%.1f(df=%d) | ResidDev=%.1f(df=%d) | AIC=%.1f | Disp=%.4f | N=%d\n",
            modelC2_wz$null.deviance, modelC2_wz$df.null, deviance(modelC2_wz),
            df.residual(modelC2_wz), AIC(modelC2_wz),
            deviance(modelC2_wz)/df.residual(modelC2_wz), nrow(data_with_zeros)))

# ============================================================
# MODEL C1-NZ: Base predictors, WITHOUT zeros
# ============================================================
modelC1_nz <- glm.nb(
  Assault ~ EP_POV150 + EP_UNEMP + EP_AGE65 + EP_NOHSDP + POP_SQMI +
    EP_AFAM + EP_HISP + EP_MALE + RUCC_2023 + offset(log(E_TOTPOP)),
  data = data_no_zeros)
cat("\n=== MODEL C1-NZ (Assault Rate, Base, WITHOUT zeros) ===\n")
print(summary(modelC1_nz))
cat("Theta:", modelC1_nz$theta, "| Dispersion:", round(deviance(modelC1_nz)/df.residual(modelC1_nz),4), "\n")
rr_C1_nz <- exp(cbind(RR = coef(modelC1_nz), confint(modelC1_nz)))
cat("Rate Ratios:\n"); print(rr_C1_nz)
cat(sprintf("Fit: NullDev=%.1f(df=%d) | ResidDev=%.1f(df=%d) | AIC=%.1f | Disp=%.4f | N=%d\n",
            modelC1_nz$null.deviance, modelC1_nz$df.null, deviance(modelC1_nz),
            df.residual(modelC1_nz), AIC(modelC1_nz),
            deviance(modelC1_nz)/df.residual(modelC1_nz), nrow(data_no_zeros)))

# ============================================================
# MODEL C2-NZ: Full model, WITHOUT zeros
# ============================================================
modelC2_nz <- glm.nb(
  Assault ~ EP_POV150 + EP_UNEMP + EP_AGE65 + EP_NOHSDP + POP_SQMI +
    EP_AFAM + EP_HISP + EP_MALE + RUCC_2023 +
    HeatIndex_2022_Mean + Green_2022 + eal_score + Total_Disasters_2005_2023 +
    offset(log(E_TOTPOP)),
  data = data_no_zeros)
cat("\n=== MODEL C2-NZ (Assault Rate, Full, WITHOUT zeros) ===\n")
print(summary(modelC2_nz))
cat("Theta:", modelC2_nz$theta, "| Dispersion:", round(deviance(modelC2_nz)/df.residual(modelC2_nz),4), "\n")
rr_C2_nz <- exp(cbind(RR = coef(modelC2_nz), confint(modelC2_nz)))
cat("Rate Ratios:\n"); print(rr_C2_nz)
cat(sprintf("Fit: NullDev=%.1f(df=%d) | ResidDev=%.1f(df=%d) | AIC=%.1f | Disp=%.4f | N=%d\n",
            modelC2_nz$null.deviance, modelC2_nz$df.null, deviance(modelC2_nz),
            df.residual(modelC2_nz), AIC(modelC2_nz),
            deviance(modelC2_nz)/df.residual(modelC2_nz), nrow(data_no_zeros)))



#Document 4


library(tidyverse)
library(MASS)
library(jsonlite)
library(readxl)
select <- dplyr::select


data <- read_csv("D:/GVA/Excel Files/Gulf_Final_Merged_MASTER.csv")
data <- data %>% mutate(POP_SQMI = E_TOTPOP / AREA_SQMI)

# Fetch % Male from Census ACS 2022
url <- "https://api.census.gov/data/2022/acs/acs5?get=NAME,B01001_001E,B01001_002E&for=county:*"
raw <- fromJSON(url)
sex_df <- as.data.frame(raw[-1,], stringsAsFactors = FALSE)
colnames(sex_df) <- c("NAME", "total_pop", "male_pop", "state", "county")
sex_df <- sex_df %>%
  mutate(
    total_pop = as.numeric(total_pop),
    male_pop  = as.numeric(male_pop),
    EP_MALE   = (male_pop / total_pop) * 100,
    FIPS      = paste0(state, county)
  ) %>%
  dplyr::select(FIPS, EP_MALE)
data <- data %>% left_join(sex_df, by = "FIPS")

rucc_new <- read_excel("C:/Desktop/Research/GVA/Excel Files/Ruralurbancontinuumcodes2023 (1).xlsx")
rucc_new_clean <- rucc_new %>%
  mutate(FIPS = str_pad(as.character(FIPS), width = 5, pad = "0")) %>%
  dplyr::select(FIPS, RUCC_2023)
data <- data %>%
  dplyr::select(-RUCC_2023) %>%
  left_join(rucc_new_clean, by = "FIPS")
cat("RUCC NAs after merge:", sum(is.na(data$RUCC_2023)), "\n")

# Replace NAs with 0 for GVA variables (true zeros)
data <- data %>%
  mutate(
    total_events = ifelse(is.na(total_events), 0, total_events),
    Homicide     = ifelse(is.na(Homicide), 0, Homicide),
    Assault      = ifelse(is.na(Assault), 0, Assault)
  )

# Zero count check
cat(sprintf("total_events zeros: %d (%.1f%%)\n",
            sum(data$total_events == 0), mean(data$total_events == 0)*100))
cat(sprintf("Homicide zeros:     %d (%.1f%%)\n",
            sum(data$Homicide == 0), mean(data$Homicide == 0)*100))
cat(sprintf("Assault zeros:      %d (%.1f%%)\n",
            sum(data$Assault == 0), mean(data$Assault == 0)*100))

# ============================================================
# COMPLETE CASE DATASET
# ============================================================
data_clean <- data %>%
  filter(
    !is.na(RUCC_2023),
    !is.na(EP_POV150), !is.na(EP_UNEMP), !is.na(EP_AGE65), !is.na(EP_NOHSDP),
    !is.na(POP_SQMI),  !is.na(EP_AFAM),  !is.na(EP_HISP),  !is.na(EP_MALE),
    !is.na(HeatIndex_2022_Mean), !is.na(Green_2022),
    !is.na(eal_score), !is.na(Total_Disasters_2005_2023),
    E_TOTPOP > 0
  )
cat("N complete cases:", nrow(data_clean), "\n")

# ============================================================
# SVI predictor block (shared across all Model 1s)
# ============================================================
svi <- "EP_POV150 + EP_UNEMP + EP_AGE65 + EP_NOHSDP + POP_SQMI + EP_AFAM + EP_HISP + EP_MALE"

# ============================================================
# FUNCTION: 4 NB models
# ============================================================
run_4_nb_models <- function(dv, offset_var = NULL, dataset, outcome_label) {
  
  off <- if (!is.null(offset_var)) paste0(" + offset(log(", offset_var, "))") else ""
  
  # Model formulas per Dr. Han's structure
  f1 <- as.formula(paste0(dv, " ~ ", svi, off))
  f2 <- as.formula(paste0(dv, " ~ ", svi, " + RUCC_2023", off))
  f3 <- as.formula(paste0(dv, " ~ ", svi, " + RUCC_2023 + Total_Disasters_2005_2023 + eal_score", off))
  f4 <- as.formula(paste0(dv, " ~ ", svi, " + RUCC_2023 + Total_Disasters_2005_2023 + eal_score + Green_2022 + HeatIndex_2022_Mean", off))
  
  formulas    <- list(f1, f2, f3, f4)
  model_names <- c(
    "Model 1: SVI Predictors",
    "Model 2: Model 1 + RUCC",
    "Model 3: Model 2 + Disaster Counts + EAL Score",
    "Model 4: Model 3 + Green Score + Heat Index"
  )
  
  results <- list()
  
  for (i in 1:4) {
    cat(sprintf("\n##########################################################\n"))
    cat(sprintf("OUTCOME: %s | %s\n", outcome_label, model_names[i]))
    cat(sprintf("##########################################################\n"))
    
    m <- glm.nb(formulas[[i]], data = dataset)
    print(summary(m))
    cat("Theta:", round(m$theta, 4),
        "| SE:", round(m$SE.theta, 4),
        "| Dispersion ratio:", round(deviance(m)/df.residual(m), 4), "\n")
    
    rr <- exp(cbind(RR = coef(m), confint(m)))
    cat("Rate Ratios (RR) and 95% CI:\n")
    print(round(rr, 7))
    
    cat(sprintf(
      "Fit: NullDev=%.1f(df=%d) | ResidDev=%.1f(df=%d) | AIC=%.1f | Disp=%.4f | N=%d\n",
      m$null.deviance, m$df.null,
      deviance(m), df.residual(m),
      AIC(m),
      deviance(m)/df.residual(m),
      nrow(dataset)
    ))
    
    results[[i]] <- m
  }
  
  return(results)
}

# ============================================================
# ALL THREE OUTCOMES
# ============================================================

cat("\n\n==========================================================\n")
cat("OUTCOME A: DV = total_events (Firearm Total Events)\n")
cat("==========================================================\n")
mod_A <- run_4_nb_models(
  dv           = "total_events",
  offset_var   = NULL,
  dataset      = data_clean,
  outcome_label= "total_events"
)

cat("\n\n==========================================================\n")
cat("OUTCOME B: DV = Homicide Rate (with population offset)\n")
cat("==========================================================\n")
mod_B <- run_4_nb_models(
  dv           = "Homicide",
  offset_var   = "E_TOTPOP",
  dataset      = data_clean,
  outcome_label= "Homicide Rate"
)

cat("\n\n==========================================================\n")
cat("OUTCOME C: DV = Assault Rate (with population offset)\n")
cat("==========================================================\n")
mod_C <- run_4_nb_models(
  dv           = "Assault",
  offset_var   = "E_TOTPOP",
  dataset      = data_clean,
  outcome_label= "Assault Rate"
)

# ============================================================
# COMPREHENSIVE MODEL FIT COMPARISON — ALL 12 MODELS
# ============================================================
cat("\n\n==========================================================\n")
cat("COMPREHENSIVE MODEL FIT COMPARISON — ALL 12 NB MODELS\n")
cat("==========================================================\n")
cat(sprintf("%-52s %10s %12s %10s %6s\n",
            "Model", "AIC", "ResidDev(df)", "Disp", "N"))
cat(strrep("-", 95), "\n")

entries <- list(
  list(label="A1: total_events  | Model 1 (SVI only)",          m=mod_A[[1]]),
  list(label="A2: total_events  | Model 2 (+RUCC)",             m=mod_A[[2]]),
  list(label="A3: total_events  | Model 3 (+Disasters+EAL)",    m=mod_A[[3]]),
  list(label="A4: total_events  | Model 4 (+Green+Heat)",       m=mod_A[[4]]),
  list(label="B1: Homicide Rate | Model 1 (SVI only)",          m=mod_B[[1]]),
  list(label="B2: Homicide Rate | Model 2 (+RUCC)",             m=mod_B[[2]]),
  list(label="B3: Homicide Rate | Model 3 (+Disasters+EAL)",    m=mod_B[[3]]),
  list(label="B4: Homicide Rate | Model 4 (+Green+Heat)",       m=mod_B[[4]]),
  list(label="C1: Assault Rate  | Model 1 (SVI only)",          m=mod_C[[1]]),
  list(label="C2: Assault Rate  | Model 2 (+RUCC)",             m=mod_C[[2]]),
  list(label="C3: Assault Rate  | Model 3 (+Disasters+EAL)",    m=mod_C[[3]]),
  list(label="C4: Assault Rate  | Model 4 (+Green+Heat)",       m=mod_C[[4]])
)

for (e in entries) {
  m <- e$m
  cat(sprintf("%-52s %10.1f %8.1f(%d) %10.4f %6d\n",
              e$label,
              AIC(m),
              deviance(m),
              df.residual(m),
              deviance(m)/df.residual(m),
              nrow(data_clean)
  ))
}






library(tidyverse)
library(MASS)
library(broom)
library(jsonlite)
library(readxl)
select <- dplyr::select

data <- read_csv("D:/GVA/Excel Files/Gulf_Final_Merged_MASTER.csv")
data <- data %>% mutate(POP_SQMI = E_TOTPOP / AREA_SQMI)

workspace_dfs <- Filter(function(x) is.data.frame(get(x)), ls())
for (nm in workspace_dfs) { df <- get(nm); if ("EP_MALE" %in% names(df)) cat("Found EP_MALE in:", nm, "\n") }
data <- data %>% left_join(sex_df %>% dplyr::select(FIPS, EP_MALE), by = "FIPS")
cat("EP_MALE NAs:", sum(is.na(data$EP_MALE)), "\n")
cat(sprintf("EP_MALE range: %.2f – %.2f (mean=%.2f)\n",
            min(data$EP_MALE, na.rm=TRUE),
            max(data$EP_MALE, na.rm=TRUE),
            mean(data$EP_MALE, na.rm=TRUE)))

# ── STEP 1: Update RUCC ──────────────────────────────────────────────────────
rucc_new <- read_excel("C:/Desktop/Research/GVA/Excel Files/Ruralurbancontinuumcodes2023 (1).xlsx")
rucc_new_clean <- rucc_new %>%
  mutate(FIPS = str_pad(as.character(FIPS), width = 5, pad = "0")) %>%
  dplyr::select(FIPS, RUCC_2023)

data <- data %>%
  dplyr::select(-RUCC_2023) %>%
  left_join(rucc_new_clean, by = "FIPS")

cat("RUCC NAs after update:", sum(is.na(data$RUCC_2023)), "\n")
cat("Urban (RUCC 1-3):", sum(data$RUCC_2023 <= 3, na.rm=TRUE), "\n")
cat("Rural  (RUCC 4-9):", sum(data$RUCC_2023 >= 4, na.rm=TRUE), "\n")

# ── STEP 2: Complete case dataset ────────────────────────────────────────────
data_clean2 <- data %>%
  filter(
    !is.na(Total_Disasters_2005_2023),
    !is.na(RUCC_2023),
    !is.na(EP_POV150), !is.na(EP_UNEMP), !is.na(EP_AGE65), !is.na(EP_NOHSDP),
    !is.na(POP_SQMI),  !is.na(EP_AFAM),  !is.na(EP_HISP),  !is.na(EP_MALE),
    !is.na(HeatIndex_2022_Mean), !is.na(Green_2022),
    E_TOTPOP > 0
  )
cat("\nN complete cases:", nrow(data_clean2), "\n")
cat(sprintf("Disaster zeros: %d (%.1f%%)\n",
            sum(data_clean2$Total_Disasters_2005_2023 == 0),
            mean(data_clean2$Total_Disasters_2005_2023 == 0) * 100))

# ── STEP 3: Overdispersion check ─────────────────────────────────────────────
pois_check <- glm(Total_Disasters_2005_2023 ~ EP_POV150 + EP_UNEMP + EP_AGE65 +
                    EP_NOHSDP + POP_SQMI + EP_AFAM + EP_HISP + EP_MALE,
                  family = "poisson", data = data_clean2)
cat("\nPoisson dispersion ratio:", round(deviance(pois_check)/df.residual(pois_check), 2),
    "— NB appropriate if > 2.0\n\n")

# ── MODEL 1: Eight SVI Predictors ────────────────────────────────────────────
model1 <- glm.nb(
  Total_Disasters_2005_2023 ~
    EP_POV150 + EP_UNEMP + EP_AGE65 + EP_NOHSDP + POP_SQMI +
    EP_AFAM + EP_HISP + EP_MALE,
  data = data_clean2)

cat("##########################################################\n")
cat("MODEL 1: SVI Predictors\n")
cat("##########################################################\n")
print(summary(model1))
cat("Theta:", round(model1$theta,4), "| SE:", round(model1$SE.theta,4),
    "| Disp:", round(deviance(model1)/df.residual(model1),4), "\n")
print(round(exp(cbind(RR=coef(model1), confint(model1))), 7))
cat(sprintf("Fit: NullDev=%.1f(df=%d) | ResidDev=%.1f(df=%d) | AIC=%.1f | Disp=%.4f | N=%d\n\n",
            model1$null.deviance, model1$df.null, deviance(model1), df.residual(model1),
            AIC(model1), deviance(model1)/df.residual(model1), nrow(data_clean2)))

# ── MODEL 2: Model 1 + RUCC ───────────────────────────────────────────────────
model2 <- glm.nb(
  Total_Disasters_2005_2023 ~
    EP_POV150 + EP_UNEMP + EP_AGE65 + EP_NOHSDP + POP_SQMI +
    EP_AFAM + EP_HISP + EP_MALE + RUCC_2023,
  data = data_clean2)

cat("##########################################################\n")
cat("MODEL 2: Model 1 + RUCC\n")
cat("##########################################################\n")
print(summary(model2))
cat("Theta:", round(model2$theta,4), "| SE:", round(model2$SE.theta,4),
    "| Disp:", round(deviance(model2)/df.residual(model2),4), "\n")
print(round(exp(cbind(RR=coef(model2), confint(model2))), 7))
cat(sprintf("Fit: NullDev=%.1f(df=%d) | ResidDev=%.1f(df=%d) | AIC=%.1f | Disp=%.4f | N=%d\n\n",
            model2$null.deviance, model2$df.null, deviance(model2), df.residual(model2),
            AIC(model2), deviance(model2)/df.residual(model2), nrow(data_clean2)))

# ── MODEL 3: Model 2 + Green + Heat ──────────────────────────────────────────
model3 <- glm.nb(
  Total_Disasters_2005_2023 ~
    EP_POV150 + EP_UNEMP + EP_AGE65 + EP_NOHSDP + POP_SQMI +
    EP_AFAM + EP_HISP + EP_MALE + RUCC_2023 +
    Green_2022 + HeatIndex_2022_Mean,
  data = data_clean2)

cat("##########################################################\n")
cat("MODEL 3: Model 2 + Green & Heat\n")
cat("##########################################################\n")
print(summary(model3))
cat("Theta:", round(model3$theta,4), "| SE:", round(model3$SE.theta,4),
    "| Disp:", round(deviance(model3)/df.residual(model3),4), "\n")
print(round(exp(cbind(RR=coef(model3), confint(model3))), 7))
cat(sprintf("Fit: NullDev=%.1f(df=%d) | ResidDev=%.1f(df=%d) | AIC=%.1f | Disp=%.4f | N=%d\n\n",
            model3$null.deviance, model3$df.null, deviance(model3), df.residual(model3),
            AIC(model3), deviance(model3)/df.residual(model3), nrow(data_clean2)))

# ── FIT COMPARISON ────────────────────────────────────────────────────────────
cat("==========================================================\n")
cat("MODEL FIT COMPARISON: All 3 NB Models\n")
cat("==========================================================\n")
cat(sprintf("%-45s %10s %14s %10s %6s\n", "Model","AIC","ResidDev(df)","Disp","N"))
cat(strrep("-", 90), "\n")
for (e in list(
  list(label="Model 1: SVI Predictors",           m=model1),
  list(label="Model 2: Model 1 + RUCC",            m=model2),
  list(label="Model 3: Model 2 + Green & Heat",    m=model3)
)) {
  cat(sprintf("%-45s %10.1f %8.1f(%d) %10.4f %6d\n",
              e$label, AIC(e$m), deviance(e$m), df.residual(e$m),
              deviance(e$m)/df.residual(e$m), nrow(data_clean2)))
}



model3_pois <- glm(
  Total_Disasters_2005_2023 ~
    EP_POV150 + EP_UNEMP + EP_AGE65 + EP_NOHSDP + POP_SQMI +
    EP_AFAM + EP_HISP + EP_MALE + RUCC_2023 +
    Green_2022 + HeatIndex_2022_Mean,
  family = "poisson", data = data_clean2)

print(summary(model3_pois))
cat("Dispersion:", round(deviance(model3_pois)/df.residual(model3_pois), 4), "\n")
print(round(exp(cbind(RR=coef(model3_pois), confint(model3_pois))), 7))