# Gulf Coast Gun Violence and Environmental Determinants Analysis Pipeline

This repository hosts a production-grade spatial epidemiology and biostatistics pipeline engineered to model, analyze, and quantify the impacts of the Social Vulnerability Index (SVI), built environments, and climate anomalies on gun violence patterns across counties in the Gulf Coast region. 

The pipeline handles end-to-end data processing, dynamic Federal API integration, overdispersion validation, and structured hierarchical multi-model regressions using Generalized Linear Models (GLMs).

---

## Core Methodological Evolution

The codebase documents a rigorous analytical pipeline developed through four foundational phases to support peer-reviewed research and executive clinical reporting:

* Phase 1: Distribution and Overdispersion Diagnostics
  Establishes baseline models testing the total event counts (total_events) and population-offset rates (Homicide) via traditional Poisson regression. It evaluates log-likelihood and variance constraints, revealing massive overdispersion ratios (far above the standard threshold of 2.0, up to 200.6 for total counts and 35.9 for rates). This quantifies the statistical necessity of pivoting from standard Poisson constraints to Negative Binomial (NB) frameworks to manage spatial-structural clustering.

* Phase 2: Structural Integrity and Structural Zero Profiling
  Restructures Rural-Urban Continuum Codes (RUCC_2023) via dynamic spatial zero-padding to ensure 100 percent data fidelity on 5-digit county FIPS codes. It handles complete-case mutations and evaluates the density of structurally true zero counts (Y = 0) across endpoints to establish clear benchmarks for model eligibility.

* Phase 3: Matrix Expansion
  Scales the pipeline to include a third critical vector: non-fatal aggressive gun violence metrics (Assault), expanding tracking across both fatal and non-fatal injury patterns.

* Phase 4: Production Hierarchical Nested Modeling
  Implements the final analytic framework requested by senior research leadership. The script programmatically executes a sequential nested modeling design across all outcomes to evaluate the variance accounted for by each additional environmental or socio-demographic block.

---

## Hierarchical Model Architecture

For each distinct clinical or behavioral gun violence endpoint, the pipeline builds and cross-compares four nested Negative Binomial models to systematically isolate structural risk factors:

* Model 1: SVI Baseline
  Includes Poverty, Unemployment, Age, Minorities, and Population Density.
* Model 2: Built Context
  Includes Model 1 Predictors plus Rural-Urban Continuum Codes (RUCC).
* Model 3: Macro-Ecology
  Includes Model 2 Predictors plus FEMA Disasters and the National Risk Index (EAL) Score.
* Model 4: Fully Adjusted Model
  Includes Model 3 Predictors plus Canopy Density (Green Score) and the Heat Index.

---

## Variables and Data Dictionary

### 1. Exposure and Covariate Blocks

* Socio-Demographic / SVI Matrices:
  * EP_POV150: Percentage of population living below 150 percent of the poverty threshold.
  * EP_UNEMP: Civilian unemployment rate.
  * EP_AGE65: Percentage of individuals aged 65 or older.
  * EP_NOHSDP: Percentage of individuals without a high school diploma.
  * POP_SQMI: Spatial population density per square mile (calculated as E_TOTPOP divided by AREA_SQMI).
  * EP_AFAM and EP_HISP: Percentages of Black/African American and Hispanic populations.
  * EP_MALE: Automated upstream population tracking of male gender concentration by percentage, ingested via live US Census ACS 5-Year API calls.

* Geographical Stratification:
  * RUCC_2023: USDA Rural-Urban Continuum Codes tracking metropolitan to deep rural axes.

* Macro-Climate and Structural Vulnerability:
  * Total_Disasters_2005_2023: Institutional FEMA Major Disaster declaration tracking.
  * eal_score: Expected Annual Loss index generated via the National Risk Index (NRI).

* Built and Micro-Environmental Context:
  * Green_2022: National canopy density / greenness index tracking urban design mitigation.
  * HeatIndex_2022_Mean: Mean heat index calculations identifying long-term localized thermal stress.

### 2. Primary Endpoints (Dependent Variables)
* total_events: Cumulative incident volume of firearm-related events.
* Homicide: Intentionally inflicted fatal firearm violence rates (modeled with a population offset).
* Assault: Non-fatal firearm-driven assault rates (modeled with a population offset).

---

## Environment Setup and Package Setup

This analytical pipeline is built to run on R version 4.0.0 or higher. Dependencies must be acquired using CRAN management utilities prior to runtime execution:

```R
install.packages(c("tidyverse", "MASS", "broom", "jsonlite", "readxl"))
