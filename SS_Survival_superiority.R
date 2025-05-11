###############################################
# Sample-Size Calculation for Superiority Log-Rank Test
# Bulky (treatment) vs Non-Bulky (control) RFS at 4 years
###############################################

# 1. Analysis time horizon (years) # This is when u evaluate OS EFS RFS at ? Yr
time_t <- 4  
#    We use 4-year RFS to derive exponential hazard rates.

# 2. Four-year RFS for each arm
#    CONTROL arm: Non-Bulky disease group: 94.4% RFS at 4 years
survival_control <- 0.944  

#    TREATMENT arm: Bulky disease group: 80.5% RFS at 4 years
survival_treatment <- 0.805  

# 3. Exponential dropout rate (per year)
#    Set to zero if dropouts are negligible or not modeled
eta <- 0  

# 4. Error rates for two-sided superiority test
alpha <- 0.05  # two-sided Type I error
beta  <- 0.20  # Type II error: 80% power
sided <- 2     # two-sided test

# 5. Observed allocation from prior study
#    Number of subjects in each arm (for ratio only)
n_treatment <- 98   # Bulky arm
n_control   <- 237  # Non-Bulky arm

# 6. Allocation ratio (treatment : control)
#    nSurv() expects ratio = n_T / n_C
ratio <- n_treatment / n_control  

# 7. Accrual (enrollment) and follow-up durations
T_enroll       <- 4  # years of uniform accrual
minfup         <- 1  # years of additional follow-up after accrual
total_duration <- T_enroll + minfup  # total study window

# 8. Uniform accrual pattern over T_enroll
#    Split into 4 one-year segments of equal length and weight
R     <- rep(1, 4)  # segment lengths (years)
gamma <- rep(1, 4)  # weights (uniform accrual)

###############################################
# Compute hazard rates and hazard ratio
###############################################

# 9. Control arm hazard rate (λ_C) from 4-year RFS
#    λ = –log(S(t)) / t under exponential assumption
control_haz_rate <- -log(survival_control) / time_t  

# 10. Treatment arm hazard rate (λ_T) likewise
treatment_haz_rate <- -log(survival_treatment) / time_t

# 11. Observed hazard ratio under alternative hypothesis
hr_alternative <- treatment_haz_rate / control_haz_rate  

# 12. Null hazard ratio (H0: no difference)
hr0 <- 1  

# 13. (Optional) Median survival of control arm for reference
MST_control <- log(2) / control_haz_rate  

###############################################
# Run sample-size calculation using gsDesign::nSurv
###############################################

# Load the gsDesign package
library(gsDesign)

# Perform the calculation
ss <- nSurv(
  # Control (Non-Bulky) hazard rate
  lambdaC = control_haz_rate,
  # Hazard ratio under H1 (Bulky vs Non-Bulky)
  hr      = hr_alternative,
  # Null hazard ratio under H0
  hr0     = hr0,
  # Dropout rate
  eta     = eta,
  # Accrual pattern: segment lengths and weights
  R       = R,
  gamma   = gamma,
  # Minimum follow-up after accrual
  minfup  = minfup,
  # Total study duration (accrual + follow-up)
  T       = total_duration,
  # Allocation ratio (Bulky:Non-Bulky)
  ratio   = ratio,
  # Type I & II error settings
  alpha   = alpha,
  beta    = beta,
  sided   = sided
)

# Print the detailed sample-size output
print(ss)
