###############################################
# Set all input parameters at the start
###############################################

# Time horizon for evaluating overall survival (OS)
time_t <- 3

# Survival proportions at time_t
survival_control <- 0.7  # Control group survival at 3 years
survival_treatment <- 0.6 # Treatment group survival at 3 years (e.g., NI margin of 10%)

# Exponential dropout rate
eta <- 0

# Hypothesized hazard ratios
hr_alternative <- 1  # For non-inferiority, H1 HR=1 means no difference
# Null hazard ratio (set from calculated hazard ratio based on survival rates)
# Will compute below after calculating hazard rates.

# Type I and Type II error
alpha <- 0.025  # One-sided Type I error
beta <- 0.2     # Type II error (80% power)

# Randomization ratio
ratio <- 1

# Enrollment and follow-up durations
T_enroll <- 2      # 2 years enrollment
minfup <- 2         # 2 years minimum follow-up after enrollment ends
total_duration <- T_enroll + minfup

# Uniform enrollment rate
gamma <- 1

###############################################
# Compute hazard rates and hazard ratio
###############################################

# Control hazard rate
control_haz_rate <- -log(survival_control) / time_t

# Treatment hazard rate
treatment_haz_rate <- -log(survival_treatment) / time_t

# Hazard ratio (treatment/control)
Hazard_Ratio <- treatment_haz_rate / control_haz_rate

# Median survival time (control)
MST_control <- log(2) / control_haz_rate

###############################################
# Run sample size calculation using nSurv
###############################################

library(gsDesign)

ss <- nSurv(
  R = T_enroll,
  gamma = gamma,
  eta = eta,
  minfup = minfup,
  T = total_duration,
  lambdaC = log(2) / MST_control, ## Hazard Rate of control / reference Treatment
  hr = hr_alternative,
  hr0 = Hazard_Ratio,
  beta = beta,
  alpha = alpha,
  sided = 1
)

print(ss)
