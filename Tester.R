library(TwoStepSDFM)

# Simulate a DGP using simFM
T <- 600 # Number of observations
N <- 100 # Number of variabes
R <- 2 # Number of factors
Sigma_epsilon <- diag(1, R) # Variance-covariance matrix of the transition errors
Lambda <- matrix(rnorm(N * R), N, R) # Factor loadings matrix
mu_xi <- rep(0, N) # Mean of the measurement error
Sigma_xi <- diag(1, N) # Variance-covariance matrix of the measurement error
Phi <- cbind(diag(0.5, R), -diag(0.25, R)) # Factor VAR coefficient matrix
P <- 2 # Order of the factor VAR process
quarterfy <- FALSE # Indicating whether or not some of the variables should be aggregated to quarterly observations (i.e., quarterfied)
corr <- FALSE # Indicating whether or not the measurement error should be internatlly cross-crossectionally correlated
beta_param <- Inf # Beta parameter governing the degree of correlation of the measurement error
m <- 0 # Ratio of monthly predictors ought to be quarterfied
seed <- 16022024 # Seed
set.seed(seed)
burn_in <- 999 # Burn-in period
rescale <- TRUE # Indicating whether the variance of the measurement error should be scaled according to the variance of the common-component 

# Draw from an exact factor model
FM_exact <- simFM(T = T, N = N, R = R, Lambda = Lambda, mu_xi = mu_xi, Sigma_xi = Sigma_xi,
                  Sigma_epsilon = Sigma_epsilon, Phi = Phi, P = P, quarterfy = FALSE, m = 0,
                  corr = corr, beta_param = beta_param, seed = seed, burn_in = burn_in, rescale = rescale)

head(t(FM_exact$F), 20) # First 20 observations of the factors
FM_exact$Sigma_xi[1:10, 1:10] # Upper-left 10x10 block of the measurement variance covariance matrix

# Draw from an approximate factor model
FM_approx <- simFM(T = T, N = N, R = R, Lambda = Lambda, mu_xi = mu_xi, Sigma_xi = Sigma_xi,
                   Sigma_epsilon = Sigma_epsilon, Phi = Phi, P = P, quarterfy = FALSE, m = 0,
                   corr = TRUE, beta_param = 1, seed = seed, burn_in = burn_in, rescale = rescale)

head(t(FM_approx$F), 20) # First 20 observations of the factors
FM_approx$Sigma_xi[1:10, 1:10] # Upper-left 10x10 block of the emasurement variance covariance matrix

# Draw from an approximate mixed frequency factor model
quarterfy <- TRUE
m <- 0.01
FM_mf <- simFM(T = 303, N = N, R = R, Lambda = Lambda, mu_xi = mu_xi, Sigma_xi = Sigma_xi,
               Sigma_epsilon = Sigma_epsilon, Phi = Phi, P = P, quarterfy = quarterfy, m = m,
               corr = TRUE, beta_param = 1, seed = 123456, burn_in = burn_in + 1, rescale = rescale)

head(t(FM_mf$X)[, 1:2]) # First 20 observations of the factors covariance matrix
# tail(t(FM_mf$X), 4) # First 20 observations of the factors covariance matrix

