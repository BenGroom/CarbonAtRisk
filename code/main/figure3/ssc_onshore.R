# Storage Security Calculator — Unregulated onshore scenario
# From Alcalde et al. (2018), adapted for standalone use in 

source("code/main/figure3/ssc_common.R")

set.seed(2100)

# Monte Carlo — ONSHORE parameter distributions
SSCMC <- function(n = 5) {
  onshore_sampler <- function() {
    list(
      CO2target        = 1.2e10,
      injectperWell    = rnorm(1, 0.75e6, 0.00415e6),
      InjectionPeriod  = 30,
      meanPlumeArea    = rlnorm(1, -0.7595, 0.1763),
      ActiveWellFreq   = rlnorm(1, -2.89, 0.7),
      SlowLeakInjector = rnorm(1, 158.5, 5.2),
      MinorBlowFreq    = runif(1, 0.062, 0.0762),
      MinorBlowout     = exp(rlnorm(1, 1.27, 0.0397)),
      MajorBlowFreq    = rnorm(1, 1.35e-4, 4.4e-5),
      CO2MajorBlowout  = exp(rlnorm(1, 2.57, 0.045)),
      KnownWellDensity = runif(1, 2.25, 2.75),
      wellUnderEst     = runif(1, 1.1, 2),
      UnPlugWells      = 0.3,
      DegradWells      = rlnorm(1, -2.89, 0.7),
      IntactHighRate   = runif(1, 0.049, 0.0595),
      IntactDegrad     = 0,
      PlugBlowoutYear  = rlnorm(1, -8.6125, 0.23),
      BlowoutWellYear  = runif(1, 1e-5, 1e-4),
      CO2largeBlowout  = rlnorm(1, 13.4, 0.35),
      CO2degraded      = runif(1, 270, 330),
      CO2intactHigh    = runif(1, 207, 253),
      CO2intactLow     = runif(1, 0.0036, 0.0044),
      NatLeakRate      = rlnorm(1, 0.693, 0.37),
      A                = rtriangle(1, 3, 53, 12),
      B                = runif(1, 0.0143, 0.5),
      res_sat          = rnorm(1, 0.58, 0.0286),
      migration        = 1
    )
  }
  SSCMC_generic(n, onshore_sampler)
}
