Copyright 2024 Brandon C. Alston
# Code for Quant Finance Algo Trading Mechanisms

## Summary of Code/Folders
### SimpleTradingBots
  - Basic ML and automatred trading bots implemented in `Python`
### Option Pricing Fundamentals
  - Derivation and foundations of the Black-Scholes-Morten model for risk-neutral asset pricing
  - The BSM PDE describes evolution of all derivatives whose payoff is a function on a single underlying asset following geometric Brownian motion and time
### Monte Carlo Methods 
  - MCM methods valuation of European and path-dependent derivatives. 
  - Various random number generators for pseudorandom, quasi-random (deterministic), Sobol, and Faure sequences. 
  - Variance reduction techniques using control variates and antithetics
    - European call option using Faure sequence for variance reduction (w/ and w/out Sobol algorithm)
    - Option pricing with anithetic, dynamic peplication
  - Hedge Control Variates
    - Monte carlo variance reduced simulation with greeks
    - Spread option valuation
  - Path-Dependent Valuation 
    - Asian option valuation with MC (geometric, control variate)
  - Brownian Bridge Technique
    - Long time horizon MC simulation
  - Jump-Diffusion Process and Constant Elasticity of Variance Diffusion Model
    - Option price valuation
### Tree Models
  - Binomial tree model for pricing European and American equity options.
  - Cox-Ross-Rubinstein Binomial Tree Model
  - Jarrow-Rudd Tree Model
  //- Two-state discrete approximation to continuous GBM
    // - The mean and variance of the binomial model match mean and variance of the lognormal distribution underlying GBM
  // - Adapted to incorporate time-varying volatility, price path-dependent options, and price derivatives depending on more than one asset with two-variable binomial trees