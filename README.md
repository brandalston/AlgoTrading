Copyright 2024 Brandon C. Alston
# Code for Quant Finance Algo Trading Mechanisms

## Summary of Code/Folders
- ### SimpleTradingBots
  - Basic ML and automatred trading bots implemented in `Python`
- ### BSM, Pricing Fundamentals
  - Derivation and foundations of the Black-Scholes-Morten model for risk-neutral asset pricing
  - The BSM PDE describes evolution of all derivatives whose payoff is a function on a single underlying asset following geometric Brownian motion and time
- ### Monte Carlo Methods 
  - MCM methods valuation of European and path-dependent derivatives. 
  - Various random number generators for pseudorandom, quasi-random (deterministic), Sobol, and Faure sequences. 
  - Variance reduction techniques using control variates and antithetics
- ### Binomial Tree
  - Binomial tree model for pricing European and American equity options. 
  - Two-state discrete approximation to continuous GBM
    - The mean and variance of the binomial model match mean and variance of the lognormal distribution underlying GBM
  - Adapted to incorporate time-varying volatility, price path-dependent options, and price derivatives depending on more than one asset with two-variable binomial trees
  - 