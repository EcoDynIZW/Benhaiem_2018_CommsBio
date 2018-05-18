# SIR-Matrix-model

The code **1_Model construction** loads the input values for the *3 epidemic periods* and builds the matrix model. It is structured as follows:

  * a) Loading input values (Multi-Event-Capture-Mark-Recapture (MECMR) parameter estimates)
  * b) Building the submatrices: Survival, State transition (demographic, social, infection) and Fecundity submatrices
  * c) Assembling the submatrices into the meta-matrix
  * d) Building the next generation matrix for the estimation of R0

This code  **2_Model analysis** presents the asymptotic analysis of the matrix model. For each epidemic period, we first load the input files containing the MECMR parameter estimates for that period and then the R script to construct the matrix model.
The code code is structured as follows:

  * a) Checking assumptions (irreducibility and ergodicity) 
  * b) Asymptotic analyses (population's growth rate, R0, stable stage distribution and reproductive values)
     * i) Pre-epidemic
     * ii) Epidemic
     * iii) Post-epidemic
     * iv) Short-term population dynamics (Fig 4)

The code **3_Model_Analysis_stochastic** presents the stochastic analysis of the matrix model. To calculate standard deviations of (λ) and R0 (Figure 1), confidence intervals for the sensitivity analysis of (λ) and R0 (Figures 2,3) and to describe changes in population size over time while accounting for parameter uncertainty (Figure 5), we used Monte Carlo iterations.

This code is structured as follows:

  * 1 Monte Carlo iterations to calcule the mean + SD of the population’s growth rate (Fig 1A), R0 (Fig 1B) and population abundance (Fig 5).
      * a) Plotting mean + SD of population’s growth rate and R0 (Fig 1)
      * b) Plotting changes in population abundance - complete model and “no rank” model (Fig 5)
  * 2 Sensitivity analysis of the population’s growth rate (λ)(Fig 2)

  * 3 Sensitivity analysis of R0 (Fig 3)
  
  Note that parameter and submatrix names may differ between main text and R codes.

