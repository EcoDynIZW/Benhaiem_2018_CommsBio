# SIR-Matrix-model

The code **1_Model construction** loads the input values for the *3 epidemic periods* and builds the matrix model. It is structured as follows:

  * a) Loading input values (Multi-Event-Capture-Mark-Recapture (MECMR) parameter estimates)
  * b) Building the submatrices: Survival, State transition (demographic, social, infection) and Fecundity submatrices
  * c) Assembling the submatrices into the meta-matrix
  * d) Building the next generation matrix for the estimation of R0


Note that parameter and submatrix names may differ between main text and R codes.

This code presents the asymptotic analysis of the matrix model. For each epidemic period, we first load the input files containing the MECMR parameter estimates for that period and then the R script to construct the matrix model.
  
The code code **2_Model analysis** is structured as follows:

  * a) Checking assumptions (irreducibility and ergodicity) 
  * b) Asymptotic analyses (population's growth rate, R0, stable stage distribution and reproductive values)
    
    i) Pre-epidemic
    ii) Epidemic
    iii) Post-epidemic
    iv) Short-term population dynamics (Fig 4)

