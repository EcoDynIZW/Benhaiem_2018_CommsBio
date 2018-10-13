
### This file has 4 functions meant to calculate the amount of uncertainty associated with the different population indicators 
# such as the population growth rate, the basic reproduction number, sensitivity values...
# normdist allows to draw the regression coefficients from a normal distribution,  
# paramdelta allows to convert the regression coefficients which are on a logit scale into a probability bounded from 0 to 1.  
# build_matrix generates for each Monte Carlo iteration the projection matrix as well as the next generation matrix necessary to caluclate the R0. It also calls the function vitalsense used to get the second order sensitivity and elasticity values for each biological parameter
# sens_elas_num is meant to calculate the sensitivity of R0 to variations in each biological parameter (demgraphic, social or infection) 


normdist<-function(value, MCiter)
{
  phi<- value[1]
  var<- value[2]
  phiMCiter<-rnorm(MCiter, phi, var) 
  return(phiMCiter)
}


paramdelta<-function(t, simu)
{
  if(t != "post-epidem2")
  {
    bioparam<-NULL
    for(z in unique(simu$parameter))
    {
      for(y in unique(simu$states))
      {
        
        Cintercept <- subset(simu, Time=="post-epidem2" & states==y & parameter== z)
        
        Ctemp <- subset(simu, Time==t & states=="all" & parameter == z)
        Totalcoeff<-rbind(Cintercept, Ctemp)
        
        if(dim(Cintercept)[1] > 0 & dim(Ctemp)[1] > 0)
        {
          av<-t(exp(Totalcoeff[1,c(1, 6: (MCiter+5))] + Totalcoeff[2,c(1, 6: (MCiter+5))])/(1+exp(Totalcoeff[1,c(1, 6: (MCiter+5))]+Totalcoeff[2,c(1, 6: (MCiter+5))])))
          colnames(av)<-Totalcoeff[1,3]
          bioparam<-cbind(bioparam, av)  
        }
        
        Cothers <- subset(simu, Time=="allperiods" & states==y & parameter== z)
        if(dim(Cothers)[1] > 0) 
        {
          otherparam<-t(exp(Cothers[1,c(1, 6: (MCiter+5))])/(1+exp(Cothers[1,c(1, 6:(MCiter+5))])))
          colnames(otherparam)<-Cothers[1,3]
          bioparam<-cbind(bioparam, otherparam)
        }
        
      }
    }
  }
  else
  {
    bioparam<-NULL
    for(z in unique(simu$parameter))
    {
      for(y in unique(simu$states))
      {
        
        Cintercept <- subset(simu, Time=="post-epidem2" & states==y & parameter== z)
        
        Totalcoeff<-Cintercept
        
        if(dim(Cintercept)[1] > 0)
        {
          av<-t(exp(Totalcoeff[1,c(1, 6: (MCiter+5))])/(1+exp(Totalcoeff[1,c(1, 6: (MCiter+5))])))
          colnames(av)<-Totalcoeff[1,3]
          bioparam<-cbind(bioparam, av)  
        }
        
        Cothers <- subset(simu, Time=="allperiods" & states==y & parameter== z)
        if(dim(Cothers)[1] > 0) 
        {
          otherparam<-t(exp(Cothers[1,c(1, 6: (MCiter+5))])/(1+exp(Cothers[1,c(1, 6:(MCiter+5))])))
          colnames(otherparam)<-Cothers[1,3]
          bioparam<-cbind(bioparam, otherparam)
        }
        
      }
    }
  }
  return (bioparam)
}



param<-data$EST 
var<-data$SE 

#########----------------------------------- STEP 1 : simulating 1000 values of each param & back-transforming them

# creates a table with all values from the Monte Carlo simulation (column) for each demographic, epidemiological, and social parameter (rows)

global<-cbind(as.numeric(param), as.numeric(var)) 
MC<-apply(global,1, function (x) normdist(x, MCiter)) 

MCtrans<-t(MC) 

simu<-cbind(data, MCtrans) 

theta<-paramdelta(period, simu) 
theta<-as.data.frame(theta) 

theta$ls <- rep(1.53,(MCiter+1))

theta$sr<-rep(0.52, (MCiter+1))
# contains the 1000 values now backtransformed, focusing only on the epidemic period 
# each column starts with the name of the parameter (C.L.S, C.L.I...)

#########----------------------------------- STEP 2 : constructing the MATRIX MODEL to get R0 

# 
#--------------------------------------------- We build our matrix with the function build_matrix
#  -------- This function creates the NEXT GENERATION MATRIX 


build_matrix <- function(theta, i, checkNodisease) {
  
  
  if(checkNodisease == TRUE)
  {
    theta[i,"C.L"] <- 1 # here we retransform it as "real" beta value
    theta[i,"C.H"] <- 1 # 
    
    theta[i,"SA&NB&B.L"] <- 1
    theta[i,"SA&NB&B.H"] <- 1
    
    theta[i,"SA&NB&B.L"] <- 1
    theta[i,"SA&NB&B.H"] <- 1
  }
  
  
  NStages<-22
  param <- list(
  phiCLS = theta[i,"C.L.S"], # cub low ranked & susceptible
  phiCLI = theta[i,"C.L.I"],# cub low ranked & infected
  
  
  phiCHS = theta[i,"C.H.S"], # cub high ranked & susceptible
  phiCHI = theta[i,"C.H.I"], # cub high ranked & infected
  
  # SUBADULT
  phiSLS = theta[i,"SA.L.S"],
  phiSLI = theta[i,"SA.I+R"],
  phiSLR = theta[i,"SA.I+R"],
  
  phiSHS = theta[i,"SA.H.S"],
  phiSHI = theta[i,"SA.I+R"],
  phiSHR = theta[i,"SA.I+R"],
  
  # BREEDER
  phiBLS = theta[i,"B"],
  phiBLI = theta[i,"B"],
  phiBLR = theta[i,"B"],
  phiBHS = theta[i,"B"],
  phiBHI = theta[i,"B"],
  phiBHR = theta[i,"B"],
  
  # NON-BREEDER
  phiNBLS = theta[i,"NB"],
  phiNBLI = theta[i,"NB"],
  phiNBLR = theta[i,"NB"],
  phiNBHS = theta[i,"NB"],
  phiNBHI = theta[i,"NB"],
  phiNBHR = theta[i,"NB"],
  
  ##########rank transition
  rLL.S = theta[i,"rLL"],  # transition probability of staying in Low social rank for S females
  rLL.I = theta[i,"rLL"], 
  rLL.R = theta[i,"rLL"], 
  
  rHH.S = theta[i,"rHH"], # transition probability of staying in High social rank for S females
  rHH.I = theta[i,"rHH"],
  rHH.R = theta[i,"rHH"], 
  
  ##########Infection 
  betaCL =  1-theta[i,"C.L"], # here we retransform it as "real" beta value
  betaCH =  1-theta[i,"C.H"], # 
  
  betaSL = 1-theta[i,"SA&NB&B.L"],
  betaSH = 1-theta[i,"SA&NB&B.H"],
  
  betaADL = 1-theta[i,"SA&NB&B.L"],
  betaADH = 1-theta[i,"SA&NB&B.H"],
  
  ####### breeding probability 
  
  bSL = theta[i,"bSA.L"],  # transition probability that subadults of Low rank become breeder
  bBL = theta[i,"bB.L"],  # transition probability that breeders of Low rank become breeder again 
  bNBL = theta[i,"bNB.L"],  # transition probability that non-breeder of Low rank become breeder
  
  bSH = theta[i,"bSA.H"],  # transition probability that subadults of High rank become breeder
  bBH = theta[i,"bB.H"],  # transition probability that breeders of High rank become breeder again 
  bNBH = theta[i,"bNB.H"], # transition probability that non-breeder of High rank become breeder
  
  
  ls = theta[i,"ls"], 
  sr = theta[i,"sr"]
 
  ) 

# the matrix Mat is the full model:    
     
Mat <-expression(0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	phiCLS*(1-betaCL)*ls* sr,	phiCLS*(1-betaCL)*ls* sr,	phiCLS*(1-betaCL)*ls* sr,	0,	0,	0,	0,	0,	0,	0,	0,	0,
              0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	phiCLI*(betaCL)*ls* sr,	phiCLI*(betaCL)*ls* sr,	phiCLI*(betaCL)*ls* sr,	0,	0,	0,	0,	0,	0,	0,	0,	0,
              0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	phiCHS*(1-betaCH)*ls* sr,	phiCHS*(1-betaCH)*ls* sr,	phiCHS*(1-betaCH)*ls* sr,	0,	0,	0,	0,	0,	0,
              0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	phiCHI*(betaCH)*ls* sr,	phiCHI*(betaCH)*ls* sr,	phiCHI*(betaCH)*ls* sr,	0,	0,	0,	0,	0,	0,
              phiSLS*(1-betaSL),	0,	 phiSLS*(1-betaSL),	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,
              phiSLI*(betaSL),	0,	 phiSLI*(betaSL),	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,
              0,	phiSLR,	0,	phiSLR,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,
              phiSHS*(1-betaSH),	0,	 phiSHS*(1-betaSH),	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,
              phiSHI*(betaSH),	0,	 phiSHI*(betaSH),	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,
              0,	phiSHR,	0,	phiSHR,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,
              0,	0,	0,	0,	phiBLS*(1-betaADL)*rLL.S*bSL,	0,	0,	phiBLS*(1-betaADL)*(1-rHH.S)*bSH,	0,	0,	phiBLS*(1-betaADL)*rLL.S*bBL,	0,	0,	phiBLS*(1-betaADL)*(1-rHH.S)*bBH,	0,	0,	phiBLS*(1-betaADL)*rLL.S*bNBL,	0,	0,	phiBLS*(1-betaADL)*(1-rHH.S)*bNBH,	0,	0,
              0,	0,	0,	0,	phiBLI*(betaADL)*rLL.I*bSL,	0,	0,	phiBLI*(betaADL)*(1-rHH.I)*bSH,	0,	0,	phiBLI*(betaADL)*rLL.I*bBL,	0,	0,	phiBLI*(betaADL)*(1-rHH.I)*bBH,	0,	0,	phiBLI*(betaADL)*rLL.I*bNBL,	0,	0,	phiBLI*(betaADL)*(1-rHH.I)*bNBH,	0,	0,
              0,	0,	0,	0,	0,	phiBLR*rLL.R*bSL,	phiBLR*rLL.R*bSL,	0,	phiBLR*(1-rHH.R)*bSH,	phiBLR*(1-rHH.R)*bSH,	0,	phiBLR*rLL.R*bBL,	phiBLR*rLL.R*bBL,	0,	phiBLR*(1-rHH.R)*bBH,	phiBLR*(1-rHH.R)*bBH,	0,	phiBLR*rLL.R*bNBL,	phiBLR*rLL.R*bNBL,	0,	phiBLR*(1-rHH.I)*bNBH,	phiBLR*(1-rHH.R)*bNBH,
              0,	0,	0,	0,	phiBHS*(1-betaADH)*(1-rLL.S)*bSL,	0,	0,	phiBHS*(1-betaADH)*(rHH.S)*bSH,	0,	0,	phiBHS*(1-betaADH)*(1-rLL.S)*bBL,	0,	0,	phiBHS*(1-betaADH)*(rHH.S)*bBH,	0,	0,	phiBHS*(1-betaADH)*(1-rLL.S)*bNBL,	0,	0,	phiBHS*(1-betaADH)*rHH.S*bNBH,	0,	0,
              0,	0,	0,	0,	phiBHI*(betaADH)*(1-rLL.I)*bSL,	0,	0,	phiBHI*(betaADH)*(rHH.I)*bSH,	0,	0,	phiBHI*(betaADH)*(1-rLL.I)*bBL,	0,	0,	phiBHI*(betaADH)*(rHH.I)*bBH,	0,	0,	phiBHI*(betaADH)*(1-rLL.I)*bNBL,	0,	0,	phiBHI*(betaADH)*rHH.I*bNBH,	0,	0,
              0,	0,	0,	0,	0,	phiBHR*(1-rLL.R)*bSL,	phiBHR*(1-rLL.R)*bSL,	0,	phiBHR*(rHH.R)*bSH,	phiBHR*(rHH.R)*bSH,	0,	phiBHR*(1-rLL.R)*bBL,	phiBHR*(1-rLL.R)*bBL,	0,	phiBHR*(rHH.R)*bBH,	phiBHR*(rHH.R)*bBH,	0,	phiBHR*(1-rLL.R)*bNBL,	phiBHR*(1-rLL.R)*bNBL,	0,	phiBHR*rHH.R*bNBH,	phiBHR*rHH.R*bNBH,
              0,	0,	0,	0,	phiNBLS*(1-betaADL)*rLL.S*bSL,	0,	0,	phiNBLS*(1-betaADL)*(1-rHH.S)*bSH,	0,	0,	phiNBLS*(1-betaADL)*rLL.S*bBL,	0,	0,	phiNBLS*(1-betaADL)*(1-rHH.S)*bBH,	0,	0,	phiNBLS*(1-betaADL)*rLL.S*bNBL,	0,	0,	phiNBLS*(1-betaADL)*(1-rHH.S)*bNBH,	0,	0,
              0,	0,	0,	0,	phiNBLI*(betaADL)*rLL.I*bSL,	0,	0,	phiNBLI*(betaADL)*(1-rHH.I)*bSH,	0,	0,	phiNBLI*(betaADL)*rLL.I*bBL,	0,	0,	phiNBLI*(betaADL)*(1-rHH.I)*bBH,	0,	0,	phiNBLI*(betaADL)*rLL.I*bNBL,	0,	0,	phiNBLI*(betaADL)*(1-rHH.I)*bNBH,	0,	0,
              0,	0,	0,	0,	0,	phiNBLR*rLL.R*bSL,	phiNBLR*rLL.R*bSL,	0,	phiNBLR*(1-rHH.R)*bSH,	phiNBLR*(1-rHH.R)*bSH,	0,	phiNBLR*rLL.R*bBL,	phiNBLR*rLL.R*bBL,	0,	phiNBLR*(1-rHH.R)*bBH,	phiNBLR*(1-rHH.R)*bBH,	0,	phiNBLR*rLL.R*bNBL,	phiNBLR*rLL.R*bNBL,	0,	phiNBLR*(1-rHH.I)*bNBH,	phiNBLR*(1-rHH.R)*bNBH,
              0,	0,	0,	0,	phiNBHS*(1-betaADH)*(1-rLL.S)*bSL,	0,	0,	phiNBHS*(1-betaADH)*(rHH.S)*bSH,	0,	0,	phiNBHS*(1-betaADH)*(1-rLL.S)*bBL,	0,	0,	phiNBHS*(1-betaADH)*(rHH.S)*bBH,	0,	0,	phiNBHS*(1-betaADH)*(1-rLL.S)*bNBL,	0,	0,	phiNBHS*(1-betaADH)*rHH.S*bNBH,	0,	0,
              0,	0,	0,	0,	phiNBHI*(betaADH)*(1-rLL.I)*bSL,	0,	0,	phiNBHI*(betaADH)*(rHH.I)*bSH,	0,	0,	phiNBHI*(betaADH)*(1-rLL.I)*bBL,	0,	0,	phiNBHI*(betaADH)*(rHH.I)*bBH,	0,	0,	phiNBHI*(betaADH)*(1-rLL.I)*bNBL,	0,	0,	phiNBHI*(betaADH)*rHH.I*bNBH,	0,	0,
              0,	0,	0,	0,	0,	phiNBHR*(1-rLL.R)*bSL,	phiNBHR*(1-rLL.R)*bSL,	0,	phiNBHR*(rHH.R)*bSH,	phiNBHR*(rHH.R)*bSH,	0,	phiNBHR*(1-rLL.R)*bBL,	phiNBHR*(1-rLL.R)*bBL,	0,	phiNBHR*(rHH.R)*bBH,	phiNBHR*(rHH.R)*bBH,	0,	phiNBHR*(1-rLL.R)*bNBL,	phiNBHR*(1-rLL.R)*bNBL,	0,	phiNBHR*rHH.R*bNBH,	phiNBHR*rHH.R*bNBH
            
  )
  


# Tr is the matrix of transitions

Tr <-expression(
      0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,
      0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,
      0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,
      0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,
      phiSLS*(1-betaSL),	0,	 phiSLS*(1-betaSL),	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,
      phiSLI*(betaSL),	0,	 phiSLI*(betaSL),	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,
      0,	phiSLR,	0,	phiSLR,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,
      phiSHS*(1-betaSH),	0,	 phiSHS*(1-betaSH),	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,
      phiSHI*(betaSH),	0,	 phiSHI*(betaSH),	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,
      0,	phiSHR,	0,	phiSHR,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,
      0,	0,	0,	0,	phiBLS*(1-betaADL)*rLL.S*bSL,	0,	0,	phiBLS*(1-betaADL)*(1-rHH.S)*bSH,	0,	0,	phiBLS*(1-betaADL)*rLL.S*bBL,	0,	0,	phiBLS*(1-betaADL)*(1-rHH.S)*bBH,	0,	0,	phiBLS*(1-betaADL)*rLL.S*bNBL,	0,	0,	phiBLS*(1-betaADL)*(1-rHH.S)*bNBH,	0,	0,
      0,	0,	0,	0,	phiBLI*(betaADL)*rLL.I*bSL,	0,	0,	phiBLI*(betaADL)*(1-rHH.I)*bSH,	0,	0,	phiBLI*(betaADL)*rLL.I*bBL,	0,	0,	phiBLI*(betaADL)*(1-rHH.I)*bBH,	0,	0,	phiBLI*(betaADL)*rLL.I*bNBL,	0,	0,	phiBLI*(betaADL)*(1-rHH.I)*bNBH,	0,	0,
      0,	0,	0,	0,	0,	phiBLR*rLL.R*bSL,	phiBLR*rLL.R*bSL,	0,	phiBLR*(1-rHH.R)*bSH,	phiBLR*(1-rHH.R)*bSH,	0,	phiBLR*rLL.R*bBL,	phiBLR*rLL.R*bBL,	0,	phiBLR*(1-rHH.R)*bBH,	phiBLR*(1-rHH.R)*bBH,	0,	phiBLR*rLL.R*bNBL,	phiBLR*rLL.R*bNBL,	0,	phiBLR*(1-rHH.I)*bNBH,	phiBLR*(1-rHH.R)*bNBH,
      0,	0,	0,	0,	phiBHS*(1-betaADH)*(1-rLL.S)*bSL,	0,	0,	phiBHS*(1-betaADH)*(rHH.S)*bSH,	0,	0,	phiBHS*(1-betaADH)*(1-rLL.S)*bBL,	0,	0,	phiBHS*(1-betaADH)*(rHH.S)*bBH,	0,	0,	phiBHS*(1-betaADH)*(1-rLL.S)*bNBL,	0,	0,	phiBHS*(1-betaADH)*rHH.S*bNBH,	0,	0,
      0,	0,	0,	0,	phiBHI*(betaADH)*(1-rLL.I)*bSL,	0,	0,	phiBHI*(betaADH)*(rHH.I)*bSH,	0,	0,	phiBHI*(betaADH)*(1-rLL.I)*bBL,	0,	0,	phiBHI*(betaADH)*(rHH.I)*bBH,	0,	0,	phiBHI*(betaADH)*(1-rLL.I)*bNBL,	0,	0,	phiBHI*(betaADH)*rHH.I*bNBH,	0,	0,
      0,	0,	0,	0,	0,	phiBHR*(1-rLL.R)*bSL,	phiBHR*(1-rLL.R)*bSL,	0,	phiBHR*(rHH.R)*bSH,	phiBHR*(rHH.R)*bSH,	0,	phiBHR*(1-rLL.R)*bBL,	phiBHR*(1-rLL.R)*bBL,	0,	phiBHR*(rHH.R)*bBH,	phiBHR*(rHH.R)*bBH,	0,	phiBHR*(1-rLL.R)*bNBL,	phiBHR*(1-rLL.R)*bNBL,	0,	phiBHR*rHH.R*bNBH,	phiBHR*rHH.R*bNBH,
      0,	0,	0,	0,	phiNBLS*(1-betaADL)*rLL.S*bSL,	0,	0,	phiNBLS*(1-betaADL)*(1-rHH.S)*bSH,	0,	0,	phiNBLS*(1-betaADL)*rLL.S*bBL,	0,	0,	phiNBLS*(1-betaADL)*(1-rHH.S)*bBH,	0,	0,	phiNBLS*(1-betaADL)*rLL.S*bNBL,	0,	0,	phiNBLS*(1-betaADL)*(1-rHH.S)*bNBH,	0,	0,
      0,	0,	0,	0,	phiNBLI*(betaADL)*rLL.I*bSL,	0,	0,	phiNBLI*(betaADL)*(1-rHH.I)*bSH,	0,	0,	phiNBLI*(betaADL)*rLL.I*bBL,	0,	0,	phiNBLI*(betaADL)*(1-rHH.I)*bBH,	0,	0,	phiNBLI*(betaADL)*rLL.I*bNBL,	0,	0,	phiNBLI*(betaADL)*(1-rHH.I)*bNBH,	0,	0,
      0,	0,	0,	0,	0,	phiNBLR*rLL.R*bSL,	phiNBLR*rLL.R*bSL,	0,	phiNBLR*(1-rHH.R)*bSH,	phiNBLR*(1-rHH.R)*bSH,	0,	phiNBLR*rLL.R*bBL,	phiNBLR*rLL.R*bBL,	0,	phiNBLR*(1-rHH.R)*bBH,	phiNBLR*(1-rHH.R)*bBH,	0,	phiNBLR*rLL.R*bNBL,	phiNBLR*rLL.R*bNBL,	0,	phiNBLR*(1-rHH.I)*bNBH,	phiNBLR*(1-rHH.R)*bNBH,
      0,	0,	0,	0,	phiNBHS*(1-betaADH)*(1-rLL.S)*bSL,	0,	0,	phiNBHS*(1-betaADH)*(rHH.S)*bSH,	0,	0,	phiNBHS*(1-betaADH)*(1-rLL.S)*bBL,	0,	0,	phiNBHS*(1-betaADH)*(rHH.S)*bBH,	0,	0,	phiNBHS*(1-betaADH)*(1-rLL.S)*bNBL,	0,	0,	phiNBHS*(1-betaADH)*rHH.S*bNBH,	0,	0,
      0,	0,	0,	0,	phiNBHI*(betaADH)*(1-rLL.I)*bSL,	0,	0,	phiNBHI*(betaADH)*(rHH.I)*bSH,	0,	0,	phiNBHI*(betaADH)*(1-rLL.I)*bBL,	0,	0,	phiNBHI*(betaADH)*(rHH.I)*bBH,	0,	0,	phiNBHI*(betaADH)*(1-rLL.I)*bNBL,	0,	0,	phiNBHI*(betaADH)*rHH.I*bNBH,	0,	0,
      0,	0,	0,	0,	0,	phiNBHR*(1-rLL.R)*bSL,	phiNBHR*(1-rLL.R)*bSL,	0,	phiNBHR*(rHH.R)*bSH,	phiNBHR*(rHH.R)*bSH,	0,	phiNBHR*(1-rLL.R)*bBL,	phiNBHR*(1-rLL.R)*bBL,	0,	phiNBHR*(rHH.R)*bBH,	phiNBHR*(rHH.R)*bBH,	0,	phiNBHR*(1-rLL.R)*bNBL,	phiNBHR*(1-rLL.R)*bNBL,	0,	phiNBHR*rHH.R*bNBH,	phiNBHR*rHH.R*bNBH
    
    )
    
    
# F is the fertility matrix:    
    
F <- expression(0, 0,           0, 0,	         0, 0,		   0, 0, 0,	        0, 0, 0,	     0, 0, 0,	          0,  0,  0,	          0, 0, 0,              0,
                            0, phiCLS*betaCL, 0, phiCHS*betaCL, 0, phiSLS*betaCL, 0, 0, phiSHS*betaCL, 0, 0, phiBLS*betaCL, 0, 0, phiBHS*betaCL, 0,  0,  phiNBLS*betaCL, 0, 0, phiNBHS*betaCL,0,
                            0, 0,  	       0, 0,  	         0, 0,		   0, 0, 0,	        0, 0, 0,	     0, 0, 0,             0,  0,  0,              0, 0, 0,              0,
                            0, phiCLS*betaCH, 0, phiCHS*betaCH, 0, phiSLS*betaCH, 0, 0, phiSHS*betaCH, 0, 0, phiBLS*betaCH, 0, 0, phiBHS*betaCH, 0,  0,  phiNBLS*betaCH, 0, 0, phiNBHS*betaCH, 0,
                            0, 0,   	       0, 0,  	         0, 0,		   0, 0, 0,	        0, 0, 0,	     0, 0, 0,             0,  0,  0,              0, 0, 0,              0,
                            0, phiCLS*betaSL, 0, phiCHS*betaSL, 0, phiSLS*betaSL, 0, 0, phiSHS*betaSL, 0, 0, phiBLS*betaSL, 0, 0, phiBHS*betaSL, 0,  0,  phiNBLS*betaSL, 0, 0, phiNBHS*betaSL, 0,
                            0, 0,             0, 0,  	         0, 0,		   0, 0, 0,	        0, 0, 0,	     0, 0, 0,             0,  0,  0,              0, 0, 0,              0,
                            0, 0,             0, 0,  	         0, 0,		   0, 0, 0,	        0, 0, 0,	     0, 0, 0,             0,  0,  0,              0, 0, 0,              0,
                            0, phiCLS*betaSH, 0, phiCHS*betaSH, 0, phiSLS*betaSH, 0, 0, phiSHS*betaSH, 0, 0, phiBLS*betaSH, 0, 0, phiBHS*betaSH, 0,  0,  phiNBLS*betaSH, 0, 0, phiNBHS*betaSH, 0,
                            0, 0,             0, 0,  	         0, 0,		   0, 0, 0,             0, 0, 0,	     0, 0, 0,             0,  0,  0,              0, 0, 0,              0,
                            0, 0,             0, 0, 	         0, 0,		   0, 0, 0,             0, 0, 0,	     0, 0, 0,             0,  0,  0,              0, 0, 0,              0,
                            0, phiCLS*betaADL, 0, phiCHS*betaADL, 0, phiSLS*betaADL, 0, 0, phiSHS*betaADL, 0, 0, phiBLS*betaADL, 0, 0, phiBHS*betaADL, 0,  0,  phiNBLS*betaADL, 0, 0, phiNBHS*betaADL, 0,
                            0, 0,             0, 0,  	         0, 0,		   0, 0, 0,             0, 0, 0,	     0, 0, 0,	       	  0,  0,  0,              0, 0, 0,              0,
                            0, 0,             0, 0,  	         0, 0,		   0, 0, 0,             0, 0, 0,	     0, 0, 0,	 	  0,  0,  0,              0, 0, 0,              0,
                            0, phiCLS*betaADH, 0, phiCHS*betaADH, 0, phiSLS*betaADH, 0, 0, phiSHS*betaADH, 0, 0, phiBLS*betaADH, 0, 0, phiBHS*betaADH, 0,  0,  phiNBLS*betaADH, 0, 0, phiNBHS*betaADH, 0,
                            0, 0,             0, 0, 	         0, 0,		   0, 0, 0,             0, 0, 0,	     0, 0, 0,             0,  0,  0,              0, 0, 0,              0,
                            0, 0,             0, 0, 	         0, 0,		   0, 0, 0,             0, 0, 0,	     0, 0, 0,             0,  0,  0,              0, 0, 0,              0,
                            0, phiCLS*betaADL, 0, phiCHS*betaADL, 0, phiSLS*betaADL, 0, 0, phiSHS*betaADL, 0, 0, phiBLS*betaADL, 0, 0, phiBHS*betaADL, 0,  0,  phiNBLS*betaADL, 0, 0, phiNBHS*betaADL, 0,
                            0, 0,             0, 0, 	         0, 0,		   0, 0, 0,             0, 0, 0,	     0, 0, 0,             0,  0,  0,              0, 0, 0,              0,
                            0, 0,             0, 0, 	         0, 0,		   0, 0, 0,             0, 0, 0,	     0, 0, 0,             0,  0,  0,              0, 0, 0,              0,
                            0, phiCLS*betaADH, 0, phiCHS*betaADH, 0, phiSLS*betaADH, 0, 0, phiSHS*betaADH, 0, 0, phiBLS*betaADH, 0, 0, phiBHS*betaADH, 0,  0,  phiNBLS*betaADH, 0, 0, phiNBHS*betaADH, 0,
                            0, 0,             0, 0, 	         0, 0,		   0, 0, 0,             0, 0, 0,	     0, 0, 0,             0,  0,  0,              0, 0, 0,              0
  )
    
  M.final<-matrix(sapply(Mat, eval, param),  nrow = (NStages), ncol = (NStages), byrow = TRUE)
  
 
  sens<-vitalsens(Mat,param)
  
   
  Fvalue<-matrix(sapply(F, eval, param),  nrow = (NStages), ncol = (NStages), byrow = TRUE)
  
  Tr.final <-  M.final
  Tr.final1 <-  M.final
  Tr.final1[1:4,] <- 0
  
  Id <- diag(22)
  
  ### --- Building N, the fundamental matrix as N = (Id-Tr.final)-1
  
  
  N <- solve((lambda(Tr.final)*Id) - Tr.final1) 
  
  
  ### --- Next generation Matrix NGM:
  
  NGM <- Fvalue %*% N  
  
  return(list(M.final, sens, NGM)) #the aim of build_matrix is to return the NGM for each value in theta 
}


#####step 2 : calculating the SD of lambda and sensitivity analysis of lambda


initmat<-build_matrix(theta, 1, checkNodisease)[[1]]
#paraminit<-initmat[4]


initpop<-stable.stage(initmat)*popsize0



############creating the population vector######
popvec<-array(0,dim=c(Tmax, MCiter))
########### matrice projection with simu MC
NStages<-22
M<-matrix(0,nrow=NStages, ncol=NStages)
Mproj <- replicate(MCiter, M, simplify=FALSE)
NGMstoch <- replicate(MCiter, M, simplify=FALSE)


sensmat<-matrix(0,nrow=42, ncol=3)

senslambda <- replicate(MCiter, sensmat, simplify=FALSE)
Rmatrix<-array(0,dim=c(NStages, NStages, MCiter))


for(i in 1: MCiter) # MCiter = 1000
{  
  allres<-build_matrix(theta, i, checkNodisease)
  Mproj[[i]] <- allres[[1]]
  senslambda[[i]] <- allres[[2]]
  NGMstoch[[i]] <- allres[[3]]
  trans <- pop.projection(Mproj[[i]],  initpop , (Tmax)) # Vector => chose an initial vector to test different scenarios  if vector="n", then project will automatically project the set of stage-biased vectors of A.
  popvec[,i] <- t(trans$pop.sizes)
}


#########----------------------------------- STEP 3 : calculating sensitivity of R0


# Next function: calculate sensitivity to parameter in position pos in vector of parameters that follows and matches the input parameters of function build_matrix 
# delta is the perturbation parameter (set to 1e-4 by default)


sens_elas_num <- function(pos, theta, delta=1e-4){
  # param char format
    # get parameters
  sensR0<-NULL
  elasR0<-NULL

  for(i in 1: MCiter) # MCiter = 1000
  {  

    # build R0 matrix
    A <- build_matrix(theta, i, checkNodisease)[[3]] # Here we apply build_matrix that returns the NGM. A is NGM
 
    # calculate growth rate ---> lambda is R0 here (the eigenvalue of NGM)
    lambda = max(Re(eigen(A)$values)) # Re primitive
  
    # get focal parameter
    c = theta[i,pos]
  
    # modify the focal parameter c by a very small amount - so for beta here it modifies proba S->S  (not S->I)
    c_new = c * (1 + delta)
    theta_new = theta
    theta_new[i,pos] = c_new
    
    # build A_new with perturbed focal parameter
    A_new <- build_matrix(theta_new, i, checkNodisease)[[3]]
    
    # calculate growth rate (new R0)
    lambda_new = max(Re(eigen(A_new)$values))
    
    # calculate sensitivity [sens = df(x)/dx = (lam.new-lam)/(c*delta)]
    sensR0 = c(sensR0, (lambda_new-lambda) / ( c * delta))
    
    # calculate elasticity [elas = sens*c/lam = (lam.new-lam)/(lam*delta)]
    elasR0 = c(elasR0, (lambda_new - lambda)/(lambda*delta))
  
  }
  
  res = list(param = colnames(theta)[pos], sensR0 = sensR0, elasR0 = elasR0)
  
  return(res) # returns a list with column names, sens, elasticity values
}
# ----> end of function sens_elas_num 

