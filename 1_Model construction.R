
Surv.C <- diag(c(phiCLS, phiCLI, phiCLR, phiCHS, phiCHI,phiCHR))
Surv.S <- diag(c(phiSLS, phiSLI, phiSLR, phiSHS, phiSHI, phiSHR)) 
Surv.B <- diag(c(phiBLS, phiBLI, phiBLR, phiBHS, phiBHI, phiBHR))
Surv.NB <- diag(c(phiNBLS, phiNBLI,phiNBLR, phiNBHS, phiNBHI, phiNBHR)) 

fB <- ls* sr * diag(c(fLS, fLI, fLR, fHS, fHI,fHR)) 

Rec.S <- diag(c(bSL, bSL, bSL, bSH, bSH,bSH))
Rec.B <- diag(c(bBL, bBL, bBL, bBH, bBH,bBH))
Rec.NB <- diag(c(bNBL, bNBL, bNBL, bNBH, bNBH,bNBH))
NRec.S <- diag(c((1-bSL), (1-bSL), (1-bSL),(1-bSH), (1-bSH),(1-bSH)))
NRec.B <- diag(c((1-bBL), (1-bBL), (1-bBL),(1-bBH), (1-bBH),(1-bBH)))
NRec.NB <- diag(c((1-bNBL), (1-bNBL), (1-bNBL),(1-bNBH), (1-bNBH),(1-bNBH)))

R <- matrix(c(
  rLL.S,       0,           0,          1-rHH.S,    0,          0,
  0,           rLL.I,       0,          0,          1-rHH.I,    0,
  0,           0,           rLL.R,      0,          0,          1-rHH.R,
  1-rLL.S,     0,           0,          rHH.S,      0,          0,
  0,           1-rLL.I,     0,          0,          rHH.I,      0,
  0,           0,           1-rLL.R,    0,          0,          rHH.R)
  
  ,byrow=T,nrow=6)

I.F <- matrix(c(  
  1,         1,         1,          0,          0,          0,
  0,         0,         0,          0,          0,          0,
  0,         0,         0,          0,          0,          0,
  0,         0,         0,          1,          1,          1,
  0,         0,         0,          0,          0,          0,
  0,         0,         0,          0,          0,          0)
  
  ,byrow=T,nrow=6)

I.C <- matrix(c(
  
  1-betaCL,  0,         0,          0,          0,          0,
  betaCL,    0,         0,          0,          0,          0,
  0,         1,         0,          0,          0,          0,
  0,         0,         0,          1-betaCH,   0,          0,
  0,         0,         0,          betaCH,     0,          0,
  0,         0,         0,          0,          1,          0)
  
  ,byrow=T,nrow=6)

I.SA <- matrix(c(
  
  1-betaSL,   0,         0,          0,          0,          0,
  betaSL,     0,         0,          0,          0,          0,
  0,         1,         1,          0,          0,          0,
  0,         0,         0,          1-betaSH,    0,          0,
  0,         0,         0,          betaSH,      0,          0,
  0,         0,         0,          0,          1,          1)
  
  ,byrow=T,nrow=6)

I.AD <- matrix(c(
  
  1-betaADL,   0,         0,          0,          0,          0,
  betaADL,     0,         0,          0,          0,          0,
  0,         1,         1,          0,          0,          0,
  0,         0,         0,          1-betaADH,    0,          0,
  0,         0,         0,          betaADH,      0,          0,
  0,         0,         0,          0,          1,          1)
  
  ,byrow=T,nrow=6)

NetF<-Surv.C%*%I.C%*%I.F%*% fB 

CS <- Surv.S%*%I.SA  
SB <- Surv.B%*%I.AD%*%R%*%Rec.S 
BB <- Surv.B%*%I.AD%*%R%*%Rec.B 
NBB <- Surv.B%*%I.AD%*%R%*%Rec.NB 
SNB <- Surv.NB%*%I.AD%*%R%*%NRec.S
BNB <- Surv.NB%*%I.AD%*%R%*%NRec.B 
NBNB<- Surv.NB%*%I.AD%*%R%*%NRec.NB 

Z <- matrix(0,6,6) 

M <-  rbind(cbind (Z,       Z,       NetF,   Z),
            cbind (CS,      Z,       Z,      Z),
            cbind (Z,       SB,      BB,     NBB),
            cbind (Z,       SNB,     BNB,    NBNB))

 
rownames(M)<-  c("CLS","CLI","CLR","CHS","CHI","CHR",
                  "SLS","SLI","SLR","SHS","SHI","SHR",
                  "BLS","BLI","BLR","BHS","BHI","BHR",
                  "NBLS","NBLI","NBLR","NBHS","NBHI","NBHR")        

colnames(M)<-  c("CLS","CLI","CLR","CHS","CHI","CHR",
                  "SLS","SLI","SLR","SHS","SHI","SHR",
                  "BLS","BLI","BLR","BHS","BHI","BHR",
                  "NBLS","NBLI","NBLR","NBHS","NBHI","NBHR")        

M.final <- M[-c(3,6), -c(3,6)] # We remove the recovered cubs as they do not exist in the data 


Tr <-  rbind(cbind (Z,       Z,       Z,      Z),
            cbind (CS,       Z,       Z,      Z),
            cbind (Z,       SB,      BB,     NBB),
            cbind (Z,       SNB,     BNB,    NBNB))

Tr.final <- Tr[-c(3,6), -c(3,6)] 

Id <- diag(22)

N <- solve((lambda(M.final)*Id) -Tr.final) # Here we include the population's growth rate

F <- matrix(c(
  
  0, 0,           0, 0,	         0, 0,		   0, 0, 0,	        0, 0, 0,	     0, 0, 0,	          0,  0,  0,	          0, 0, 0,              0,
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
      
  ),byrow=T,nrow=22)


NGM <- F %*% N 

rownames(NGM)<-  c("CLS","CLI","CHS","CHI",
                 "SLS","SLI","SLR","SHS","SHI","SHR",
                 "BLS","BLI","BLR","BHS","BHI","BHR",
                 "NBLS","NBLI","NBLR","NBHS","NBHI","NBHR")        

colnames(NGM)<-  c("CLS","CLI","CHS","CHI",
                 "SLS","SLI","SLR","SHS","SHI","SHR",
                 "BLS","BLI","BLR","BHS","BHI","BHR",
                 "NBLS","NBLI","NBLR","NBHS","NBHI","NBHR")





