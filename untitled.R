W
WG
W both
W steve
W imac

E
EG

#Exponential Model
          eCov <- function(x) { ## function to optimize
            x1 <- x[1] #lambda
            x2 <- x[2] # beta
            imatr$exb=exp(imatr$adv..bn.*x2)
            imatr$A=cumsum(imatr$exb)
            imatr$Pt<-1-exp(-x1*imatr$A)
            imatr$dP[2:nrow(imatr)]=diff(imatr$Pt)
            imatr$dP[1]=imatr$Pt[1]
            imatr$LL=imatr$sales*log(imatr$dP)
            -sum(imatr$LL, (150452-sum(imatr[,3]))*log(1-imatr$Pt[40]))
          }
#Weibull Model
          wCov <- function(x) { ## function to optimize
            x1 <- x[1] #lambda
            x2 <- x[2] # beta
            x3 <- x[3] # c
            imatr$exb=exp(imatr$adv..bn.*x2)
            imatr$A=cumsum(((imatr$t^x3) - (c(0, imatr$t[1:39]))^x3)*imatr$exb)
            imatr$Pt<-1-exp(-x1*imatr$A)
            imatr$dP[2:nrow(imatr)]=diff(imatr$Pt)
            imatr$dP[1]=imatr$Pt[1]
            imatr$LL=imatr$sales*log(imatr$dP)
            -sum(imatr$LL, (150452-sum(imatr[,3]))*log(1-imatr$Pt[40]))
          }
#Exponential, iMac Covar
        if (modelType == 1 && length(covariates) == 1 && covariates[1] == 'ni') {
          e2Cov <- function(x) { ## function to optimize
            x1 <- x[1] #lambda
            x2 <- x[2] # beta
            x4 <- x[3] #b_iMac
            imatr$exb=exp(imatr$adv..bn.*x2 + imatr$iMac*x4)
            imatr$A=imatr$exb)
            imatr$Pt<-1-exp(-x1*imatr$A)
            imatr$dP[2:nrow(imatr)]=diff(imatr$Pt)
            imatr$dP[1]=imatr$Pt[1]
            imatr$LL=imatr$sales*log(imatr$dP)
            -sum(imatr$LL, (150452-sum(imatr[,3]))*log(1-imatr$Pt[40]))
          }