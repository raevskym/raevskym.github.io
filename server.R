library(shiny)

#Data Source: Statista https://proxy.library.upenn.edu:5546/statistics/263444/sales-of-apple-mac-computers-since-first-quarter-2006/
iData <- read.csv("iData.csv")
imatr <- iData
imatr[,3] <- as.numeric(gsub(",","",imatr[,3]))
range.ex <- function(x, fac=1.2) { xrg <- range(x);  m <- mean(xrg);  (xrg - m)*fac + m }

tripsPlot <- function(startQuarter, modelType, covariates, mixture, benchmark, num) {
  
  #adjustment for user-inputed 2016 advertisements
  imatr[41:44,4] <- num
  
        if (modelType == 1 && !mixture) {
          #Exponential Model
          eCov <- function(x) { ## function to optimize
            x1 <- x[1] #lambda
            x2 <- x[2] # beta
            imatr$exb=exp(imatr$adv..bn.*x2)
            imatr$A=cumsum(imatr$exb)
            imatr$Pt<-1-exp(-x1*imatr$A)
            imatr$dP[2:nrow(imatr)]=diff(imatr$Pt)
            imatr$dP[1]=imatr$Pt[1]
            imatr$LL <- c((imatr$sales[1:40])*(log(imatr$dP)[1:40]), (150452-sum(imatr[1:40,3]))*log(1-imatr$Pt[40]), NA, NA, NA)
            -sum(imatr$LL[1:41])
          }
          initial_guess=c(.005, 2)
          x <- optim(initial_guess, eCov)$par
          x1 <- x[1] #lambda
          x2 <- x[2] # beta
          imatr$exb=exp(imatr$adv..bn.*x2)
          imatr$A=cumsum(imatr$exb)
          imatr$Pt<-1-exp(-x1*imatr$A)
          imatr$dP[2:nrow(imatr)]=diff(imatr$Pt)
          imatr$dP[1]=imatr$Pt[1]
          imatr$LL <- c((imatr$sales[1:40])*(log(imatr$dP)[1:40]), (150452-sum(imatr[1:40,3]))*log(1-imatr$Pt[40]), NA, NA, NA)
          LL <- -sum(imatr$LL[1:41])
          imatr$exp <- imatr$dP * 150542
          MAPE <- mean(abs((imatr$sales[1:40]-imatr$exp[1:40])/imatr$sales[1:40]))
        
          z <- cbind(lambda = x1, b_ads = x2, LL = -LL, MAPE = MAPE)
          rownames(z) <- "parameters"
          
          
          }
  
        if (modelType == 1 && mixture) {
          #Exponential Gamma EG
          EG <- function(x) { ## function to optimize
            x1 <- x[1] #r
            x2 <- x[2] #alpha
            x3 <- x[3] #b
            imatr$exb=exp(imatr$adv..bn.*x3)
            imatr$A=cumsum(imatr$exb)
            imatr$Pt<-1-((x2/(x2+imatr$A)))^x1
            imatr$dP[2:nrow(imatr)]=diff(imatr$Pt)
            imatr$dP[1]=imatr$Pt[1]
            imatr$LL <- c((imatr$sales[1:40])*(log(imatr$dP)[1:40]), (150452-sum(imatr[1:40,3]))*log(1-imatr$Pt[40]), NA, NA, NA)
            -sum(imatr$LL[1:41])
          }
          initial_guess=c(300, 50000, 1.5)
          optim(initial_guess, EG)
          x1 <- optim(initial_guess, EG)$par[1] #r
          x2 <- optim(initial_guess, EG)$par[2] #alpha
          x3 <- optim(initial_guess, EG)$par[3] #b
          imatr$exb=exp(imatr$adv..bn.*x3)
          imatr$A=cumsum(imatr$exb)
          imatr$Pt<-1-((x2/(x2+imatr$A)))^x1
          imatr$dP[2:nrow(imatr)]=diff(imatr$Pt)
          imatr$dP[1]=imatr$Pt[1]
          imatr$LL <- c((imatr$sales[1:40])*(log(imatr$dP)[1:40]), (150452-sum(imatr[1:40,3]))*log(1-imatr$Pt[40]), NA, NA, NA)
          LL <- -sum(imatr$LL[1:41])
          imatr$exp <- imatr$dP * 150542
          MAPE <- mean(abs((imatr$sales[1:40]-imatr$exp[1:40])/imatr$sales[1:40]))
        
          z <- cbind(r = x1, alpha = x2, b_ads = x3, LL = -LL, MAPE = MAPE)
          rownames(z) <- "parameters"
          
          
          }
              
  #Exponential, iMac Covar
            if (modelType == 1 && length(covariates) == 1 && covariates[1] == 'ni') {
              e2Cov <- function(x) { ## function to optimize
                x1 <- x[1] #lambda
                x2 <- x[2] # beta
                x4 <- x[3] #b_iMac
                imatr$exb=exp(imatr$adv..bn.*x2 + imatr$iMac*x4)
                imatr$A=cumsum(imatr$exb)
                imatr$Pt<-1-exp(-x1*imatr$A)
                imatr$dP[2:nrow(imatr)]=diff(imatr$Pt)
                imatr$dP[1]=imatr$Pt[1]
                imatr$LL <- c((imatr$sales[1:40])*(log(imatr$dP)[1:40]), (150452-sum(imatr[1:40,3]))*log(1-imatr$Pt[40]), NA, NA, NA)
                -sum(imatr$LL[1:41])
              }
                initial_guess=c(.005, 1.5, 0.5)
                x <- optim(initial_guess, e2Cov)$par
                x1 <- x[1] # lambda
                x2 <- x[2] # beta
                x4 <- x[3] # b_iMac
                imatr$exb=exp(imatr$adv..bn.*x2 + imatr$iMac*x4)
                imatr$A=cumsum(imatr$exb)
                imatr$Pt<-1-exp(-x1*imatr$A)
                imatr$dP[2:nrow(imatr)]=diff(imatr$Pt)
                imatr$dP[1]=imatr$Pt[1]
                imatr$LL <- c((imatr$sales[1:40])*(log(imatr$dP)[1:40]), (150452-sum(imatr[1:40,3]))*log(1-imatr$Pt[40]), NA, NA, NA)
                LL <- -sum(imatr$LL[1:41])
                imatr$exp <- imatr$dP * 150542
                MAPE <- mean(abs((imatr$sales[1:40]-imatr$exp[1:40])/imatr$sales[1:40]))
               
                 
                z <- cbind(lambda = x1, b_ads = x2, b_iMac = x4, LL = -LL, MAPE = MAPE)
                rownames(z) <- "parameters"
                
                
                }
      
  #Exponential, Steve Covar
  if (modelType == 1 && length(covariates) == 1 && covariates[1] == 'jd') {
    e2CovJobs <- function(x) { ## function to optimize
      x1 <- x[1] #lambda
      x2 <- x[2] # beta
      x5 <- x[3] #b_steve
      x6 <- x[4] #gamma
      x7 <- x[5] #delta
      imatr$sJobs <- c(array(0, c(1,23)), 1-x6*(1-exp(-x7*abs(imatr$t-24)))[24:44])
      imatr$exb=exp(imatr$adv..bn.*x2 + imatr$sJobs*x5)
      imatr$A=cumsum(imatr$exb)
      imatr$Pt<-1-exp(-x1*imatr$A)
      imatr$dP[2:nrow(imatr)]=diff(imatr$Pt)
      imatr$dP[1]=imatr$Pt[1]
      imatr$LL <- c((imatr$sales[1:40])*(log(imatr$dP)[1:40]), (150452-sum(imatr[1:40,3]))*log(1-imatr$Pt[40]), NA, NA, NA)
      -sum(imatr$LL[1:41])
    }
    initial_guess=c(0.005, 1.5, 0.5, 0.2, 14)
    x <- optim(initial_guess, e2CovJobs)$par
    x1 <- x[1] #lambda
    x2 <- x[2] # beta
    x5 <- x[3] #b_steve
    x6 <- x[4] #gamma
    x7 <- x[5] #delta
    imatr$sJobs <- c(array(0, c(1,23)), 1-x6*(1-exp(-x7*abs(imatr$t-24)))[24:44])
    imatr$exb=exp(imatr$adv..bn.*x2 + imatr$sJobs*x5)
    imatr$A=cumsum(imatr$exb)
    imatr$Pt<-1-exp(-x1*imatr$A)
    imatr$dP[2:nrow(imatr)]=diff(imatr$Pt)
    imatr$dP[1]=imatr$Pt[1]
    imatr$LL <- c((imatr$sales[1:40])*(log(imatr$dP)[1:40]), (150452-sum(imatr[1:40,3]))*log(1-imatr$Pt[40]), NA, NA, NA)
    LL <- -sum(imatr$LL[1:41])
    imatr$exp <- imatr$dP * 150542
    MAPE <- mean(abs((imatr$sales[1:40]-imatr$exp[1:40])/imatr$sales[1:40]))
 
    
    z <- cbind(lambda = x1, b_ads = x2, b_jobs = x5, gamma = x6, delta = x7, LL = -LL, MAPE = MAPE)
    rownames(z) <- "parameters"
    
    
    }
  
  #Exponential, Both Covar
  if (modelType == 1 && length(covariates) == 2) {
      e3Cov <- function(x) { ## function to optimize
        x1 <- x[1] #lambda
        x2 <- x[2] # beta
        x4 <- x[3] #b_iMac
        x5 <- x[4] #b_steve
        x6 <- x[5] #gamma
        x7 <- x[6] #delta
        imatr$sJobs <- c(array(0, c(1,23)), 1-x6*(1-exp(-x7*abs(imatr$t-24)))[24:44])
        imatr$exb=exp(imatr$adv..bn.*x2 + imatr$iMac*x4 + imatr$sJobs*x5)
        imatr$A=cumsum(imatr$exb)
        imatr$Pt<-1-exp(-x1*imatr$A)
        imatr$dP[2:nrow(imatr)]=diff(imatr$Pt)
        imatr$dP[1]=imatr$Pt[1]
        imatr$LL <- c((imatr$sales[1:40])*(log(imatr$dP)[1:40]), (150452-sum(imatr[1:40,3]))*log(1-imatr$Pt[40]), NA, NA, NA)
        -sum(imatr$LL[1:41])
      }
      initial_guess=c(0.005, 1.5, 0.5, 0.5, 0.2, 14)
      optim(initial_guess, e3Cov)
      x1 <- optim(initial_guess, e3Cov)$par[1] #lamda
      x2 <- optim(initial_guess, e3Cov)$par[2] #b_ads
      x4 <- optim(initial_guess, e3Cov)$par[3] #b_iMac
      x5 <- optim(initial_guess, e3Cov)$par[4] #b_steve
      x6 <- optim(initial_guess, e3Cov)$par[5] #gamma
      x7 <- optim(initial_guess, e3Cov)$par[6] #delta
      imatr$sJobs <- c(array(0, c(1,23)), 1-x6*(1-exp(-x7*abs(imatr$t-24)))[24:44])
      imatr$exb=exp(imatr$adv..bn.*x2 + imatr$iMac*x4 + imatr$sJobs*x5)
      imatr$A=cumsum(imatr$exb)
      imatr$Pt<-1-exp(-x1*imatr$A)
      imatr$dP[2:nrow(imatr)]=diff(imatr$Pt)
      imatr$dP[1]=imatr$Pt[1]
      imatr$LL <- c((imatr$sales[1:40])*(log(imatr$dP)[1:40]), (150452-sum(imatr[1:40,3]))*log(1-imatr$Pt[40]), NA, NA, NA)
      LL <- -sum(imatr$LL[1:41])
      imatr$exp <- imatr$dP * 150542
      MAPE <- mean(abs((imatr$sales[1:40]-imatr$exp[1:40])/imatr$sales[1:40]))
      
      z <- cbind(lambda = x1, b_ads = x2, b_iMac = x4, b_jobs = x5, gamma = x6, delta = x7, LL = -LL, MAPE = MAPE)
      rownames(z) <- "parameters"
      
      
      }
  
  #Weibull Model
        if (modelType == 2 && length(covariates) == 0 && !mixture) {
          wCov <- function(x) { ## function to optimize
            x1 <- x[1] #lambda
            x2 <- x[2] # beta
            x3 <- x[3] # c
            imatr$exb=exp(imatr$adv..bn.*x2)
            imatr$A=cumsum(((imatr$t^x3) - (c(0, imatr$t[1:43]))^x3)*imatr$exb)
            imatr$Pt<-1-exp(-x1*imatr$A)
            imatr$dP[2:nrow(imatr)]=diff(imatr$Pt)
            imatr$dP[1]=imatr$Pt[1]
            imatr$LL <- c((imatr$sales[1:40])*(log(imatr$dP)[1:40]), (150452-sum(imatr[1:40,3]))*log(1-imatr$Pt[40]), NA, NA, NA)
            -sum(imatr$LL[1:41])
          }
          initial_guess=c(.001, 1.5, 1.5)
          x <- optim(initial_guess, wCov)$par
          x1 <- x[1] #lambda
          x2 <- x[2] # beta
          x3 <- x[3] # c
          imatr$exb=exp(imatr$adv..bn.*x2)
          imatr$A=cumsum(((imatr$t^x3) - (c(0, imatr$t[1:43]))^x3)*imatr$exb)
          imatr$Pt<-1-exp(-x1*imatr$A)
          imatr$dP[2:nrow(imatr)]=diff(imatr$Pt)
          imatr$dP[1]=imatr$Pt[1]
          imatr$LL <- c((imatr$sales[1:40])*(log(imatr$dP)[1:40]), (150452-sum(imatr[1:40,3]))*log(1-imatr$Pt[40]), NA, NA, NA)
          LL <- -sum(imatr$LL[1:41])
          imatr$exp <- imatr$dP * 150542
          MAPE <- mean(abs((imatr$sales[1:40]-imatr$exp[1:40])/imatr$sales[1:40]))
         
          z <- cbind(lambda = x1, b_ads = x2, c = x3, LL = -LL, MAPE = MAPE)
          rownames(z) <- "parameters"
          
          
           }
  
  
    if (modelType == 2 && mixture) {
          #WG Weibull Gamma
          WG <- function(x) { ## function to optimize
            x1 <- x[1] # r
            x2 <- x[2] # beta
            x3 <- x[3] # c
            x4 <- x[4] #alpha
            imatr$exb=exp(imatr$adv..bn.*x2)
            imatr$A=cumsum(((imatr$t^x3) - (c(0, imatr$t[1:43]))^x3)*imatr$exb)
            imatr$Pt<-1-((x4/(x4+imatr$A)))^x1
            imatr$dP[2:nrow(imatr)]=diff(imatr$Pt)
            imatr$dP[1]=imatr$Pt[1]
            imatr$LL <- c((imatr$sales[1:40])*(log(imatr$dP)[1:40]), (150452-sum(imatr[1:40,3]))*log(1-imatr$Pt[40]), NA, NA, NA)
            -sum(imatr$LL[1:41])
          }
          initial_guess=c(80, 1.5, 1.5, 30000)
          x <- optim(initial_guess, WG)$par
          x1 <- x[1] #r
          x2 <- x[2] # beta
          x3 <- x[3] # c
          x4 <- x[4] #alpha
          imatr$exb=exp(imatr$adv..bn.*x2)
          imatr$A=cumsum(((imatr$t^x3) - (c(0, imatr$t[1:43]))^x3)*imatr$exb)
          imatr$Pt<-1-((x4/(x4+imatr$A)))^x1
          imatr$dP[2:nrow(imatr)]=diff(imatr$Pt)
          imatr$dP[1]=imatr$Pt[1]
          imatr$LL <- c((imatr$sales[1:40])*(log(imatr$dP)[1:40]), (150452-sum(imatr[1:40,3]))*log(1-imatr$Pt[40]), NA, NA, NA)
          LL <- -sum(imatr$LL[1:41])
          imatr$exp <- imatr$dP * 150542
          MAPE <- mean(abs((imatr$sales[1:40]-imatr$exp[1:40])/imatr$sales[1:40]))
         
          
          z <- cbind(r = x1, alpha = x4, b_ads = x2, c = x3, LL = -LL, MAPE = MAPE)
          rownames(z) <- "parameters"
          
          
           }
          
        #Weibull, iMac Covar
        if (modelType == 2 && length(covariates) == 1 && covariates[1] == 'ni') {
          w2Cov <- function(x) { ## function to optimize
            x1 <- x[1] #lambda
            x2 <- x[2] # beta
            x3 <- x[3] # c
            x4 <- x[4] #b_iMac
            imatr$exb=exp(imatr$adv..bn.*x2 + imatr$iMac*x4)
            imatr$A=cumsum(((imatr$t^x3) - (c(0, imatr$t[1:43]))^x3)*imatr$exb)
            imatr$Pt<-1-exp(-x1*imatr$A)
            imatr$dP[2:nrow(imatr)]=diff(imatr$Pt)
            imatr$dP[1]=imatr$Pt[1]
            imatr$LL <- c((imatr$sales[1:40])*(log(imatr$dP)[1:40]), (150452-sum(imatr[1:40,3]))*log(1-imatr$Pt[40]), NA, NA, NA)
            -sum(imatr$LL[1:41])
          }
          initial_guess=c(0.005, 1.5, 1.5, 0.5)
          x <- optim(initial_guess, w2Cov)$par
          x1 <- x[1] #lamda
          x2 <- x[2] #b_ads
          x3 <- x[3] #c
          x4 <- x[4] #b_iMac
          imatr$exb=exp(imatr$adv..bn.*x2 + imatr$iMac*x4)
          imatr$A=cumsum(((imatr$t^x3) - (c(0, imatr$t[1:43]))^x3)*imatr$exb)
          imatr$Pt<-1-exp(-x1*imatr$A)
          imatr$dP[2:nrow(imatr)]=diff(imatr$Pt)
          imatr$dP[1]=imatr$Pt[1]
          imatr$LL <- c((imatr$sales[1:40])*(log(imatr$dP)[1:40]), (150452-sum(imatr[1:40,3]))*log(1-imatr$Pt[40]), NA, NA, NA)
          LL <- -sum(imatr$LL[1:41])
          imatr$exp <- imatr$dP * 150542
          MAPE <- mean(abs((imatr$sales[1:40]-imatr$exp[1:40])/imatr$sales[1:40]))
        
          z <- cbind(lambda = x1, b_ads = x2, c = x3, b_iMac = x4, LL = -LL, MAPE = MAPE)
          rownames(z) <- "parameters"
          
          
        }
  #Weibull, Steve Covar
  if (modelType == 2 && length(covariates) == 1 && covariates[1] == 'jd') {
    w2CovJobs <- function(x) { ## function to optimize
      x1 <- x[1] #lambda
      x2 <- x[2] # beta
      x3 <- x[3] # c
      x5 <- x[4] #b_steve
      x6 <- x[5] #gamma
      x7 <- x[6] #delta
      imatr$sJobs <- c(array(0, c(1,23)), 1-x6*(1-exp(-x7*abs(imatr$t-24)))[24:44])
      imatr$exb=exp(imatr$adv..bn.*x2 + imatr$sJobs*x5)
      imatr$A=cumsum(((imatr$t^x3) - (c(0, imatr$t[1:43]))^x3)*imatr$exb)
      imatr$Pt<-1-exp(-x1*imatr$A)
      imatr$dP[2:nrow(imatr)]=diff(imatr$Pt)
      imatr$dP[1]=imatr$Pt[1]
      imatr$LL <- c((imatr$sales[1:40])*(log(imatr$dP)[1:40]), (150452-sum(imatr[1:40,3]))*log(1-imatr$Pt[40]), NA, NA, NA)
      -sum(imatr$LL[1:41])
    }
    initial_guess=c(0.005, 1.5, 1.5, 0.5, 0.2, 14)
    optim(initial_guess, w2CovJobs)
    x1 <- optim(initial_guess, w2CovJobs)$par[1] #lamda
    x2 <- optim(initial_guess, w2CovJobs)$par[2] #b_ads
    x3 <- optim(initial_guess, w2CovJobs)$par[3] #c
    x5 <- optim(initial_guess, w2CovJobs)$par[4] #b_steve
    x6 <- optim(initial_guess, w2CovJobs)$par[5] #gamma
    x7 <- optim(initial_guess, w2CovJobs)$par[6] #delta
    imatr$sJobs <- c(array(0, c(1,23)), 1-x6*(1-exp(-x7*abs(imatr$t-24)))[24:44])
    imatr$exb=exp(imatr$adv..bn.*x2 + imatr$sJobs*x5)
    imatr$A=cumsum(((imatr$t^x3) - (c(0, imatr$t[1:43]))^x3)*imatr$exb)
    imatr$Pt<-1-exp(-x1*imatr$A)
    imatr$dP[2:nrow(imatr)]=diff(imatr$Pt)
    imatr$dP[1]=imatr$Pt[1]
    imatr$LL <- c((imatr$sales[1:40])*(log(imatr$dP)[1:40]), (150452-sum(imatr[1:40,3]))*log(1-imatr$Pt[40]), NA, NA, NA)
    LL <- -sum(imatr$LL[1:41])
    imatr$exp <- imatr$dP * 150542
    MAPE <- mean(abs((imatr$sales[1:40]-imatr$exp[1:40])/imatr$sales[1:40]))
    
  
    z <- cbind(lambda = x1, b_ads = x2, c = x3, b_jobs = x5, gamma = x6, delta = x7, LL = -LL, MAPE = MAPE)
    rownames(z) <- "parameters"
   
    
  }
  
  #Weibull, iMac & Steve Covars
  if (modelType == 2 && length(covariates) == 2) {
    w3Cov <- function(x) { ## function to optimize
      x1 <- x[1] #lambda
      x2 <- x[2] # beta
      x3 <- x[3] # c
      x4 <- x[4] #b_iMac
      x5 <- x[5] #b_steve
      x6 <- x[6] #gamma
      x7 <- x[7] #delta
      imatr$sJobs <- c(array(0, c(1,23)), 1-x6*(1-exp(-x7*abs(imatr$t-24)))[24:44])
      imatr$exb=exp(imatr$adv..bn.*x2 + imatr$iMac*x4 + imatr$sJobs*x5)
      imatr$A=cumsum(((imatr$t^x3) - (c(0, imatr$t[1:43]))^x3)*imatr$exb)
      imatr$Pt<-1-exp(-x1*imatr$A)
      imatr$dP[2:nrow(imatr)]=diff(imatr$Pt)
      imatr$dP[1]=imatr$Pt[1]
      imatr$LL <- c((imatr$sales[1:40])*(log(imatr$dP)[1:40]), (150452-sum(imatr[1:40,3]))*log(1-imatr$Pt[40]), NA, NA, NA)
      -sum(imatr$LL[1:41])
    }
    initial_guess=c(0.005, 1.5, 1.5, 0.5, 0.2, 0.2, 13.5)
    optim(initial_guess, w3Cov)
    x1 <- optim(initial_guess, w3Cov)$par[1] #lamda
    x2 <- optim(initial_guess, w3Cov)$par[2] #b_ads
    x3 <- optim(initial_guess, w3Cov)$par[3] #c
    x4 <- optim(initial_guess, w3Cov)$par[4] #b_iMac
    x5 <- optim(initial_guess, w3Cov)$par[5] #b_steve
    x6 <- optim(initial_guess, w3Cov)$par[6] #gamma
    x7 <- optim(initial_guess, w3Cov)$par[7] #delta
    imatr$sJobs <- c(array(0, c(1,23)), 1-x6*(1-exp(-x7*abs(imatr$t-24)))[24:44])
    imatr$exb=exp(imatr$adv..bn.*x2 + imatr$iMac*x4 + imatr$sJobs*x5)
    imatr$A=cumsum(((imatr$t^x3) - (c(0, imatr$t[1:43]))^x3)*imatr$exb)
    imatr$Pt<-1-exp(-x1*imatr$A)
    imatr$dP[2:nrow(imatr)]=diff(imatr$Pt)
    imatr$dP[1]=imatr$Pt[1]
    imatr$LL <- c((imatr$sales[1:40])*(log(imatr$dP)[1:40]), (150452-sum(imatr[1:40,3]))*log(1-imatr$Pt[40]), NA, NA, NA)
    LL <- -sum(imatr$LL[1:41])
    imatr$exp <- imatr$dP * 150542
    MAPE <- mean(abs((imatr$sales[1:40]-imatr$exp[1:40])/imatr$sales[1:40]))
    
    
    z <- cbind(lambda = x1, b_ads = x2, c = x3, b_iMac = x4, b_jobs = x5, gamma = x6, delta = x7, LL = -LL, MAPE = MAPE)
    rownames(z) <- "parameters"
  }
  
  
  plotMin <- min((imatr$dP[startQuarter[1]:startQuarter[2]] * 150452), (imatr[startQuarter[1]:(startQuarter[2]-1),3]), (imatr[startQuarter[1]:startQuarter[2],4]*5000))
  plotMax <- max((imatr$dP[startQuarter[1]:startQuarter[2]] * 150452), (imatr[startQuarter[1]:(startQuarter[2]-1),3]), (imatr[startQuarter[1]:startQuarter[2],4]*5000))
  
  taxiPlot <- plot (imatr[1:43,1], imatr$dP[1:43] * 150452, type="l", col="blue", xlab="Quarter", xlim = c(startQuarter[1], startQuarter[2]), ylab="iMac Sales (in 1,000s)", pch=16, cex=1.2,
                    lwd=2, ylim = c(0, plotMax))
  
  lines (imatr[1:43,1], imatr[1:43,3], pch=16, cex=1.2, lwd=2, col="red")
  lines (imatr[1:40,1], imatr[1:40,4]*5000, pch=16, cex=1.2, lwd=2, lty=2)
  segments(40, 1.8*5000, 41, num*5000, pch=16, cex=1.2, lwd=2, lty=2)
  segments(41, num*5000, 43, num*5000, pch=16, cex=1.2, lwd=2, lty=2)
  lines (within(data.frame(c(1:43)), ypred <- predict(lm((imatr$dP[1:43] * 150452)~poly(c(1:43),benchmark,raw=TRUE)), data.frame(c(1:43)))), col = "grey", pch=16, cex=1.2, lwd=2, lty=2)
  legend (startQuarter[1], plotMax, c("Projected Sales (model)", "Actual Sales", "Ad Expenditures (covariate)","Polynomial Regression (benchmark)"),
            lty = c(1,1,2,2), lwd = c(2,2,2,2),
            col = c("blue","red","black","grey"))
    
  z
}

description <- function(startQuarter, modelType, covariates, mixture) {
  #Weibull, iMac & Steve Covars
  if (mixture) {
    "Adding a mixing distribution has no impact on this particular dataset, implying relatively homogenous purchase behaviors among the customers. To visualize, try removing both extra covariates, and toggling the gamma mixture."
  }
  else {
    if (modelType == 1) {
      "Exponential distributions are essentially Weibull distributions that don't account for duration dependence, c, which measures how the time since most recent purchase impacts repurchase likelihood. Thus for any c significantly greater than 1, Weibull will have better fit (lower MAPE)."
    }
    else {
      if (modelType == 2 && length(covariates) == 2) {
        "With a Mean Average Percent Error (MAPE) of less than .13, the Weibull (with all 3 covariates) demonstrates better fit than all other probability models and traditional regression methods.

Note: At time of analysis (Q4 2016), Apple disclosed Q1-Q3 2016 sales but had not disclosed 2016 ad expenditures. The model initially assumes that ad expenditures in 2016 equal those in 2015 ($1.8bn); Apple could build the most accurate model by inputting their true 2016 ad exenditures."
      }
        else {
        "Unlike tradional regression methods (curve fitting), probability models predict sales by analyzing purchase behaviors of individual customers. Lambda measures the average customer's purchase propensity in a given period, betas measure the covariance between purchase likelihood and their given co-variable, Log Likelihood is the parameter we seek to maximize, and MAPE is a goodness-of-fit measure for comparing models."
        }
      }
    }
}

  
shinyServer(
  function(input, output) {
    output$plot <- renderPlot({tripsPlot(input$startQuarter, input$modelType, input$covariates, input$mixture, input$benchmark, input$num)})
    output$table <- renderTable({tripsPlot(input$startQuarter, input$modelType, input$covariates, input$mixture, input$benchmark, input$num)}, digits = 4)
    output$description <- renderText({description(input$startQuarter, input$modelType, input$covariates, input$mixture)})
    }
  )