library(shiny)

#Data Source, m3 sales: 
#Data Source, oil: http://www.macrotrends.net/1369/crude-oil-price-history-chart

salesmatrix <- read.csv("teslasales.csv")[1:16,]
range.ex <- function(x, fac=1.2) { xrg <- range(x);  m <- mean(xrg);  (xrg - m)*fac + m }

TeslaPlot <- function(startQuarter, modelType, covariates, benchmark) {
#total number of Model S sales, N = 129377
  
        if (modelType == 1) {
          #Exponential Model
          eCov <- function(x) { ## function to optimize
            x1 <- x[1] #lambda
            x2 <- x[2] # beta
            salesmatrix$exb=exp(salesmatrix[,7]*x2)
            salesmatrix$A=cumsum(salesmatrix$exb)
            salesmatrix$Pt<-1-exp(-x1*salesmatrix$A)
            salesmatrix$dP[2:nrow(salesmatrix)]=diff(salesmatrix$Pt)
            salesmatrix$dP[1]=salesmatrix$Pt[1]
            salesmatrix$LL <- c(salesmatrix$sales*log(salesmatrix$dP))
            -sum(salesmatrix$LL)
          }
          initial_guess=c(2.5, -4)
          res <- optim(initial_guess, eCov)
          x1 <- res$par[1] #lambda
          x2 <- res$par[2] # beta
          salesmatrix$exb=exp(salesmatrix[,7]*x2)
          salesmatrix$A=cumsum(salesmatrix$exb)
          salesmatrix$Pt<-1-exp(-x1*salesmatrix$A)
          salesmatrix$dP[2:nrow(salesmatrix)]=diff(salesmatrix$Pt)
          salesmatrix$dP[1]=salesmatrix$Pt[1]
          print(salesmatrix$dP[1:16] * 129377)
          salesmatrix$LL <- c(salesmatrix$sales*log(salesmatrix$dP))
          LL <- -sum(salesmatrix$LL)
          salesmatrix$exp <- salesmatrix$dP * 129377
          MAPE <- mean(abs((salesmatrix$sales[1:16]-salesmatrix$exp[1:16])/salesmatrix$sales[1:16]))
        
          z <- cbind(lambda = x1, b_oil = x2, LL = -LL, MAPE = MAPE)
          rownames(z) <- "parameters"
        }

  #Weibull
        if (modelType == 2) {
          wCov <- function(x) { ## function to optimize
            x1 <- x[1] #lambda
            x2 <- x[2] #beta
            x3 <- x[3] #c
            salesmatrix$exb=exp(salesmatrix[,7]*x2)
            salesmatrix$A=cumsum((salesmatrix$t^x3)*salesmatrix$exb)
            salesmatrix$Pt<-1-exp(-x1*salesmatrix$A)
            salesmatrix$dP[2:nrow(salesmatrix)]=diff(salesmatrix$Pt)
            salesmatrix$dP[1]=salesmatrix$Pt[1]
            salesmatrix$LL <- c(salesmatrix$sales*log(salesmatrix$dP))
            -sum(salesmatrix$LL)
          }
          initial_guess=c(.001, -2, 2)
          x <- optim(initial_guess, wCov)$par
          x1 <- x[1] #lambda
          x2 <- x[2] # beta
          x3 <- x[3] # c
          salesmatrix$exb=exp(salesmatrix[,7]*x2)
          salesmatrix$A=cumsum((salesmatrix$t^x3)*salesmatrix$exb)
          salesmatrix$Pt<-1-exp(-x1*salesmatrix$A)
          salesmatrix$dP[2:nrow(salesmatrix)]=diff(salesmatrix$Pt)
          salesmatrix$dP[1]=salesmatrix$Pt[1]
          salesmatrix$LL <- c(salesmatrix$sales*log(salesmatrix$dP))
          LL <- -sum(salesmatrix$LL)
          salesmatrix$exp <- salesmatrix$dP * 129377
          MAPE <- mean(abs((salesmatrix$sales[1:16]-salesmatrix$exp[1:16])/salesmatrix$sales[1:16]))
         
          z <- cbind(lambda = x1, b_oil = x2, c = x3, LL = -LL, MAPE = MAPE)
          rownames(z) <- "parameters"
          
          
        }
  
  #Weibull, Model X covar
  if (modelType == 3) {
    w2CovJobs <- function(x) { ## function to optimize
      x1 <- x[1] #lambda
      x2 <- x[2] # beta
      x3 <- x[3] # c
      x5 <- x[4] #b_mx
      x6 <- x[5] #gamma
      x7 <- x[6] #delta
      salesmatrix$mx <- c(array(0, c(1,12)), 1-x6*(1-exp(-x7*abs(salesmatrix$t-13)))[13:16])
      salesmatrix$exb=exp(salesmatrix[,7]*x2 + salesmatrix$mx*x5)
      salesmatrix$A=cumsum((salesmatrix$t^x3)*salesmatrix$exb)
      salesmatrix$Pt<-1-exp(-x1*salesmatrix$A)
      salesmatrix$dP[2:nrow(salesmatrix)]=diff(salesmatrix$Pt)
      salesmatrix$dP[1]=salesmatrix$Pt[1]
      salesmatrix$LL <- c(salesmatrix$sales*log(salesmatrix$dP))
      -sum(salesmatrix$LL)
    }
    initial_guess=c(.01, -.5, 2, 0.3, -83, .03)
    optim(initial_guess, w2CovJobs)
    x1 <- optim(initial_guess, w2CovJobs)$par[1] #lamda
    x2 <- optim(initial_guess, w2CovJobs)$par[2] #b_ads
    x3 <- optim(initial_guess, w2CovJobs)$par[3] #c
    x5 <- optim(initial_guess, w2CovJobs)$par[4] #b_mx
    x6 <- optim(initial_guess, w2CovJobs)$par[5] #gamma
    x7 <- optim(initial_guess, w2CovJobs)$par[6] #delta
    salesmatrix$mx <- c(array(0, c(1,12)), 1-x6*(1-exp(-x7*abs(salesmatrix$t-13)))[13:16])
    salesmatrix$exb=exp(salesmatrix[,7]*x2 + salesmatrix$mx*x5)
    salesmatrix$A=cumsum(salesmatrix$t^x3*salesmatrix$exb)
    salesmatrix$Pt<-1-exp(-x1*salesmatrix$A)
    salesmatrix$dP[2:nrow(salesmatrix)]=diff(salesmatrix$Pt)
    salesmatrix$dP[1]=salesmatrix$Pt[1]
    salesmatrix$LL <- c(salesmatrix$sales*log(salesmatrix$dP))
    LL <- -sum(salesmatrix$LL)
    salesmatrix$exp <- salesmatrix$dP * 129377
    MAPE <- mean(abs((salesmatrix$sales[1:16]-salesmatrix$exp[1:16])/salesmatrix$sales[1:16]))
    
    
    z <- cbind(lambda = x1, b_oil = x2, c = x3, b_mX = x5, gamma = x6, delta = x7, LL = -LL, MAPE = MAPE)
    rownames(z) <- "parameters"
    
    
  }
  
  plotMin <- min((salesmatrix$dP[startQuarter[1]:startQuarter[2]] * 129377), (salesmatrix[startQuarter[1]:(startQuarter[2]-1),5]), (salesmatrix[startQuarter[1]:startQuarter[2],6]*100))
  plotMax <- max((salesmatrix$dP[startQuarter[1]:startQuarter[2]] * 129377), (salesmatrix[startQuarter[1]:(startQuarter[2]-1),5]), (salesmatrix[startQuarter[1]:startQuarter[2],6]*100))
  
  salesPlot <- plot (c(1:16), salesmatrix$dP[1:16] * 129377, type="l", col="blue", xlab="Quarter", xlim = c(startQuarter[1], startQuarter[2]), ylab="Number of Model S sold", pch=16, cex=1.2,
                    lwd=2, ylim = c(0, plotMax))
  print (salesmatrix$dP[1:16] * 129377)
  lines (c(1:16), salesmatrix[1:16,5], pch=16, cex=1.2, lwd=2, col="red")
  lines (c(1:16), salesmatrix[1:16,6]*100, pch=16, cex=1.2, lwd=2, lty=2)
  legend (startQuarter[1], plotMax, c("Projected Sales (model)", "Actual Sales", "Oil Price (covariate)"),
            lty = c(1,1,2), lwd = c(2,2,2),
            col = c("blue","red","black"))
    
  z
}


  
shinyServer(
  function(input, output) {
    output$plot <- renderPlot({TeslaPlot(input$startQuarter, input$modelType, input$covariates, input$benchmark)})
    output$table <- renderTable({TeslaPlot(input$startQuarter, input$modelType, input$covariates, input$benchmark)}, digits = 4)
    }
  )