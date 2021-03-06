#########################################################
# Mean Variance Optimisation -  Modern Portfolio Theory #
# Creating efficient frontiers for optimial allocations #
#########################################################
# Marcus Williamson 18/06/16 - Orig: Economist at Large #
#########################################################

library(stockPortfolio) # Base package for retrieving returns
library(ggplot2) # Used to graph efficient frontier
library(reshape2) # Used to manipulate the data
library(quadprog) #Needed for solve.QP, solving quadratic problems

setwd("~/Documents/Repos/AlternativeReturn/MVO") # Set directory

returns <- read.csv("Data/ar_monthly_returns_80.csv", header = TRUE) # Loan in our monthly returns

# Remove date and convert to numeric data frame
returns <- returns[,-1]
returns[, 1:4] <- sapply(returns[, 1:4], as.numeric)

#### Efficient Frontier function ###############################################
eff.frontier <- function (returns, short="no", max.allocation=NULL,
                          risk.premium.up=.5, risk.increment=.005){
  # return argument should be a m x n matrix with one column per security
  # short argument is whether short-selling is allowed; default is no (short
  # selling prohibited) max.allocation is the maximum % allowed for any one
  # security (reduces concentration) risk.premium.up is the upper limit of the
  # risk premium modeled (see for loop below) and risk.increment is the
  # increment (by) value used in the for loop
  
  covariance <- cov(returns)
  #print(covariance)
  n <- ncol(covariance)
  
  # Create initial Amat and bvec assuming only equality constraint
  # (short-selling is allowed, no allocation constraints)
  Amat <- matrix (1, nrow=n)
  bvec <- 1
  meq <- 1
  
  # Then modify the Amat and bvec if short-selling is prohibited
  if(short=="no"){
    Amat <- cbind(1, diag(n))
    bvec <- c(bvec, rep(0, n))
  }
  
  # And modify Amat and bvec if a max allocation (concentration) is specified
  if(!is.null(max.allocation)){
    if(max.allocation > 1 | max.allocation <0){
      stop("max.allocation must be greater than 0 and less than 1")
    }
    if(max.allocation * n < 1){
      stop("Need to set max.allocation higher; not enough assets to add to 1")
    }
    Amat <- cbind(Amat, -diag(n))
    bvec <- c(bvec, rep(-max.allocation, n))
  }
  
  # Calculate the number of loops
  loops <- risk.premium.up / risk.increment + 1
  loop <- 1
  
  # Initialize a matrix to contain allocation and statistics
  # This is not necessary, but speeds up processing and uses less memory
  effd <- matrix(nrow=loops, ncol=n+3)
  # Now I need to give the matrix column names
  colnames(eff) <- c(colnames(returns), "Std.Dev", "Exp.Return", "sharpe")
  
  # Loop through the quadratic program solver
  for (i in seq(from=0, to=risk.premium.up, by=risk.increment)){
    dvec <- colMeans(returns) * i # This moves the solution along the EF
    sol <- solve.QP(2*covariance, dvec=dvec, Amat=Amat, bvec=bvec, meq=meq)
    eff[loop,"Std.Dev"] <- sqrt(sum(sol$solution*colSums((covariance*sol$solution))))
    eff[loop,"Exp.Return"] <- as.numeric(as.character((sol$solution %*% colMeans(returns))))
    eff[loop,"sharpe"] <- eff[loop,"Exp.Return"] / eff[loop,"Std.Dev"]
    eff[loop,1:n] <- sol$solution
    loop <- loop+1
  }
  
  return(as.data.frame(eff))
}
#### End of Efficient Frontier Function ######################################

# Create Simulations over a rolling time period

rows <- seq(1,ceiling(0.3*nrow(returns)))
cols <- ncol(returns)+3
optimalpoints <- matrix(nrow=ceiling(0.7*nrow(returns)), ncol=cols)
colnames(optimalpoints) <- c(colnames(returns), "Std.Dev", "Exp.Return", "sharpe")


for (i in seq(1,nrow(returns)*0.7,1)){
  
  returnsample <- returns[rows,]
  eff <- eff.frontier(returns=returnsample, short = "no", max.allocation =0.5,
                      risk.premium.up=0.5, risk.increment =.001) # calc frontier
  
  eff.optimal.point <- eff[eff$sharpe==max(eff$sharpe),] # find optimal point
 
   #something is not working here cannot use "optimalpoints[i,1:8] <- eff.optimal.point"...
  optimalpoints[i,1] <- eff.optimal.point[,1]  
  optimalpoints[i,2] <- eff.optimal.point[,2]
  optimalpoints[i,3] <- eff.optimal.point[,3]
  optimalpoints[i,4] <- eff.optimal.point[,4]
  optimalpoints[i,5] <- eff.optimal.point[,5]
  optimalpoints[i,6] <- eff.optimal.point[,6]
  optimalpoints[i,7] <- eff.optimal.point[,7]
  optimalpoints[i,8] <- eff.optimal.point[,8]
  
  rows = rows + 1
}

# Run the eff.frontier function based on no short and 50% alloc. restrictions
eff <- eff.frontier(returns=returns, short="no", max.allocation=0.5,
                    risk.premium.up=0.5, risk.increment=.001)

# Find the optimal portfolios
eff.optimal.point <- eff[eff$sharpe==max(eff$sharpe),]
eff.tenpercent <- eff[eff$sharpe==max(eff$sharpe[eff$Exp.Return>((1.1^(1/12))-1)]),]

# Add simulations to eff data

eff <- rbind(eff,optimalpoints) #bind them

# graph efficient frontier
# Start with color scheme
ealred <- "#7D110C"
ealtan <- "#CDC4B6"
eallighttan <- "#F7F6F0"
ealdark <- "#423C30"

#### Plot it
ggplot(eff, aes(x=Std.Dev, y=Exp.Return)) + geom_point(alpha=.1, color=ealdark) +

  #Plotting the optimal portfolio first
  geom_point(data=eff.optimal.point, aes(x=Std.Dev, y=Exp.Return, label=sharpe),
             color=ealred, size=5) +
  annotate(geom="text", x=eff.optimal.point$Std.Dev,
           y=eff.optimal.point$Exp.Return,
           label=paste("Optimal Portfolio\nRisk: ",
                       round(((1+eff.optimal.point$Std.Dev)^(12)-1)*100, digits=2),"%\nReturn: ",
                       round(((1+eff.optimal.point$Exp.Return)^(12)-1)*100, digits=2),"%\nSharpe: ",
                       round(eff.optimal.point$sharpe*sqrt(12), digits=4), sep=""),
           hjust=-0.2, vjust=-0.2) +
  
  #Plotting the optimised 10% return portfolio
  geom_point(data=eff.tenpercent, aes(x=Std.Dev, y=Exp.Return, label=sharpe),
             color=ealdark, size=5) +
  annotate(geom="text", x=eff.tenpercent$Std.Dev,
           y=eff.tenpercent$Exp.Return,
           label=paste("10% Return Optimised\nRisk: ",
                       round(((1+eff.tenpercent$Std.Dev)^(12)-1)*100, digits=2),"%\nReturn: ",
                       round(((1+eff.tenpercent$Exp.Return)^(12)-1)*100, digits=2),"%\nSharpe: ",
                       round(eff.tenpercent$sharpe*sqrt(12), digits=4), sep=""),
           hjust=-0.2, vjust=-0.2) +
  #End plotting
  
  ggtitle("Efficient Frontier and Optimal Portfolio by Sharpe Ratio") +
  labs(x="Monthly Risk (standard deviation of portfolio)", y="Monthly Return") +
  theme(panel.background=element_rect(fill=eallighttan),
        text=element_text(color=ealdark),
        plot.title=element_text(size=24, color=ealred))

# Ouput our optimal allocations and stats
meanoptimal <- colMeans(optimalpoints) # mean of simulations
meanoptimal
eff.optimal.point
eff.tenpercent

meanoptimal <- colMeans(optimalpoints)
