---
title: "Portfolio Optimisation with alternative asset classes using MVO"
author: "Marcus Williamson"
date: "18 June 2016"
output: pdf_document
---

In order to construct the fund's portfolio efficiently over at [AlternativeReturn](www.alternativereturn.com), for this optimisation [Modern Portfolio Theory](https://en.wikipedia.org/wiki/Modern_portfolio_theory) was called upon. 

This piece is roughly broken down into three main areas:

+ Data collection
+ Run data through optimisation script
+ Review the output

### Data Collection

Monthly returns had to be collected for the assets that were to be invested in. These returns would be fed into the optimisation script.

The assets in question during this optimisation were:

+ Crowdfunded Debt (Commercial, Personal, Hybrid & Property based)
+ Crowdfunded Property
+ Whisky 
+ Gold
+ Stocks & Bonds

For __Crowdfunded Debt__ I used a [CDS measure](http://www.alternativereturn.com/#!bad-debt-rates/tlttb) and scaled it to match the historic default rates of the lending platforms, then backed out from this the returns that would have been seen if you had held these loans back to 2005 with reinvestment. There is an extent of guesswork here to generate some of the returns.

For __Crowdfunded Property__ I used the [Halifax House Price Index](http://www.halifax.co.uk/house-price-index/) data, then added on rental yields to to reach the rates seen on the lending platforms currently, this was again simulated historically to achieve a monthly returns series from 2005.

For __Whisky__ I used a [composite index](https://www.whiskystats.net/wwi/) of all Whisky auction prices which have taken place since 2005, I have then adjustd to include the frictional monthly storage fee, we are left with an conservative average annual return as compared to research of the platform which I use to buy & sell Whisky.

For __Gold__ I used the [price of Gold](https://www.bullionvault.com/gold-price-chart.do) in GBP since 2005, nothing clever here!

For __Stocks & Bonds__ I used the Vanguard LifeStrategy Growth Acc [fund](https://www.quandl.com/data/YAHOO/FUND_VASGX-VASGX-Vanguard-Lifestrategy-Growth-Fd) which represents an 80% equity 20% Bonds mix.

We end up with something like this:

```{r, echo=FALSE}
returns <- read.csv("Data/ar_monthly_returns_80.csv",header=TRUE)
head(returns)
```


### Optimisation Script

Now turning to the script to optimise, I must give credit to the [Economist at Large](http://economistatlarge.com/portfolio-theory) for his work here and R script!

Getting setup:
```{r}
library(stockPortfolio) # Base package for retrieving returns
library(ggplot2) # Used to graph efficient frontier
library(reshape2) # Used to manipulate the data
library(quadprog) #Needed for solve.QP, solving quadratic problems

setwd("~/Documents/Repos/AlternativeReturn/MVO") # Set directory

returns <- read.csv("Data/ar_monthly_returns_80.csv", header = TRUE) # Loan in our monthly returns

# Remove date and convert to numeric data frame
returns <- returns[,-1]
returns[, 1:4] <- sapply(returns[, 1:4], as.numeric)
```

The efficient frontier algorithm:
```{r}
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
   eff <- matrix(nrow=loops, ncol=n+3)
   # Now I need to give the matrix column names
   colnames(eff) <- c(colnames(returns), "Std.Dev", "Exp.Return", "sharpe")
   
   # Loop through the quadratic program solver
   for (i in seq(from=0, to=risk.premium.up, by=risk.increment)){
     dvec <- colMeans(returns) * i # This moves the solution along the EF
     sol <- solve.QP(covariance, dvec=dvec, Amat=Amat, bvec=bvec, meq=meq)
     eff[loop,"Std.Dev"] <- sqrt(sum(sol$solution*colSums((covariance*sol$solution))))
     eff[loop,"Exp.Return"] <- as.numeric(as.character((sol$solution %*% colMeans(returns))))
     eff[loop,"sharpe"] <- eff[loop,"Exp.Return"] / eff[loop,"Std.Dev"]
     eff[loop,1:n] <- sol$solution
     loop <- loop+1
   }
   
   return(as.data.frame(eff))
 }
```

A couple of things to note here, you can find the full breakdown and method on [Wikipedia](https://en.wikipedia.org/wiki/Modern_portfolio_theory), however in terms of the key variables:
+ Covariance - The covariance matrix calculated from returns data
+ dvec - Vector of average returns for each security, to find minimum portfolio variance we set these to zero
+ Amat - A matrix of constraints, such as the weights must sum to 1, no negative weights (if no short selling) and limiting the max weight of any security (if its been set)
+ bvec - This is a vector of values which forms a constraint in the minimisation problem, it is set to zero by default
+ meq - This tell the solve.QP function which columns of the Amat matrix to treat as equality constraints, for our example this is only 1, where all weights must sum to 1, the other remaining constraints are inequalities

We go on now to use this function we have create:
```{r}
eff <- eff.frontier(returns=returns, short="no", max.allocation=0.5,
                    risk.premium.up=0.5, risk.increment=.001)
head(eff)
```

We can see that for each allocation for the securities there are risk, return and Sharp ratio metrics calculated, these are monthly figures.

We go on to extract the optimal portfolios from the frontier generated on the conditions of:

+ Greatest Sharpe ratio
+ Greatest Sharpe ratio with returns over 10%

```{r}
eff.optimal.point <- eff[eff$sharpe==max(eff$sharpe),]
eff.tenpercent <- eff[eff$sharpe==max(eff$sharpe[eff$Exp.Return>((1.1^(1/12))-1)]),]
```

We now have all we need to review our optimal portfolios!


### Reviewing Findings

First we shall use ggplot to generate a seductive chart:
```{r}
# graph efficient frontier
# Start with color scheme
 ealred <- "#7D110C"
 ealtan <- "#CDC4B6"
 eallighttan <- "#F7F6F0"
 ealdark <- "#423C30"
 
 # Plot it
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

   ggtitle("Efficient Frontier and Optimal Portfolios") +
   labs(x="Monthly Risk (standard deviation of portfolio)", y="Monthly Return") +
   theme(panel.background=element_rect(fill=eallighttan),
         text=element_text(color=ealdark),
         plot.title=element_text(size=24, color=ealred))
```

We see the chart above shows us the efficient frontier, we seek to almost always remain on this frontier as this is a simplified view of the best return obtainable for a unit of risk. 

We have the optimal portfolio plotted in red here and target return optimal portfolio plotted in black, now we shall review the corresponding suggested allocations.

For the most optimal portfolio (with highest Sharpe Ratio) we find the following allocations suggested:
```{r,echo=FALSE}
eff.optimal.point 
```

For a portfolio yielding over 10% per year with the highest Sharpe ratio we find the following suggested allocations:
```{r, echo=FALSE}
eff.tenpercent
```

__There is a couple of things to note here:__

+ These figures are monthly stats, on the prior chart I suitably scaled them up to annual, the Sharpe ratio scaling was done with sqrt(12) which is open to debate...
+ This optimisation will always look to max out Debt, this is because its incredibly stable (Sharpe ratio > 5), so it has been capped at 50%, I will look to revisit this in the future
+ This is a very simple optimisation, it does not take into account many factors which would dicatate a final portoflio allocation, but it is a good starting point!

I am yet to see this thinking or approach replicated elsewhere on the internet, it may be because its a waste of time in the eyes of some - but investment returns dont lie...

I believe MPT provides a great framework to analyse securities / assets from a multidimentional perspective which includes covariance which is really something to not be overlooked.

__I hope you found the above useful, please get in contact with [myself](https://uk.linkedin.com/in/williamsonmarcus), if you have any questions or suggestions__