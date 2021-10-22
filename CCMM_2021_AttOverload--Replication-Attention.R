################################################################################
# Replication file for Cattaneo, Cheung, Ma, and Masatlioglu (2021)
# 20-Oct-2021
################################################################################

# This replication file provides R code for revealed attention for 2000
#   simulations. Effective sample size is 200.

################################################################################
# Load package
################################################################################

# The R code provided relies on the package "ramchoice", which is available on
#   CRAN. Run the following if it is not already installed.
#
#                 install.packages("ramchoice")

rm(list=ls(all=TRUE))
library("ramchoice")

################################################################################
# set seed
################################################################################
set.seed(42)

################################################################################
# varying parameters
################################################################################
n <- 200 # the effective sample size

varsig <- 2 # the parameter in the logit attention rule specification

# list of preferences
pref <- matrix(c(1, 2, 3, 4, 5, 6), ncol=6, byrow=T)

# choice problems
S <- matrix(c(1, 1, 1, 0, 0, 0,
              1, 1, 1, 0, 0, 1,
              1, 1, 1, 0, 1, 0,
              1, 1, 1, 1, 0, 0,
              1, 1, 1, 0, 1, 1,
              1, 1, 1, 1, 0, 1,
              1, 1, 1, 1, 1, 0,
              1, 1, 1, 1, 1, 1), ncol=6, byrow=TRUE)

repe   <- 2000 # number of Monte Carlo repetitions
repeSS <- 2000 # number of simulations for critical values

################################################################################
# simulation
################################################################################

### Change the following to TRUE for full simulation
### This may take 2-5 hours
if (FALSE) {
# three alternatives * 8 lower bounds * 8 upper bounds
Result <- matrix(NA, nrow=repe, ncol=3*nrow(S)*2)

ptm <- proc.time()
for (i in 1:repe) {
    # generate data
    menu <- choice <- matrix(0, nrow=0, ncol=6)
    for (j in 6:2) {
      temp <- logitSimu(n = n, uSize = 6, mSize = j, a = varsig)
      menu <- rbind(menu, temp$menu); choice <- rbind(choice, temp$choice)
    }

    for (j in 1:3) { # three alternatives
      AOM <- TRUE; RAM <- FALSE; attBinary <- 1
      temp <- revealAtte(menu = menu, choice = choice,
                         alternative = j, S = S,
                         lower = T, upper = T,
                         pref = pref)
      Result[i, (j-1)*nrow(S)*2           + (1:nrow(S))] <- temp$lowerBound[1, ]
      Result[i, (j-1)*nrow(S)*2 + nrow(S) + (1:nrow(S))] <- temp$upperBound[1, ]
    }
}
proc.time() - ptm

colMeans(Result)
write.table(Result, file="CCMM_2021_AttOverload--Replication-Attention.txt", sep=",", row.names=FALSE, col.names=FALSE)
}

### Simulation results are stored in the file, which can be used to generate the bounds on the attention frequencies
# lower bounds (95th percentile across simulations)
matrix(apply(read.table("CCMM_2021_AttOverload--Replication-Attention.txt", sep=",")[, c(1:8, 16+(1:8), 32+(1:8))], MARGIN=2, FUN=function(x) quantile(x, 0.95, na.rm=T)), nrow=3, byrow=T)
# upper bounds (5th percentile across simulations)
matrix(apply(read.table("CCMM_2021_AttOverload--Replication-Attention.txt", sep=",")[, c(8+(1:8), 24+(1:8), 40+(1:8))], MARGIN=2, FUN=function(x) quantile(x, 0.05, na.rm=T)), nrow=3, byrow=T)
