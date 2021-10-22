################################################################################
# Replication file for Cattaneo, Cheung, Ma, and Masatlioglu (2021)
# 20-Oct-2021
################################################################################

# This replication file provides R code to test against 4 preferences for 2000
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

# additional parameters, AOM, RAM, Att-at-binary
jobPars <- matrix(c(1, 0, 1,
                    0, 1, 1,
                    1, 1, 1,
                    1, 1, 0.9,
                    1, 1, 0.8,
                    1, 1, 0.7,
                    1, 1, 0.6), ncol=3, byrow=TRUE)

# list of preferences
prefList <- matrix(c(1, 2, 3, 4, 5, 6,
                     2, 3, 4, 5, 6, 1,
                     1, 2, 6, 5, 4, 3,
                     1, 6, 5, 4, 3, 2
                     ), ncol=6, byrow=T)

repe   <- 2000 # number of Monte Carlo repetitions
repeSS <- 2000 # number of simulations for critical values

################################################################################
# simulation
################################################################################

### Change the following to TRUE for full simulation
### This may take 10-15 hours
if (FALSE) {
  Result <- matrix(NA, nrow=repe, ncol=4*nrow(jobPars))
  
  ptm <- proc.time()
  for (i in 1:repe) {
    # generate data
    menu <- choice <- matrix(0, nrow=0, ncol=6)
    for (j in 6:2) {
      temp <- logitSimu(n = n, uSize = 6, mSize = j, a = varsig)
      menu <- rbind(menu, temp$menu); choice <- rbind(choice, temp$choice)
    }
    
    for (j in 1:nrow(jobPars)) {
      AOM <- TRUE; RAM <- FALSE; attBinary <- 1
      temp <- revealPref(menu = menu, choice = choice, pref_list = prefList, method = "2MS",
                         nCritSimu = repeSS,
                         BARatio2MS = 0.1, BARatio2UB = 0.1, MNRatioGMS = NULL,
                         AOM = (jobPars[j, 1] == 1),
                         RAM = (jobPars[j, 2] == 1),
                         limDataCorr = TRUE,
                         attBinary = jobPars[j, 3])
      Result[i, (j-1)*4 + (1:4)] <- (temp$critVal$MS[, 2] < temp$Tstat) * 1
    }
  }
  proc.time() - ptm
  
  colMeans(Result)
  write.table(Result, file="CCMM_2021_AttOverload--Replication-Preference.txt", sep=",", row.names=FALSE, col.names=FALSE)
}

### Simulation results are stored in the file, which can be used to generate the rejection probabilities
matrix(colMeans(read.table("CCMM_2021_AttOverload--Replication-Preference.txt", sep=",")), nrow=4)

