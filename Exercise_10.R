### Exercise_10

setwd("/Users/Samantha/Documents/5_5th-Year/Biocomputing/Exercises/Exercise_10/ICB2019_Exercise10/")

## In this exercise, you will use a population model to investigate evolution of drug
## resistance in tumors. 

## Imagine a cancer cell in a tumor that spontaneously exhibited a mutation that confers
## drug resistance. The mutation does not have any positive or begative effects on growth
## rate of that sub-population when the cancer drug is absent. However, when the cancer
## drug is present, the mutant sub-population grows at 50% of its growth rate in the absence
## of the drug and the non-mutant sub-population declines rapidly. 

# The model we will use to represent the growth of the two sub-populations is this:

  # N(t+1) = Nt + rN(Nt)(1-(Nt+Mt)/K)
  # M(t+1) = Mt + rM(Mt)(1-(Nt+Mt)/K)

# Assume in the absence of the cancer drug that the cells grow at a rate of 0.1 per 
# day (rN = rM = 0.1) and the carrying capacity (K) of the tumor is one million cells.
# The mutation of a single cell occured early in the tumor growth and when it occured
# there were 100 total cells in the tumor. Drug treatment of non-mutant cells results
# in a negative growth rate of -0.1.

## Generate a script that simulates growth of the two-sub populations in the tumor to
## equilibrium followed by drug treatment. Plot your results using a line graph. 

  # Set parameters

  timesteps = 1000 

## Simulation for Cell Population Growth
  
tumorSim <- function(N0=99,M0 = 1, rN= 0.1, rM = 0.1, K=1000000, timesteps=1000){ # Parameters from problem statement)
  Nt = numeric(length=timesteps) # Creates empty vectors for normal and mutant population
  Nt[1]=N0 # First value in the vector is the initial value 
  Mt = numeric(length=timesteps)
  Mt[1]=M0
  
  for(t in 1:(timesteps-1)){
    if (t >= 160){ # The population reaches carrying capacity around timestep 125; drug treatment starts here
      rN=-0.1 # Change growth rates after drug treatment
      rM= 0.1*.5
      Nt[t+1]= Nt[t] + rN*Nt[t]*(1-((Nt[t]+Mt[t])/K))
      Mt[t+1]= Mt[t] + rM*Mt[t]*(1-((Nt[t]+Mt[t])/K))
    }else { # Growth before drug treatment
      Nt[t+1]= Nt[t] + rN*Nt[t]*(1-((Nt[t]+Mt[t])/K))
      Mt[t+1]= Mt[t] + rM*Mt[t]*(1-((Nt[t]+Mt[t])/K))
    }
  }
  popOut <- data.frame(time=1:timesteps) # Create dataframe to store population data
  popOut$Nt <- Nt # Normal cells
  popOut$Mt <- Mt # Mutant cells
  return(popOut)
}

popOut <- tumorSim() # Run simulation
 
library(ggplot2) # Load ggplot tool


## Plot the two population numbers over time 

ggplot() + geom_line(popOut, mapping = aes(x = time, y = Nt), color = 'blue') +
  geom_line(popOut, mapping = aes(x = time, y = Mt), color = 'red') + 
  xlab("Timesteps") + ylab("Number of Cells") + theme_bw()


ggsave("tumorSimPlot.jpg", last_plot()) # Save the plot as a JPEG
