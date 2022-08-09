# CS
#### R-code accompanying Chapter 6: testing for central symmetry

### Required packages
library(lessR)
library(copula)
library(QBAsyDist)
library(doParallel)
library(quantmod)
library(LaplacesDemon)
library(metRology)
library(QBAsyDist)
library(geometry)
library(mrfDepth)
library(mvtnorm)
library(MTS)
library(ggplot2)
library(gridExtra)
library(simukde)




## files:
Data Examples with tests.R: Tests for central symmetry on two datasets
Plots data examples.R: Plots of results data examples; requires data examples with tests to have run as this generates the results
simulation symmetry tests: Contains the simulation study used for the first three settings
simulation 4 symmetry tests: Contains the simulation study used in setting four
funtions.R: Contains all functions used in the previous scripts

skew normal copula: Sampler, density and fitter for a skew-normal copula (all in pseudo observations), also contained in functions.R
skew t copula: Sampler, density and fitter for a skew-normal copula (all in pseudo observations), also contained in functions.R
plot skew normal copula.R: Plots the pseudo-observations for 15 different setting in a skew-normal copula, used for illustrative purposes

GLRT LCQBA.R: example generalized likelihood ratio test for central symmetry in a linear combination of QBA-distributions model
GLRT SNC.R: example of GLRT for central symmetry in skew-normal copula for fit with assumed marginal symmetry 
GLRT SPC: example of GLRT for central symmetry in skew-normal copula for fit without assumed marginal symmetry
