GPflex

R package for estimating the unknown parameters of a Gaussian process with flexible correlation structures.

I. Introduction

The current version (v0.1.0) implements GP with three correlation structures: hybrid power-exponential, full power-exponential and hybrid matern structure. Please copy the following codes to R console to see the help files for inputs to the functions.

######BEGIN

>library(devtools)

>devtools::install_github("nelsonch/GPflex")

>library(GPflex)

##1 Hybrid power-exponential

>help(runBMCMC.power.hybrid)

##2 Full power-exponential

>help(runBMCMC.power.full)

##3 Hybrid matern

>help(runBMCMC.matern.hybrid)
######END

DEMO
A simple DEMO with three correlation structures (Power_Hybrid, Power_Full and Matern_Hybrid) is available at http://htmlpreview.github.io/?https://github.com/nelsonch/GPflex/blob/master/DEMO.html

To cite: Chen, H., Loeppky, J.L., and Welch, W.J. (2017), Flexible Correlation Structure for Accurate Prediction and Uncertainty Quantification in Bayesian Gaussian Process Emulation of a Computer Model, SIAM/ASA Journal on Uncertainty Quantification, 5 (1), 598-620.

Feel free to email the author (Hao (Nelson) Chen) should you have any questions at hao.chen@stat.ubc.ca
