#=
                       ------ MAIN INTERFACE -----
This code is a companion to the paper Surrogate Monte Carlo (SMC) -arxiv 2021. This is the main
code used to generate SMC data. However the default setting here are selected
for the code to run withing minutes (<10 min).
Authors: Andrei Da Silva & A.Christian Silva
e-mail: csilva@idatafactory.com
version on 10/2020
=#

using Random
using CSV
using DataFrames
using StatsBase

include("CDFJulia.jl")  # random numbers generation
include("utilities.jl") # opt functions
##
# Parameters
cgoal = 1e-3  # opt goal
nlag  = 50    # max lag for autocorrelation
n2lag = nlag  # max lag for autocorrelation of abs()
maxN  = 1e6   # max numb of terms in loop
##
# Data: downloaded from yahoo finance
sp = CSV.read("spyData.csv", DataFrame)
##
# Algorithm
# Step 1: Generate random numbers that have the same unconditional probability distribution
Random.seed!(28)
rs = getrnd(sp.R, length(sp.R))
#Step 2: Calculate the components of the optimization function
aR  = autocor(sp.R, 1:nlag)
bR  = autocor(abs.(sp.R), 1:n2lag)
cR  = crosscor(sp.R, abs.(sp.R), -nlag:nlag)
bR2 = autocor((sp.R).^2, 1:n2lag)
# Step 3: Calculate the optimization function. This is the starting point of the optimization
cmin = efun6(aR, bR, cR, bR2, rs, nlag, n2lag)
# Step 4: Perform optimization. Random permutation of rs should create final time series
# which minimizes the objective function and it is consequenty dynamically similar to the original
# time series
# ra: artificial time series
# E:  optimization function value
@time ra,E = annelOpt(aR,bR,cR,bR2,cmin,rs,nlag,n2lag,cgoal,maxN)
##
# How good is the opt ?
# Example:
#(1) Compare the autocorrelation of the abs(ra) with the "real" autocorrelation
B = autocor(abs.(ra), 1:n2lag)
plot(B)
plot!(bR)
#(2) Compare the crosscorrelation of the ra and abs(ra) with the measured crosscorrelation
C  = crosscor(ra, abs.(ra), -nlag:nlag)
plot(C)
plot!(cR)
