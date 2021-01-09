#=
MAIN INTERFACE
This code is a companion to the paper Surrogate Monte Carlo (SMC). It implements the main result of the article

=#

using Random
using CSV
using DataFrames
using StatsBase

include("CDFJulia.jl")
include("utilities.jl")
##
# Parameters
cgoal = 1e-3
nlag  = 50
n2lag = 5*nlag
##
# Data
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
# Step 3: Calculate the optimization function using white noise random variables. This is the starting point of the optimization
cmin = efun6(aR, bR, cR, bR2, rs, nlag, n2lag)
# Step 4: perform optimization: random permutation of rs should create final time series which minimizes the objective function
# and consequenty dynamically similar to the original SP
ra,E = annelOpt(aR,bR,cR,bR2,cmin,rs,nlag,n2lag,cgoal)

# How good is the opt ?
# Example:
B = autocor(abs.(ra), 1:n2lag)
plot(B)
plot!(bR)
