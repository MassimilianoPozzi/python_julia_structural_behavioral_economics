### A Pluto.jl notebook ###
# v0.14.0

using Markdown
using InteractiveUtils

# ╔═╡ 030f5512-2a66-11ec-1a7a-21059e6facac
begin

# Import packages used

using Random
using DataFrames
using CSV
using StatFiles
using Statistics
	
end

# ╔═╡ af9b5960-2a65-11ec-3015-53185d92610c
# Replication of DellaVigna and Pope, 2018, "What Motivates Effort? Evidence and Expert Forecasts"

# The code in this notebook performs the estimates to replicate columns (1) (3) panel A and columns (1) (4) panel B of Table 5 using minimum distance.

# ╔═╡ 03372860-2a66-11ec-2b31-298b2c9afc9d
# Authors: 
# Massimiliano Pozzi (Bocconi University, pozzi.massimiliano@studbocconi.it)
# Salvatore Nunnari (Bocconi University, https://snunnari.github.io, salvatore.nunnari@unibocconi.it)

# ╔═╡ 7c88c030-55f8-11ec-0159-99a5b81bf850
# This notebook works with both Neptune and Pluto; the "begin" and "end" tags in each cell are necessary for Pluto notebooks but do not play any role in Neptune 

# ╔═╡ 7dcdea60-55f8-11ec-2fdc-1d683f29a01c
# This notebook was tested with the following packages versions:

# Pozzi: julia 1.6.1, Pluto 0.14.7, Neptune 0.14.0, DataFrames 1.2.2, CSV 0.9.5

# ╔═╡ 03008800-2a66-11ec-0c80-1767685eb4fa
begin

# 1. Data Cleaning and Data Preparation

# We import the relevant dataset containing data on the number of buttonpresses in the different treatments and for different piece rates wage that the participants received when completing the task. We then compute the means for the treatments that are used to compute the minimum distance estimates.

dt = DataFrame(load(raw"../input/mturk_clean_data_short.dta"))

emp_moments = round.(combine(groupby(dt, :treatment), :buttonpresses => mean)[2], digits=0)
	
end;

# ╔═╡ 944bb54e-55f9-11ec-32b3-c953374ecd8f
# We set up the Bootstrap procedure for the standard errors. The following cell creates 2000 new samples and computes the mean for the relevant treatments. 

# ╔═╡ 030039e0-2a66-11ec-3766-ad9593bffbd3
begin

# Treatments are as follows:
# 1.1: benchmark specification with piece rate of 0.01  
# 1.2: benchmark specification with piece rate of 0.10 
# 1.3: benchmark specification with piece rate of 0.00 
# 3.1: social preferences (charity) with piece rate of 0.01 
# 3.2: social preferences (charity) with piece rate of 0.10  
# 10 : social preferences (gift exchange) bonus of 40 cents (independently of nr buttonpresses)
# 4.1: time discounting with extra 0.01 paid two weeks later
# 4.2: time discounting with extra 0.01 paid four weeks later

# define the vector containing the treatments (used for the loop) and the vectors that will store the rounded mean buttonpresses in our 2000 new samples. 
# E11 is a 1x1000 vector that will contain the average buttonpresses in each of our 1000 samples for treatment 1.1 etc.


treatment = ["1.1", "1.2", "1.3", "3.1", "3.2", "10", "4.1", "4.2"]
E11, E12, E13, E31, E32, E10, E41, E42 = Int32[], Int32[], Int32[], Int32[], Int32[], Int32[], Int32[], Int32[]

# We first get a smaller dataframe containing only observations for a specific treatment, we then resample the observations and compute the rounded mean of buttonpresses, save the result and then pass onto the next treatment. we do this 1000 times.

for i=1:1000                 # we want 1000 new samples
    for treat in treatment   # we want to compute the mean for all treatments
        db = copy(dt)        # copy the original dataset
        db = db[db[:,2] .== treat,:] # keep only observation for treatment 'treat'
        bootsample = db[vec(rand(1:size(db[!,1],1),1,size(db[!,1],1))),1] # resample the dataset by taking random rows from the first column
		if treat == "1.1" append!(E11, round(mean(bootsample))) end
        if treat == "1.2" append!(E12, round(mean(bootsample))) end
        if treat == "1.3" append!(E13, round(mean(bootsample))) end
        if treat == "3.1" append!(E31, round(mean(bootsample))) end
        if treat == "3.2" append!(E32, round(mean(bootsample))) end
        if treat == "10"  append!(E10, round(mean(bootsample))) end
        if treat == "4.1" append!(E41, round(mean(bootsample))) end
        if treat == "4.2" append!(E42, round(mean(bootsample))) end
	end
end
	
end

# ╔═╡ ced0482e-2a66-11ec-06ca-0b961f35a057
# For a more detailed explanation of the model we refer the reader to the python notebook

# Define the function to estimate the parameters using minimum distance:

# There are two specifications: exponential cost function and power cost function
# The parameters are found by imposing the theoretical means found from the agent's problem being equal to the empirical means found in the data.
# E11 up to E42 are the relevant empirical moments
# Specification can either be 'Exp' or 'Power'
# P is a vector containing the different piece-rates
# To simplify the computations the authors take the log of the variables to estimate gamma k and s.

function mindisest(E11,E12,E13,E31,E32,E10,E41,E42,specification)
	
    P = [0, 0.01, 0.1]
	
	# gamma, k and s are found by solving the system of N equations in N unknowns
	
    if specification == "Exp"	
        log_k = (log(P[3]) - log(P[2])*(E12)/(E11))/(1 - (E12)/(E11))
        log_gamma = log((log(P[2]) - log_k)/(E11))
        log_s = exp(log_gamma)*(E13) + log_k	
	end
	
    if specification == "Power"	
        log_k = (log(P[3]) - log(P[2])*log(E12)/log(E11))/(1 - log(E12)/log(E11))
        log_gamma = log((log(P[2]) - log_k)/log(E11))
        log_s = exp(log_gamma)*log(E13) + log_k	
	end
	
    k = exp(log_k)
    g = exp(log_gamma)
    s = exp(log_s)
	
    if specification == "Exp"
        EG31, EG32, EG10, EG41, EG42 = exp.(E31*g), exp.(E32*g), exp.(E10*g), exp.(E41*g), exp.(E42*g)
		
        alpha = 100/9*k*(EG32-EG31)
        a = 100*k*EG31-100*s-alpha
        s_ge = k*EG10 - s
        delta = sqrt((k*EG42-s)/(k*EG41-s))
        beta  = 100*(k*EG41-s)/(delta^2)
	end
	
    if specification == "Power"
        alpha = 100/9*k*(E32^g-E31^g)
        a = 100*k*E31^g-100*s-alpha
        s_ge = k*E10^g - s
        delta = sqrt((k*E42^g-s)/(k*E41^g-s))
        beta  = 100*(k*E41^g-s)/(delta^2)
	end
	
    return k, g, s, alpha, a, s_ge, beta, delta
end

# ╔═╡ 9ac87bc0-5608-11ec-3f39-3713d508fdfb
begin

# Table 5 minimum distance estimates: columns (1) (3) panel A and columns (1) (4) panel B
    
Table5Exp = mindisest(emp_moments[1],emp_moments[2],emp_moments[3],emp_moments[7],emp_moments[8],emp_moments[5],emp_moments[9],emp_moments[10],"Exp")

Table5Power = mindisest(emp_moments[1],emp_moments[2],emp_moments[3],emp_moments[7],emp_moments[8],emp_moments[5],emp_moments[9],emp_moments[10],"Power")
	
end;

# ╔═╡ 4b1369f0-2a70-11ec-2cf3-e514a4899e06
begin

# Store mean and standard error of estimates for the exponential cost function specification. Our function returns an array containing 1000 vectors of 8 elements each. The estimate for k is in position 1 in each 1000 vectors etc

estimatesExp = mindisest.(E11,E12,E13,E31,E32,E10,E41,E42,"Exp")
k_exp_mean = mean([item[1] for item in estimatesExp])
k_exp_sd = std([item[1] for item in estimatesExp])
g_exp_mean = mean([item[2] for item in estimatesExp])
g_exp_sd = std([item[2] for item in estimatesExp])
s_exp_mean = mean([item[3] for item in estimatesExp])
s_exp_sd = std([item[3] for item in estimatesExp])
alpha_exp_mean = mean([item[4] for item in estimatesExp])
alpha_exp_sd = std([item[4] for item in estimatesExp])
a_exp_mean = mean([item[5] for item in estimatesExp])
a_exp_sd = std([item[5] for item in estimatesExp])
s_ge_exp_mean = mean([item[6] for item in estimatesExp])
s_ge_exp_sd = std([item[6] for item in estimatesExp])
beta_exp_mean = mean([item[7] for item in estimatesExp])
beta_exp_sd = std([item[7] for item in estimatesExp])
delta_exp_mean = mean([item[8] for item in estimatesExp])
delta_exp_sd = std([item[8] for item in estimatesExp])
	
end;

# ╔═╡ 2e15f490-2a6f-11ec-3bec-0f631d9fbd70
begin

# Store mean and standard error of estimates for the power cost function specification

estimatesPower = mindisest.(E11,E12,E13,E31,E32,E10,E41,E42,"Power")
k_power_mean = mean([item[1] for item in estimatesPower])
k_power_sd = std([item[1] for item in estimatesPower])
g_power_mean = mean([item[2] for item in estimatesPower])
g_power_sd = std([item[2] for item in estimatesPower])
s_power_mean = mean([item[3] for item in estimatesPower])
s_power_sd = std([item[3] for item in estimatesPower])

# A couple of elements in the arrays take value "NaN". We need to remove them before computing the mean and std

alpha_power = [item[4] for item in estimatesPower]
alpha_power_mean = mean(alpha_power[(!isnan).(alpha_power)])
alpha_power_sd = std(alpha_power[(!isnan).(alpha_power)])
a_power = [item[5] for item in estimatesPower]
a_power_mean = mean(a_power[(!isnan).(a_power)])
a_power_sd = std(a_power[(!isnan).(a_power)])
s_ge_power = [item[6] for item in estimatesPower]
s_ge_power_mean = mean(s_ge_power[(!isnan).(s_ge_power)])
s_ge_power_sd = std(s_ge_power[(!isnan).(s_ge_power)])
beta_power = [item[7] for item in estimatesPower]
beta_power_mean = mean(beta_power[(!isnan).(beta_power)])
beta_power_sd = std(beta_power[(!isnan).(beta_power)])
delta_power = [item[8] for item in estimatesPower]
delta_power_mean = mean(delta_power[(!isnan).(delta_power)])
delta_power_sd = std(delta_power[(!isnan).(delta_power)])
	
end;

# ╔═╡ c2d506c0-2a6f-11ec-0871-fdb55e56de67
begin

# To obtain confidence intervals in panel B table 5. CI are the 2.5% and 97.5% quantiles of the distribution of our parameters vectors. We obtain CI only for alpha, a, s_ge, beta, delta as in the paper

CI_Exp, CI_Power = [], []
for ci=4:8
    a  = sort!([item[ci] for item in estimatesExp]) # sort from smallest to biggest number
    append!(CI_Exp, [quantile!(a,0.025) ,quantile!(a,0.975)])
end
vecest = [alpha_power, a_power, s_ge_power, beta_power, delta_power]
for ci in vecest
    a  = sort!(ci[(!isnan).(ci)]) 
    append!(CI_Power, [quantile!(a,0.025) ,quantile!(a,0.975)])
end
	
end

# ╔═╡ a9cedae0-2a7c-11ec-3b2f-a3577f3decf7
begin

# Print our results for Table 5 (GMM estimates)
	
# Standard errors are different since the seed we used for the bootstrap procedure is different from the one used by the authors since random generation across softwares/languages is not easily replicated (each software uses its own algorithm)

sd_exp = [k_exp_sd,g_exp_sd,s_exp_sd,alpha_exp_sd,a_exp_sd,s_ge_exp_sd,beta_exp_sd,delta_exp_sd]

sd_power = [k_power_sd,g_power_sd,s_power_sd,alpha_power_sd,a_power_sd,s_ge_power_sd,beta_power_sd,delta_power_sd]

params_name = ["Level k of cost of effort", "Curvature γ of cost function","Intrinsic motivation s","Social preferences α", "Warm glow coefficient a","Gift exchange Δs", "Present bias β","(Weekly) discount factor δ"]
  
Table5Results = DataFrame()
Table5Results.params = params_name
	
# Table5Power is saved as a NTuple{8,Float64} but DataFrame only accepts vectors as columns. Using Collect transforms the NTuple into a vector
	
Table5Results.estimatesPower = collect(Table5Power)
Table5Results.sdPower = sd_power
Table5Results.estimatesExp = collect(Table5Exp)
Table5Results.sdExp = sd_exp
Table5Results
	
end

# ╔═╡ 88d81640-5609-11ec-1234-23840357ca64
# Save Table5Results as a csv file

CSV.write(raw"../output/table5GMM_julia.csv", Table5Results)

# ╔═╡ Cell order:
# ╠═af9b5960-2a65-11ec-3015-53185d92610c
# ╠═03372860-2a66-11ec-2b31-298b2c9afc9d
# ╠═7c88c030-55f8-11ec-0159-99a5b81bf850
# ╠═7dcdea60-55f8-11ec-2fdc-1d683f29a01c
# ╠═030f5512-2a66-11ec-1a7a-21059e6facac
# ╠═03008800-2a66-11ec-0c80-1767685eb4fa
# ╠═944bb54e-55f9-11ec-32b3-c953374ecd8f
# ╠═030039e0-2a66-11ec-3766-ad9593bffbd3
# ╠═ced0482e-2a66-11ec-06ca-0b961f35a057
# ╠═9ac87bc0-5608-11ec-3f39-3713d508fdfb
# ╠═4b1369f0-2a70-11ec-2cf3-e514a4899e06
# ╠═2e15f490-2a6f-11ec-3bec-0f631d9fbd70
# ╠═c2d506c0-2a6f-11ec-0871-fdb55e56de67
# ╠═a9cedae0-2a7c-11ec-3b2f-a3577f3decf7
# ╠═88d81640-5609-11ec-1234-23840357ca64
