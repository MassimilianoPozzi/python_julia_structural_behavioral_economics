### A Pluto.jl notebook ###
# v0.14.0

using Markdown
using InteractiveUtils

# ╔═╡ 086c95a0-279f-11ec-3378-6924b1d7777e
begin

# Import packages used

using Distributed
using DelimitedFiles
using DataFrames
using Optim
using CSV
using LsqFit
using ForwardDiff
using LinearAlgebra
using StatFiles
	
end

# ╔═╡ aa3922f0-279e-11ec-33b6-6bb73cf7c16a
# Replication of DellaVigna and Pope, 2018, "What Motivates Effort? Evidence and Expert Forecasts"

# The code in this notebook estimates the parameters of column (2) (4) panel A and columns (3) (6) panel B of "Table 5: Estimates of behavioural parameters I: Mturkers actual effort and expert beliefs" and table 6 panel A: "Estimates of reference-dependent parameters: Mturker actual effort and expert beliefs" using Non-Linear-Leat-Squares

# ╔═╡ 08a0ec12-279f-11ec-12a3-f15a088873f4
# Authors: 
# Massimiliano Pozzi (Bocconi University, pozzi.massimiliano@studbocconi.it)
# Salvatore Nunnari (Bocconi University, https://snunnari.github.io, salvatore.nunnari@unibocconi.it)

# ╔═╡ 6783e0c0-5904-11ec-2a41-5784eec17590
# This notebook works with both Neptune and Pluto; the "begin" and "end" tags in each cell are necessary for Pluto notebooks but do not play any role in Neptune 

# ╔═╡ 903cd5d0-5904-11ec-0a95-75a9f8821a44
# This notebook was tested with the following packages versions:

# Pozzi: julia 1.6.1, Pluto 0.14.7, Neptune 0.14.0, DataFrames 1.2.2, CSV 0.9.5

# ╔═╡ 086c6e90-279f-11ec-2ce7-ef05700c68b3
begin

# 1. Data Cleaning and Data Preparation

# We import the relevant dataset containing data on the number of buttonpresses in the different treatments and for different piece rates wage that the participants received to complete the task. We then create a series of variables that are needed for estimation.

dt = DataFrame(load(raw"../input/mturk_clean_data_short.dta"))
	
# Create new variables needed for estimation:

# Create payoff per 100 button presses (p)

dt.payoff_per_100 = 0.0
dt[(dt[:treatment] .== "1.1"), :payoff_per_100] = 0.01
dt[(dt[:treatment] .== "1.2"), :payoff_per_100] = 0.1
dt[(dt[:treatment] .== "1.3"), :payoff_per_100] = 0.0
dt[(dt[:treatment] .== "2"),   :payoff_per_100] = 0.001
dt[(dt[:treatment] .== "1.4"), :payoff_per_100] = 0.04
dt[(dt[:treatment] .== "4.1"), :payoff_per_100] = 0.01
dt[(dt[:treatment] .== "4.2"), :payoff_per_100] = 0.01
dt[(dt[:treatment] .== "6.2"), :payoff_per_100] = 0.02
dt[(dt[:treatment] .== "6.1"), :payoff_per_100] = 1.0

# (ALPHA) create payoff per 100 to charity and dummy charity

dt.payoff_charity_per_100 = 0.0
dt[(dt[:treatment] .== "3.1"), :payoff_charity_per_100] = 0.01
dt[(dt[:treatment] .== "3.2"), :payoff_charity_per_100] = 0.1
dt.dummy_charity = 0
dt[(dt[:treatment] .== "3.1"), :dummy_charity] = 1
dt[(dt[:treatment] .== "3.2"), :dummy_charity] = 1

# (BETA/DELTA) create payoff per 100 delayed by 2 weeks and dummy delay

dt.delay_wks = 0
dt[(dt[:treatment] .== "4.1"), :delay_wks] = 2
dt[(dt[:treatment] .== "4.2"), :delay_wks] = 4
dt.delay_dummy = 0
dt[(dt[:treatment] .== "4.1"), :delay_dummy] = 1
dt[(dt[:treatment] .== "4.2"), :delay_dummy] = 1

# probability weights to back out curvature and dummy

dt.prob = 1.0
dt[(dt[:treatment] .== "6.2"), :prob] = 0.5
dt[(dt[:treatment] .== "6.1"), :prob] = 0.01
dt.weight_dummy = 0
dt[(dt[:treatment] .== "6.1"), :weight_dummy] = 1

# dummy for gift exchange

dt.gift_dummy = 0
dt[(dt[:treatment] .== "10"), :gift_dummy] = 1

# generating effort and log effort. authors round buttonpressed to nearest 100 value. If 0 set it to 25. julia rounds 50 to 0, while stata to 100. by adding a small value to the column 'buttonpresses' we avoid this mismatch
	
dt.buttonpresses = dt.buttonpresses .+ 0.1 
dt.buttonpresses_nearest_100 = round.(dt.buttonpresses, digits = -2)
dt[(dt[:buttonpresses_nearest_100] .== 0.0), :buttonpresses_nearest_100] = 25
dt.logbuttonpresses_nearest_100 = log.(dt.buttonpresses_nearest_100)	
	
end;

# ╔═╡ 082f8ca0-279f-11ec-07c3-838f8c01603e
# 2. Define the Model and the estimation technique

# For a more detailed explanation of the model we refer the reader to our python notebook or to section 2 in the paper.

# The model is one of costly effort, where an agent needs to choose the optimal effort (in this case the number of buttons pressed in a 10 minute session) to solve a simple tradeoff problem between disutility of effort and consumption utility derived from the consequent payment. On top of this simple problem the authors use 18 different treatments to examine the effects of standard monetary incentives, behavioral factors like social preferences and reference dependence and non-monetary incentives. 

# ╔═╡ 531e546e-27a0-11ec-091f-ab9fb76190b2
begin

# 3. Point estimates and standard errors

# Define the benchmark sample by creating dummies equal to one if in treatment 1.1, 1.2, 1.3 

dt.t11 = dt.treatment .== "1.1"
dt.t12 = dt.treatment .== "1.2"
dt.t13 = dt.treatment .== "1.3"
dt.dummy1 = dt.t11 .+ dt.t12 .+ dt.t13

# Set the initial values for the optimization procedure and scalers for k and s in the exp cost function case

gamma_init_exp, k_init_exp, s_init_exp =  0.015645717, 1.69443, 3.69198
st_values_exp = [gamma_init_exp, k_init_exp, s_init_exp]
k_scaler_exp, s_scaler_exp = 1e+16,1e+6
	
end;

# ╔═╡ 7b307e00-27a7-11ec-2070-9f6d80983e32
begin

# Since there are many different specifications (5 columns for the power cost function and 5 for the exponential cost function) we preferred to write each function separately instead of writing a single function with many if statements. Hopefully this will make each specification more clear.
	
# Estimate procedure for s, k, gamma in benchmark case with Exp cost function
	
# The non-linear least square procedure will choose the parameters that minimize the sum of squared residuals (y-f(x,θ))^2. Here we are computing f(x,θ) while y is the number of rounded buttonpresses. 
	
# pay100 is the column we created containing the piece rate for different treatments
# g, k, s are the parameters to estimate (our θ vector). g stands for gamma.

function benchmarkExp(pay100, params)
	
	g, k, s = params
	
    check1= k / k_scaler_exp  # 'first' component to compute f(x,θ). We call                                           it check1 since it will enter a log, so we                                             need to be careful with its value being > 0
   	check2= s / s_scaler_exp .+ pay100 
	
    f_x = (-1 /g * log.(check1) .+1 /g * log.(check2))
	
    return f_x
end

# Find the solution to the problem by non-linear least squares

sol = curve_fit(benchmarkExp,
				dt[(dt[:dummy1] .== 1.0), :payoff_per_100],
				dt[(dt[:dummy1] .== 1.0), :buttonpresses_nearest_100],
				st_values_exp,
		        autodiff=:forwarddiff) 
be54 = sol.param     # Point estimates
se54 = stderror(sol) # Standard errors. They are the square root of the diagonal                              elements of the estimate for the var-cov matrix
end;

# ╔═╡ d153ea60-590c-11ec-09a4-d7d6bffb1d5a
begin

# Estimate procedure for s, k, gamma in benchmark case with power cost function

gamma_init_power, k_init_power, s_init_power =  19.8117987, 1.66306e-10, 7.74996
st_values_power = [gamma_init_power, k_init_power, s_init_power]
k_scaler_power, s_scaler_power = 1e+57, 1e+6

# Define f(x,θ) in the power case

function benchmarkPower(pay100, params)
	
	g, k, s = params
	
    check1= max(k, 1e-50) / k_scaler_power 
   	check2= max(s, 1e-10) / s_scaler_power .+ pay100 
	
    f_x = (-1 /g .* log.(check1) .+1 /g .* log.(check2))
	
    return f_x
end

# Find the solution to the problem by non-linear least squares

sol52 = curve_fit(benchmarkPower,
				dt[(dt[:dummy1] .== 1.0), :payoff_per_100],
				dt[(dt[:dummy1] .== 1.0), :logbuttonpresses_nearest_100],
				st_values_power,
		        autodiff=:forwarddiff) 
bp52 = sol52.param     
sp52 = stderror(sol52) 
	
end;

# ╔═╡ 05e3eb00-7527-11ec-3380-cdcafa241aaa
# We find slightly different estimates compared to the ones obtained by the authors, but we see also see that different minimization algorithms implemented with the same programming language (julia) result in different estimates and/or values of the objective function, so it is not surprising that there are small discrepancies between our estimates and the authors' estimates (the authors use the Gauss-Newton minimization algorithm implemented in Stata). At the same time, the differences are small in absolute value (and limited to the NLSS estimation method, there are no discrepancies when using GMM) and the estimated values of k and s are always statistically indistinguishable from 0. More importantly, the economic implications of the estimated parameters and the qualitative conclusions on what motivates effort in the experiment are unaffected by the choice of programming language and minimization algorithm. We report the results we obtained with the curve_fit function since this also returns an estimate for the variance-covariance matrix for the parameters

# ╔═╡ 66dd6fb0-2b32-11ec-2a0c-5593db218da8
begin

# Allnoweight Exp. Create dummies for this specification

dt.t31 = dt.treatment .== "3.1"
dt.t32 = dt.treatment .== "3.2"
dt.t41 = dt.treatment .== "4.1"
dt.t42 = dt.treatment .== "4.2"
dt.t10 = dt.treatment .== "10"
dt.samplenw = dt.dummy1 .+ dt.t31 .+ dt.t32 .+ dt.t41 .+ dt.t42 .+ dt.t10

# Define the initial guesses for all parameters

alpha_init, a_init, beta_init, delta_init, gift_init = 0.003, 0.13, 1.16, 0.75, 5e-6
stvale_spec = [alpha_init, a_init, gift_init, beta_init, delta_init]
st_valuesnoweight_exp = vcat(st_values_exp,stvale_spec)

end;

# ╔═╡ 27a92130-2c64-11ec-2c22-c55b0094de30
begin

# Define the f(x,θ) to estimate all parameters but the probability weight in the exp case

# xdata is the vector containing the explanatory variables x:

# gd is gift dummy
# dd is delay dummy
# dw is delay weeks
# paychar is pay in charity treatment
# dc is dummy charity

# parameters:

# g, k, s are the same parameters from before
# alpha is the pure altruism coefficient
# a is the warm glow coefficient
# gift is the gift exchange coefficient Δs
# beta is the present bias paramater
# delta is the (weekly) discount factor
	
function noweightExp(pay100, params)
	
	# We know that passing variables from outside the function instead of using them         as function arguments is bad practice and it slows down the code, but curve_fit       will return an error if the first argument for the independent variable is a           vector of vectors. The problem lies in the fact that the package checks for           quality of the data using isinf() on the independent variables, but fails for a       vector of vectors since the code does not check element-wise. Writing isinf.()         in the package code and then check for any true should solve the problem.
		
    gd = dt[(dt[:samplenw] .== 1.0), :gift_dummy]
    dd = dt[(dt[:samplenw] .== 1.0), :delay_dummy]
    dw = dt[(dt[:samplenw] .== 1.0), :delay_wks]
    paychar = dt[(dt[:samplenw] .== 1.0), :payoff_charity_per_100]
    dc = dt[(dt[:samplenw] .== 1.0), :dummy_charity]
	g, k, s, alpha, a, gift, beta, delta = params

    check1 = k/k_scaler_exp  
   	check2 = s/s_scaler_exp .+ gift.*0.4.*gd .+ (beta.^dd).*(delta.^dw).*pay100 .+ alpha.*paychar .+ a.*0.01.*dc 
	
    f_x = (-1/g .* log.(check1) .+ 1/g .* log.(check2))
	
    return f_x
end

# Find the solution to the problem by non-linear least squares

sol56 = curve_fit(noweightExp,
				dt[(dt[:samplenw] .== 1.0), :payoff_per_100],
				dt[(dt[:samplenw] .== 1.0), :buttonpresses_nearest_100],
				st_valuesnoweight_exp,
		        autodiff=:forwarddiff) 
be56 = sol56.param     
se56 = stderror(sol56) 
	
end;

# ╔═╡ 592baf2e-599f-11ec-3f47-1782a76c828d
begin

# Define the f(x,θ) to estimate all parameters but the probability weight in the power case
	
function noweightPower(pay100, params)
	
    gd = dt[(dt[:samplenw] .== 1.0), :gift_dummy]
    dd = dt[(dt[:samplenw] .== 1.0), :delay_dummy]
    dw = dt[(dt[:samplenw] .== 1.0), :delay_wks]
    paychar = dt[(dt[:samplenw] .== 1.0), :payoff_charity_per_100]
    dc = dt[(dt[:samplenw] .== 1.0), :dummy_charity]
	g, k, s, alpha, a, gift, beta, delta = params

    check1 = max(k,1e-50)/k_scaler_power  
   	check2 = max(s,1e-10)/s_scaler_power .+ gift.*0.4.*gd .+ (beta.^dd).*(delta.^dw).*pay100 .+ alpha.*paychar .+ a.*0.01.*dc 
	
    f_x = (-1/g .* log.(check1) .+ 1/g .* log.(check2))
	
    return f_x
end

# Find the solution to the problem by non-linear least squares

st_valuesnoweight_power = vcat(st_values_power, stvale_spec)
sol53 = curve_fit(noweightPower,
				  dt[(dt[:samplenw] .== 1.0), :payoff_per_100],
				  dt[(dt[:samplenw] .== 1.0), :logbuttonpresses_nearest_100],
				  st_valuesnoweight_power,
		          autodiff=:forwarddiff) 
bp53 = sol53.param     
sp53 = stderror(sol53) 
	
end;

# ╔═╡ eff34620-59a0-11ec-37b4-f372da908fac
begin

# Create the arrays that will be the columns for table 5. Create the dataframe for table 5 and save it as a csv file
	
# Point estimates for power case do not coincide precisely as explained above. Standard errors do not coincide precisely because of the differences in the point estimates and because we leave here non-robust standard errors provided by curve_fit. To see an implementation of the formula for robust standard errors please refer to the python or julia notebooks for table_1 of augenblick-rabin or table_1 of bruhin-fehr-schunk. The formula is the same as in the cited notebooks without considering the clustering at the individual level.

params_name = ["Curvature γ of cost function", "Level k of cost of effort", "Intrinsic motivation s","Social preferences α","Warm glow coefficient a","Gift exchange Δs", "Present bias β","(Weekly) discount factor δ"]

be5 = [be54[1], be54[2]/1e+16, be54[3]/1e+6, be56[4], be56[5], be56[6], be56[7], be56[8]]
se5 = [se54[1], se54[2]/1e+16, se54[3]/1e+6, se56[4], se56[5], se56[6], se56[7], se56[8]]
bp5 = [bp52[1], bp52[2]/1e+57, bp52[3]/1e+6, bp53[4], bp53[5], bp53[6], bp53[7], bp53[8]]
sp5 = [sp52[1], sp52[2]/1e+57, sp52[3]/1e+6, sp53[4], sp53[5], sp53[6], sp53[7], sp53[8]]
t5 = DataFrame()
t5.parameters_name = params_name
t5.power_est=bp5
t5.power_se=sp5
t5.exp_est=be5
t5.exp_se=se5
CSV.write(raw"../output/table5NLS_julia.csv", t5)
t5	
	
end

# ╔═╡ 8e8624e0-59a3-11ec-1934-6bf58626ec8a
# We now look at panel A of Table 6

# ╔═╡ 8ca6c9f0-2c66-11ec-1cb8-bf7d2915a69d
begin

# Create the sample used for Table 6 panel A

dt.t61 = dt.treatment .== "6.1"
dt.t62 = dt.treatment .== "6.2"
dt.samplepr = dt.dummy1 .+ dt.t61 .+ dt.t62

end;

# ╔═╡ 2f9046b0-2c6b-11ec-3d18-51475bf3a66c
begin

# Define f(x,θ) for the exponential cost function. Here we assume curvature of utility over piece rate = 1, (Column 4)
# wd is the weight_dummy
# prob is the prob_dummy
# g, k and s are the same parameters as before
# p_weight is the probability weighting coefficient under the assumption of linear value function in this case 
# curv is the curvature of the value function. Here curv = 1

function probweight4Exp(pay100, params)
	
	wd = dt[(dt[:samplepr] .== 1.0), :weight_dummy]
    prob = dt[(dt[:samplepr] .== 1.0), :prob]
	g, k, s, p_weight = params
	
    check1= k / k_scaler_exp  
   	check2= s / s_scaler_exp .+ p_weight.^wd.*prob.*pay100
	
    f_x = (-1 /g * log.(check1) .+1 /g * log.(check2))
	
    return f_x
end

# Find the solution to the problem by non-linear least squares

prob_weight_init = [0.2]
st_valuesprobweight_exp = vcat(st_values_exp, prob_weight_init)

sol64 = curve_fit(probweight4Exp,
				dt[(dt[:samplepr] .== 1.0), :payoff_per_100],
				dt[(dt[:samplepr] .== 1.0), :buttonpresses_nearest_100],
				st_valuesprobweight_exp,
		        autodiff=:forwarddiff) 
be64 = sol64.param     
se64 = stderror(sol64) 

# we rescale the variables	
	
be64 = [be64[1];be64[2]/1e+16;be64[3]/1e+6;be64[4]]
se64 = [se64[1];se64[2]/1e+16;se64[3]/1e+6;se64[4]]
	
end;

# ╔═╡ 3dd7ee00-59a5-11ec-021d-7fb70238f5a9
begin

# Define f(x,θ). Here we assume curvature of utility over piece rate = 0.88, Column (5)

function probweight5Exp(pay100, params)
	
	wd = dt[(dt[:samplepr] .== 1.0), :weight_dummy]
    prob = dt[(dt[:samplepr] .== 1.0), :prob]
	g, k, s, p_weight = params
	
    check1= k / k_scaler_exp  
   	check2= s / s_scaler_exp .+ p_weight.^wd.*prob.*pay100.^0.88
	
    f_x = (-1 /g * log.(check1) .+1 /g * log.(check2))
	
    return f_x
end

# Find the solution to the problem by non-linear least squares

sol65 = curve_fit(probweight5Exp,
				dt[(dt[:samplepr] .== 1.0), :payoff_per_100],
				dt[(dt[:samplepr] .== 1.0), :buttonpresses_nearest_100],
				st_valuesprobweight_exp,
		        autodiff=:forwarddiff) 
be65 = sol65.param     
se65 = stderror(sol65) 

be65 = [be65[1];be65[2]/1e+16;be65[3]/1e+6;be65[4]]
se65 = [se65[1];se65[2]/1e+16;se65[3]/1e+6;se65[4]]

# Define f(x,θ). Here we assume curvature of utility over piece rate = 0.88, Column (5)

function probweight6Exp(pay100, params)
	
	wd = dt[(dt[:samplepr] .== 1.0), :weight_dummy]
    prob = dt[(dt[:samplepr] .== 1.0), :prob]
	g, k, s, p_weight, curv = params
	
    check1= k / k_scaler_exp  
   	check2= s / s_scaler_exp .+ p_weight.^wd.*prob.*pay100.^curv
	
    f_x = (-1 /g * log.(check1) .+1 /g * log.(check2))
	
    return f_x
end

# Find the solution to the problem by non-linear least squares

curv_init = [0.5]
st_valuesprobweight6_exp = vcat(st_valuesprobweight_exp, curv_init)
	
sol66 = curve_fit(probweight6Exp,
				dt[(dt[:samplepr] .== 1.0), :payoff_per_100],
				dt[(dt[:samplepr] .== 1.0), :buttonpresses_nearest_100],
				st_valuesprobweight6_exp,
		        autodiff=:forwarddiff) 
be66 = sol66.param     
se66 = stderror(sol66)
	
be66 = [be66[1];be66[2]/1e+16;be66[3]/1e+6;be66[4:5]]
se66 = [se66[1];se66[2]/1e+16;se66[3]/1e+6;se66[4:5]]
	
end;

# ╔═╡ b5b94c80-2c6d-11ec-0ab0-71387355858d
begin

# We do the same for the power cost function specification

# For column 4

function probweight4Power(pay100, params)
	
	wd = dt[(dt[:samplepr] .== 1.0), :weight_dummy]
    prob = dt[(dt[:samplepr] .== 1.0), :prob]
	g, k, s, p_weight = params
	
    check1= max(k,1e-50) / k_scaler_power  
   	check2= max(s,1e-10) / s_scaler_power .+ p_weight.^wd.*prob.*pay100
	
    f_x = (-1 /g * log.(check1) .+1 /g * log.(check2))
	
    return f_x
end

# Find the solution to the problem by non-linear least squares

st_valuesprobweight_power = vcat(st_values_power, prob_weight_init)

sol61 = curve_fit(probweight4Power,
				dt[(dt[:samplepr] .== 1.0), :payoff_per_100],
				dt[(dt[:samplepr] .== 1.0), :logbuttonpresses_nearest_100],
				st_valuesprobweight_power,
		        autodiff=:forwarddiff) 
bp61 = sol61.param     
sp61 = stderror(sol61) 
	
bp61 = [bp61[1];bp61[2]/1e+57;bp61[3]/1e+6;bp61[4]]
sp61 = [sp61[1];sp61[2]/1e+57;sp61[3]/1e+6;sp61[4]]

# For column 5

function probweight5Power(pay100, params)
	
	wd = dt[(dt[:samplepr] .== 1.0), :weight_dummy]
    prob = dt[(dt[:samplepr] .== 1.0), :prob]
	g, k, s, p_weight = params
	
    check1= max(k,1e-50) / k_scaler_power  
   	check2= max(s,1e-10) / s_scaler_power .+ p_weight.^wd.*prob.*pay100.^0.88
	
    f_x = (-1 /g * log.(check1) .+1 /g * log.(check2))
	
    return f_x
end

# Find the solution to the problem by non-linear least squares

sol62 = curve_fit(probweight5Power,
				dt[(dt[:samplepr] .== 1.0), :payoff_per_100],
				dt[(dt[:samplepr] .== 1.0), :logbuttonpresses_nearest_100],
				st_valuesprobweight_power,
		        autodiff=:forwarddiff) 
bp62 = sol62.param     
sp62 = stderror(sol62) 
	
bp62 = [bp62[1];bp62[2]/1e+57;bp62[3]/1e+6;bp62[4]]
sp62 = [sp62[1];sp62[2]/1e+57;sp62[3]/1e+6;sp62[4]]

# For column 6

function probweight6Power(pay100, params)
	
	wd = dt[(dt[:samplepr] .== 1.0), :weight_dummy]
    prob = dt[(dt[:samplepr] .== 1.0), :prob]
	g, k, s, p_weight, curv = params
	
    check1= max(k,1e-50) / k_scaler_power  
   	check2= max(s,1e-10) / s_scaler_power .+ p_weight.^wd.*prob.*pay100.^curv
	
    f_x = (-1 /g * log.(check1) .+1 /g * log.(check2))
	
    return f_x
end

# Find the solution to the problem by non-linear least squares

st_valuesprobweight6_power = vcat(st_valuesprobweight_power, curv_init)
	
sol63 = curve_fit(probweight6Power,
				dt[(dt[:samplepr] .== 1.0), :payoff_per_100],
				dt[(dt[:samplepr] .== 1.0), :logbuttonpresses_nearest_100],
				st_valuesprobweight6_power,
		        autodiff=:forwarddiff) 
bp63 = sol63.param     
sp63 = stderror(sol63)
	
bp63 = [bp63[1];bp63[2]/1e+57;bp63[3]/1e+6;bp63[4:5]]
sp63 = [sp63[1];sp63[2]/1e+57;sp63[3]/1e+6;sp63[4:5]]
	
# columns must be of the same length
abp61 = append!(bp61,1)
asp61 = append!(sp61,0)
abp62 = append!(bp62,0.88)
asp62 = append!(sp62,0)
abe64 = append!(be64,1)
ase64 = append!(se64,0)
abe65 = append!(be65,0.88)
ase65 = append!(se65,0)

end;

# ╔═╡ d7f36f90-2c6f-11ec-155a-db4f5d543604
begin

# Create the dataframe relative to table 6 and save it as a csv file
	
# There are some differences in the standard errors since we leave here non robust standard errors provided by curve_fit. Point estimates for the power cost function are again a little different from the authors', while they are the same for the exponential cost function

t6 = DataFrame()

pnames = ["Curvature γ of cost function", "Level k of cost of effort", "Intrinsic motivation s", "Probability weighting π (1%) (in %)", "Curvature of utility over piece rate"]

t6.parameters=pnames
t6.p_est1 = abp61 
t6.p_se1  = asp61
t6.p_est2 = abp62
t6.p_se2  = asp62
t6.p_est3 = bp63
t6.p_se3  = sp63
t6.e_est4 = abe64
t6.e_se4  = ase64
t6.e_est5 = abe65
t6.e_se5  = ase65
t6.e_est6 = be66
t6.e_se6  = se66
	
CSV.write(raw"../output/table6_julia.csv", t6)
t6
	
end

# ╔═╡ Cell order:
# ╠═aa3922f0-279e-11ec-33b6-6bb73cf7c16a
# ╠═08a0ec12-279f-11ec-12a3-f15a088873f4
# ╠═6783e0c0-5904-11ec-2a41-5784eec17590
# ╠═903cd5d0-5904-11ec-0a95-75a9f8821a44
# ╠═086c95a0-279f-11ec-3378-6924b1d7777e
# ╠═086c6e90-279f-11ec-2ce7-ef05700c68b3
# ╠═082f8ca0-279f-11ec-07c3-838f8c01603e
# ╠═531e546e-27a0-11ec-091f-ab9fb76190b2
# ╠═7b307e00-27a7-11ec-2070-9f6d80983e32
# ╠═d153ea60-590c-11ec-09a4-d7d6bffb1d5a
# ╠═05e3eb00-7527-11ec-3380-cdcafa241aaa
# ╠═66dd6fb0-2b32-11ec-2a0c-5593db218da8
# ╠═27a92130-2c64-11ec-2c22-c55b0094de30
# ╠═592baf2e-599f-11ec-3f47-1782a76c828d
# ╠═eff34620-59a0-11ec-37b4-f372da908fac
# ╠═8e8624e0-59a3-11ec-1934-6bf58626ec8a
# ╠═8ca6c9f0-2c66-11ec-1cb8-bf7d2915a69d
# ╠═2f9046b0-2c6b-11ec-3d18-51475bf3a66c
# ╠═3dd7ee00-59a5-11ec-021d-7fb70238f5a9
# ╠═b5b94c80-2c6d-11ec-0ab0-71387355858d
# ╠═d7f36f90-2c6f-11ec-155a-db4f5d543604
