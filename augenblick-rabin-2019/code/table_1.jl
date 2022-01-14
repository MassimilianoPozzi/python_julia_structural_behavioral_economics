### A Pluto.jl notebook ###
# v0.14.0

using Markdown
using InteractiveUtils

# ╔═╡ 3bf57e20-267b-11ec-35eb-a983675c3502
begin 

# Import packages used

using Distributions
using DataFrames
using Optim
using CSV
using DelimitedFiles
using StatFiles
using Random
using Statistics
using ForwardDiff	
using LinearAlgebra
using FiniteDiff
	
end

# ╔═╡ 0c4e6562-267b-11ec-2ca3-bbbd9a335a0d
# Augenblick and Rabin, 2019, "An Experiment on Time Preference and Misprediction in Unpleasant Tasks", Table 1

# ╔═╡ 3c1a1d20-267b-11ec-03dd-2b946e7960e3
# Authors: 

# Massimiliano Pozzi (Bocconi University, pozzi.massimiliano@studbocconi.it)
# Salvatore Nunnari (Bocconi University, https://snunnari.github.io, salvatore.nunnari@unibocconi.it)

# The code in this Jupyter notebook performs the aggregate estimates to replicate column 1 of Table 1

# This notebook was tested with the following packages versions:

# Pozzi: julia 1.6.1, Pluto 0.14.7, Neptune 0.14.0, Distributions 0.25.17, DataFrames 1.2.2, Optim 1.4.1, CSV 0.9.5, ForwardDiff 0.10.20

# ╔═╡ c0908a20-4b8e-11ec-374c-9f6e5b07e54d
# This notebook works with both Neptune and Pluto; the "begin" and "end" tags in each cell are necessary for Pluto notebooks but do not play any role in Neptune 

# ╔═╡ 01bac3e0-267c-11ec-1a92-c5ff23e97f0e
# To import a new package (es. StatFiles) use this cell. Change the name inside the brackets to change Pkg

# import Pkg; Pkg.add("StatFiles")

# ╔═╡ 0f4aa740-4b8f-11ec-1749-df8a9c57903d
begin

# 1. Data Cleaning and Data Preparation

# We import the dataset containing the choices of all 100 individuals who participated to the experiment. To guarantee consistency with the authors' results, we then construct the primary sample used for the aggregate estimates. This sample consists of 72 individuals whose individual parameter estimates converged in less than 200 iterations when using the authors' Stata algorithm. In particular, we run the "03MergeIndMLEAndConstructMainSample.do" file provided by the authors. This script creates a file named "ind_to_keep.csv" which contains the identifiers of the individuals to keep.

dt = DataFrame(load(raw"../input/decisions_data.dta"))             # full sample
data_keep = readdlm(raw"../input/ind_to_keep.csv",',',header=true) # import csv with ID of subjects to keep
	
	
# drop subjects whose IDs are not listed in the data_drop dataframe (28 individuals)
# this is the primary sample for the aggregate estimates (72 individuals)
	
	
dt = dt[[i for i=1:length(dt.wid) if dt.wid[i] in(data_keep[1][:,1])],:]
	
# We remove observations when a bonus was offered and create the following dummy variables that will be useful for estimation: pb is equal to one if the subject completed 10 mandatory tasks on subject-day (this is used to estimate the projection bias parameter &alpha;); ind_effort10 and ind_effort110 are equal to one if, respectively, the subject completed 10 or 110 tasks (and they are used for the Tobit correction when computing the likelihood).

dt = dt[dt.bonusoffered .!=1,:]  # remove observations when a bonus was offered
dt.pb = dt.workdone1 ./ 10       # pb dummy variable. workdone1 can either be 10 or 0, so dividing the variable by 10 creates our dummy

dt.ind_effort10  = dt.effort .== 10    # ind_effort10 dummy
dt.ind_effort110 = dt.effort .== 110   # ind_effort110 dummy	
	
end;

# ╔═╡ 2079e14e-47f9-11ec-1620-4b0af35d1fe0
# 2. Define the Model and the Likelihood (Section 3 in the Paper)

# For a more detailed explanation of the model please refer to the python notebook

# Define the criterion function to minimize. This function computes the negative log likelihood of observing the choices made by the agents. We first compute the predicted choice that is found by solving the agents' utility maximization problem. We then assume that the observed effort is distributed as a Gaussian with mean the predicted choice and a given value of sigma (to be estimated). We then work out the likelihood of observing the effort found in the data, take the log and the negative of it to define our objective function.

# parameters:

# beta is the present bias parameter
# betahat is the perceived present bias parameter
# delta is the usual time-discounting parameter
# gamma and phi are the two parameters controlling the cost of effort function
# alpha is the projection bias parameter
# sigma is the standard deviation of the normal error term ϵ

# args:

# netdistance is the difference between the payment date T and the work time t
# wage is the amount paid per task in a certain session
# today is a dummy variable equal to one if the decision involves the choice of work today
# prediction is a dummy variable equal to one if the decision involves the choice of work in the future
# pb is a dummy equal to one if the subject completed 10 mandatory tasks on subject-day
# effort is the number of tasks completed by a subject in a session. It can range from a minimum of 10 to a maximum of 110
# ind_effort10 is a dummy equal to one if the subject's effort was equal to 10
# ind_effort110 is a dummy equal to one if the subject's effort was equal to 110

function negloglike(params, args)
	
    beta, betahat, delta, gamma, phi, alpha, sigma = params
    netdistance,wage,today,prediction,pb,effort,ind_effort10,ind_effort110 = args
    
	# compute the predicted choice
    predchoice=((phi.*(delta.^netdistance).*(beta.^today).*(betahat.^prediction).*wage).^(1 ./(gamma.-1))).-pb.*alpha
    
	# compute the vector of individual observations' likelihood, a 1x8049 vector
	
    prob = pdf.(Normal(0,sigma),predchoice.-effort) 
    
	# find index when effort is 10 ot 110 to apply Tobit correction
	len = length(effort)
    index_10  = [i for i=1:len if effort[i]==10]
    index_110 = [i for i=1:len if effort[i]==110]   
    
	# Apply Tobit correction
    for i in index_10
        prob[i] = 1 - cdf.(Normal(0,1),(predchoice[i]-effort[i])/sigma)
	end
        
    for i in index_110
        prob[i] = cdf.(Normal(0,1),(predchoice[i]-effort[i])/sigma)
	end
       
	# find index when prob is equal to 0 or 1 and change it. this is needed to avoid problems when taking logs
    index_p0 = [i for i=1:len if prob[i]==0]
    index_p1 = [i for i=1:len if prob[i]==1]
    
	# change prob when 0 or 1
    for i in index_p0
        prob[i] = 1E-4
	end
        
    for i in index_p1
        prob[i] = 1 - 1E-4
	end
    
    negll = - sum(log.(prob)) # negative log likelihood
    
    return negll
	
end

# ╔═╡ d3378530-47fa-11ec-31b8-53005e5a1992
begin

# Define the initial guesses. They are the same as the ones used by the authors in their do.file

beta_init, betahat_init, delta_init, gamma_init, phi_init, alpha_init, sigma_init = 0.8, 1.0, 1.0, 2.0, 500.0, 7.0, 40.0

params_initPb = [beta_init, betahat_init, delta_init, gamma_init, phi_init, alpha_init, sigma_init]

# define the arguments for the criterion function

mle_args = [dt.netdistance,dt.wage,dt.today,dt.prediction,dt.pb,dt.effort,dt.ind_effort10,dt.ind_effort110]

# The function to minimize in Optim can take only one argument, so we create this function g that includes already the arguments needed to compute the negative log likelihood

function g(params)
	nll = negloglike(params, mle_args)
	return nll
end
	
end

# ╔═╡ 70860f00-47fb-11ec-35df-c5ccd76d49c8
begin
	
# Compute the estimates and the value of the log likelihood
	
sol = optimize(g,params_initPb, method = NelderMead())
restable_1 = Optim.minimizer(sol) # point estimates
logvalue = -Optim.minimum(sol)    # value of the log likelihood
		
end;

# ╔═╡ bce2cac0-4b96-11ec-35bc-033c65df8502
begin

# Compute the individual cluster robust standard errors. 
# For a more detailed explanation on how to compute cluster robust standard errors we refer the readero to the python notebook or to David A. Freedman, 2006, "On The So-Called 'Huber Sandwich Estimator' and 'Robust Standard Errors'", (The American Statistician 60:4, 299-302)
	
# The formula is the following:

# adj * (inv_hessian * G * inv_hessian). 

# adj is a correction for the degree of freedoms and the number of clusters. inv_hessian is the inverse of the hessian of the negative log-likelihood evaluated at the minimum. G is is a 5x5 matrix of gradient contributions. It is computed as the sum over all clusters of the outer product of the sum of the single observation gradient evaluated in the minimum for all observations of an individual cluster.

# Define an auxiliary function.
# This function takes as input a dataset and a scalar and returns a single observation dataframe

function getobs(dataset, obs)
    dataset_ind = dataset[obs, :]
    return dataset_ind
end

# Define the function that computes the matrix of individual gradient contribution 

function gradcontr(dataset, parameters)
    G = zeros(length(parameters), length(parameters))
	vsingle_grad = []
		
    for obs=1:length(dataset.wid)        # loop over all observations
        dt_ind = getobs(dataset, obs)    # dataframe for observation 'obs'
			
		# define the arguments that are needed to compute the negative log likelihood
	    mle_arg_ind=[[dt_ind.netdistance],[dt_ind.wage],[dt_ind.today],[dt_ind.prediction], [dt_ind.pb],[dt_ind.effort],[dt_ind.ind_effort10],[dt_ind.ind_effort110]]

		# gradient of the negloglike function evaluated at the minimum

         single_grad = ForwardDiff.gradient(x->negloglike(x, mle_arg_ind), parameters) 
		 
		append!(vsingle_grad, [single_grad])
			
	end
		
	dg = DataFrame()
	dg.wid = dataset.wid
	dg.gradient = vsingle_grad
	
	for wid in unique(dataset.wid)                        # loop over individuals' IDs
		ind_grad = sum(dg[(dg[:wid] .== wid), :gradient]) # sum gradients element wise
		G .+= (ind_grad*ind_grad')
	end
			
	return G
end
	
end

# ╔═╡ 4a11de80-4b98-11ec-046a-5bfdf5943448
begin

# Compute the individual cluster robust standard errors

# Compute the hessian

hessian = ForwardDiff.hessian(g,restable_1)  # numerical hessian
hess_inv = inv(hessian)                      # inverse of the hessian

# Compute the matrix of gradient contribution

grad_contribution = gradcontr(dt, restable_1)

# Compute the adjustment for degree of freedoms and number of clusters

adj = (length(dt.wid)-1)/(length(dt.wid)-length(restable_1)) * length(unique(dt.wid))/(length(unique(dt.wid))-1)

# var-cov matrix of our estimates

varcov_estimates = adj .*(hess_inv * grad_contribution * hess_inv) 

# individual cluster robust standard errors

se_cluster = sqrt.(diag((varcov_estimates)))
	
end;

# ╔═╡ 5691e8e0-4ba6-11ec-07c5-b90ec99a3b44
begin

# Hypothesis testing

# We now do some hypothesis testing on the parameters we obtained. We compute the z-test statistics and the corresponding p-values to check if beta, betahat or delta are statistically different from one. We then compute the p-value of a z-test to check if the parameter for projection bias is statistically different from zero.
	
# Compute the z-test statistics and the corresponding p-values to check if beta, betahat, delta are statistically different from one

zvalues_1 = (restable_1[1:3].-1)./se_cluster[1:3] # the first three elements are for beta (position 1), betahat (position 2) and delta (position 3)
pvalues_1 = 2*(1 .-cdf.(Normal(0,1),abs.(zvalues_1)))

# Now compute the z-test statistics and the corresponding p-value for H0: alpha different from 0

zvalue_a = restable_1[6]/se_cluster[6]
pvalue_a = 2*(1 .-cdf.(Normal(0,1),zvalue_a))
	
end;

# ╔═╡ 2882c3c0-268f-11ec-2c70-6d123ee6cc45
begin

# Create a Dataframe with the results and save it as a csv file

parameters_name = ["Present Bias β",
                   "Naive Pres. Bias β_h",
                   "Discount Factor δ",
                   "Cost Curvature γ",
                   "Cost Slope ϕ",
                   "Proj Task Reduction α",
                   "Sd of error term σ"]

table_1 = DataFrame()
table_1.parameters = parameters_name
table_1.estimates = round.(restable_1,digits=3)
table_1.standarderr = round.(se_cluster,digits=3)

CSV.write(raw"../output/table1_julia.csv", table_1)
table_1
	
end

# ╔═╡ 2c3d28b0-4ba7-11ec-17ea-8b4d0d602c2c
begin

# in this cell we show the nr. of observations, the nr. of participants, the value of the log likelihood and the p-values of the tests
	
t1 = DataFrame()
t1.description = ["Observations","Participants","Log Likelihood","H0(β=1)","H0(β_h=1)","H0(α=0)","H0(δ=1)"]
t1.values = [length(dt.wid),length(unique(dt.wid)),-Optim.minimum(sol),round(pvalues_1[1],digits=3),round(pvalues_1[2],digits=2),round(pvalue_a,digits=3),round(pvalues_1[3],digits=2)]
t1
	
end

# ╔═╡ Cell order:
# ╠═0c4e6562-267b-11ec-2ca3-bbbd9a335a0d
# ╠═3c1a1d20-267b-11ec-03dd-2b946e7960e3
# ╠═c0908a20-4b8e-11ec-374c-9f6e5b07e54d
# ╠═01bac3e0-267c-11ec-1a92-c5ff23e97f0e
# ╠═3bf57e20-267b-11ec-35eb-a983675c3502
# ╠═0f4aa740-4b8f-11ec-1749-df8a9c57903d
# ╠═2079e14e-47f9-11ec-1620-4b0af35d1fe0
# ╠═d3378530-47fa-11ec-31b8-53005e5a1992
# ╠═70860f00-47fb-11ec-35df-c5ccd76d49c8
# ╠═bce2cac0-4b96-11ec-35bc-033c65df8502
# ╠═4a11de80-4b98-11ec-046a-5bfdf5943448
# ╠═5691e8e0-4ba6-11ec-07c5-b90ec99a3b44
# ╠═2882c3c0-268f-11ec-2c70-6d123ee6cc45
# ╠═2c3d28b0-4ba7-11ec-17ea-8b4d0d602c2c
