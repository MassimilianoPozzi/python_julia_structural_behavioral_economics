### A Pluto.jl notebook ###
# v0.14.0

using Markdown
using InteractiveUtils

# ╔═╡ a0736b60-3c7d-11ec-1b71-27fd38cbf6e4
begin
	
# Import the necessary packages

using Distributions
using DelimitedFiles
using DataFrames
using Optim
using CSV
using Random
using Distributions
using Statistics
using LinearAlgebra
using ForwardDiff
	
end

# ╔═╡ 7b7526d0-3c9d-11ec-3575-37a66c520386
# Replication of Bruhin, Fehr, and Schunk, 2019, "Many Faces of Human Sociality: Uncovering the Distribution and Stability of Social Preferences"

# ╔═╡ 9d2d5b40-6023-11ec-0e3b-e5fe5c8d527b
# The code in this notebook performs the finite mixture estimates as in table 2

# ╔═╡ 7b36490e-3c9d-11ec-1ea6-dbf0082362dc
# Authors: 
# Massimiliano Pozzi (Bocconi University, pozzi.massimiliano@studbocconi.it)
# Salvatore Nunnari (Bocconi University, https://snunnari.github.io, salvatore.nunnari@unibocconi.it)

# ╔═╡ 600cd6e0-4ea4-11ec-2071-13c721a5b90e
# This notebook works with both Neptune and Pluto; the "begin" and "end" tags in each cell are necessary for Pluto notebooks but do not play any role in Neptune 

# ╔═╡ 61e5b590-4ea4-11ec-1bf8-a900f11ecb18
# This notebook was tested with the following packages versions:

# Pozzi: julia 1.6.1, Pluto 0.14.7, Neptune 0.14.0, Distributions 0.25.17, DataFrames 1.2.2, Optim 1.4.1, CSV 0.9.5, FiniteDiff 2.8.1, ForwardDiff 0.10.20

# ╔═╡ e0965590-3c7d-11ec-13b6-a7d15db3fd9f
begin

# 1. Data Cleaning and Preparation

# We import the relevant datasets containing the data on the 39 dictator games and 78 reciprocity games in Session 1 and Session 2. As the authors, we remove from these datasets those subjects that behaved very inconsistenly throughout the games. These subjects are identified from the individual estimates that we do not run in this notebook. The file "dropped_subjects_section4paragraph2.csv" contains the ID of these subjects.
	
# Load the data for Session 1

data, header = readdlm(raw"../input/choices_exp1.csv",',',header=true)
dt1 = DataFrame(data, vec(header))
	
# Import data with ID of subjects to drop

data_drop = readdlm(raw"../input/dropped_subjects_section4paragraph2.csv",',',header=false)

# Drop the session 1 subjects whose IDs are listed in the data_drop dataframe (14 individuals)
	
dt1 = dt1[[i for i=1:length(dt1.sid) if !(dt1.sid[i] in(data_drop))],:]
	
sort!(dt1,:sid) # individuals were not in order. We sort the dataframe so that the first 117 obs are for individual 1 and so on.
	
# Import session 2 data, remove subjects to drop and sort the dataframe as we did for session 1 

data2, header2 = readdlm(raw"../input/choices_exp2.csv",',',header=true)
dt2 = DataFrame(data2, vec(header2))
dt2 = dt2[[i for i=1:length(dt2.sid) if !(dt2.sid[i] in(data_drop))],:]
sort!(dt2,:sid)	

# import the aggregate estimates for session 1 and session 2 saved in the table1_julia.csv in the output folder. These will be used as a starting point for the EM-algorithm
	
res_s1=readdlm(raw"../output/table1_julia.csv",',',header=true)[1][:,2]
res_s1=convert(Array{Float64,1}, [res_s1[1:4];log(res_s1[5])])	
	
res_s2=readdlm(raw"../output/table1_julia.csv",',',header=true)[1][:,6]
res_s2=convert(Array{Float64,1}, [res_s2[1:4];log(res_s2[5])])	
	
end

# ╔═╡ d8928282-3c85-11ec-00ca-c5215c731b33
begin

# Create indicators_x and indicators_y. These are the columns s, r, q, v

indicators_x = hcat(dt1.s_x,dt1.r_x,dt1.q,dt1.v)
indicators_y = hcat(dt1.s_y,dt1.r_y,dt1.q,dt1.v)

indicators_x2 = hcat(dt2.s_x,dt2.r_x,dt2.q,dt2.v)
indicators_y2 = hcat(dt2.s_y,dt2.r_y,dt2.q,dt2.v)
	
end;

# ╔═╡ 47094410-4ea6-11ec-2865-47e35e532cf1
# For a more detailed explanation on the model we refer the reader to the python notebook.

# Below, we explain each step of the algoritm, but for a more rigorous discussion please refer to Dempster, A.P, N. M. Laird and D. B. Rubin, 1977, "Maximum Likelihood from Incomplete Data via the EM Algorithm", Journal of the Royal Statistical Society, Series B (Methodological), 39(1): 1&ndash;38 and to McLachlan, Geoffrey, Sharon Lee, and Suren Rathnayake, 2019, "Finite Mixture Models", Annual Review of Statistics and Its Application, 6: 355&ndash;378.

# Define the component density function. This would be f(θ_k, σ_k; X, Y, C_i) above for a specific set of parameters (θ_k, σ_k).
# v is the vector of parameters (θ_k, σ_k)
# y is the choice of the player
# self_x is the payoff when choosing x (left)
# self_y is the payoff when choosing y (right)
# indicators are the ones explained above
# This function returns a 1x160 vector with the likelihood contribution for each individual (nr. of unique individual is 160)

function singlecomp(v,y,self_x, other_x, self_y,other_y,indicators_x,indicators_y)
	
	beta  = v[1:(length(v)-1)] # vector of parameters alpha,beta,gamma,delta
    sigma = exp(v[length(v)])  # sigma. exp to force it positive during search
	
    lli = (indicators_x * beta)  # We get the matrix with rows [αs,βr,γq,δv]
    rli = (indicators_y * beta)
	
	uleft  = sigma.*( (1 .-lli) .*self_x .+lli .*other_x)
	uright = sigma.*( (1 .-rli) .*self_y .+rli .*other_y)
	
	probs  = (exp.(uleft)./(exp.(uleft).+exp.(uright))).^y .* (exp.(uright)./(exp.(uleft).+exp.(uright))).^(1 .-y)
	
	indprob = reshape(probs, 117, 160) # 117x160 matrix. each column is the likelihood of the observations for a single individual
	
	indlike = prod(indprob, dims=1) # 1x160 vector containing the likelihood contribution for a single individual. We took the product of a column
	
	# in theory it should not happen that a probability is equal to zero, but for some combinations of parameters it does. If prob = 0 we set it equal to a small value 
	
	index_p0 = [i for i=1:length(indlike) if indlike[i]==0.0]
	
	for i in index_p0
        indlike[i] = 1E-8
	end
	
	return indlike
	
end	

# ╔═╡ ca59c720-47a5-11ec-3cf6-f1b17e96c3e9
# This function returns a 160x3 matrix. Column one contains the individual likelihood contribution using the parameters for type one. Column two contains the individual likelihood contribution for type two etc. Using the notation from above a generic element i in the first column is what we called f(θ_1, σ_1; X, Y, C_i). Remember that we assumed there are three types (K=3).
# vall is a 1x15 vector that contains the parameters (θ,σ) for component 1,2 and 3, what we called (θ_1,σ_1,θ_2,σ_2,θ_3,σ_3).


function compdens(vall,y,self_x, other_x,self_y,other_y,indicators_x,indicators_y)
	
	theta1 = vall[1:5]    # parameters for first component  α1, β1, γ1, δ1, σ1
	theta2 = vall[6:10]   # parameters for second component α2, β2, γ2, δ2, σ2
	theta3 = vall[11:15]  # parameters for third component  α3, β3, γ3, δ3, σ3
	
	comp1 = singlecomp(theta1,y,self_x, other_x, self_y,other_y,indicators_x,indicators_y)
	
	comp2 = singlecomp(theta2,y,self_x, other_x, self_y,other_y,indicators_x,indicators_y)
	
	comp3 = singlecomp(theta3,y,self_x, other_x, self_y,other_y,indicators_x,indicators_y)
	
	comp = hcat(comp1',comp2',comp3') # 160x3 matrix explained above
	
	return comp

end

# ╔═╡ b77e0c80-4780-11ec-2a25-cbc0a9d7ca2f
# This function returns a 160x3 matrix containing the posterior probability of an individual to be of type one (first column), two (second column) and three (third column). Element in the first row first column contains the ex-post probability of individual 1 being of type 1, element in the first row second column contains the ex-post probability of individual 1 being of type 2 etc.
# mix is a 1x3 vector containing the shares of type one, two and three in the population, what we called (π_1, π_2, π_3)

function estep(vall,y,self_x, other_x, self_y,other_y,indicators_x,indicators_y,mix)
	
	comp = compdens(vall,y,self_x, other_x,self_y,other_y,indicators_x,indicators_y)
	
# m is a 160x3 matrix. first colunm is the weight for the first component, second column the weight for the second component etc
	
	m = ones(160,3) .* mix'
	
	numerator = m .* comp # 160x3 matrix. each element is π_k f(θ_k,σ_k;X,Y,C_i)
	
	denominator = ones(160,3) .* (sum(numerator, dims=2)) # 160x3 matrix. elements in the same rows are equal. They are the denominator for individual 1 up to 160
	
# tau is our 160x3 matrix of ex-post probabilities for individual i belonging to type k
	
	tau = ones(160,3) .* (numerator ./ denominator)

	return tau
	
end

# ╔═╡ 930a4d20-477e-11ec-38a8-d57693e5c84e
# This function returns the negative of the sum of the individual observation log likelihood weighted by the ex-post probabilities for individuak i to belong to type k. This is the function to minimize in the M-step

function maxloglike(vall,y,self_x, other_x, self_y,other_y,indicators_x,indicators_y, tau)
	
	comp = compdens(vall,y,self_x, other_x,self_y,other_y,indicators_x,indicators_y)
	
	mll = -sum(tau .* log.(comp))
	
	return mll 
	
end

# ╔═╡ 93006210-477e-11ec-09d7-192bb2e62b61
# Define the M-step function. This function minimizes maxloglike and returns two vectors. mix is a 1x3 vector containing the shares of the types in the population, namely (π_1, π_2, π_3) that will be used in the E-step. vall_res is a 1x15 vector containing the new sets of parameters (θ_1,σ_1,θ_2,σ_2,θ_3,σ_3) that will be used in the E-step

function mstep(vall,y,self_x, other_x, self_y,other_y,indicators_x,indicators_y, tau)
	
	ret = optimize(x->maxloglike(x,y,self_x, other_x, self_y,other_y,indicators_x,indicators_y, tau),
			vall,method=NelderMead(),iterations = 10000)
	
	mix = mean(tau, dims=1)
	vall_res = Optim.minimizer(ret)
	
	return mix, vall_res 
	
end

# ╔═╡ 0cf5a4f0-47c4-11ec-2f4f-b7efee6740c6
# This function returns the finite mixture negative log likelihood we want to minimize. This is the sum over individuals of the negative log of what we called L(Ψ;X,Y,C_i). psi is a 1x17 vector containing the 5 parameters for each type (5x3) and the two mixtures. It's what we defined as Ψ.

function totallike(psi,y,self_x, other_x, self_y,other_y,indicators_x,indicators_y)
	
	vall = psi[1:end-2] 
	mix  = psi[end-1:end] 
	mix  = [mix; 1-sum(mix)]
	
	comp = compdens(vall,y,self_x, other_x,self_y,other_y,indicators_x,indicators_y)
	m = ones(160,3) .* mix'
	
	tll = -sum(log.(sum(m .* comp, dims=2)))
	
	return tll
	
end

# ╔═╡ 92c13630-477e-11ec-13b8-3b6951af272c
# The function estimates runs the EM-algoritm and returns a series of objects.
# theta_v is the 1x15 vector containing the estimates on the parameters for each type (θ_1,σ_1, θ_2,σ_2, θ_3,σ_3)
# mix is a 1x3 vector containing the types' shares (π_1, π_2, π_3)
# the third object returned is the value of the log likelihood obtained
# tau is a 160x3 matrix of ex-post probabilities
# psi is a 1x17 vector containing theta_v and the first two mixtures

# Note that the prints show up in the console, not in the notebook

function estimate(y,self_x, other_x, self_y,other_y,indicators_x,indicators_y,maxiter,diffcrit)

# Starting conditions. random starting guesses close to the pooled ones. dimension is 1x15. 5 parameters for each type preferences (3 in our case)
	
	theta_v = repeat(res_s1,3) .* rand(Uniform(0.8,1.2), 15) 
	
# Starting guesses for the mixtures
	
	m0 = rand(Uniform(0.5,2.0), 3)
	mix = m0/sum(m0)
	
# Initial EM control flow parameters
	
	iter  = 1 # nr iterations
	diff  = 1 # difference
	llold = 0 # old log like
	llnew = 1 # new log like

	while iter < maxiter+1 && abs(diff) > diffcrit 
	
		println("EM type maximization: ", iter, "/", maxiter)
		println("\nWorking on E-step. \n")
		
		tau_em = estep(theta_v,y,self_x, other_x, self_y,other_y,indicators_x,indicators_y,mix)
		
		println("Working on M-step. \n")
			
		m_step = mstep(theta_v,y,self_x, other_x, self_y,other_y,indicators_x,indicators_y, tau_em)
		
		mix = m_step[1][:]
		theta_v = m_step[2]
		psi = [m_step[2];mix[1:2]]
		
		println(" sum of mixtures:", sum(mix), "\n")
		
# Adjust EM control flow parameters	
		
		iter =  iter + 1
		llold = llnew
		llnew = totallike(psi,y,self_x, other_x, self_y,other_y,indicators_x,indicators_y)
		diff = llnew - llold
		
		println("Achieved Log Likelihood:", round(llnew,digits=4),"\n \n")
	end
	
	println("End of EM-steps \n")
	
# Maximize the likelihood directly. 
		
	res = optimize(x->totallike(x,y,self_x, other_x, self_y,other_y,indicators_x,indicators_y),
			psi,method=NelderMead(),iterations = 1000)

	println("Achieved final Log Likelihood:", round(Optim.minimum(res),digits=4),"\n")
	
	psi = Optim.minimizer(res)
	theta_v = psi[1:end-2]
	mix = psi[end-1:end]
	mix = [mix;1-sum(mix)]
	
# Retrieve the taus
	
	tau = estep(theta_v,y,self_x, other_x, self_y,other_y,indicators_x,indicators_y,mix)
	
	return theta_v, mix, Optim.minimum(res), tau, psi
end

# ╔═╡ 1d4d9680-47cb-11ec-021f-5511d2dd1dc0
begin

# This cell computes the estimates. Rarely the last minimization of the objective function returns a domain error when too close to the minimum. If that happens just restart the cell. 
# Note that prints will be shown in the terminal and not in the notebook

solmixture = estimate(dt1.choice_x,dt1.self_x,dt1.other_x,dt1.self_y,dt1.other_y,indicators_x,indicators_y,15,5e-5)
	
end;

# ╔═╡ 77c1e520-4c3b-11ec-3db8-a58c85ca0e8f
begin

# the results we want

theta_v = solmixture[1]           # parameters for the three types
mixture = solmixture[2]           # the mixture probabilities
likelihood = solmixture[3]        # value of the likelihood
posterior_prob = solmixture[4]    # individual posterior prob. of being of a type
psi = solmixture[5] 
results_sm1 = round.([psi[1:4];exp(psi[5]);psi[6:9];exp(psi[10]);psi[11:14];exp(psi[15]);psi[16:17];1-sum(psi[16:17])],digits=3) 
	
end;

# ╔═╡ 2c418020-4c44-11ec-24f8-13b7ff1ffcdf
begin

# For a more detailed explanation on how to compute cluster robust standard errors we refer the reader to the python notebook or to David A. Freedman, 2006, "On The So-Called 'Huber Sandwich Estimator' and 'Robust Standard Errors'", The American Statistician, 60:4, 299-302).
	
# Compute the individual cluster robust standard errors. 
# The formula is the following:

# adj * (inv_hessian * grad_con * inv_hessian). 

# adj is a correction for the degree of freedoms and the number of clusters. inv_hessian is the inverse of the hessian of the negative log-likelihood evaluated at the minimum. grad_con is is a 5x5 matrix of gradient contributions. It is computed as the sum over all clusters of the outer product of the sum of the single observation gradient evaluated in the minimum for all observations of an individual cluster.


# Function to compute mixture negative log likelihood. Same as totallike but here mix is the full vector, so that the hessian is a 18x18 matrix

function totallike_mix(psi,y,self_x, other_x, self_y,other_y,indicators_x,indicators_y)
	
	# psi contains the 15 parameters, 5 for each type and the three mixture prob.
	
	vall = psi[1:15]   # 1x15 vector containing parameters for each type
	mix  = psi[16:end] # three mixture prob

	comp = compdens(vall,y,self_x, other_x,self_y,other_y,indicators_x,indicators_y)
	m = ones(160,3) .* mix'
	
	tll = -sum(log.(sum(m .* comp, dims=2)))
	
	return tll
	
end

# compute the inverse of the hessian

hessian_ms1 = ForwardDiff.hessian(x->totallike_mix(x,dt1.choice_x,dt1.self_x,dt1.other_x,dt1.self_y,dt1.other_y,indicators_x,indicators_y),[theta_v;mixture])

inv_hess_ms1 = inv(hessian_ms1)
	
end;

# ╔═╡ f2b8a210-4df7-11ec-1676-6dd04ccc9a65
begin

# Define an auxiliary function.
# This function takes as input a dataset and a scalar and returns a dataframe containing only observations for individual whose ID is equal to that scalar

function getind(dataset,sid)
    dataset_ind = dataset[dataset.sid .== sid,:]
    return dataset_ind
end

function ind_singlecomp(v,y,self_x, other_x, self_y,other_y,indicators_x,indicators_y)
	
	beta  = v[1:(length(v)-1)] # vector of parameters alpha,beta,gamma,delta
    sigma = exp(v[length(v)])  # sigma. exp to force it positive during search
	
    lli = (indicators_x * beta)  # We get the matrix with rows [αs,βr,γq,δv]
    rli = (indicators_y * beta)
	
	uleft  = sigma .* ( (1 .- lli) .* self_x .+ lli .* other_x)
	uright = sigma .* ( (1 .- rli) .* self_y .+ rli .* other_y)
	
	probs  = (exp.(uleft)./(exp.(uleft).+exp.(uright))).^y .* (exp.(uright)./(exp.(uleft).+exp.(uright))).^(1 .-y)
	
	indlike = prod(probs)
		
	if indlike==0.0 indlike=1e-8 end
	
	return indlike
	
end	

function ind_compdens(vall,y,self_x, other_x,self_y,other_y,indicators_x,indicators_y)
	
	# vall is a 1x15 vector that contains the parameters θ for component 1,2 and 3	
	
	theta1 = vall[1:5]    # parameters for first component  α1, β1, γ1, δ1, σ1
	theta2 = vall[6:10]   # parameters for second component α2, β2, γ2, δ2, σ2
	theta3 = vall[11:15]  # parameters for third component  α3, β3, γ3, δ3, σ3
	
	comp1 = ind_singlecomp(theta1,y,self_x, other_x, self_y,other_y,indicators_x,indicators_y)
	
	comp2 = ind_singlecomp(theta2,y,self_x, other_x, self_y,other_y,indicators_x,indicators_y)
	
	comp3 = ind_singlecomp(theta3,y,self_x, other_x, self_y,other_y,indicators_x,indicators_y)
	
# comp is a 1x3 vector. The first column is component1 (where each element in column is the likelihood of individual i when we use as parameters the ones for the first component). Same for column 2 and 3.
	
	comp = [comp1,comp2,comp3]

end

function ind_totallike(psi,y,self_x, other_x, self_y,other_y,indicators_x,indicators_y)
	
	# psi contains the 15 parameters, 5 for each type and the two mixture prob.
	
	vall = psi[1:15] # 1x15 vector containing parameters for each type
	mix  = psi[16:end] 
	comp = ind_compdens(vall,y,self_x, other_x,self_y,other_y,indicators_x,indicators_y)
	
	tll = -log.(sum(mix .* comp))
	
	return tll
	
end
	
end

# ╔═╡ 985c0c70-4c40-11ec-0c1a-a395646d2ef7
# Define the function that computes the matrix of individual gradient contribution 

function grad_contr(dataset, parameters)
	G_sm1 = zeros(18,18)
	for sid in unique(dataset.sid)
		dt1_ind = getind(dataset,sid)
	 ind_grad = ForwardDiff.gradient(x->ind_totallike(x,dt1_ind.choice_x,dt1_ind.self_x,dt1_ind.other_x,dt1_ind.self_y,dt1_ind.other_y,hcat(dt1_ind.s_x,dt1_ind.r_x,dt1_ind.q,dt1_ind.v),hcat(dt1_ind.s_y,dt1_ind.r_y,dt1_ind.q,dt1_ind.v)), parameters)
	G_sm1 .+= (ind_grad*ind_grad')
	end
	return G_sm1
end

# ╔═╡ aeef235e-4c5d-11ec-0a38-fbe40474af45
begin

# Do the same for session 2 

solmixture2 = estimate(dt2.choice_x,dt2.self_x,dt2.other_x,dt2.self_y,dt2.other_y,indicators_x2,indicators_y2,15,5e-3);

# the results we want

theta_v2 = solmixture2[1]           # parameters for the three types
mixture2 = solmixture2[2]           # the mixture probabilities
likelihood2 = solmixture2[3]        # value of the likelihood
posterior_prob2 = solmixture2[4]    # individual posterior prob. of being of a type
psi2 = solmixture2[5] 
results_sm2 = round.([psi2[1:4];exp(psi2[5]);psi2[6:9];exp(psi2[10]);psi2[11:14];exp(psi2[15]);psi2[16:17];1-sum(psi2[16:17])],digits=3) 
	
end;

# ╔═╡ e85bb2fe-4c6e-11ec-3e51-3929725157f8
begin

# compute the inverse of the hessian

hessian_ms2 = ForwardDiff.hessian(x->totallike_mix(x,dt2.choice_x,dt2.self_x,dt2.other_x,dt2.self_y,dt2.other_y,indicators_x2,indicators_y2),[theta_v2;mixture2])

inv_hess_ms2 = inv(hessian_ms2)

# compute gradient contribution G

grad_contr_ms2 = grad_contr(dt2, [theta_v2;mixture2])

# degree of freedoms and number of clusters adjustment

adj2 = (length(dt2.self_x)-1) / (length(dt2.self_x)-length([theta_v;mixture])) * (160/159)

# Sandwich formula for the clustered standard errors

varcov_ms2 = adj2 .* (inv_hess_ms2 * grad_contr_ms2 * inv_hess_ms2)
se_sm2 = sqrt.(diag(varcov_ms2))
se_sm2 = [se_sm2[1:4]; sqrt(exp(psi2[5])^2*se_sm2[5]^2);se_sm2[6:9]; sqrt(exp(psi2[10])^2*se_sm2[10]^2);se_sm2[11:14]; sqrt(exp(psi2[15])^2*se_sm2[15]^2);se_sm2[16:18]]

# sqrt(exp(res_s1[5])^2*se_s1[5]^2) is the formula to go from standard error exp(sigma) to standard error of sigma

se_sm2 = round.(se_sm2, digits=3)
	
end;

# ╔═╡ 5de4a510-4eb7-11ec-0a34-41d9b6a79e59
begin 
	
# Hypothesis Testing

# We now do some hypothesis testing on the parameters we obtained in session 1 and session 2. We compute the z-test statistics and the corresponding p-values to check if each parameter we obtained is statistically different from zero.

# Compute the z-test statistics and the corresponding p-values to check if the parameters are statistically different from zero

# Session 1

zvalues_s1 = abs.(results_sm1 ./ se_sm1)
pvalues_s1 = 2*(1 .-cdf.(Normal(0,1),zvalues_s1))

# Session 2

zvalues_s2 = abs.(results_sm2 ./ se_sm2)
pvalues_s2 = 2*(1 .-cdf.(Normal(0,1),zvalues_s2))

end;

# ╔═╡ b5b26ef0-4dfc-11ec-3539-b70db077382f
begin

# print the results for session 1. Note that the parameters will not be necessarily in this order: strongly altruistic, moderately altruistic, behindness averse. It depends on the starting conditions of the mixtures and on how the algorithm moves.
	
# There are some discrepancies between the standard errors we obtain and those in Table 2 of the paper. The differences in the SEs for the point estimates of Session 1 are due to small differences in the point estimates themselves which are likely due to the use of a different minimization algorithm (we use Nelder-Mead, which is gradient-free, while the authors write down derivatives in closed form and use BFGS). The differences in the SEs for the mixture probabilities of Sessions 1 and 2 are probably due to the fact that our code takes a shortcut and computes the gradient of a subject's log likelihood instead of computing the sum of the gradients of each single choice's log likelihood for the same subject (which is the theoretically correct way of computing the SEs; see, for example, the discussion in the paper cited above). The reason for taking this shortcut is that the theoretically correct object is complex to compute, given the nature of the log likelihood function, and that the shortcut actually returns the same SEs as the authors' for the other parameters. 

param_names = ["π: Type's shares in the population",
			   "α: Weight on other's payoff when behind",
	           "β: Weight on other's payoff when ahead",
	           "γ: Measure of positive reciprocity",
	           "δ: Measure of negative reciprocity",
			   "σ: Choice sensitivity"]

table21 = DataFrame()
table21.param_name = param_names
table21.str_altr_est = [results_sm1[end];results_sm1[11:15]]
table21.str_altr_se = [se_sm1[end];se_sm1[11:15]]
table21.pval_str = round.([pvalues_s1[end];pvalues_s1[11:15]],digits=3)
table21.mod_altr_est = [results_sm1[end-2];results_sm1[1:5]]
table21.mod_altr_se = [se_sm1[end-2];se_sm1[1:5]]
table21.pval_mod = round.([pvalues_s1[end-2];pvalues_s1[1:5]],digits=3)
table21.beh_aver_est = [results_sm1[end-1];results_sm1[6:10]]
table21.beh_aver_se = [se_sm1[end-1];se_sm1[6:10]]
table21.pval_beh = round.([pvalues_s1[end-1];pvalues_s1[6:10]],digits=3)
table21
	
end

# ╔═╡ b7f90390-4dfc-11ec-2c4d-b3dedc1ed52d
begin

# print the results for session 2

table22 = DataFrame()
table22.param_name = param_names
table22.str_altr_est = [results_sm2[end];results_sm2[11:15]]
table22.str_altr_se = [se_sm2[end];se_sm2[11:15]]
table22.pval_str = round.([pvalues_s2[end];pvalues_s2[11:15]],digits=3)
table22.mod_altr_est = [results_sm2[end-2];results_sm2[1:5]]
table22.mod_altr_se = [se_sm2[end-2];se_sm2[1:5]]
table22.pval_mod = round.([pvalues_s2[end-2];pvalues_s2[1:5]],digits=3)
table22.beh_aver_est = [results_sm2[end-1];results_sm2[6:10]]
table22.beh_aver_se = [se_sm2[end-1];se_sm2[6:10]]
table22.pval_beh = round.([pvalues_s2[end-1];pvalues_s2[6:10]],digits=3)
table22
	
end

# ╔═╡ 06309850-4dff-11ec-04ba-e3b11d017c68
begin

# save the results as a csv file

table2 = DataFrame()
table2.param_name = param_names
table2.str_altr_est_s1 = [results_sm1[end];results_sm1[11:15]]
table2.str_altr_se_s1 = [se_sm1[end];se_sm1[11:15]]
table2.pval_str_s1 = round.([pvalues_s1[end];pvalues_s1[11:15]],digits=3)
table2.mod_altr_est_s1 = [results_sm1[end-2];results_sm1[1:5]]
table2.mod_altr_se_s1 = [se_sm1[end-2];se_sm1[1:5]]
table2.pval_mod_s1 = round.([pvalues_s1[end-2];pvalues_s1[1:5]],digits=3)
table2.beh_aver_est_s1 = [results_sm1[end-1];results_sm1[6:10]]
table2.beh_aver_se_s1 = [se_sm1[end-1];se_sm1[6:10]]
table2.pval_beh_s1 = round.([pvalues_s1[end-1];pvalues_s1[6:10]],digits=3)
table2.str_altr_est_s2 = [results_sm2[end-1];results_sm2[6:10]]
table2.str_altr_se_s2 = [se_sm2[end-1];se_sm2[6:10]]
table2.pval_str_s2 = round.([pvalues_s2[end-1];pvalues_s2[6:10]],digits=3)
table2.mod_altr_est_s2 = [results_sm2[end];results_sm2[11:15]]
table2.mod_altr_se_s2 = [se_sm2[end];se_sm2[11:15]]
table2.pval_mod_s2 = round.([pvalues_s2[end];pvalues_s2[11:15]],digits=3)
table2.beh_aver_est_s2 = [results_sm2[end-2];results_sm2[1:5]]
table2.beh_aver_se_s2 = [se_sm2[end-2];se_sm2[1:5]]
table2.pval_beh_s2 = round.([pvalues_s2[end-2];pvalues_s2[1:5]],digits=3)

CSV.write(raw"../output/table2_julia.csv", table2)
	
end

# ╔═╡ Cell order:
# ╠═7b7526d0-3c9d-11ec-3575-37a66c520386
# ╠═9d2d5b40-6023-11ec-0e3b-e5fe5c8d527b
# ╠═7b36490e-3c9d-11ec-1ea6-dbf0082362dc
# ╠═600cd6e0-4ea4-11ec-2071-13c721a5b90e
# ╠═61e5b590-4ea4-11ec-1bf8-a900f11ecb18
# ╠═a0736b60-3c7d-11ec-1b71-27fd38cbf6e4
# ╠═e0965590-3c7d-11ec-13b6-a7d15db3fd9f
# ╠═d8928282-3c85-11ec-00ca-c5215c731b33
# ╠═47094410-4ea6-11ec-2865-47e35e532cf1
# ╠═ca59c720-47a5-11ec-3cf6-f1b17e96c3e9
# ╠═b77e0c80-4780-11ec-2a25-cbc0a9d7ca2f
# ╠═930a4d20-477e-11ec-38a8-d57693e5c84e
# ╠═93006210-477e-11ec-09d7-192bb2e62b61
# ╠═0cf5a4f0-47c4-11ec-2f4f-b7efee6740c6
# ╠═92c13630-477e-11ec-13b8-3b6951af272c
# ╠═1d4d9680-47cb-11ec-021f-5511d2dd1dc0
# ╠═77c1e520-4c3b-11ec-3db8-a58c85ca0e8f
# ╠═2c418020-4c44-11ec-24f8-13b7ff1ffcdf
# ╠═f2b8a210-4df7-11ec-1676-6dd04ccc9a65
# ╠═985c0c70-4c40-11ec-0c1a-a395646d2ef7
# ╠═aeef235e-4c5d-11ec-0a38-fbe40474af45
# ╠═e85bb2fe-4c6e-11ec-3e51-3929725157f8
# ╠═5de4a510-4eb7-11ec-0a34-41d9b6a79e59
# ╠═b5b26ef0-4dfc-11ec-3539-b70db077382f
# ╠═b7f90390-4dfc-11ec-2c4d-b3dedc1ed52d
# ╠═06309850-4dff-11ec-04ba-e3b11d017c68
