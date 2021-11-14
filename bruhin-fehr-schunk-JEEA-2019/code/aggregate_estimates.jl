### A Pluto.jl notebook ###
# v0.14.0

using Markdown
using InteractiveUtils

# ╔═╡ a0736b60-3c7d-11ec-1b71-27fd38cbf6e4
# Import packages used

using Distributions
using DelimitedFiles
using DataFrames
using Optim
using CSV
using Random
using Distributions
using Statistics
using FiniteDiff
using LinearAlgebra
using ForwardDiff

# ╔═╡ 7b7526d0-3c9d-11ec-3575-37a66c520386
# The code in this notebook performs the aggregate estimates (table 1) for the model in the paper 'The Many Faces of Human Sociality: Uncovering the Distribution and Stability 
# of Social Preferences' by Adrian Bruhin, Ernst Fehr and Daniel Schunk (JEEA-2019).

# ╔═╡ 7b36490e-3c9d-11ec-1ea6-dbf0082362dc
# Authors: 
# Massimiliano Pozzi (pozzi.massimiliano@studbocconi.it)
# Salvatore Nunnari (salvatore.nunnari@unibocconi.it)

# ╔═╡ e0965590-3c7d-11ec-13b6-a7d15db3fd9f
# 1. Data cleaning and data preparation

# We import the relevant dataset containing the data on the 39 dictator games and 78 reciprocity games in Session 1. We then remove from the dataset those subjects that 
# behaved very inconsistenly throughout the games. These subjects are identified from the individual estimates that we do not run in this notebook. 
# dropped_subjects_section4paragraph2.csv contains the id of these subjects.

# Load the data on Session 1
data, header = readdlm(raw"../input/choices_exp1.csv",',',header=true)
dt1 = DataFrame(data, vec(header))

# Load dropped_subjects_section4paragraph2.csv 
data_drop = readdlm(raw"../input/dropped_subjects_section4paragraph2.csv",',',header=false)

# drop the subjects whose ids are listed in the data_drop dataframe (14 individuals)
for i=1:length(data_drop)
	global dt1=dt1[dt1.sid .!= data_drop[i],:] # need global else variables inside the loop can’t access variables outside the loop
end;

# ╔═╡ d8928282-3c85-11ec-00ca-c5215c731b33
# Create indicators_x and indicators_y. These are the columns s, r, q, v

indicators_x = hcat(dt1.s_x,dt1.r_x,dt1.q,dt1.v)
indicators_y = hcat(dt1.s_y,dt1.r_y,dt1.q,dt1.v);

# ╔═╡ 56a26400-422c-11ec-3599-77a4d3d0e180
# Define the function to minimize. This is the negative of the log likelihood of observing our data given the parameters of the model. For more information refer to equation
# (4) and (5) in the paper

# v is the vector of parameters (θ,σ)
# y is the choice of the player
# self_x is the payoff when choosing x (left)
# self_y is the payoff when choosing y (right)
# indicators are the ones explained above

function loglikeM(v,y,self_x, other_x, self_y,other_y,indicators_x,indicators_y)
	
	beta  = v[1:(length(v)-1)]  # vector of parameters alpha,beta,gamma,delta
    sigma = exp(v[length(v)])   # choice sensitivity σ. We take the exp since it can only take positive values
    lli = (indicators_x * beta) # we obtain a 1x18720 vector. Each element is (αs+βr+γq+δv) for the single game
    rli = (indicators_y * beta)
	uleft  = sigma.*( (1 .-lli) .*self_x .+lli .*other_x) # utility when choosing allocation X (left)
	uright = sigma.*( (1 .-rli) .*self_y .+rli .*other_y) # utility when choosing 	allocation Y (right)
	# probs is a vector 1x18720 containing the likelihood of observing the data of a single game
	probs  = (exp.(uleft)./(exp.(uleft).+exp.(uright))).^y .* (exp.(uright)./
		(exp.(uleft).+exp.(uright))).^(1 .-y)
	nll = -sum(log.(probs)) # negative log-likelihood
	return nll
	end

# ╔═╡ dfcceb10-3c7d-11ec-2664-9d36301add75
# initialize random starting guesses. Same random initial guesses as in the R code

beta_init = rand(Uniform(0.01,0.02),4) # close to zero for α, β, γ, δ
sigma_init = log.(rand(Uniform(0.05,0.8),1)/mean([mean(dt1.self_x),mean(dt1.other_x),mean(dt1.self_y),mean(dt1.other_y)])) # log since in function we take exp
v0 = [beta_init;sigma_init]

# ╔═╡ a3f76eb0-3d51-11ec-20e7-0d245610572b
# Wrap our function. Optim.optimize can only take one argument

function ll(v)
	fll = loglikeM(v,dt1.choice_x,dt1.self_x,dt1.other_x,dt1.self_y,dt1.other_y,indicators_x,indicators_y)
	return fll
end

# ╔═╡ 1baa58c0-4494-11ec-1794-6b4d74f6812b
# This cell computes the aggregate estimate for session 1 using the BFGS algortihm. We do not provide the analytical gradient since the one computed numerically by the 
# algorithm is precise enough for this problem

sol = optimize(ll,v0,method=BFGS(),iterations = 10000,store_trace=true, extended_trace=true)
res_s1 = Optim.minimizer(sol) # estimates
results_s1 = [res_s1[1:4];exp(res_s1[5])] # Need to take exp for sigma
vloglike_s1 = Optim.minimum(sol) # Value of loglikelihood evaluated in the minimum

# ╔═╡ f576544e-3d55-11ec-009e-5ddcc502108b
# Compute the individual cluster robust standard errors. The formula is the following:

# adj * (inv_hessian * grad_con * inv_hessian). 

# adj is a correction for the degree of freedoms and the number of clusters.
# inv_hessian is the inverse of the hessian of the negative log-likelihood evaluated at the minimum.
# grad_con is is a 5x5 matrix of gradient contributions. It is computed as the sum over all clusters of the outer product of the sum of the single observation gradient 
# evaluated in the minimum for all observations of an individual cluster.


# Define the function that computes the matrix of individual gradient contribution 

function congradloglike(v,y,self_x,other_x,self_y,other_y,indicators_x,indicators_y,clusters) 
	
	# parameters. Same as before
	beta  = v[1:(length(v)-1)] 
    sigma = exp(v[length(v)])
	
	# This part is the same as before
    lli = (indicators_x * beta)  
    rli = (indicators_y * beta)
	utl  = sigma.*( (1 .-lli) .*self_x .+lli .*other_x)
	utr  = sigma.*( (1 .-rli) .*self_y .+rli .*other_y)
	probs  = (exp.(utl)./(exp.(utl).+exp.(utr))).^y .* (exp.(utr)./(exp.(utl).+exp.(utr))).^(1 .-y)
	
	# Compute the analytical gradient
	probsm = ones(size(indicators_x,1), size(indicators_x,2)) .* probs
	u = ones(size(indicators_x,1), size(indicators_x,2)) .*exp.(utl) 
	w = ones(size(indicators_x,1), size(indicators_x,2)) .*(exp.(utl) .+ exp.(utr))
	up = sigma * indicators_x .* (ones(size(indicators_x,1), size(indicators_x,2)) .* (other_x .-self_x)) .* u
	wp = up .+ sigma .* indicators_y .* (ones(size(indicators_x,1), size(indicators_x,2)) .*((other_y.-self_y).*exp.(utr)))

	bgradi = (-1).^(1 .-y).*((up.*w .- u.*wp)./w.^2)./probsm
	up2 = utl .* u[:,1]
	wp2 = up2 .+ utr .* exp.(utr)
	ggradi = ((-1).^(1 .-y) .* (((up2.*w[:,1] - u[:,1].*wp2)./w[:,1].^2)./probs))

	gradi = [bgradi ggradi] # 18720x5 matrix of partial derivatives for each obs. 
	
	cl = unique(clusters)
	j = length(cl) # nr of clusters is length of unique individual = 160
	k = size(gradi,2) # nr of columns in gradi = 5
	sandwich = zeros(k,k)
	for i=1:j # sum columns by columns for individual i
		sel = clusters .== cl[i] # vector. 1 if same cluster else 0
		sel = findall(>(0), sel) # indices of 1
		gradsel = sum(gradi[sel,:],dims=1)
		sandwich .+= (gradsel' * gradsel) # sum of outer product
	end
	return sandwich
end

# ╔═╡ f558e140-3d55-11ec-075d-9d31f531f525
# Compute the cluster robust standard errors. 

# This would be the code to use the inverse of the hessian stored in the BFGS algorithm. This is less precise than the one computed with ForwardDiff, so we don't use it but we 
# leave the code since it could be interesting to know.

# inverse_hessian_s1 = sol.trace[end].metadata["~inv(H)"] 


# Compute the hessian numerically using the ForwardDiff package 
inv_hess_num = inv(ForwardDiff.hessian(ll,res_s1)) 

# degree of freedom and cluster adjustment
adj = (length(dt1.self_x)-1) / (length(dt1.self_x)-length(res_s1)) * (160/159)

grad_contribution = congradloglike(res_s1,dt1.choice_x,dt1.self_x,dt1.other_x,dt1.self_y,dt1.other_y,indicators_x,indicators_y,dt1.sid)

# Sandwich formula for the clustered standard errors

varcov_s1 = adj .* (inv_hess_num * grad_contribution * inv_hess_num)
se_s1 = sqrt.(diag(varcov_s1))
se_s1 = [se_s1[1:4]; sqrt(exp(res_s1[5])^2*se_s1[5]^2)]

# sqrt(exp(res_s1[5])^2*se_s1[5]^2) is the formula to go from standard error log(sigma) to standard error of sigma. We use the Delta method

# ╔═╡ f970b080-3cb2-11ec-0e09-65acfacc9837
# Compute aggregate estimates for session 2

# ╔═╡ 9e5d90a0-3cb2-11ec-2efc-8189c078aee7
# Load the data for Session 2 and remove the dropped subjects as before

data, header = readdlm(raw"../input/choices_exp2.csv",',',header=true)
dt2 = DataFrame(data, vec(header))
for i=1:length(data_drop)
	global dt2=dt2[dt2.sid .!= data_drop[i],:] # need global else variables inside the loop can’t access variables outside the loop
end;

# ╔═╡ e6ab8b50-3cb2-11ec-0b39-e5424b4fdbd2
# Create indicators_x and indicators_y. These are the columns s, r, q, v

indicators_x2 = hcat(dt2.s_x,dt2.r_x,dt2.q,dt2.v)
indicators_y2 = hcat(dt2.s_y,dt2.r_y,dt2.q,dt2.v);

# ╔═╡ d554ea40-3cb2-11ec-1642-99c3aadfd3e4
# This cell computes the aggregate estimate for session 2

function ll2(v)
	fll = loglikeM(v,dt2.choice_x,dt2.self_x,dt2.other_x,dt2.self_y,dt2.other_y,indicators_x2,indicators_y2)
	return fll
end

sol = optimize(ll2,v0,method=BFGS(),iterations = 10000,store_trace=true, extended_trace=true)
res_s2 = Optim.minimizer(sol)
results_s2 = [res_s2[1:4];exp(res_s2[5])]
vloglike_s2 = Optim.minimum(sol)

# ╔═╡ 586df880-3cd7-11ec-348c-6133b9ddba73
# Compute the cluster robust standard errors. 

inv_hessian_num2 = inv(ForwardDiff.hessian(ll2,res_s2))

adj = (length(dt2.self_x)-1) / (length(dt2.self_x)-length(res_s2)) * (160/159) # degree of freedom adjustment

grad_contribution = congradloglike(res_s2,dt2.choice_x,dt2.self_x,dt2.other_x,dt2.self_y,dt2.other_y,indicators_x2,indicators_y2,dt2.sid)

# Sandwich formula for the clustered standard errors

varcov_s2 = adj .* (inv_hessian_num2 * grad_contribution * inv_hessian_num2)
se_s2 = sqrt.(diag(varcov_s2))
se_s2 = [se_s2[1:4]; sqrt(exp(res_s2[5])^2*se_s2[5]^2)]

# ╔═╡ 802c89c0-3d60-11ec-1d69-0dcf94880bf3
# Save the results in a csv file

table1 = DataFrame()
table1.est_s1 = results_s1
table1.se_s1 = se_s1
table1.est_s2 = results_s2
table1.se_s2 = se_s2
CSV.write(raw"../output/table1_julia.csv", table1)

# ╔═╡ e071b690-3c7d-11ec-0911-bd86ef5668aa
# We write here the same negative log-likelihood function as before written using a loop over all observations. This is slower but could be interesting to see

function loglike_dummy_ind(v,y,self_x,other_x,self_y,other_y,indicators_x,indicators_y)
	beta  = v[1:(length(v)-1)] # vector of parameters alpha,beta,gamma,delta
    sigma = exp(v[length(v)])  # sigma. exp to force it positive during search
	len = length(dt1.self_x)
	probs = Vector{Float64}(undef, len)
	for j=1:len
    lli = (indicators_x[j,1:4] .* beta)  # We get the vector [αs,βr,γq,δv]
    rli = (indicators_y[j,1:4] .* beta)
    uleft  = (sigma*((1-lli[1]-lli[2]-lli[3]-lli[4])*self_x[j]+ (lli[1]+lli[2]+lli[3]+lli[4]) *other_x[j])) # Utility when choosing left
    uright = (sigma*((1-rli[1]-rli[2]-rli[3]-rli[4])*self_y[j]+ (rli[1]+rli[2]+rli[3]+rli[4]) *other_y[j])) # Utility when choosing right
    probability = log((exp(uleft)/(exp(uleft)+exp(uright)))^y[j] * (exp(uright)/(exp(uleft)+exp(uright)))^(1-y[j])) # probability. See equation (4) in the paper
	probs[j]=probability
	end
    nsumloglike = -sum((probs)) # negative log likelihood
    return nsumloglike
end

# ╔═╡ Cell order:
# ╠═7b7526d0-3c9d-11ec-3575-37a66c520386
# ╠═7b36490e-3c9d-11ec-1ea6-dbf0082362dc
# ╠═a0736b60-3c7d-11ec-1b71-27fd38cbf6e4
# ╠═e0965590-3c7d-11ec-13b6-a7d15db3fd9f
# ╠═d8928282-3c85-11ec-00ca-c5215c731b33
# ╠═56a26400-422c-11ec-3599-77a4d3d0e180
# ╠═dfcceb10-3c7d-11ec-2664-9d36301add75
# ╠═a3f76eb0-3d51-11ec-20e7-0d245610572b
# ╠═1baa58c0-4494-11ec-1794-6b4d74f6812b
# ╠═f576544e-3d55-11ec-009e-5ddcc502108b
# ╠═f558e140-3d55-11ec-075d-9d31f531f525
# ╠═f970b080-3cb2-11ec-0e09-65acfacc9837
# ╠═9e5d90a0-3cb2-11ec-2efc-8189c078aee7
# ╠═e6ab8b50-3cb2-11ec-0b39-e5424b4fdbd2
# ╠═d554ea40-3cb2-11ec-1642-99c3aadfd3e4
# ╠═586df880-3cd7-11ec-348c-6133b9ddba73
# ╠═802c89c0-3d60-11ec-1d69-0dcf94880bf3
# ╠═e071b690-3c7d-11ec-0911-bd86ef5668aa
