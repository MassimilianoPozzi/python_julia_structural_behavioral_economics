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
using Statistics
using FiniteDiff
using LinearAlgebra
using ForwardDiff
using JuMP 
import Ipopt
	
end

# ╔═╡ 7b7526d0-3c9d-11ec-3575-37a66c520386
# Replication of Bruhin, Fehr, and Schunk, 2019, "Many Faces of Human Sociality: Uncovering the Distribution and Stability of Social Preferences"

# ╔═╡ 7b36490e-3c9d-11ec-1ea6-dbf0082362dc
# Authors: 
# Massimiliano Pozzi (Bocconi University, pozzi.massimiliano@studbocconi.it)
# Salvatore Nunnari (Bocconi University, https://snunnari.github.io, salvatore.nunnari@unibocconi.it)

# ╔═╡ c7fb5b79-7d87-49c5-9459-a360e1057715
# This notebook works with both Neptune and Pluto; the "begin" and "end" tags in each cell are necessary for Pluto notebooks but do not play any role in Neptune 

# ╔═╡ 7d517760-4913-11ec-24fb-03122e4dda9a
# This notebook was tested with the following packages versions:

# Pozzi: julia 1.6.1, Pluto 0.14.7, Neptune 0.14.0, Distributions 0.25.17, DataFrames 1.2.2, Optim 1.4.1, CSV 0.9.5, FiniteDiff 2.8.1, ForwardDiff 0.10.20, JuMP 0.21.10


# ╔═╡ e0965590-3c7d-11ec-13b6-a7d15db3fd9f
begin

# 1. Data Cleaning and Preparation

# We import the relevant datasets containing the data on the 39 dictator games and 78 reciprocity games in Session 1 and Session 2. As the authors, we remove from these datasets those subjects that behaved very inconsistenly throughout the games. These subjects are identified from the individual estimates that we do not run in this notebook. The file "dropped_subjects_section4paragraph2.csv" contains the ID of these subjects.

# Import data for session 1 and create a DataFrame
data, header = readdlm(raw"../input/choices_exp1.csv",',',header=true)
dt1 = DataFrame(data, vec(header))

# Import data with ID of subjects to drop
data_drop = readdlm(raw"../input/dropped_subjects_section4paragraph2.csv",',',header=false)
	
# Drop the session 1 subjects whose IDs are listed in the data_drop dataframe (14 individuals)

dt1 = dt1[[i for i=1:length(dt1.sid) if !(dt1.sid[i] in(data_drop))],:]
	
end

# ╔═╡ d8928282-3c85-11ec-00ca-c5215c731b33
begin

# Create indicators_x and indicators_y. These are the columns s, r, q, v

indicators_x = hcat(dt1.s_x,dt1.r_x,dt1.q,dt1.v)
indicators_y = hcat(dt1.s_y,dt1.r_y,dt1.q,dt1.v)
	
end;

# ╔═╡ 56a26400-422c-11ec-3599-77a4d3d0e180
# Define the function to minimize. This is the negative of the log likelihood of observing our data given the parameters of the model.

# For more information we refer the reader to equation (4) and (5) in the paper or to the python notebook in this repository

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
	
	# probs is a 1x18720 vector containing the likelihood of observing the data of a single game
	probs  = (exp.(uleft)./(exp.(uleft).+exp.(uright))).^y .* (exp.(uright)./
		(exp.(uleft).+exp.(uright))).^(1 .-y)
	nll = -sum(log.(probs)) # negative log-likelihood
	return nll
	end

# ╔═╡ dfcceb10-3c7d-11ec-2664-9d36301add75
begin

# initialize random starting guesses. Same random initial guesses as in the authors' R code

beta_init = rand(Uniform(0.01,0.02),4) # close to zero for α, β, γ, δ
sigma_init = log.(rand(Uniform(0.05,0.8),1)/mean([mean(dt1.self_x),mean(dt1.other_x),mean(dt1.self_y),mean(dt1.other_y)])) # log since in function we take exp
v0 = [beta_init;sigma_init]
	
end

# ╔═╡ a3f76eb0-3d51-11ec-20e7-0d245610572b
# Wrap our function. Optim.optimize can only take one argument

function ll(v)
	fll = loglikeM(v,dt1.choice_x,dt1.self_x,dt1.other_x,dt1.self_y,dt1.other_y,indicators_x,indicators_y)
	return fll
end

# ╔═╡ 1baa58c0-4494-11ec-1794-6b4d74f6812b
begin

# This cell computes the aggregate estimate for session 1 using the BFGS algortihm. We do not provide the analytical gradient since the one computed by the algorithm is precise enough for this problem

sol = optimize(ll,v0,method=BFGS(),iterations = 10000,store_trace=true, extended_trace=true, autodiff=:forward)
res_s1 = Optim.minimizer(sol) # estimates
results_s1 = [res_s1[1:4];exp(res_s1[5])] # Need to take exp for sigma
vloglike_s1 = -Optim.minimum(sol) # Value of loglikelihood evaluated in the minimum

end

# ╔═╡ f576544e-3d55-11ec-009e-5ddcc502108b
# Compute the individual cluster robust standard errors. 

# For a more detailed explanation on how to compute cluster robust standard errors we refer the reader to the python notebook or to David A. Freedman, 2006, "On The So-Called 'Huber Sandwich Estimator' and 'Robust Standard Errors'", (The American Statistician 60:4, 299-302)

# The formula is the following:

# adj * (inv_hessian * grad_con * inv_hessian). 

# adj is a correction for the degree of freedoms and the number of clusters. inv_hessian is the inverse of the hessian of the negative log-likelihood evaluated at the minimum. grad_con is is a 5x5 matrix of gradient contributions. It is computed as the sum over all clusters of the outer product of the sum of the single observation gradient evaluated in the minimum for all observations of an individual cluster.


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
begin

# Compute the cluster robust standard errors. 

# This is the code to use the inverse of the hessian stored in the BFGS algorithm. We use the one computed using the Pkg ForwardDiff but they should be the same since we used the option autodiff=:forward 

inverse_hessian_s1 = sol.trace[end].metadata["~inv(H)"] 

# Compute the hessian using the ForwardDiff package 
	
inv_hess_num = inv(ForwardDiff.hessian(ll,res_s1)) 

# degree of freedom and cluster adjustment
	
adj = (length(dt1.self_x)-1) / (length(dt1.self_x)-length(res_s1)) * (160/159)

grad_contribution = congradloglike(res_s1,dt1.choice_x,dt1.self_x,dt1.other_x,dt1.self_y,dt1.other_y,indicators_x,indicators_y,dt1.sid)

# Sandwich formula for the clustered standard errors

varcov_s1 = adj .* (inv_hess_num * grad_contribution * inv_hess_num)
se_s1 = sqrt.(diag(varcov_s1))
se_s1 = [se_s1[1:4]; sqrt(exp(res_s1[5])^2*se_s1[5]^2)]

# sqrt(exp(res_s1[5])^2*se_s1[5]^2) is the formula to go from standard error log(sigma) to standard error of sigma. We use the Delta method
	
end

# ╔═╡ f970b080-3cb2-11ec-0e09-65acfacc9837
# Compute aggregate estimates for session 2

# ╔═╡ 9e5d90a0-3cb2-11ec-2efc-8189c078aee7
begin

# Import session 2 data and remove subjects to drop as we did for session 1 above

data2, header2 = readdlm(raw"../input/choices_exp2.csv",',',header=true)
dt2 = DataFrame(data2, vec(header2))
dt2 = dt2[[i for i=1:length(dt2.sid) if !(dt2.sid[i] in(data_drop))],:]
	
end

# ╔═╡ e6ab8b50-3cb2-11ec-0b39-e5424b4fdbd2
begin

# Create indicators_x and indicators_y. These are the columns s, r, q, v

indicators_x2 = hcat(dt2.s_x,dt2.r_x,dt2.q,dt2.v)
indicators_y2 = hcat(dt2.s_y,dt2.r_y,dt2.q,dt2.v)

end;

# ╔═╡ d554ea40-3cb2-11ec-1642-99c3aadfd3e4
begin

# This cell computes the aggregate estimate for session 2

function ll2(v)
	fll = loglikeM(v,dt2.choice_x,dt2.self_x,dt2.other_x,dt2.self_y,dt2.other_y,indicators_x2,indicators_y2)
	return fll
end

sol2 = optimize(ll2,v0,method=BFGS(),iterations = 10000,store_trace=true, extended_trace=true, autodiff=:forward)
res_s2 = Optim.minimizer(sol2)
results_s2 = [res_s2[1:4];exp(res_s2[5])]
vloglike_s2 = -Optim.minimum(sol2)
	
end

# ╔═╡ 586df880-3cd7-11ec-348c-6133b9ddba73
begin

# Compute the cluster robust standard errors. 

inv_hessian_num2 = inv(ForwardDiff.hessian(ll2,res_s2))

adj2 = (length(dt2.self_x)-1) / (length(dt2.self_x)-length(res_s2)) * (160/159) # degree of freedom adjustment

grad_contribution2 = congradloglike(res_s2,dt2.choice_x,dt2.self_x,dt2.other_x,dt2.self_y,dt2.other_y,indicators_x2,indicators_y2,dt2.sid)

# Sandwich formula for the clustered standard errors

varcov_s2 = adj2 .* (inv_hessian_num2 * grad_contribution2 * inv_hessian_num2)
se_s2 = sqrt.(diag(varcov_s2))
se_s2 = [se_s2[1:4]; sqrt(exp(res_s2[5])^2*se_s2[5]^2)]
	
end

# ╔═╡ c32cedfe-daf7-4b0b-81bc-e8028a050294
begin

# Hypothesis Testing

# We now do some hypothesis testing on the parameters we obtained in session 1 and session 2. We first compute the z-test statistics and the corresponding p-values to check if each parameter we obtained is statistically different from zero. We then compute the p-value of a z-test to check if the parameters we obtained in session 1 are statistically different from the parameters we obtained in session 2.

# Compute the z-test statistics and the corresponding p-values to check if the parameters are statistically different from zero

# Session 1

zvalues_s1 = results_s1 ./ se_s1
pvalues_s1 = 2*(1 .-cdf.(Normal(0,1),zvalues_s1))

# Session 2

zvalues_s2 = results_s2 ./ se_s2
pvalues_s2 = 2*(1 .-cdf.(Normal(0,1),zvalues_s2))

# Check if parameters obtained in session 1 are statistically different from parameters in session 2

# First we need the variance for the parameters in session 1 and session2

var_s1 = se_s1.^2
var_s2 = se_s2.^2

# Now we compute the p-values of the z-test statistics 

zvalues_s1s2 = abs.(results_s1.-results_s2) ./ sqrt.(var_s1.+var_s2)
pvalues_s1s2 = 2*(1 .-cdf.(Normal(0,1),zvalues_s1s2))
	
end

# ╔═╡ 802c89c0-3d60-11ec-1d69-0dcf94880bf3
begin

# Create a new table and save the results in a csv file. We round results up to the third decimal

parameters_name = ["α: Weight on other's payoff when behind",
                   "β: Weight on other's payoff when ahead",
                   "γ: Measure of positive reciprocity",
                   "δ: Measure of negative reciprocity",
                   "σ: Choice sensitivity"]

table1 = DataFrame()
table1.parameters = parameters_name
table1.estimates_s1 = round.(results_s1, digits=3)
table1.standarderr_s1 = round.(se_s1, digits=3)
table1.zstat_s1 = round.(zvalues_s1, digits=3)
table1.pval_s1 = round.(pvalues_s1, digits=3)
table1.estimates_s2 = round.(results_s2, digits=3)
table1.standarderr_s2 = round.(se_s2, digits=3)
table1.zstat_s2 = round.(zvalues_s2, digits=3)
table1.pval_s2 = round.(pvalues_s2, digits=3)
table1.pval_s1s2 = round.(pvalues_s1s2, digits=3)

CSV.write(raw"../output/table1_julia.csv", table1)
	
table1
	
end

# ╔═╡ 4ee3a2da-1a91-450b-97d8-4108fa08944e
begin

# We now compute the point estimates for session 1 using JuMP.

# Define a function that takes the variables separately and returns the negative log likelihood. This simplifies the way we can define the variables in the JuMP model later on

function lljump(alpha,beta,gamma,delta,sigma)
		v = [alpha,beta,gamma,delta,sigma]
		flljump = loglikeM(v,dt1.choice_x,dt1.self_x,dt1.other_x,dt1.self_y,dt1.other_y,indicators_x,indicators_y)
	return flljump	
end
	
# Define a new model using JuMP. We use the Ipopt optimizer

model = Model(Ipopt.Optimizer; bridge_constraints = false)
	
# You need to clear the julia variables, else you will not be able to run this cell a second time. This is a known problem when using Pluto 
	
alpha,beta,gamma,delta,sigma= nothing, nothing, nothing, nothing, nothing
	
# Define the variables and the starting values. Our problem is unconstrained.

@variable(model, alpha, start = v0[1])
@variable(model, beta,  start = v0[2])
@variable(model, gamma, start = v0[3])
@variable(model, delta, start = v0[4])
@variable(model, sigma, start = v0[5])
@NLobjective(
        model,
        Min,
        lljump(alpha,beta,gamma,delta,sigma) # the function to minimize
    )
	
optimize!(model)
solution_summary(model)
	
end

# ╔═╡ bbe59950-4915-11ec-289a-c97076f93ca1
begin

# In this cell we show that the point estimates for Session 1 are the same when using Optim or JuMP

JuMP_s1 = [round(value(alpha), digits=3); round(value(beta), digits=3); round(value(gamma), digits=3); round(value(delta), digits=3); round(exp(value(sigma)), digits=3)]

table_diff = DataFrame()
table_diff.parameters = parameters_name
table_diff.Optim_s1 = round.(results_s1, digits=3)
table_diff.JuMP_s1  = JuMP_s1
table_diff
	
end

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
# ╠═c7fb5b79-7d87-49c5-9459-a360e1057715
# ╠═7d517760-4913-11ec-24fb-03122e4dda9a
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
# ╠═c32cedfe-daf7-4b0b-81bc-e8028a050294
# ╠═802c89c0-3d60-11ec-1d69-0dcf94880bf3
# ╠═4ee3a2da-1a91-450b-97d8-4108fa08944e
# ╠═bbe59950-4915-11ec-289a-c97076f93ca1
# ╠═e071b690-3c7d-11ec-0911-bd86ef5668aa
