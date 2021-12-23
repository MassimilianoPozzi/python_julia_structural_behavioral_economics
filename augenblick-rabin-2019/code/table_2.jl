### A Pluto.jl notebook ###
# v0.14.0

using Markdown
using InteractiveUtils

# ╔═╡ feb1b280-268f-11ec-1574-e5037552b641
begin

# Import packages used

using Distributions
using DataFrames
using Optim
using CSV
using DelimitedFiles
using StatFiles
using Statistics
using HypothesisTests	
using Distributions

end

# ╔═╡ d08aaa10-268f-11ec-067f-4322e82ba779
# Augenblick and Rabin, 2019, "An Experiment on Time Preference and Misprediction in Unpleasant Tasks", Table 2

# ╔═╡ fee9d980-268f-11ec-3901-3705290b7c24
# Authors: 

# Massimiliano Pozzi (Bocconi University, pozzi.massimiliano@studbocconi.it)
# Salvatore Nunnari (Bocconi University, https://snunnari.github.io, salvatore.nunnari@unibocconi.it)

# The code in this Jupyter notebook performs the individual estimates to replicate columns 1, 2, 3, 4 of Table 2

# This notebook was tested with the following packages versions:

# Pozzi: julia 1.6.1, Pluto 0.14.7, Neptune 0.14.0, Distributions 0.25.17, DataFrames 1.2.2, Optim 1.4.1, CSV 0.9.5

# ╔═╡ f8c1b410-51b4-11ec-1d65-11df2949bf23
# This notebook works with both Neptune and Pluto; the "begin" and "end" tags in each cell are necessary for Pluto notebooks but do not play any role in Neptune 

# ╔═╡ 0de27c80-51b5-11ec-05e9-a5b37469b0c7
# To import a new package (es. StatFiles) use this cell. Change the name inside the brackets to change Pkg

# import Pkg; Pkg.add("StatFiles")

# ╔═╡ fe9ddc60-268f-11ec-112a-59ea860b2c81
begin

# 1. Data Cleaning and Data Preparation

# We import the dataset containing observations for all 100 individuals. We then remove observations when a bonus was offered and create dummy variables that will be useful for estimation: pb is equal to one if the subject completed 10 mandatory tasks on subject-day (this is used to estimate the projection bias parameter α), ind_effort10 and ind_effort110 are equal to one if the subject completed 10 or 110 tasks respectively (and they are used for the Tobit correction when computing the likelihood)

dt = DataFrame(load(raw"../input/decisions_data.dta")) # full sample
dt = dt[dt.bonusoffered .!=1,:]  # remove observations when a bonus was offered

# We remove observations when a bonus was offered and create the following dummy variables that will be useful for estimation: pb is equal to one if the subject completed 10 mandatory tasks on subject-day (this is used to estimate the projection bias parameter α); ind_effort10 and ind_effort110 are equal to one if, respectively, the subject completed 10 or 110 tasks (and they are used for the Tobit correction when computing the likelihood)

dt.pb = dt.workdone1 ./ 10 # pb dummy variable. workdone1 can either be 10 or 0, so                                dividing the variable by 10 creates our dummy

dt.ind_effort10  = dt.effort .== 10    # ind_effort10 dummy
dt.ind_effort110 = dt.effort .== 110   # ind_effort110 dummy	
	
end;

# ╔═╡ 74270f10-51b5-11ec-2edd-af650c9bec6b
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
	
    netdistance, wage, today, prediction, pb, effort, ind_effort10, ind_effort110, spec = args
    
# Column 4 specification has 7 parameters, the others have 6 parameters. Use if statement to specify the correct params vector	

	if spec!=4 
		beta, betahat, delta, gamma, phi, sigma  = params 
	else
		beta, betahat, delta, gamma, phi, alpha, sigma = params
	end
	
	# compute the predicted choice. To avoid domain error we use max for beta,               betahat, delta, phi
	
	predchoice=((max(phi,1e-6).*(max(delta,1e-6).^netdistance).*(max(beta,1e-6).^today).*(max(betahat,1e-6).^prediction).*wage).^(1 ./(gamma.-1))) 
	
	if spec==4 predchoice = predchoice .- pb.*alpha end
    
	# compute the vector of individual observations' likelihood, a 1x8049 vector
	
    prob = pdf.(Normal(0,max(sigma,0.1)),predchoice.-effort) 
    
	# find index when effort is 10 ot 110 to apply Tobit correction
	len = length(effort)
    index_10  = [i for i=1:len if effort[i]==10]
    index_110 = [i for i=1:len if effort[i]==110]   
    
	# Apply Tobit correction
    for i in index_10
        prob[i] = 1 - cdf.(Normal(0,1),(predchoice[i]-effort[i])/max(sigma,0.1))
	end
        
    for i in index_110
        prob[i] = cdf.(Normal(0,1),(predchoice[i]-effort[i])/max(sigma,0.1))
	end
       
	# find index when prob is equal to 0 or 1 and change it. this is needed to avoid         problems when taking logs
    index_p0 = [i for i=1:len if prob[i]==0]
    index_p1 = [i for i=1:len if prob[i]==1]
    
	# change prob when 0 or 1
    for i in index_p0
        prob[i] = 1e-6
	end
        
    for i in index_p1
        prob[i] = 1 - 1e-6
	end
    
    negll = - sum(log.(prob)) # negative log likelihood
    
    return negll
	
end

# ╔═╡ ce679cb0-51b5-11ec-39ba-3f65d7b1ce33
# Estimation

# We define the function that computes the individual estimates for column 1 up to 4 of table 2. There are four specifications to consider: the first column uses the entire sample, the second only when the date of decision is less than 4, the third one when the date decision is from time 4 onwards and the fourth one uses the full sample but estimates also the projection bias parameter.


# Define an auxiliary function.
# This function takes as input a dataset, a scalar and a specification argument and returns the relevant individual-level dataframe containing only observations for individual whose ID is equal to that scalar for a given specification

function getind(dataset, wid, spec)
    dataset_ind = dataset[dataset.wid .== wid,:]
	# second column specification
	if spec == 2 dataset_ind = dataset_ind[dataset_ind.decisiondatenum .< 4,:]  end 
	# third column specification    
	if spec == 3 dataset_ind = dataset_ind[dataset_ind.decisiondatenum .>= 4,:] end 
    return dataset_ind

end			

# ╔═╡ a5d0925e-5535-11ec-06ca-85bfd5b5081b
# Define the function that estimates the model

function estimate(dataset, spec)
    
# Define the vectors containing the initial guesses and the vectors that will contain our results. res_params will be a 1x100 vector, each element is a 1x6 or 1x7 vector containing the estimates for the parameters. res_like will be a 1x100 vector containing the log likelihood obtained, while res_wid will be a 1x100 vector containing the wid of the individual for which we find the estimates. res_conv is a 1x100 vector of 1 (converged) or 0 (not converges)
    
    beta_init, betahat_init, delta_init, gamma_init, phi_init, alpha_init, sigma_init = 1.0, 1, 1, 2, 250, 3, 50
	
	res_param = []
    res_like  = Vector{Float64}(undef, 100)
    res_wid   = Vector{Float64}(undef, 100)
	res_conv  = Vector{Float64}(undef, 100)
	
# compute the individual estimates and save the results
    
    for wid in unique(dataset.wid)  # loop over individuals' ID
        
        params_init = [beta_init, betahat_init, delta_init, gamma_init, phi_init,                            alpha_init, sigma_init]
        params_nopb = [beta_init, betahat_init, delta_init, gamma_init, phi_init,                            sigma_init]
        
        dt_ind = getind(dataset, wid, spec) # get the relevant individual                                                           dataframe that depends also on the                                                     specification
        
# args needed to compute the individual likelihood
		
        args_ind=[dt_ind.netdistance, dt_ind.wage, dt_ind.today, dt_ind.prediction, dt_ind.pb, dt_ind.effort, dt_ind.ind_effort10, dt_ind.ind_effort110, spec]
        
# We run each estimation twice starting from the previous argmin to see if we were stuck on a local minima. This is usually not necessary, but also not expensive in terms of time so it can be worth it to do it
                
        if spec != 4 # first three columns do not have the projection bias parameter
                
            j = 0
            while (j<3)
                global sol = optimize(x->negloglike(x,args_ind),
				   params_nopb,
				   method = NelderMead(), 
				   iterations = 2500)
				   params_nopb = Optim.minimizer(sol)
                j = j+1
			end
        else
            j = 0
            while (j<3)
                global sol = optimize(x->negloglike(x,args_ind),
				   params_init,
				   method = NelderMead(), 
				   iterations = 4000) # projetion bias needs more iterations
				   params_init = Optim.minimizer(sol)
                j = j+1
			end
		end

    append!(res_param, [Optim.minimizer(sol)])
    res_like[round(Int32, wid)]  = -Optim.minimum(sol)
    res_wid[round(Int32, wid)]   = wid
	res_conv[round(Int32, wid)]  = Optim.converged(sol)
			
	end # end for wid loop
	
    return res_param, res_like, res_conv
																				
end

# ╔═╡ f2177020-5536-11ec-0db6-1b2a59aa1c0e
begin

# This cell computes the estimates for column (1) and stores the results.

res1 = estimate(dt,1)
param1 = res1[1]
like1  = res1[2]
conv1  = res1[3]
	
end;

# ╔═╡ 99348ad0-5543-11ec-274d-1d59ea3e02b0
begin

# This cell computes the estimates for column (2) and stores the results.

res2 = estimate(dt,2)
param2 = res2[1]
like2  = res2[2]
conv2  = res2[3]
	
end;

# ╔═╡ a49523d0-5543-11ec-3eec-fbe419cf66af
begin

# This cell computes the estimates for column (3) and stores the results.

res3 = estimate(dt,3)
param3 = res3[1]
like3  = res3[2]
conv3  = res3[3]
	
end;

# ╔═╡ a3f7d260-5543-11ec-3bd8-3d1836a9524e
begin

# This cell computes the estimates for column (4) and stores the results.

res4 = estimate(dt,4)
param4 = res4[1]
like4  = res4[2]
conv4  = res4[3]
	
end;

# ╔═╡ 7ccc8cbe-5544-11ec-3b68-31c23aaacd87
begin

# unpack the results:

# The params vectors we obtained are a list of lists where the first element is a list containing the estimates for the first individual in the dataset, the second element is a list containing the estimates for the second individual in the dataset, etc. 

# For column 1:
b1  = [item[1] for item in param1]; bh1 = [item[2] for item in param1]
d1  = [item[3] for item in param1]; g1  = [item[4] for item in param1]
p1  = [item[5] for item in param1]; s1  = [item[6] for item in param1]

# For column 2:
b2  = [item[1] for item in param2]; bh2 = [item[2] for item in param2]
d2  = [item[3] for item in param2]; g2  = [item[4] for item in param2]
p2  = [item[5] for item in param2]; s2  = [item[6] for item in param2]

# For column 3:
b3  = [item[1] for item in param3]; bh3 = [item[2] for item in param3]
d3  = [item[3] for item in param3]; g3  = [item[4] for item in param3]
p3  = [item[5] for item in param3]; s3  = [item[6] for item in param3]

# For column 4:
b4  = [item[1] for item in param4]; bh4 = [item[2] for item in param4]
d4  = [item[3] for item in param4]; g4  = [item[4] for item in param4]
p4  = [item[5] for item in param4]; s4  = [item[6] for item in param4]
a4 =  [item[7] for item in param4]
	
end;

# ╔═╡ 8617c2c0-55d2-11ec-007f-315cd7023c9f
# 4. Print and Save Estimation Results

# From the initial 100 subjects, the authors drop those subjects whose estimates their Stata algorithm was unable to find, took longer than 200 iterations to converge, or were outliers (identified by a grubbs test). Since we are using a different programming language and a different minimization algorithm (which, in our case, is gradient-free), we managed to find estimates for all 100 individuals. Below, we show results for both the subset of individuals considered by the authors (which we identified by running their "03MergeIndMLEAndConstructMainSample.do" file) and the whole sample.

# First, we show results for the subset of individuals considered by the authors.

# ╔═╡ 861774a0-55d2-11ec-23d0-2feb0b638bc2
begin

# Load the dataset containing the wid of the individuals to keep in the different specifications. We created this csv file by running their "03MergeIndMLEAndConstructMainSample" 

# 4 columns, each containing individuals to keep in each specification
ind_keep = readdlm(raw"../input/ind_to_keep.csv",',',header=true)

# keep parameters only for individuals in ind_keep. First create a dataframe containing the estimates for column 1, 2, 3, 4 and then remove individuals to drop

# For column 1:
mc1 = DataFrame()
mc1.wid=unique(dt.wid); mc1.b=b1; mc1.bh=bh1; mc1.d=d1; mc1.g=g1;
dc1 = mc1[[i for i=1:length(mc1.wid) if mc1.wid[i] in(ind_keep[1][:,1])],:]

# For column 2:
mc2 = DataFrame()
mc2.wid=unique(dt.wid); mc2.b=b2; mc2.bh=bh2; mc2.d=d2; mc2.g=g2;
dc2 = mc2[[i for i=1:length(mc2.wid) if mc2.wid[i] in(ind_keep[1][:,2])],:]

# For column 3:
mc3 = DataFrame()
mc3.wid=unique(dt.wid); mc3.b=b3; mc3.bh=bh3; mc3.d=d3; mc3.g=g3;
dc3 = mc3[[i for i=1:length(mc3.wid) if mc3.wid[i] in(ind_keep[1][:,3])],:]

# For column 4:
mc4 = DataFrame()
mc4.wid=unique(dt.wid); mc4.b=b4; mc4.bh=bh4; mc4.d=d4; mc4.g=g4; mc4.a=a4
dc4 = mc4[[i for i=1:length(mc4.wid) if mc4.wid[i] in(ind_keep[1][:,4])],:]

end;

# ╔═╡ a49dbdd0-55dc-11ec-1e48-21d25efc1b43
begin

# Print the results using the authors' subset of individuals

pnames1 = ["mean(β)", "median(β)", "sd(β)", "mean(β_h)", "median(β_h)", "sd(β_h)",               "mean(δ)", "median(δ)", "sd(δ)", "mean(γ)", "median(γ)", "sd(γ)","mean(α)",           "median(α)", "sd(α)", "P[β]<1", "P[β_h]<1", "r(β,β_h)","p-value r(β,β_h)",
          "Observations"]

at1 = DataFrame()
at1.parameters = pnames1
at1.Primary_Estimation = round.([mean(dc1.b);median(dc1.b);std(dc1.b);mean(dc1.bh);median(dc1.bh);std(dc1.bh);mean(dc1.d);median(dc1.d);std(dc1.d);mean(dc1.g);median(dc1.g);std(dc1.g);0.0;0.0;0.0;mean(dc1.b.<1);mean(dc1.bh.<1);cor(dc1.b,dc1.bh);pvalue(CorrelationTest(dc1.b,dc1.bh));length(dc1.wid)],digits=3)
at1
at1.Early_Decisions = round.([mean(dc2.b);median(dc2.b);std(dc2.b);mean(dc2.bh);median(dc2.bh);std(dc2.bh);mean(dc2.d);median(dc2.d);std(dc2.d);mean(dc2.g);median(dc2.g);std(dc2.g);0.0;0.0;0.0;mean(dc2.b.<1);mean(dc2.bh.<1);cor(dc2.b,dc2.bh);pvalue(CorrelationTest(dc2.b,dc2.bh));length(dc2.wid)],digits=3)
at1
at1.Later_Decisions = round.([mean(dc3.b);median(dc3.b);std(dc3.b);mean(dc3.bh);median(dc3.bh);std(dc3.bh);mean(dc3.d);median(dc3.d);std(dc3.d);mean(dc3.g);median(dc3.g);std(dc3.g);0.0;0.0;0.0;mean(dc3.b.<1);mean(dc3.bh.<1);cor(dc3.b,dc3.bh);pvalue(CorrelationTest(dc3.b,dc3.bh));length(dc3.wid)],digits=3)
at1
at1.Proj_Bias = round.([mean(dc4.b);median(dc4.b);std(dc4.b);mean(dc4.bh);median(dc4.bh);std(dc4.bh);mean(dc4.d);median(dc4.d);std(dc4.d);mean(dc4.g);median(dc4.g);std(dc4.g);mean(dc4.a);median(dc4.a);std(dc4.a);mean(dc4.b.<1);mean(dc4.bh.<1);cor(dc4.b,dc4.bh);pvalue(CorrelationTest(dc4.b,dc4.bh));length(dc4.wid)],digits=3)
at1
	
end

# ╔═╡ 9087d8e0-55e0-11ec-35df-21166bf459e5
begin

# We now look at our estimates. We only remove outliers using a grubbs test with alpha at 1% since we found estimates for all individuals. There is no package in Julia that computes the grubbs test for us, so we need to write a function ourselves:

function grubbs_test(y::Vector{Float64}, alpha)
	
	# compute test statistics 
	x = copy(y)
	avg_x = mean(x)
	abs_val_minus_avg = abs.(x .- avg_x)
	max_dev = maximum(abs_val_minus_avg)
	std_x = std(x)
	G = max_dev / std_x
	index = findall(x->x==max_dev, abs_val_minus_avg)
	
	# compute critical value
	
	t_dist = quantile.(TDist(length(x)-2), alpha / (2*length(x)))
	numerator = (length(x)-1) * sqrt(t_dist^2)
	denominator = sqrt(length(x)) * sqrt(length(x)-2+t_dist^2)
	critical_value = numerator / denominator
	
	if G > critical_value deleteat!(x,index) end   # remove if outlier
	
	return x
	
end

# We now write a function for the Generalized ESD test for outliers. We basically reiterate the function grubbs_test until no outliers are detected

function esd_test(x::Vector{Float64}, alpha)

	diff = 1
	
	while diff != 0	
		oldlength = length(x)
		x = grubbs_test(x,alpha)
		newlength = length(x)
		diff = oldlength - newlength	
	end
	
	return x
	
end
	
end

# ╔═╡ 08059440-55f3-11ec-0d98-438f14370058
begin

# Apply ESD test

beta1 = esd_test(b1, 0.01);  betahat1 = esd_test(bh1, 0.01)
delta1 = esd_test(d1, 0.01); gamma1 = esd_test(g1, 0.01)
beta2 = esd_test(b2, 0.01);  betahat2 = esd_test(bh2, 0.01)
delta2 = esd_test(d2, 0.01); gamma2 = esd_test(g2, 0.01)
beta3 = esd_test(b3, 0.01);  betahat3 = esd_test(bh3, 0.01)
delta3 = esd_test(d3, 0.01); gamma3 = esd_test(g3, 0.01)
beta4 = esd_test(b4, 0.01);  betahat4 = esd_test(bh4, 0.01)
delta4 = esd_test(d4, 0.01); gamma4 = esd_test(g4, 0.01)
alpha4 = esd_test(a4, 0.01)

# we keep only those individuals for which none of the parameters are removed by the grubbs test

# For column 1: 

mdc1 = mc1[[i for i=1:length(mc1.wid) if mc1.b[i] in(beta1)],:]
mdc1 = mdc1[[i for i=1:length(mdc1.wid) if mdc1.bh[i] in(betahat1)],:]
mdc1 = mdc1[[i for i=1:length(mdc1.wid) if mdc1.d[i] in(delta1)],:]
mdc1 = mdc1[[i for i=1:length(mdc1.wid) if mdc1.g[i] in(gamma1)],:]

# For column 2:

mdc2 = mc2[[i for i=1:length(mc2.wid) if mc2.b[i] in(beta2)],:]
mdc2 = mdc2[[i for i=1:length(mdc2.wid) if mdc2.bh[i] in(betahat2)],:]
mdc2 = mdc2[[i for i=1:length(mdc2.wid) if mdc2.d[i] in(delta2)],:]
mdc2 = mdc2[[i for i=1:length(mdc2.wid) if mdc2.g[i] in(gamma2)],:]

# For column 3:

mdc3 = mc3[[i for i=1:length(mc3.wid) if mc3.b[i] in(beta3)],:]
mdc3 = mdc3[[i for i=1:length(mdc3.wid) if mdc3.bh[i] in(betahat3)],:]
mdc3 = mdc3[[i for i=1:length(mdc3.wid) if mdc3.d[i] in(delta3)],:]
mdc3 = mdc3[[i for i=1:length(mdc3.wid) if mdc3.g[i] in(gamma3)],:]

# For column 4:

mdc4 = mc4[[i for i=1:length(mc4.wid) if mc4.b[i] in(beta4)],:]
mdc4 = mdc4[[i for i=1:length(mdc4.wid) if mdc4.bh[i] in(betahat4)],:]
mdc4 = mdc4[[i for i=1:length(mdc4.wid) if mdc4.d[i] in(delta4)],:]
mdc4 = mdc4[[i for i=1:length(mdc4.wid) if mdc4.g[i] in(gamma4)],:]
mdc4 = mdc4[[i for i=1:length(mdc4.wid) if mdc4.a[i] in(alpha4)],:]
	
end;

# ╔═╡ c4bcd1a0-55f5-11ec-357b-87eef6d8494e
begin

# Print our results

pnames2 = ["mean(β)", "median(β)", "sd(β)", "mean(β_h)", "median(β_h)", "sd(β_h)",               "mean(δ)", "median(δ)", "sd(δ)", "mean(γ)", "median(γ)", "sd(γ)","mean(α)",           "median(α)", "sd(α)", "P[β]<1", "P[β_h]<1", "r(β,β_h)","p-value r(β,β_h)",
          "Observations"]

at2 = DataFrame()
at2.parameters = pnames2
at2.Primary_Estimation = round.([mean(mdc1.b);median(mdc1.b);std(mdc1.b);mean(mdc1.bh);median(mdc1.bh);std(mdc1.bh);mean(mdc1.d);median(mdc1.d);std(mdc1.d);mean(mdc1.g);median(mdc1.g);std(mdc1.g);0.0;0.0;0.0;mean(mdc1.b.<1);mean(mdc1.bh.<1);cor(mdc1.b,mdc1.bh);pvalue(CorrelationTest(mdc1.b,mdc1.bh));length(mdc1.wid)],digits=3)
at2.Early_Decisions = round.([mean(mdc2.b);median(mdc2.b);std(mdc2.b);mean(mdc2.bh);median(mdc2.bh);std(mdc2.bh);mean(mdc2.d);median(mdc2.d);std(mdc2.d);mean(mdc2.g);median(mdc2.g);std(mdc2.g);0.0;0.0;0.0;mean(mdc2.b.<1);mean(mdc2.bh.<1);cor(mdc2.b,mdc2.bh);pvalue(CorrelationTest(mdc2.b,mdc2.bh));length(mdc2.wid)],digits=3)
at2.Later_Decisions = round.([mean(mdc3.b);median(mdc3.b);std(mdc3.b);mean(mdc3.bh);median(mdc3.bh);std(mdc3.bh);mean(mdc3.d);median(mdc3.d);std(mdc3.d);mean(mdc3.g);median(mdc3.g);std(mdc3.g);0.0;0.0;0.0;mean(mdc3.b.<1);mean(mdc3.bh.<1);cor(mdc3.b,mdc3.bh);pvalue(CorrelationTest(mdc3.b,mdc3.bh));length(mdc3.wid)],digits=3)
at2.Proj_Bias = round.([mean(mdc4.b);median(mdc4.b);std(mdc4.b);mean(mdc4.bh);median(mdc4.bh);std(mdc4.bh);mean(mdc4.d);median(mdc4.d);std(mdc4.d);mean(mdc4.g);median(mdc4.g);std(mdc4.g);mean(mdc4.a);median(mdc4.a);std(mdc4.a);mean(mdc4.b.<1);mean(mdc4.bh.<1);cor(mdc4.b,mdc4.bh);pvalue(CorrelationTest(mdc4.b,mdc4.bh));length(mdc4.wid)],digits=3)
at2
	
end

# ╔═╡ c072cbde-55f5-11ec-1a99-c17e1667c4e2
begin

# Save our raw results 

table2_julia = DataFrame()
table2_julia.beta1 = b1
table2_julia.betahat1 = bh1
table2_julia.delta1 = d1
table2_julia.gamma1 = g1
table2_julia.phi1 = p1
table2_julia.sigma1 = s1
table2_julia.beta2 = b2
table2_julia.betahat2 = bh2
table2_julia.delta2 = d2
table2_julia.gamma2 = g2
table2_julia.phi2 = p2
table2_julia.sigma2 = s2
table2_julia.beta3 = b3
table2_julia.betahat3 = bh3
table2_julia.delta3 = d3
table2_julia.gamma3 = g3
table2_julia.phi3 = p3
table2_julia.sigma3 = s3
table2_julia.beta4 = b4
table2_julia.betahat4 = bh4
table2_julia.delta4 = d4
table2_julia.gamma4 = g4
table2_julia.phi4 = p4
table2_julia.sigma4 = s4
table2_julia.alpha4 = a4

CSV.write(raw"../output/table2_julia.csv", table2_julia)
	
end

# ╔═╡ fe89b822-268f-11ec-1759-67648ab10a25


# ╔═╡ fe75baf0-268f-11ec-3b55-15fbd2599c0f


# ╔═╡ 4c4dd50e-2757-11ec-3b42-f76583bd9178


# ╔═╡ e5515240-2758-11ec-0aa4-c1dba4b2b512


# ╔═╡ 62c2f7fe-2759-11ec-0cd8-33f81746ff23


# ╔═╡ ca7ac9a2-2759-11ec-1411-afe4da3e4d7d


# ╔═╡ Cell order:
# ╠═d08aaa10-268f-11ec-067f-4322e82ba779
# ╠═fee9d980-268f-11ec-3901-3705290b7c24
# ╠═f8c1b410-51b4-11ec-1d65-11df2949bf23
# ╠═0de27c80-51b5-11ec-05e9-a5b37469b0c7
# ╠═feb1b280-268f-11ec-1574-e5037552b641
# ╠═fe9ddc60-268f-11ec-112a-59ea860b2c81
# ╠═74270f10-51b5-11ec-2edd-af650c9bec6b
# ╠═ce679cb0-51b5-11ec-39ba-3f65d7b1ce33
# ╠═a5d0925e-5535-11ec-06ca-85bfd5b5081b
# ╠═f2177020-5536-11ec-0db6-1b2a59aa1c0e
# ╠═99348ad0-5543-11ec-274d-1d59ea3e02b0
# ╠═a49523d0-5543-11ec-3eec-fbe419cf66af
# ╠═a3f7d260-5543-11ec-3bd8-3d1836a9524e
# ╠═7ccc8cbe-5544-11ec-3b68-31c23aaacd87
# ╠═8617c2c0-55d2-11ec-007f-315cd7023c9f
# ╠═861774a0-55d2-11ec-23d0-2feb0b638bc2
# ╠═a49dbdd0-55dc-11ec-1e48-21d25efc1b43
# ╠═9087d8e0-55e0-11ec-35df-21166bf459e5
# ╠═08059440-55f3-11ec-0d98-438f14370058
# ╠═c4bcd1a0-55f5-11ec-357b-87eef6d8494e
# ╠═c072cbde-55f5-11ec-1a99-c17e1667c4e2
# ╟─fe89b822-268f-11ec-1759-67648ab10a25
# ╟─fe75baf0-268f-11ec-3b55-15fbd2599c0f
# ╟─4c4dd50e-2757-11ec-3b42-f76583bd9178
# ╟─e5515240-2758-11ec-0aa4-c1dba4b2b512
# ╟─62c2f7fe-2759-11ec-0cd8-33f81746ff23
# ╟─ca7ac9a2-2759-11ec-1411-afe4da3e4d7d
