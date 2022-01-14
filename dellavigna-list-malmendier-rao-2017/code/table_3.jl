### A Pluto.jl notebook ###
# v0.14.0

using Markdown
using InteractiveUtils

# ╔═╡ 15e16bb0-3560-11ec-238e-e3018fb17cf9
begin

# Import packages used

using Distributed
using DelimitedFiles
using DataFrames
using Optim
using CSV
using LinearAlgebra
using Statistics
using Distributions
using MAT
using ForwardDiff
using SparseArrays
	
end

# ╔═╡ 9c0d53a0-38c3-11ec-0d4d-1d549264f0f0
# Replication of DellaVigna, List, Malmendier and Rao, 2017, "Voting to Tell Others"
# The code in this notebook performs the benchmark estimates with heterogeneous auxiliary parameters (Table 3, column 1)

# ╔═╡ 92cfa5a2-5ce9-11ec-0e3e-61fe154ce557
# Authors: 
# Massimiliano Pozzi (Bocconi University, pozzi.massimiliano@studbocconi.it)
# Salvatore Nunnari (Bocconi University, https://snunnari.github.io, salvatore.nunnari@unibocconi.it)

# ╔═╡ b34d84a0-5ce9-11ec-1316-39f597588bef
# This notebook works with both Neptune and Pluto; the "begin" and "end" tags in each cell are necessary for Pluto notebooks but do not play any role in Neptune 

# ╔═╡ 73767a00-5db9-11ec-1bef-b9e6f6e711cc
# This notebook was tested with the following packages versions:

# Pozzi: julia 1.6.1, Pluto 0.14.7, Neptune 0.14.0, Distributions 0.25.17, DataFrames 1.2.2, Optim 1.4.1, CSV 0.9.5, MAT 0.10.1, ForwardDiff 0.10.20

# ╔═╡ c6a68930-3588-11ec-08e1-7d0922b9816b
begin

# 1. Data Cleaning and Data Preparation	
	
# We import the relevant 100 empirical moments and variance-covariance matrix. In addition to those moments we manually add the baseline turnout and its standard error. For the benchmark specification we then use 101 moments.

file = matopen(raw"../input/Moments.mat")
emp_moments = vec(read(file, "Moments"))
emp_moments = [emp_moments; 0.6000] # empirical moments. We manually add baseline                                           turnout
	
# We now import the variance-covariance matrix of the empirical moments. sparse to use blockdiag in the next lines

VCcontrol = sparse(read(file, "VCcontrol")) 
	
# We need to add an elements on the diagonal of VCcontrol. Matrix to go back from sparse
	
VCcontrol = Matrix(blockdiag(VCcontrol,sparse(diagm([0.0109^2])))) 

# W is a diagonal weighting matrix. This is the inverse of the diagonal of VCcontrol. We use 101 moments in the benchmark estimates, so this is a 101x101 diagonal matrix

W = inv(Diagonal(Diagonal(VCcontrol)))
	
end;

# ╔═╡ c6731d20-3588-11ec-3127-51520038c40c
# 2. Briefly define the model and the estimation strategy

# The estimation method is simulated minimum distance, where we impose simulated moments being equal to the empirical moments observed in the data. These simulated moments come from the fact that some of the parameters in the model are assumed to be stochastic and heterogeneous. What we will be minimizing is the weighted sum of squared distances between simulated and empirical moments.

# For more information about the model please refer to the python notebook or section 2 and 5 in the paper

# This function computes the relevant 101 moments by simulating N individuals. 
# rand_set is a vector of vectors of random draws for the simulated parameters

function voteSimEndogenousVoting_vary(parameters, rand_set)

# Parameters: 
# h0  = baseline probability of being at home
# r   = prob of seeing flyer
# eta = elasticity of response to sorting in and out of the house
# s   = how much you like doing a survey
# S   = social pressure cost of doing a survey
# sv  = value of saying you voted (appearing as a voter)
# sn  = value of saying you didn't vote (appearing as a non-voter)
# rho = correlation between sv and sn (assumed to be 0 in the benchmark estimates)
# eps = all other reasons to vote
	
N=5.4    # Number of times asked if you voted in the Congressional elections
N_P=10.1 # Times asked in presidential
	
# Auxiliary parameters. These parameters vary between voters and non-voters. sub v for voters and sub nv for non-voters
	
h0_v	=	parameters[1]
h0_nv	=	parameters[2]
r_v	    =	parameters[3]
r_nv	=	parameters[4]
eta_v	=	parameters[5]
eta_nv	=	parameters[6]
mu_s_v	=	parameters[7]
mu_s_nv	=	parameters[8]
sigma_s_v	=	parameters[9]
sigma_s_nv	=	parameters[10]
S_svy_v	    =	parameters[11]
S_svy_nv	=	parameters[12]
timeval_v	=	parameters[13]
timeval_nv	=	parameters[14]

# Other parameters. Here they don't vary between voters and non-voters
mu_sv	=	parameters[15]
mu_sn	=	parameters[16]
sigma_svn = parameters[17]  # sigma_sv = sigma_sn in benchmark estimates 
L = parameters[18]
mu_eps = parameters[19]
sigma_eps = parameters[20] 
rho	= 0                     # Set to zero in benchmark
	
# Draws from a standard normal. The dimension of these vectors is equal to the number of simulated individuals
	
rand_set1 = rand_set[1] # used to compute s
rand_set2 = rand_set[2] # used to compute s_v
rand_set3 = rand_set[3] # used to compute s_n
rand_set4 = rand_set[4] # used to compute eps

# Simulate epsilons ("other" reasons to vote)
	
eps = mu_eps .+ sigma_eps.*rand_set2

# Simulate sv and sn
	
sv = mu_sv .+ sigma_svn.*rand_set3
sn = mu_sn .+ sigma_svn.*rand_set4

	
# !!!Vote or not!!! Does our simulated individual vote or not?
	
	
#Net expected utility of voting (relative to not voting)
	
sigVal = (max.(sv,sn.-L) .- max.(sn,sv.-L))# net sig utility voter - non-voter (1 ask)
sigVal_x_N = N.*sigVal                     # asked N times
utilityVoting = sigVal_x_N .+ eps          # net utility
voted = utilityVoting .>0                  # Vote if net utility > 0

# Return a list of indices (e.g., 1, 3, 4) identifying voters and non-voters
	
voterIndex    = findall(==(1), voted)
nonvoterIndex = findall(==(0), voted)

# Share who turn out in the control group (no GOTV intervention)
	
Turnout_control = mean(voted) 

# Make vectors with voter and non-voter parameters based on whether voted in the control experiment
	
h0 = voted.*h0_v .+ (1 .-voted).*h0_nv
r = voted.*r_v .+ (1 .-voted).*r_nv
eta = voted.*eta_v .+ (1 .-voted).*eta_nv
S_svy = voted.*S_svy_v .+ (1 .-voted).*S_svy_nv
timeval = voted.*timeval_v .+ (1 .-voted).*timeval_nv

# GOTV (use parameters as assigned by turnout in the control). Assumes intervention = expect to be asked h0 times more

# net utility of voting if seen the GOTV flyer
sigVal_x_N_GOTV = (N.+h0).*sigVal
utilityVoting_GOTV = sigVal_x_N_GOTV .+ eps 
voted_GOTV = utilityVoting_GOTV .>0

# Share who turn out with the GOTV intervention assume everyone sees the flyer, counts as N+h0 times asked
Turnout_GOTV = mean(voted_GOTV)

# PRESIDENTIAL
sigVal_x_N_P = N_P.*sigVal; # asked N_P times
sigVal_x_N_P_GOTV = (N_P.+h0).*sigVal # asked N_P times

utilityVoting_P = sigVal_x_N_P .+ eps
utilityVoting_P_GOTV = sigVal_x_N_P_GOTV .+ eps 

voted_P = utilityVoting_P .>0 
voted_P_GOTV = utilityVoting_P_GOTV .>0 

# Turnout 
Turnout_P_control = mean(voted_P)
Turnout_P_GOTV = mean(voted_P_GOTV)

# Draw random s, Utilities of doing the survey for voters/non-voters 

# Simulate utility of doing a 0d10m survey 
s = voted.*mu_s_v .+ voted.*sigma_s_v.*rand_set1 .+ (1 .-voted).*mu_s_nv .+ (1 .-voted).*sigma_s_nv.*rand_set1

# Values of survey incentives (relative to 0d10m) XdYm = X dollars and Y min
D_0d10m_v = 0
D_0d5m_v = timeval_v*5/60
D_10d5m_v = 10+timeval_v*5/60
D_10d10m_v = 10

D_0d10m_nv = 0
D_0d5m_nv = timeval_nv*5/60
D_10d5m_nv = 10+timeval_nv*5/60
D_10d10m_nv = 10

# Extra incentive if say "not vote". 5m survey: +1m + $5. 10m survey: -8m
I_5d1m_v = 5 .-timeval_v*1/60
I_8m_v = timeval_v*8/60

I_5d1m_nv = 5 .-timeval_nv*1/60
I_8m_nv = timeval_nv*8/60

# Lying if asked

# UtilVotingQuestion = utility you get from being asked one time (max of lie or not lie)
wouldLieIfAsked = voted.*(sn.-L .>sv) .+ (1 .-voted).*(sv.-L .>sn)
utilVotingQuestion = voted.*max.(sn.-L,sv) .+ (1 .-voted).*max.(sv.-L,sn)

# Response to incentives to say "not vote"
wouldLieIfAsked_5d1m = voted.*(sn.-L.+I_5d1m_v .>sv) .+ (1 .-voted).*(sv.-L.>sn.+I_5d1m_nv)
wouldLieIfAsked_8m = voted.*(sn.-L.+I_8m_v .>sv) .+ (1 .-voted).*(sv.-L .>sn.+I_8m_nv)

	
# !!!!! Compute Moments !!!!!

	
# Utility from doing survey

# NF= no flyer, F = flyer, FV= flyer + voting, OO = opt-out, OOV = opt-out + voting

# Anticipated utility from doing survey and voting survey (VF or I)
util_svyOnly_0d5m = s .+ voted.*D_0d5m_v .+ (1 .-voted).*D_0d5m_nv
util_svyPlusVotingQues_0d5m = util_svyOnly_0d5m .+ utilVotingQuestion

util_svyOnly_10d10m = s .+ voted.*D_10d10m_v .+ (1 .-voted).*D_10d10m_nv
util_svyPlusVotingQues_10d10m = util_svyOnly_10d10m .+ utilVotingQuestion

util_svyOnly_10d5m = s .+ voted.*D_10d5m_v .+ (1 .-voted).*D_10d5m_nv
util_svyPlusVotingQues_10d5m = util_svyOnly_10d5m .+ utilVotingQuestion

# If asked, do survey if greater than the social pressure cost. NI = not informed that survey is about voting, I=informed (VF or I)
doesSvyIfAsked_NI_0d5m = util_svyOnly_0d5m .> -S_svy
doesSvyIfAsked_I_0d5m = util_svyPlusVotingQues_0d5m .> -S_svy

doesSvyIfAsked_NI_10d10m = util_svyOnly_10d10m .> -S_svy
doesSvyIfAsked_I_10d10m = util_svyPlusVotingQues_10d10m .> -S_svy

doesSvyIfAsked_NI_10d5m = util_svyOnly_10d5m .> -S_svy
doesSvyIfAsked_I_10d5m = util_svyPlusVotingQues_10d5m .> -S_svy


# Anticipated utility given you are asked to do the survey
anticipatedUtil_Svy_NI_0d5m = max.(util_svyOnly_0d5m,-S_svy) 
anticipatedUtil_Svy_I_0d5m = max.(util_svyPlusVotingQues_0d5m,-S_svy)

anticipatedUtil_Svy_NI_10d10m = max.(util_svyOnly_10d10m,-S_svy)
anticipatedUtil_Svy_I_10d10m = max.(util_svyPlusVotingQues_10d10m,-S_svy)

anticipatedUtil_Svy_NI_10d5m = max.(util_svyOnly_10d5m,-S_svy)
anticipatedUtil_Svy_I_10d5m = max.(util_svyPlusVotingQues_10d5m,-S_svy)


# Opt-out if anticipated utility is negative
optsOutIfSees_OO_0d5m  = anticipatedUtil_Svy_NI_0d5m .< 0
optsOutIfSees_OOV_0d5m = anticipatedUtil_Svy_I_0d5m .< 0

optsOutIfSees_OO_10d10m  = anticipatedUtil_Svy_NI_10d10m .< 0
optsOutIfSees_OOV_10d10m = anticipatedUtil_Svy_I_10d10m .< 0

optsOutIfSees_OO_10d5m  = anticipatedUtil_Svy_NI_10d5m .< 0
optsOutIfSees_OOV_10d5m = anticipatedUtil_Svy_I_10d5m .< 0


# Choosing probability of being at home is bounded between 0 and 1

hStar_F_0d5m = max.(0,min.(1,h0.+eta.*anticipatedUtil_Svy_NI_0d5m))
hStar_FV_0d5m = max.(0,min.(1,h0.+eta.*anticipatedUtil_Svy_I_0d5m))

hStar_F_10d10m = max.(0,min.(1,h0.+eta.*anticipatedUtil_Svy_NI_10d10m))
hStar_FV_10d10m = max.(0,min.(1,h0.+eta.*anticipatedUtil_Svy_I_10d10m))

hStar_F_10d5m = max.(0,min.(1,h0+eta.*anticipatedUtil_Svy_NI_10d5m))
hStar_FV_10d5m = max.(0,min.(1,h0+eta.*anticipatedUtil_Svy_I_10d5m))

# Separate Voters and Nonvoters

# Split into separate vectors of voters and non-voter vectors note: they will be of different length

# Voters
hStar_F_0d5m_v = hStar_F_0d5m[voterIndex]
hStar_FV_0d5m_v = hStar_FV_0d5m[voterIndex]
doesSvyIfAsked_NI_0d5m_v = doesSvyIfAsked_NI_0d5m[voterIndex]
doesSvyIfAsked_I_0d5m_v = doesSvyIfAsked_I_0d5m[voterIndex]
optsOutIfSees_OO_0d5m_v=optsOutIfSees_OO_0d5m[voterIndex]
optsOutIfSees_OOV_0d5m_v=optsOutIfSees_OOV_0d5m[voterIndex]

hStar_F_10d10m_v = hStar_F_10d10m[voterIndex]
hStar_FV_10d10m_v = hStar_FV_10d10m[voterIndex]
doesSvyIfAsked_NI_10d10m_v = doesSvyIfAsked_NI_10d10m[voterIndex]
doesSvyIfAsked_I_10d10m_v = doesSvyIfAsked_I_10d10m[voterIndex]
optsOutIfSees_OO_10d10m_v=optsOutIfSees_OO_10d10m[voterIndex]
optsOutIfSees_OOV_10d10m_v=optsOutIfSees_OOV_10d10m[voterIndex]

hStar_F_10d5m_v = hStar_F_10d5m[voterIndex]
hStar_FV_10d5m_v = hStar_FV_10d5m[voterIndex]
doesSvyIfAsked_NI_10d5m_v = doesSvyIfAsked_NI_10d5m[voterIndex]
doesSvyIfAsked_I_10d5m_v = doesSvyIfAsked_I_10d5m[voterIndex]
optsOutIfSees_OO_10d5m_v=optsOutIfSees_OO_10d5m[voterIndex]
optsOutIfSees_OOV_10d5m_v=optsOutIfSees_OOV_10d5m[voterIndex]

wouldLieIfAsked_v = wouldLieIfAsked[voterIndex]
wouldLieIfAsked_5d1m_v = wouldLieIfAsked_5d1m[voterIndex]
wouldLieIfAsked_8m_v=wouldLieIfAsked_8m[voterIndex]

# Non-voters
hStar_F_0d5m_nv = hStar_F_0d5m[nonvoterIndex]
hStar_FV_0d5m_nv = hStar_FV_0d5m[nonvoterIndex]
doesSvyIfAsked_NI_0d5m_nv = doesSvyIfAsked_NI_0d5m[nonvoterIndex]
doesSvyIfAsked_I_0d5m_nv = doesSvyIfAsked_I_0d5m[nonvoterIndex]
optsOutIfSees_OO_0d5m_nv=optsOutIfSees_OO_0d5m[nonvoterIndex]
optsOutIfSees_OOV_0d5m_nv=optsOutIfSees_OOV_0d5m[nonvoterIndex]

hStar_F_10d10m_nv = hStar_F_10d10m[nonvoterIndex]
hStar_FV_10d10m_nv = hStar_FV_10d10m[nonvoterIndex]
doesSvyIfAsked_NI_10d10m_nv = doesSvyIfAsked_NI_10d10m[nonvoterIndex]
doesSvyIfAsked_I_10d10m_nv = doesSvyIfAsked_I_10d10m[nonvoterIndex]
optsOutIfSees_OO_10d10m_nv=optsOutIfSees_OO_10d10m[nonvoterIndex]
optsOutIfSees_OOV_10d10m_nv=optsOutIfSees_OOV_10d10m[nonvoterIndex]

hStar_F_10d5m_nv = hStar_F_10d5m[nonvoterIndex]
hStar_FV_10d5m_nv = hStar_FV_10d5m[nonvoterIndex]
doesSvyIfAsked_NI_10d5m_nv = doesSvyIfAsked_NI_10d5m[nonvoterIndex]
doesSvyIfAsked_I_10d5m_nv = doesSvyIfAsked_I_10d5m[nonvoterIndex]
optsOutIfSees_OO_10d5m_nv=optsOutIfSees_OO_10d5m[nonvoterIndex]
optsOutIfSees_OOV_10d5m_nv=optsOutIfSees_OOV_10d5m[nonvoterIndex]

wouldLieIfAsked_nv = wouldLieIfAsked[nonvoterIndex]
wouldLieIfAsked_5d1m_nv = wouldLieIfAsked_5d1m[nonvoterIndex]
wouldLieIfAsked_8m_nv=wouldLieIfAsked_8m[nonvoterIndex]

	
# !!!!! Disaggregated Moments !!!!!

	
# Voters

# PH = probability of being at home
PH_NF_0d5m_v = h0_v;
PH_F_0d5m_v = (1 .-r_v)*h0_v .+ r_v.*mean(hStar_F_0d5m_v);
PH_FV_0d5m_v = (1 .-r_v)*h0_v .+ r_v.*mean(hStar_FV_0d5m_v);
PH_OO_0d5m_v = (1 .-r_v)*h0_v .+ r_v.*mean((1 .-optsOutIfSees_OO_0d5m_v).*hStar_F_0d5m_v); 
PH_OOV_0d5m_v = (1 .-r_v)*h0_v .+ r_v.*mean((1 .-optsOutIfSees_OOV_0d5m_v).*hStar_FV_0d5m_v);

PH_NF_10d10m_v = h0_v;
PH_F_10d10m_v = (1 .-r_v)*h0_v .+ r_v.*mean(hStar_F_10d10m_v);
PH_FV_10d10m_v = (1 .-r_v)*h0_v .+ r_v.*mean(hStar_FV_10d10m_v);
PH_OO_10d10m_v = (1 .-r_v)*h0_v .+ r_v.*mean((1 .-optsOutIfSees_OO_10d10m_v).*hStar_F_10d10m_v); 
PH_OOV_10d10m_v = (1 .-r_v)*h0_v .+ r_v.*mean((1 .-optsOutIfSees_OOV_10d10m_v).*hStar_FV_10d10m_v);

PH_NF_10d5m_v = h0_v;
PH_F_10d5m_v = (1 .-r_v)*h0_v .+ r_v.*mean(hStar_F_10d5m_v);
PH_FV_10d5m_v = (1 .-r_v)*h0_v .+ r_v.*mean(hStar_FV_10d5m_v);
PH_OO_10d5m_v = (1 .-r_v)*h0_v .+ r_v.*mean((1 .-optsOutIfSees_OO_10d5m_v).*hStar_F_10d5m_v); 
PH_OOV_10d5m_v = (1 .-r_v)*h0_v .+ r_v.*mean((1 .-optsOutIfSees_OOV_10d5m_v).*hStar_FV_10d5m_v);


# PSV=unconditional prob of doing the survey (not cond on opening door). PSV < PH mechanically

# 0d5m
PSV_NF_NI_0d5m_v = h0_v*mean(doesSvyIfAsked_NI_0d5m_v);
PSV_NF_I_0d5m_v = h0_v*mean(doesSvyIfAsked_I_0d5m_v);

PSV_F_NI_0d5m_v = (1 .-r_v).*PSV_NF_NI_0d5m_v .+ r_v.*mean(hStar_F_0d5m_v.*doesSvyIfAsked_NI_0d5m_v);
PSV_F_I_0d5m_v = (1 .-r_v)*PSV_NF_I_0d5m_v .+ r_v.*mean(hStar_F_0d5m_v.*doesSvyIfAsked_I_0d5m_v);

PSV_FV_NI_0d5m_v = (1 .-r_v)*PSV_NF_NI_0d5m_v .+ r_v.*mean(hStar_FV_0d5m_v.*doesSvyIfAsked_I_0d5m_v);
PSV_FV_I_0d5m_v = (1 .-r_v)*PSV_NF_I_0d5m_v .+ r_v.*mean(hStar_FV_0d5m_v.*doesSvyIfAsked_I_0d5m_v);

PSV_OO_NI_0d5m_v = (1 .-r_v)*PSV_NF_NI_0d5m_v .+ r_v*mean((1 .-optsOutIfSees_OO_0d5m_v).*hStar_F_0d5m_v.*doesSvyIfAsked_NI_0d5m_v);
PSV_OO_I_0d5m_v = (1 .-r_v)*PSV_NF_I_0d5m_v .+ r_v*mean((1 .-optsOutIfSees_OO_0d5m_v).*hStar_F_0d5m_v.*doesSvyIfAsked_I_0d5m_v);

PSV_OOV_NI_0d5m_v = (1 .-r_v)*PSV_NF_NI_0d5m_v .+ r_v*mean((1 .-optsOutIfSees_OOV_0d5m_v).*hStar_FV_0d5m_v.*doesSvyIfAsked_I_0d5m_v);
PSV_OOV_I_0d5m_v = (1 .-r_v)*PSV_NF_I_0d5m_v .+ r_v*mean((1 .-optsOutIfSees_OOV_0d5m_v).*hStar_FV_0d5m_v.*doesSvyIfAsked_I_0d5m_v);

# 10d10m
PSV_NF_NI_10d10m_v = h0_v*mean(doesSvyIfAsked_NI_10d10m_v);
PSV_NF_I_10d10m_v = h0_v*mean(doesSvyIfAsked_I_10d10m_v);
 
PSV_F_NI_10d10m_v = (1-r_v)*PSV_NF_NI_10d10m_v + r_v*mean(hStar_F_10d10m_v.*doesSvyIfAsked_NI_10d10m_v);
PSV_F_I_10d10m_v = (1-r_v)*PSV_NF_I_10d10m_v + r_v*mean(hStar_F_10d10m_v.*doesSvyIfAsked_I_10d10m_v);
 
PSV_FV_NI_10d10m_v = (1-r_v)*PSV_NF_NI_10d10m_v + r_v*mean(hStar_FV_10d10m_v.*doesSvyIfAsked_I_10d10m_v);
PSV_FV_I_10d10m_v = (1-r_v)*PSV_NF_I_10d10m_v + r_v*mean(hStar_FV_10d10m_v.*doesSvyIfAsked_I_10d10m_v);
 
PSV_OO_NI_10d10m_v = (1-r_v)*PSV_NF_NI_10d10m_v + r_v*mean((1 .-optsOutIfSees_OO_10d10m_v).*hStar_F_10d10m_v.*doesSvyIfAsked_NI_10d10m_v);
PSV_OO_I_10d10m_v = (1-r_v)*PSV_NF_I_10d10m_v + r_v*mean((1 .-optsOutIfSees_OO_10d10m_v).*hStar_F_10d10m_v.*doesSvyIfAsked_I_10d10m_v);
 
PSV_OOV_NI_10d10m_v = (1-r_v)*PSV_NF_NI_10d10m_v + r_v*mean((1 .-optsOutIfSees_OOV_10d10m_v).*hStar_FV_10d10m_v.*doesSvyIfAsked_I_10d10m_v);
PSV_OOV_I_10d10m_v = (1-r_v)*PSV_NF_I_10d10m_v + r_v*mean((1 .-optsOutIfSees_OOV_10d10m_v).*hStar_FV_10d10m_v.*doesSvyIfAsked_I_10d10m_v);

# 10d5m
PSV_NF_NI_10d5m_v = h0_v*mean(doesSvyIfAsked_NI_10d5m_v);
PSV_NF_I_10d5m_v = h0_v*mean(doesSvyIfAsked_I_10d5m_v);
 
PSV_F_NI_10d5m_v = (1-r_v)*PSV_NF_NI_10d5m_v + r_v*mean(hStar_F_10d5m_v.*doesSvyIfAsked_NI_10d5m_v);
PSV_F_I_10d5m_v = (1-r_v)*PSV_NF_I_10d5m_v + r_v*mean(hStar_F_10d5m_v.*doesSvyIfAsked_I_10d5m_v);
 
PSV_FV_NI_10d5m_v = (1-r_v)*PSV_NF_NI_10d5m_v + r_v*mean(hStar_FV_10d5m_v.*doesSvyIfAsked_I_10d5m_v);
PSV_FV_I_10d5m_v = (1-r_v)*PSV_NF_I_10d5m_v + r_v*mean(hStar_FV_10d5m_v.*doesSvyIfAsked_I_10d5m_v);
 
PSV_OO_NI_10d5m_v = (1-r_v)*PSV_NF_NI_10d5m_v + r_v*mean((1 .-optsOutIfSees_OO_10d5m_v).*hStar_F_10d5m_v.*doesSvyIfAsked_NI_10d5m_v);
PSV_OO_I_10d5m_v = (1-r_v)*PSV_NF_I_10d5m_v + r_v*mean((1 .-optsOutIfSees_OO_10d5m_v).*hStar_F_10d5m_v.*doesSvyIfAsked_I_10d5m_v);
 
PSV_OOV_NI_10d5m_v = (1-r_v)*PSV_NF_NI_10d5m_v + r_v*mean((1 .-optsOutIfSees_OOV_10d5m_v).*hStar_FV_10d5m_v.*doesSvyIfAsked_I_10d5m_v);
PSV_OOV_I_10d5m_v = (1-r_v)*PSV_NF_I_10d5m_v + r_v*mean((1 .-optsOutIfSees_OOV_10d5m_v).*hStar_FV_10d5m_v.*doesSvyIfAsked_I_10d5m_v);


# POO=prob of opting out (not conditional on seeing flyer). Scaled by baseline likelihood of being at home
POO_OO_0d5m_v =  h0_v*r_v*mean(optsOutIfSees_OO_0d5m_v);
POO_OOV_0d5m_v = h0_v*r_v*mean(optsOutIfSees_OOV_0d5m_v);

POO_OO_10d10m_v =  h0_v*r_v*mean(optsOutIfSees_OO_10d10m_v);
POO_OOV_10d10m_v = h0_v*r_v*mean(optsOutIfSees_OOV_10d10m_v);

POO_OO_10d5m_v =  h0_v*r_v*mean(optsOutIfSees_OO_10d5m_v);
POO_OOV_10d5m_v = h0_v*r_v*mean(optsOutIfSees_OOV_10d5m_v);


# Empirical moments are total lying in treatments / total doing survey in treatments

# PSVL = unconditional percent who do survey and lie. No flyer treatment only, simplifies later code

# PL=cond on agreeing to do the survey, did you lie? Incentive to lie is a surprise later (doesn't affect PH or PSV)

# 0d5m, 5d1m incentive
PSVL_NF_NI_0d5m_v = mean(h0_v.*doesSvyIfAsked_NI_0d5m_v.*wouldLieIfAsked_v);
PSVL_NF_I_0d5m_v = mean(h0_v.*doesSvyIfAsked_I_0d5m_v.*wouldLieIfAsked_v);
PSVL_NF_NI_0d5m_5d1m_v = mean(h0_v.*doesSvyIfAsked_NI_0d5m_v.*wouldLieIfAsked_5d1m_v);
PSVL_NF_I_0d5m_5d1m_v = mean(h0_v.*doesSvyIfAsked_I_0d5m_v.*wouldLieIfAsked_5d1m_v);

PSVL_F_NI_0d5m_v = (1-r_v)*PSVL_NF_NI_0d5m_v + r_v*mean(hStar_F_0d5m_v.*doesSvyIfAsked_NI_0d5m_v.*wouldLieIfAsked_v);
PSVL_F_I_0d5m_v = (1-r_v)*PSVL_NF_I_0d5m_v + r_v*mean(hStar_F_0d5m_v.*doesSvyIfAsked_I_0d5m_v.*wouldLieIfAsked_v);
PSVL_F_NI_0d5m_5d1m_v = (1-r_v)*PSVL_NF_NI_0d5m_5d1m_v + r_v*mean(hStar_F_0d5m_v.*doesSvyIfAsked_NI_0d5m_v.*wouldLieIfAsked_5d1m_v);
PSVL_F_I_0d5m_5d1m_v = (1-r_v)*PSVL_NF_I_0d5m_5d1m_v + r_v*mean(hStar_F_0d5m_v.*doesSvyIfAsked_I_0d5m_v.*wouldLieIfAsked_5d1m_v);

PSVL_FV_NI_0d5m_v = (1-r_v)*PSVL_NF_NI_0d5m_v + r_v*mean(hStar_FV_0d5m_v.*doesSvyIfAsked_I_0d5m_v.*wouldLieIfAsked_v);
PSVL_FV_I_0d5m_v = (1-r_v)*PSVL_NF_I_0d5m_v + r_v*mean(hStar_FV_0d5m_v.*doesSvyIfAsked_I_0d5m_v.*wouldLieIfAsked_v);
PSVL_FV_NI_0d5m_5d1m_v = (1-r_v)*PSVL_NF_NI_0d5m_5d1m_v + r_v*mean(hStar_FV_0d5m_v.*doesSvyIfAsked_I_0d5m_v.*wouldLieIfAsked_5d1m_v);
PSVL_FV_I_0d5m_5d1m_v = (1-r_v)*PSVL_NF_I_0d5m_5d1m_v + r_v*mean(hStar_FV_0d5m_v.*doesSvyIfAsked_I_0d5m_v.*wouldLieIfAsked_5d1m_v);

PSVL_OO_NI_0d5m_v = (1-r_v)*PSVL_NF_NI_0d5m_v + r_v*mean((1 .-optsOutIfSees_OO_0d5m_v).*hStar_F_0d5m_v.*doesSvyIfAsked_NI_0d5m_v.*wouldLieIfAsked_v);
PSVL_OO_I_0d5m_v = (1-r_v)*PSVL_NF_I_0d5m_v + r_v*mean((1 .-optsOutIfSees_OO_0d5m_v).*hStar_F_0d5m_v.*doesSvyIfAsked_I_0d5m_v.*wouldLieIfAsked_v);
PSVL_OO_NI_0d5m_5d1m_v = (1-r_v)*PSVL_NF_NI_0d5m_5d1m_v + r_v*mean((1 .-optsOutIfSees_OO_0d5m_v).*hStar_F_0d5m_v.*doesSvyIfAsked_NI_0d5m_v.*wouldLieIfAsked_5d1m_v);
PSVL_OO_I_0d5m_5d1m_v = (1-r_v)*PSVL_NF_I_0d5m_5d1m_v + r_v*mean((1 .-optsOutIfSees_OO_0d5m_v).*hStar_F_0d5m_v.*doesSvyIfAsked_I_0d5m_v.*wouldLieIfAsked_5d1m_v);

PSVL_OOV_NI_0d5m_v = (1-r_v)*PSVL_NF_NI_0d5m_v + r_v*mean((1 .-optsOutIfSees_OOV_0d5m_v).*hStar_FV_0d5m_v.*doesSvyIfAsked_I_0d5m_v.*wouldLieIfAsked_v);
PSVL_OOV_I_0d5m_v = (1-r_v)*PSVL_NF_I_0d5m_v + r_v*mean((1 .-optsOutIfSees_OOV_0d5m_v).*hStar_FV_0d5m_v.*doesSvyIfAsked_I_0d5m_v.*wouldLieIfAsked_v);
PSVL_OOV_NI_0d5m_5d1m_v = (1-r_v)*PSVL_NF_NI_0d5m_5d1m_v + r_v*mean((1 .-optsOutIfSees_OOV_0d5m_v).*hStar_FV_0d5m_v.*doesSvyIfAsked_I_0d5m_v.*wouldLieIfAsked_5d1m_v);
PSVL_OOV_I_0d5m_5d1m_v = (1-r_v)*PSVL_NF_I_0d5m_5d1m_v + r_v*mean((1 .-optsOutIfSees_OOV_0d5m_v).*hStar_FV_0d5m_v.*doesSvyIfAsked_I_0d5m_v.*wouldLieIfAsked_5d1m_v);


# 10d10m, 8m incentive
PSVL_NF_NI_10d10m_v = mean(h0_v.*doesSvyIfAsked_NI_10d10m_v.*wouldLieIfAsked_v);
PSVL_NF_I_10d10m_v = mean(h0_v.*doesSvyIfAsked_I_10d10m_v.*wouldLieIfAsked_v);
PSVL_NF_NI_10d10m_8m_v = mean(h0_v.*doesSvyIfAsked_NI_10d10m_v.*wouldLieIfAsked_8m_v);
PSVL_NF_I_10d10m_8m_v = mean(h0_v.*doesSvyIfAsked_I_10d10m_v.*wouldLieIfAsked_8m_v);

PSVL_F_NI_10d10m_v = (1-r_v)*PSVL_NF_NI_10d10m_v + r_v*mean(hStar_F_10d10m_v.*doesSvyIfAsked_NI_10d10m_v.*wouldLieIfAsked_v);
PSVL_F_I_10d10m_v = (1-r_v)*PSVL_NF_I_10d10m_v + r_v*mean(hStar_F_10d10m_v.*doesSvyIfAsked_I_10d10m_v.*wouldLieIfAsked_v);
PSVL_F_NI_10d10m_8m_v = (1-r_v)*PSVL_NF_NI_10d10m_8m_v + r_v*mean(hStar_F_10d10m_v.*doesSvyIfAsked_NI_10d10m_v.*wouldLieIfAsked_8m_v);
PSVL_F_I_10d10m_8m_v = (1-r_v)*PSVL_NF_I_10d10m_8m_v + r_v*mean(hStar_F_10d10m_v.*doesSvyIfAsked_I_10d10m_v.*wouldLieIfAsked_8m_v);
 
PSVL_FV_NI_10d10m_v = (1-r_v)*PSVL_NF_NI_10d10m_v + r_v*mean(hStar_FV_10d10m_v.*doesSvyIfAsked_I_10d10m_v.*wouldLieIfAsked_v);
PSVL_FV_I_10d10m_v = (1-r_v)*PSVL_NF_I_10d10m_v + r_v*mean(hStar_FV_10d10m_v.*doesSvyIfAsked_I_10d10m_v.*wouldLieIfAsked_v);
PSVL_FV_NI_10d10m_8m_v = (1-r_v)*PSVL_NF_NI_10d10m_8m_v + r_v*mean(hStar_FV_10d10m_v.*doesSvyIfAsked_I_10d10m_v.*wouldLieIfAsked_8m_v);
PSVL_FV_I_10d10m_8m_v = (1-r_v)*PSVL_NF_I_10d10m_8m_v + r_v*mean(hStar_FV_10d10m_v.*doesSvyIfAsked_I_10d10m_v.*wouldLieIfAsked_8m_v);
 
PSVL_OO_NI_10d10m_v = (1-r_v)*PSVL_NF_NI_10d10m_v + r_v*mean((1 .-optsOutIfSees_OO_10d10m_v).*hStar_F_10d10m_v.*doesSvyIfAsked_NI_10d10m_v.*wouldLieIfAsked_v);
PSVL_OO_I_10d10m_v = (1-r_v)*PSVL_NF_I_10d10m_v + r_v*mean((1 .-optsOutIfSees_OO_10d10m_v).*hStar_F_10d10m_v.*doesSvyIfAsked_I_10d10m_v.*wouldLieIfAsked_v);
PSVL_OO_NI_10d10m_8m_v = (1-r_v)*PSVL_NF_NI_10d10m_8m_v + r_v*mean((1 .-optsOutIfSees_OO_10d10m_v).*hStar_F_10d10m_v.*doesSvyIfAsked_NI_10d10m_v.*wouldLieIfAsked_8m_v);
PSVL_OO_I_10d10m_8m_v = (1-r_v)*PSVL_NF_I_10d10m_8m_v + r_v*mean((1 .-optsOutIfSees_OO_10d10m_v).*hStar_F_10d10m_v.*doesSvyIfAsked_I_10d10m_v.*wouldLieIfAsked_8m_v);
 
PSVL_OOV_NI_10d10m_v = (1-r_v)*PSVL_NF_NI_10d10m_v + r_v*mean((1 .-optsOutIfSees_OOV_10d10m_v).*hStar_FV_10d10m_v.*doesSvyIfAsked_I_10d10m_v.*wouldLieIfAsked_v);
PSVL_OOV_I_10d10m_v = (1-r_v)*PSVL_NF_I_10d10m_v + r_v*mean((1 .-optsOutIfSees_OOV_10d10m_v).*hStar_FV_10d10m_v.*doesSvyIfAsked_I_10d10m_v.*wouldLieIfAsked_v);
PSVL_OOV_NI_10d10m_8m_v = (1-r_v)*PSVL_NF_NI_10d10m_8m_v + r_v*mean((1 .-optsOutIfSees_OOV_10d10m_v).*hStar_FV_10d10m_v.*doesSvyIfAsked_I_10d10m_v.*wouldLieIfAsked_8m_v);
PSVL_OOV_I_10d10m_8m_v = (1-r_v)*PSVL_NF_I_10d10m_8m_v + r_v*mean((1 .-optsOutIfSees_OOV_10d10m_v).*hStar_FV_10d10m_v.*doesSvyIfAsked_I_10d10m_v.*wouldLieIfAsked_8m_v);

# 10d5m, 5d1m incentive
PSVL_NF_NI_10d5m_v = mean(h0_v.*doesSvyIfAsked_NI_10d5m_v.*wouldLieIfAsked_v);
PSVL_NF_I_10d5m_v = mean(h0_v.*doesSvyIfAsked_I_10d5m_v.*wouldLieIfAsked_v);
PSVL_NF_NI_10d5m_5d1m_v = mean(h0_v.*doesSvyIfAsked_NI_10d5m_v.*wouldLieIfAsked_5d1m_v);
PSVL_NF_I_10d5m_5d1m_v = mean(h0_v.*doesSvyIfAsked_I_10d5m_v.*wouldLieIfAsked_5d1m_v);
 
PSVL_F_NI_10d5m_v = (1-r_v)*PSVL_NF_NI_10d5m_v + r_v*mean(hStar_F_10d5m_v.*doesSvyIfAsked_NI_10d5m_v.*wouldLieIfAsked_v);
PSVL_F_I_10d5m_v = (1-r_v)*PSVL_NF_I_10d5m_v + r_v*mean(hStar_F_10d5m_v.*doesSvyIfAsked_I_10d5m_v.*wouldLieIfAsked_v);
PSVL_F_NI_10d5m_5d1m_v = (1-r_v)*PSVL_NF_NI_10d5m_5d1m_v + r_v*mean(hStar_F_10d5m_v.*doesSvyIfAsked_NI_10d5m_v.*wouldLieIfAsked_5d1m_v);
PSVL_F_I_10d5m_5d1m_v = (1-r_v)*PSVL_NF_I_10d5m_5d1m_v + r_v*mean(hStar_F_10d5m_v.*doesSvyIfAsked_I_10d5m_v.*wouldLieIfAsked_5d1m_v);
 
PSVL_FV_NI_10d5m_v = (1-r_v)*PSVL_NF_NI_10d5m_v + r_v*mean(hStar_FV_10d5m_v.*doesSvyIfAsked_I_10d5m_v.*wouldLieIfAsked_v);
PSVL_FV_I_10d5m_v = (1-r_v)*PSVL_NF_I_10d5m_v + r_v*mean(hStar_FV_10d5m_v.*doesSvyIfAsked_I_10d5m_v.*wouldLieIfAsked_v);
PSVL_FV_NI_10d5m_5d1m_v = (1-r_v)*PSVL_NF_NI_10d5m_5d1m_v + r_v*mean(hStar_FV_10d5m_v.*doesSvyIfAsked_I_10d5m_v.*wouldLieIfAsked_5d1m_v);
PSVL_FV_I_10d5m_5d1m_v = (1-r_v)*PSVL_NF_I_10d5m_5d1m_v + r_v*mean(hStar_FV_10d5m_v.*doesSvyIfAsked_I_10d5m_v.*wouldLieIfAsked_5d1m_v);
 
PSVL_OO_NI_10d5m_v = (1-r_v)*PSVL_NF_NI_10d5m_v + r_v*mean((1 .-optsOutIfSees_OO_10d5m_v).*hStar_F_10d5m_v.*doesSvyIfAsked_NI_10d5m_v.*wouldLieIfAsked_v);
PSVL_OO_I_10d5m_v = (1-r_v)*PSVL_NF_I_10d5m_v + r_v*mean((1 .-optsOutIfSees_OO_10d5m_v).*hStar_F_10d5m_v.*doesSvyIfAsked_I_10d5m_v.*wouldLieIfAsked_v);
PSVL_OO_NI_10d5m_5d1m_v = (1-r_v)*PSVL_NF_NI_10d5m_5d1m_v + r_v*mean((1 .-optsOutIfSees_OO_10d5m_v).*hStar_F_10d5m_v.*doesSvyIfAsked_NI_10d5m_v.*wouldLieIfAsked_5d1m_v);
PSVL_OO_I_10d5m_5d1m_v = (1-r_v)*PSVL_NF_I_10d5m_5d1m_v + r_v*mean((1 .-optsOutIfSees_OO_10d5m_v).*hStar_F_10d5m_v.*doesSvyIfAsked_I_10d5m_v.*wouldLieIfAsked_5d1m_v);
 
PSVL_OOV_NI_10d5m_v = (1-r_v)*PSVL_NF_NI_10d5m_v + r_v*mean((1 .-optsOutIfSees_OOV_10d5m_v).*hStar_FV_10d5m_v.*doesSvyIfAsked_I_10d5m_v.*wouldLieIfAsked_v);
PSVL_OOV_I_10d5m_v = (1-r_v)*PSVL_NF_I_10d5m_v + r_v*mean((1 .-optsOutIfSees_OOV_10d5m_v).*hStar_FV_10d5m_v.*doesSvyIfAsked_I_10d5m_v.*wouldLieIfAsked_v);
PSVL_OOV_NI_10d5m_5d1m_v = (1-r_v)*PSVL_NF_NI_10d5m_5d1m_v + r_v*mean((1 .-optsOutIfSees_OOV_10d5m_v).*hStar_FV_10d5m_v.*doesSvyIfAsked_I_10d5m_v.*wouldLieIfAsked_5d1m_v);
PSVL_OOV_I_10d5m_5d1m_v = (1-r_v)*PSVL_NF_I_10d5m_5d1m_v + r_v*mean((1 .-optsOutIfSees_OOV_10d5m_v).*hStar_FV_10d5m_v.*doesSvyIfAsked_I_10d5m_v.*wouldLieIfAsked_5d1m_v);


# !!!Non-Voters!!!
	
	
# PH = probability of being at home
PH_NF_0d5m_nv = h0_nv;
PH_F_0d5m_nv = (1-r_nv)*h0_nv + r_nv.*mean(hStar_F_0d5m_nv);
PH_FV_0d5m_nv = (1-r_nv)*h0_nv + r_nv.*mean(hStar_FV_0d5m_nv);
PH_OO_0d5m_nv = (1-r_nv)*h0_nv + r_nv.*mean((1 .-optsOutIfSees_OO_0d5m_nv).*hStar_F_0d5m_nv); 
PH_OOV_0d5m_nv = (1-r_nv)*h0_nv + r_nv.*mean((1 .-optsOutIfSees_OOV_0d5m_nv).*hStar_FV_0d5m_nv);
 
PH_NF_10d10m_nv = h0_nv;
PH_F_10d10m_nv = (1-r_nv)*h0_nv + r_nv.*mean(hStar_F_10d10m_nv);
PH_FV_10d10m_nv = (1-r_nv)*h0_nv + r_nv.*mean(hStar_FV_10d10m_nv);
PH_OO_10d10m_nv = (1-r_nv)*h0_nv + r_nv.*mean((1 .-optsOutIfSees_OO_10d10m_nv).*hStar_F_10d10m_nv); 
PH_OOV_10d10m_nv = (1-r_nv)*h0_nv + r_nv.*mean((1 .-optsOutIfSees_OOV_10d10m_nv).*hStar_FV_10d10m_nv);
 
PH_NF_10d5m_nv = h0_nv;
PH_F_10d5m_nv = (1-r_nv)*h0_nv + r_nv.*mean(hStar_F_10d5m_nv);
PH_FV_10d5m_nv = (1-r_nv)*h0_nv + r_nv.*mean(hStar_FV_10d5m_nv);
PH_OO_10d5m_nv = (1-r_nv)*h0_nv + r_nv.*mean((1 .-optsOutIfSees_OO_10d5m_nv).*hStar_F_10d5m_nv); 
PH_OOV_10d5m_nv = (1-r_nv)*h0_nv + r_nv.*mean((1 .-optsOutIfSees_OOV_10d5m_nv).*hStar_FV_10d5m_nv);
 
 
# PSV=unconditional prob of doing the survey (not cond on opening door). PSV < PH mechanically
 
# 0d5m
PSV_NF_NI_0d5m_nv = h0_nv*mean(doesSvyIfAsked_NI_0d5m_nv);
PSV_NF_I_0d5m_nv = h0_nv*mean(doesSvyIfAsked_I_0d5m_nv);
 
PSV_F_NI_0d5m_nv = (1-r_nv)*PSV_NF_NI_0d5m_nv + r_nv*mean(hStar_F_0d5m_nv.*doesSvyIfAsked_NI_0d5m_nv);
PSV_F_I_0d5m_nv = (1-r_nv)*PSV_NF_I_0d5m_nv + r_nv*mean(hStar_F_0d5m_nv.*doesSvyIfAsked_I_0d5m_nv);
 
PSV_FV_NI_0d5m_nv = (1-r_nv)*PSV_NF_NI_0d5m_nv + r_nv*mean(hStar_FV_0d5m_nv.*doesSvyIfAsked_I_0d5m_nv);
PSV_FV_I_0d5m_nv = (1-r_nv)*PSV_NF_I_0d5m_nv + r_nv*mean(hStar_FV_0d5m_nv.*doesSvyIfAsked_I_0d5m_nv);
 
PSV_OO_NI_0d5m_nv = (1-r_nv)*PSV_NF_NI_0d5m_nv + r_nv*mean((1 .-optsOutIfSees_OO_0d5m_nv).*hStar_F_0d5m_nv.*doesSvyIfAsked_NI_0d5m_nv);
PSV_OO_I_0d5m_nv = (1-r_nv)*PSV_NF_I_0d5m_nv + r_nv*mean((1 .-optsOutIfSees_OO_0d5m_nv).*hStar_F_0d5m_nv.*doesSvyIfAsked_I_0d5m_nv);
 
PSV_OOV_NI_0d5m_nv = (1-r_nv)*PSV_NF_NI_0d5m_nv + r_nv*mean((1 .-optsOutIfSees_OOV_0d5m_nv).*hStar_FV_0d5m_nv.*doesSvyIfAsked_I_0d5m_nv);
PSV_OOV_I_0d5m_nv = (1-r_nv)*PSV_NF_I_0d5m_nv + r_nv*mean((1 .-optsOutIfSees_OOV_0d5m_nv).*hStar_FV_0d5m_nv.*doesSvyIfAsked_I_0d5m_nv);
 
# 10d10m
PSV_NF_NI_10d10m_nv = h0_nv*mean(doesSvyIfAsked_NI_10d10m_nv);
PSV_NF_I_10d10m_nv = h0_nv*mean(doesSvyIfAsked_I_10d10m_nv);
 
PSV_F_NI_10d10m_nv = (1-r_nv)*PSV_NF_NI_10d10m_nv + r_nv*mean(hStar_F_10d10m_nv.*doesSvyIfAsked_NI_10d10m_nv);
PSV_F_I_10d10m_nv = (1-r_nv)*PSV_NF_I_10d10m_nv + r_nv*mean(hStar_F_10d10m_nv.*doesSvyIfAsked_I_10d10m_nv);
 
PSV_FV_NI_10d10m_nv = (1-r_nv)*PSV_NF_NI_10d10m_nv + r_nv*mean(hStar_FV_10d10m_nv.*doesSvyIfAsked_I_10d10m_nv);
PSV_FV_I_10d10m_nv = (1-r_nv)*PSV_NF_I_10d10m_nv + r_nv*mean(hStar_FV_10d10m_nv.*doesSvyIfAsked_I_10d10m_nv);
 
PSV_OO_NI_10d10m_nv = (1-r_nv)*PSV_NF_NI_10d10m_nv + r_nv*mean((1 .-optsOutIfSees_OO_10d10m_nv).*hStar_F_10d10m_nv.*doesSvyIfAsked_NI_10d10m_nv);
PSV_OO_I_10d10m_nv = (1-r_nv)*PSV_NF_I_10d10m_nv + r_nv*mean((1 .-optsOutIfSees_OO_10d10m_nv).*hStar_F_10d10m_nv.*doesSvyIfAsked_I_10d10m_nv);
 
PSV_OOV_NI_10d10m_nv = (1-r_nv)*PSV_NF_NI_10d10m_nv + r_nv*mean((1 .-optsOutIfSees_OOV_10d10m_nv).*hStar_FV_10d10m_nv.*doesSvyIfAsked_I_10d10m_nv);
PSV_OOV_I_10d10m_nv = (1-r_nv)*PSV_NF_I_10d10m_nv + r_nv*mean((1 .-optsOutIfSees_OOV_10d10m_nv).*hStar_FV_10d10m_nv.*doesSvyIfAsked_I_10d10m_nv);
 
# 10d5m
PSV_NF_NI_10d5m_nv = h0_nv*mean(doesSvyIfAsked_NI_10d5m_nv);
PSV_NF_I_10d5m_nv = h0_nv*mean(doesSvyIfAsked_I_10d5m_nv);
 
PSV_F_NI_10d5m_nv = (1-r_nv)*PSV_NF_NI_10d5m_nv + r_nv*mean(hStar_F_10d5m_nv.*doesSvyIfAsked_NI_10d5m_nv);
PSV_F_I_10d5m_nv = (1-r_nv)*PSV_NF_I_10d5m_nv + r_nv*mean(hStar_F_10d5m_nv.*doesSvyIfAsked_I_10d5m_nv);
 
PSV_FV_NI_10d5m_nv = (1-r_nv)*PSV_NF_NI_10d5m_nv + r_nv*mean(hStar_FV_10d5m_nv.*doesSvyIfAsked_I_10d5m_nv);
PSV_FV_I_10d5m_nv = (1-r_nv)*PSV_NF_I_10d5m_nv + r_nv*mean(hStar_FV_10d5m_nv.*doesSvyIfAsked_I_10d5m_nv);
 
PSV_OO_NI_10d5m_nv = (1-r_nv)*PSV_NF_NI_10d5m_nv + r_nv*mean((1 .-optsOutIfSees_OO_10d5m_nv).*hStar_F_10d5m_nv.*doesSvyIfAsked_NI_10d5m_nv);
PSV_OO_I_10d5m_nv = (1-r_nv)*PSV_NF_I_10d5m_nv + r_nv*mean((1 .-optsOutIfSees_OO_10d5m_nv).*hStar_F_10d5m_nv.*doesSvyIfAsked_I_10d5m_nv);
 
PSV_OOV_NI_10d5m_nv = (1-r_nv)*PSV_NF_NI_10d5m_nv + r_nv*mean((1 .-optsOutIfSees_OOV_10d5m_nv).*hStar_FV_10d5m_nv.*doesSvyIfAsked_I_10d5m_nv);
PSV_OOV_I_10d5m_nv = (1-r_nv)*PSV_NF_I_10d5m_nv + r_nv*mean((1 .-optsOutIfSees_OOV_10d5m_nv).*hStar_FV_10d5m_nv.*doesSvyIfAsked_I_10d5m_nv);
 
 
# POO=prob of opting out (not conditional on seeing flyer). Scaled by baseline likelihood of being at home
POO_OO_0d5m_nv =  h0_nv*r_nv*mean(optsOutIfSees_OO_0d5m_nv);
POO_OOV_0d5m_nv = h0_nv*r_nv*mean(optsOutIfSees_OOV_0d5m_nv);
 
POO_OO_10d10m_nv =  h0_nv*r_nv*mean(optsOutIfSees_OO_10d10m_nv);
POO_OOV_10d10m_nv = h0_nv*r_nv*mean(optsOutIfSees_OOV_10d10m_nv);
 
POO_OO_10d5m_nv =  h0_nv*r_nv*mean(optsOutIfSees_OO_10d5m_nv);
POO_OOV_10d5m_nv = h0_nv*r_nv*mean(optsOutIfSees_OOV_10d5m_nv);

 
# PSVL = unconditional percent who do survey and lie No flyer treatment only, simplifies later code
 
# PL=cond on agreeing to do the survey, did you lie? Incentive to lie is a surprise later (doesn't affect PH or PSV)
 
# 0d5m, 5d1m incentive
PSVL_NF_NI_0d5m_nv = mean(h0_nv.*doesSvyIfAsked_NI_0d5m_nv.*wouldLieIfAsked_nv);
PSVL_NF_I_0d5m_nv = mean(h0_nv.*doesSvyIfAsked_I_0d5m_nv.*wouldLieIfAsked_nv);
PSVL_NF_NI_0d5m_5d1m_nv = mean(h0_nv.*doesSvyIfAsked_NI_0d5m_nv.*wouldLieIfAsked_5d1m_nv);
PSVL_NF_I_0d5m_5d1m_nv = mean(h0_nv.*doesSvyIfAsked_I_0d5m_nv.*wouldLieIfAsked_5d1m_nv);

PSVL_F_NI_0d5m_nv = (1-r_nv)*PSVL_NF_NI_0d5m_nv + r_nv*mean(hStar_F_0d5m_nv.*doesSvyIfAsked_NI_0d5m_nv.*wouldLieIfAsked_nv);
PSVL_F_I_0d5m_nv = (1-r_nv)*PSVL_NF_I_0d5m_nv + r_nv*mean(hStar_F_0d5m_nv.*doesSvyIfAsked_I_0d5m_nv.*wouldLieIfAsked_nv);
PSVL_F_NI_0d5m_5d1m_nv = (1-r_nv)*PSVL_NF_NI_0d5m_5d1m_nv + r_nv*mean(hStar_F_0d5m_nv.*doesSvyIfAsked_NI_0d5m_nv.*wouldLieIfAsked_5d1m_nv);
PSVL_F_I_0d5m_5d1m_nv = (1-r_nv)*PSVL_NF_I_0d5m_5d1m_nv + r_nv*mean(hStar_F_0d5m_nv.*doesSvyIfAsked_I_0d5m_nv.*wouldLieIfAsked_5d1m_nv);
 
PSVL_FV_NI_0d5m_nv = (1-r_nv)*PSVL_NF_NI_0d5m_nv + r_nv*mean(hStar_FV_0d5m_nv.*doesSvyIfAsked_I_0d5m_nv.*wouldLieIfAsked_nv);
PSVL_FV_I_0d5m_nv = (1-r_nv)*PSVL_NF_I_0d5m_nv + r_nv*mean(hStar_FV_0d5m_nv.*doesSvyIfAsked_I_0d5m_nv.*wouldLieIfAsked_nv);
PSVL_FV_NI_0d5m_5d1m_nv = (1-r_nv)*PSVL_NF_NI_0d5m_5d1m_nv + r_nv*mean(hStar_FV_0d5m_nv.*doesSvyIfAsked_I_0d5m_nv.*wouldLieIfAsked_5d1m_nv);
PSVL_FV_I_0d5m_5d1m_nv = (1-r_nv)*PSVL_NF_I_0d5m_5d1m_nv + r_nv*mean(hStar_FV_0d5m_nv.*doesSvyIfAsked_I_0d5m_nv.*wouldLieIfAsked_5d1m_nv);
 
PSVL_OO_NI_0d5m_nv = (1-r_nv)*PSVL_NF_NI_0d5m_nv + r_nv*mean((1 .-optsOutIfSees_OO_0d5m_nv).*hStar_F_0d5m_nv.*doesSvyIfAsked_NI_0d5m_nv.*wouldLieIfAsked_nv);
PSVL_OO_I_0d5m_nv = (1-r_nv)*PSVL_NF_I_0d5m_nv + r_nv*mean((1 .-optsOutIfSees_OO_0d5m_nv).*hStar_F_0d5m_nv.*doesSvyIfAsked_I_0d5m_nv.*wouldLieIfAsked_nv);
PSVL_OO_NI_0d5m_5d1m_nv = (1-r_nv)*PSVL_NF_NI_0d5m_5d1m_nv + r_nv*mean((1 .-optsOutIfSees_OO_0d5m_nv).*hStar_F_0d5m_nv.*doesSvyIfAsked_NI_0d5m_nv.*wouldLieIfAsked_5d1m_nv); 
PSVL_OO_I_0d5m_5d1m_nv = (1-r_nv)*PSVL_NF_I_0d5m_5d1m_nv + r_nv*mean((1 .-optsOutIfSees_OO_0d5m_nv).*hStar_F_0d5m_nv.*doesSvyIfAsked_I_0d5m_nv.*wouldLieIfAsked_5d1m_nv);
 
PSVL_OOV_NI_0d5m_nv = (1-r_nv)*PSVL_NF_NI_0d5m_nv + r_nv*mean((1 .-optsOutIfSees_OOV_0d5m_nv).*hStar_FV_0d5m_nv.*doesSvyIfAsked_I_0d5m_nv.*wouldLieIfAsked_nv);
PSVL_OOV_I_0d5m_nv = (1-r_nv)*PSVL_NF_I_0d5m_nv + r_nv*mean((1 .-optsOutIfSees_OOV_0d5m_nv).*hStar_FV_0d5m_nv.*doesSvyIfAsked_I_0d5m_nv.*wouldLieIfAsked_nv);
PSVL_OOV_NI_0d5m_5d1m_nv = (1-r_nv)*PSVL_NF_NI_0d5m_5d1m_nv + r_nv*mean((1 .-optsOutIfSees_OOV_0d5m_nv).*hStar_FV_0d5m_nv.*doesSvyIfAsked_I_0d5m_nv.*wouldLieIfAsked_5d1m_nv);
PSVL_OOV_I_0d5m_5d1m_nv = (1-r_nv)*PSVL_NF_I_0d5m_5d1m_nv + r_nv*mean((1 .-optsOutIfSees_OOV_0d5m_nv).*hStar_FV_0d5m_nv.*doesSvyIfAsked_I_0d5m_nv.*wouldLieIfAsked_5d1m_nv);
 
# 10d10m, 8m incentive
PSVL_NF_NI_10d10m_nv = mean(h0_nv.*doesSvyIfAsked_NI_10d10m_nv.*wouldLieIfAsked_nv);
PSVL_NF_I_10d10m_nv = mean(h0_nv.*doesSvyIfAsked_I_10d10m_nv.*wouldLieIfAsked_nv);
PSVL_NF_NI_10d10m_8m_nv = mean(h0_nv.*doesSvyIfAsked_NI_10d10m_nv.*wouldLieIfAsked_8m_nv);
PSVL_NF_I_10d10m_8m_nv = mean(h0_nv.*doesSvyIfAsked_I_10d10m_nv.*wouldLieIfAsked_8m_nv);
 
PSVL_F_NI_10d10m_nv = (1-r_nv)*PSVL_NF_NI_10d10m_nv + r_nv*mean(hStar_F_10d10m_nv.*doesSvyIfAsked_NI_10d10m_nv.*wouldLieIfAsked_nv);
PSVL_F_I_10d10m_nv = (1-r_nv)*PSVL_NF_I_10d10m_nv + r_nv*mean(hStar_F_10d10m_nv.*doesSvyIfAsked_I_10d10m_nv.*wouldLieIfAsked_nv);
PSVL_F_NI_10d10m_8m_nv = (1-r_nv)*PSVL_NF_NI_10d10m_8m_nv + r_nv*mean(hStar_F_10d10m_nv.*doesSvyIfAsked_NI_10d10m_nv.*wouldLieIfAsked_8m_nv);
PSVL_F_I_10d10m_8m_nv = (1-r_nv)*PSVL_NF_I_10d10m_8m_nv + r_nv*mean(hStar_F_10d10m_nv.*doesSvyIfAsked_I_10d10m_nv.*wouldLieIfAsked_8m_nv);
 
PSVL_FV_NI_10d10m_nv = (1-r_nv)*PSVL_NF_NI_10d10m_nv + r_nv*mean(hStar_FV_10d10m_nv.*doesSvyIfAsked_I_10d10m_nv.*wouldLieIfAsked_nv);
PSVL_FV_I_10d10m_nv = (1-r_nv)*PSVL_NF_I_10d10m_nv + r_nv*mean(hStar_FV_10d10m_nv.*doesSvyIfAsked_I_10d10m_nv.*wouldLieIfAsked_nv);
PSVL_FV_NI_10d10m_8m_nv = (1-r_nv)*PSVL_NF_NI_10d10m_8m_nv + r_nv*mean(hStar_FV_10d10m_nv.*doesSvyIfAsked_I_10d10m_nv.*wouldLieIfAsked_8m_nv);
PSVL_FV_I_10d10m_8m_nv = (1-r_nv)*PSVL_NF_I_10d10m_8m_nv + r_nv*mean(hStar_FV_10d10m_nv.*doesSvyIfAsked_I_10d10m_nv.*wouldLieIfAsked_8m_nv);
 
PSVL_OO_NI_10d10m_nv = (1-r_nv)*PSVL_NF_NI_10d10m_nv + r_nv*mean((1 .-optsOutIfSees_OO_10d10m_nv).*hStar_F_10d10m_nv.*doesSvyIfAsked_NI_10d10m_nv.*wouldLieIfAsked_nv);
PSVL_OO_I_10d10m_nv = (1-r_nv)*PSVL_NF_I_10d10m_nv + r_nv*mean((1 .-optsOutIfSees_OO_10d10m_nv).*hStar_F_10d10m_nv.*doesSvyIfAsked_I_10d10m_nv.*wouldLieIfAsked_nv);
PSVL_OO_NI_10d10m_8m_nv = (1-r_nv)*PSVL_NF_NI_10d10m_8m_nv + r_nv*mean((1 .-optsOutIfSees_OO_10d10m_nv).*hStar_F_10d10m_nv.*doesSvyIfAsked_NI_10d10m_nv.*wouldLieIfAsked_8m_nv);
PSVL_OO_I_10d10m_8m_nv = (1-r_nv)*PSVL_NF_I_10d10m_8m_nv + r_nv*mean((1 .-optsOutIfSees_OO_10d10m_nv).*hStar_F_10d10m_nv.*doesSvyIfAsked_I_10d10m_nv.*wouldLieIfAsked_8m_nv);
 
PSVL_OOV_NI_10d10m_nv = (1-r_nv)*PSVL_NF_NI_10d10m_nv + r_nv*mean((1 .-optsOutIfSees_OOV_10d10m_nv).*hStar_FV_10d10m_nv.*doesSvyIfAsked_I_10d10m_nv.*wouldLieIfAsked_nv);
PSVL_OOV_I_10d10m_nv = (1-r_nv)*PSVL_NF_I_10d10m_nv + r_nv*mean((1 .-optsOutIfSees_OOV_10d10m_nv).*hStar_FV_10d10m_nv.*doesSvyIfAsked_I_10d10m_nv.*wouldLieIfAsked_nv);
PSVL_OOV_NI_10d10m_8m_nv = (1-r_nv)*PSVL_NF_NI_10d10m_8m_nv + r_nv*mean((1 .-optsOutIfSees_OOV_10d10m_nv).*hStar_FV_10d10m_nv.*doesSvyIfAsked_I_10d10m_nv.*wouldLieIfAsked_8m_nv);
PSVL_OOV_I_10d10m_8m_nv = (1-r_nv)*PSVL_NF_I_10d10m_8m_nv + r_nv*mean((1 .-optsOutIfSees_OOV_10d10m_nv).*hStar_FV_10d10m_nv.*doesSvyIfAsked_I_10d10m_nv.*wouldLieIfAsked_8m_nv);
 
# 10d5m, 5d1m incentive
PSVL_NF_NI_10d5m_nv = mean(h0_nv.*doesSvyIfAsked_NI_10d5m_nv.*wouldLieIfAsked_nv);
PSVL_NF_I_10d5m_nv = mean(h0_nv.*doesSvyIfAsked_I_10d5m_nv.*wouldLieIfAsked_nv);
PSVL_NF_NI_10d5m_5d1m_nv = mean(h0_nv.*doesSvyIfAsked_NI_10d5m_nv.*wouldLieIfAsked_5d1m_nv);
PSVL_NF_I_10d5m_5d1m_nv = mean(h0_nv.*doesSvyIfAsked_I_10d5m_nv.*wouldLieIfAsked_5d1m_nv);
 
PSVL_F_NI_10d5m_nv = (1-r_nv)*PSVL_NF_NI_10d5m_nv + r_nv*mean(hStar_F_10d5m_nv.*doesSvyIfAsked_NI_10d5m_nv.*wouldLieIfAsked_nv);
PSVL_F_I_10d5m_nv = (1-r_nv)*PSVL_NF_I_10d5m_nv + r_nv*mean(hStar_F_10d5m_nv.*doesSvyIfAsked_I_10d5m_nv.*wouldLieIfAsked_nv);
PSVL_F_NI_10d5m_5d1m_nv = (1-r_nv)*PSVL_NF_NI_10d5m_5d1m_nv + r_nv*mean(hStar_F_10d5m_nv.*doesSvyIfAsked_NI_10d5m_nv.*wouldLieIfAsked_5d1m_nv);
PSVL_F_I_10d5m_5d1m_nv = (1-r_nv)*PSVL_NF_I_10d5m_5d1m_nv + r_nv*mean(hStar_F_10d5m_nv.*doesSvyIfAsked_I_10d5m_nv.*wouldLieIfAsked_5d1m_nv);
 
PSVL_FV_NI_10d5m_nv = (1-r_nv)*PSVL_NF_NI_10d5m_nv + r_nv*mean(hStar_FV_10d5m_nv.*doesSvyIfAsked_I_10d5m_nv.*wouldLieIfAsked_nv);
PSVL_FV_I_10d5m_nv = (1-r_nv)*PSVL_NF_I_10d5m_nv + r_nv*mean(hStar_FV_10d5m_nv.*doesSvyIfAsked_I_10d5m_nv.*wouldLieIfAsked_nv);
PSVL_FV_NI_10d5m_5d1m_nv = (1-r_nv)*PSVL_NF_NI_10d5m_5d1m_nv + r_nv*mean(hStar_FV_10d5m_nv.*doesSvyIfAsked_I_10d5m_nv.*wouldLieIfAsked_5d1m_nv);
PSVL_FV_I_10d5m_5d1m_nv = (1-r_nv)*PSVL_NF_I_10d5m_5d1m_nv + r_nv*mean(hStar_FV_10d5m_nv.*doesSvyIfAsked_I_10d5m_nv.*wouldLieIfAsked_5d1m_nv);
 
PSVL_OO_NI_10d5m_nv = (1-r_nv)*PSVL_NF_NI_10d5m_nv + r_nv*mean((1 .-optsOutIfSees_OO_10d5m_nv).*hStar_F_10d5m_nv.*doesSvyIfAsked_NI_10d5m_nv.*wouldLieIfAsked_nv);
PSVL_OO_I_10d5m_nv = (1-r_nv)*PSVL_NF_I_10d5m_nv + r_nv*mean((1 .-optsOutIfSees_OO_10d5m_nv).*hStar_F_10d5m_nv.*doesSvyIfAsked_I_10d5m_nv.*wouldLieIfAsked_nv);
PSVL_OO_NI_10d5m_5d1m_nv = (1-r_nv)*PSVL_NF_NI_10d5m_5d1m_nv + r_nv*mean((1 .-optsOutIfSees_OO_10d5m_nv).*hStar_F_10d5m_nv.*doesSvyIfAsked_NI_10d5m_nv.*wouldLieIfAsked_5d1m_nv);
PSVL_OO_I_10d5m_5d1m_nv = (1-r_nv)*PSVL_NF_I_10d5m_5d1m_nv + r_nv*mean((1 .-optsOutIfSees_OO_10d5m_nv).*hStar_F_10d5m_nv.*doesSvyIfAsked_I_10d5m_nv.*wouldLieIfAsked_5d1m_nv);
 
PSVL_OOV_NI_10d5m_nv = (1-r_nv)*PSVL_NF_NI_10d5m_nv + r_nv*mean((1 .-optsOutIfSees_OOV_10d5m_nv).*hStar_FV_10d5m_nv.*doesSvyIfAsked_I_10d5m_nv.*wouldLieIfAsked_nv);
PSVL_OOV_I_10d5m_nv = (1-r_nv)*PSVL_NF_I_10d5m_nv + r_nv*mean((1 .-optsOutIfSees_OOV_10d5m_nv).*hStar_FV_10d5m_nv.*doesSvyIfAsked_I_10d5m_nv.*wouldLieIfAsked_nv);
PSVL_OOV_NI_10d5m_5d1m_nv = (1-r_nv)*PSVL_NF_NI_10d5m_5d1m_nv + r_nv*mean((1 .-optsOutIfSees_OOV_10d5m_nv).*hStar_FV_10d5m_nv.*doesSvyIfAsked_I_10d5m_nv.*wouldLieIfAsked_5d1m_nv);
PSVL_OOV_I_10d5m_5d1m_nv = (1-r_nv)*PSVL_NF_I_10d5m_5d1m_nv + r_nv*mean((1 .-optsOutIfSees_OOV_10d5m_nv).*hStar_FV_10d5m_nv.*doesSvyIfAsked_I_10d5m_nv.*wouldLieIfAsked_5d1m_nv);
 
	
# 30 PH moments (5 tx * 3 length * 2 v/nv)
    sm = Vector{Real}(undef, 101)
    sm[1]	=		PH_NF_0d5m_v		
    sm[2]	=		PH_NF_10d10m_v		
    sm[3]	=		PH_NF_10d5m_v		
    sm[4]	=		PH_F_0d5m_v	
    sm[5]	=		PH_F_10d10m_v		
    sm[6]	=		PH_F_10d5m_v		
    sm[7]	=		PH_FV_0d5m_v		
    sm[8]	=		PH_FV_10d10m_v		
    sm[9]	=		PH_FV_10d5m_v		
    sm[10]	=		PH_OO_0d5m_v		
    sm[11]	=		PH_OO_10d10m_v		
    sm[12]	=		PH_OO_10d5m_v		
    sm[13]	=		PH_OOV_0d5m_v		
    sm[14]	=		PH_OOV_10d10m_v		
    sm[15]	=		PH_OOV_10d5m_v		
    sm[16]	=		PH_NF_0d5m_nv		
    sm[17]	=		PH_NF_10d10m_nv		
    sm[18]	=		PH_NF_10d5m_nv		
    sm[19]	=		PH_F_0d5m_nv		
    sm[20]	=		PH_F_10d10m_nv		
    sm[21]	=		PH_F_10d5m_nv		
    sm[22]	=		PH_FV_0d5m_nv		
    sm[23]	=		PH_FV_10d10m_nv		
    sm[24]	=		PH_FV_10d5m_nv		
    sm[25]	=		PH_OO_0d5m_nv		
    sm[26]	=		PH_OO_10d10m_nv		
    sm[27]	=		PH_OO_10d5m_nv		
    sm[28]	=		PH_OOV_0d5m_nv		
    sm[29]	=		PH_OOV_10d10m_nv		
    sm[30]	=		PH_OOV_10d5m_nv	

# 30 PSV moments (5 tx * 3 length * 2 v/nv). Taking 50/50 average across NI and I treatments
    sm[31]  =   mean([  PSV_NF_I_0d5m_v,     PSV_NF_NI_0d5m_v    ])
    sm[32]  =   mean([  PSV_NF_I_10d10m_v,   PSV_NF_NI_10d10m_v  ])
    sm[33]  =   mean([  PSV_NF_I_10d5m_v,    PSV_NF_NI_10d5m_v   ]) 
    sm[34]  =   mean([  PSV_F_I_0d5m_v,      PSV_F_NI_0d5m_v     ]) 
    sm[35]  =   mean([  PSV_F_I_10d10m_v,    PSV_F_NI_10d10m_v   ]) 
    sm[36]  =   mean([  PSV_F_I_10d5m_v,     PSV_F_NI_10d5m_v    ]) 
    sm[37]  =   mean([  PSV_FV_I_0d5m_v,     PSV_FV_NI_0d5m_v    ]) 
    sm[38]  =   mean([  PSV_FV_I_10d10m_v,   PSV_FV_NI_10d10m_v  ]) 
    sm[39]  =   mean([  PSV_FV_I_10d5m_v,    PSV_FV_NI_10d5m_v   ]) 
    sm[40] =   mean([  PSV_OO_I_0d5m_v,     PSV_OO_NI_0d5m_v    ]) 
    sm[41] =   mean([  PSV_OO_I_10d10m_v,   PSV_OO_NI_10d10m_v  ]) 
    sm[42] =   mean([  PSV_OO_I_10d5m_v,    PSV_OO_NI_10d5m_v   ]) 
    sm[43] =   mean([  PSV_OOV_I_0d5m_v,    PSV_OOV_NI_0d5m_v   ]) 
    sm[44] =   mean([  PSV_OOV_I_10d10m_v,  PSV_OOV_NI_10d10m_v ]) 
    sm[45] =   mean([  PSV_OOV_I_10d5m_v,   PSV_OOV_NI_10d5m_v  ]) 
    sm[46] =   mean([  PSV_NF_I_0d5m_nv,    PSV_NF_NI_0d5m_nv   ]) 
    sm[47] =   mean([  PSV_NF_I_10d10m_nv,  PSV_NF_NI_10d10m_nv ]) 
    sm[48] =   mean([  PSV_NF_I_10d5m_nv,   PSV_NF_NI_10d5m_nv  ]) 
    sm[49] =   mean([  PSV_F_I_0d5m_nv,     PSV_F_NI_0d5m_nv    ]) 
    sm[50] =   mean([  PSV_F_I_10d10m_nv,   PSV_F_NI_10d10m_nv  ]) 
    sm[51] =   mean([  PSV_F_I_10d5m_nv,    PSV_F_NI_10d5m_nv   ]) 
    sm[52] =   mean([  PSV_FV_I_0d5m_nv,    PSV_FV_NI_0d5m_nv   ]) 
    sm[53] =   mean([  PSV_FV_I_10d10m_nv,  PSV_FV_NI_10d10m_nv ]) 
    sm[54] =   mean([  PSV_FV_I_10d5m_nv,   PSV_FV_NI_10d5m_nv  ]) 
    sm[55] =   mean([  PSV_OO_I_0d5m_nv,    PSV_OO_NI_0d5m_nv   ]) 
    sm[56] =   mean([  PSV_OO_I_10d10m_nv,  PSV_OO_NI_10d10m_nv ]) 
    sm[57] =   mean([  PSV_OO_I_10d5m_nv,   PSV_OO_NI_10d5m_nv  ]) 
    sm[58] =   mean([  PSV_OOV_I_0d5m_nv,   PSV_OOV_NI_0d5m_nv  ]) 
    sm[59] =   mean([  PSV_OOV_I_10d10m_nv, PSV_OOV_NI_10d10m_nv]) 
    sm[60] =   mean([  PSV_OOV_I_10d5m_nv,  PSV_OOV_NI_10d5m_nv ]) 	

# 12 POO moments (2 tx * 3 length * 2 v/nv)
    sm[61]	=		POO_OO_0d5m_v		
    sm[62]	=		POO_OO_10d10m_v	
    sm[63]	=		POO_OO_10d5m_v		
    sm[64]	=		POO_OOV_0d5m_v		
    sm[65]	=		POO_OOV_10d10m_v		
    sm[66]	=		POO_OOV_10d5m_v			
    sm[67]	=		POO_OO_0d5m_nv			
    sm[68]	=		POO_OO_10d10m_nv		
    sm[69]	=		POO_OO_10d5m_nv		
    sm[70]	=		POO_OOV_0d5m_nv		
    sm[71]	=		POO_OOV_10d10m_nv		
    sm[72]	=		POO_OOV_10d5m_nv

# 20 PSV by info moments (5 tx * 2 I/NI * 2 v/nv)
    sm[73]   =   mean([  PSV_NF_NI_0d5m_v,    PSV_NF_NI_10d10m_v,      PSV_NF_NI_10d5m_v   ])
    sm[74]   =   mean([  PSV_NF_I_0d5m_v,     PSV_NF_I_10d10m_v,       PSV_NF_I_10d5m_v  ])
    sm[75]   =   mean([  PSV_F_NI_0d5m_v,     PSV_F_NI_10d10m_v,       PSV_F_NI_10d5m_v    ])
    sm[76]   =   mean([  PSV_F_I_0d5m_v,      PSV_F_I_10d10m_v,       PSV_F_I_10d5m_v   ])
    sm[77]   =   mean([  PSV_FV_NI_0d5m_v,    PSV_FV_NI_10d10m_v,      PSV_FV_NI_10d5m_v   ])
    sm[78]   =   mean([  PSV_FV_I_0d5m_v,     PSV_FV_I_10d10m_v,       PSV_FV_I_10d5m_v  ])
    sm[79]   =   mean([  PSV_OO_NI_0d5m_v,    PSV_OO_NI_10d10m_v,      PSV_OO_NI_10d5m_v   ])
    sm[80]   =   mean([  PSV_OO_I_0d5m_v,     PSV_OO_I_10d10m_v,       PSV_OO_I_10d5m_v  ])
    sm[81]   =   mean([  PSV_OOV_NI_0d5m_v,   PSV_OOV_NI_10d10m_v,     PSV_OOV_NI_10d5m_v  ])
    sm[82]   =  mean([  PSV_OOV_I_0d5m_v,    PSV_OOV_I_10d10m_v,      PSV_OOV_I_10d5m_v ])
    sm[83]  =   mean([  PSV_NF_NI_0d5m_nv,   PSV_NF_NI_10d10m_nv,     PSV_NF_NI_10d5m_nv  ])
    sm[84]  =   mean([  PSV_NF_I_0d5m_nv,    PSV_NF_I_10d10m_nv,      PSV_NF_I_10d5m_nv ])
    sm[85]  =   mean([  PSV_F_NI_0d5m_nv,    PSV_F_NI_10d10m_nv,      PSV_F_NI_10d5m_nv   ])
    sm[86]  =   mean([  PSV_F_I_0d5m_nv,     PSV_F_I_10d10m_nv,       PSV_F_I_10d5m_nv  ])
    sm[87]  =   mean([  PSV_FV_NI_0d5m_nv,   PSV_FV_NI_10d10m_nv,     PSV_FV_NI_10d5m_nv  ])
    sm[88]  =   mean([  PSV_FV_I_0d5m_nv,    PSV_FV_I_10d10m_nv,      PSV_FV_I_10d5m_nv ])
    sm[89]  =   mean([  PSV_OO_NI_0d5m_nv,   PSV_OO_NI_10d10m_nv,     PSV_OO_NI_10d5m_nv  ])
    sm[90]  =   mean([  PSV_OO_I_0d5m_nv,    PSV_OO_I_10d10m_nv,      PSV_OO_I_10d5m_nv ])
    sm[91]  =   mean([  PSV_OOV_NI_0d5m_nv,  PSV_OOV_NI_10d10m_nv,    PSV_OOV_NI_10d5m_nv ])
    sm[92]  =   mean([  PSV_OOV_I_0d5m_nv,   PSV_OOV_I_10d10m_nv,     PSV_OOV_I_10d5m_nv ])

# 8 PL moments (1 tx * 2 10m/5m * 2 incentives)

# Empirical moments are sum of people lying in relevant tx divided by the sum of people answering the survey in relevant tx.
    sm[93] = mean([PSVL_NF_NI_0d5m_v, PSVL_NF_I_0d5m_v, PSVL_NF_NI_10d5m_v, PSVL_NF_I_10d5m_v,
                    PSVL_F_NI_0d5m_v, PSVL_F_I_0d5m_v, PSVL_F_NI_10d5m_v, PSVL_F_I_10d5m_v,
                    PSVL_FV_NI_0d5m_v, PSVL_FV_I_0d5m_v, PSVL_FV_NI_10d5m_v, PSVL_FV_I_10d5m_v,
                    PSVL_OO_NI_0d5m_v, PSVL_OO_I_0d5m_v, PSVL_OO_NI_10d5m_v, PSVL_OO_I_10d5m_v,
                    PSVL_OOV_NI_0d5m_v, PSVL_OOV_I_0d5m_v, PSVL_OOV_NI_10d5m_v, PSVL_OOV_I_10d5m_v])/mean([PSV_NF_NI_0d5m_v,
                    PSV_NF_I_0d5m_v, PSV_NF_NI_10d5m_v, PSV_NF_I_10d5m_v,
                    PSV_F_NI_0d5m_v, PSV_F_I_0d5m_v, PSV_F_NI_10d5m_v, PSV_F_I_10d5m_v,
                    PSV_FV_NI_0d5m_v, PSV_FV_I_0d5m_v, PSV_FV_NI_10d5m_v, PSV_FV_I_10d5m_v,
                    PSV_OO_NI_0d5m_v, PSV_OO_I_0d5m_v,  PSV_OO_NI_10d5m_v, PSV_OO_I_10d5m_v,
                    PSV_OOV_NI_0d5m_v, PSV_OOV_I_0d5m_v, PSV_OOV_NI_10d5m_v, PSV_OOV_I_10d5m_v])

    sm[94] = mean([PSVL_NF_NI_0d5m_5d1m_v, PSVL_NF_I_0d5m_5d1m_v, PSVL_NF_NI_10d5m_5d1m_v, PSVL_NF_I_10d5m_5d1m_v,
                    PSVL_F_NI_0d5m_5d1m_v, PSVL_F_I_0d5m_5d1m_v, PSVL_F_NI_10d5m_5d1m_v, PSVL_F_I_10d5m_5d1m_v,
                    PSVL_FV_NI_0d5m_5d1m_v, PSVL_FV_I_0d5m_5d1m_v, PSVL_FV_NI_10d5m_5d1m_v, PSVL_FV_I_10d5m_5d1m_v,
                    PSVL_OO_NI_0d5m_5d1m_v, PSVL_OO_I_0d5m_5d1m_v, PSVL_OO_NI_10d5m_5d1m_v, PSVL_OO_I_10d5m_5d1m_v,
                    PSVL_OOV_NI_0d5m_5d1m_v, PSVL_OOV_I_0d5m_5d1m_v, PSVL_OOV_NI_10d5m_5d1m_v, PSVL_OOV_I_10d5m_5d1m_v])/ mean([PSV_NF_NI_0d5m_v,
                    PSV_NF_I_0d5m_v, PSV_NF_NI_10d5m_v, PSV_NF_I_10d5m_v,
                    PSV_F_NI_0d5m_v, PSV_F_I_0d5m_v, PSV_F_NI_10d5m_v, PSV_F_I_10d5m_v,
                    PSV_FV_NI_0d5m_v, PSV_FV_I_0d5m_v, PSV_FV_NI_10d5m_v, PSV_FV_I_10d5m_v,
                    PSV_OO_NI_0d5m_v, PSV_OO_I_0d5m_v, PSV_OO_NI_10d5m_v, PSV_OO_I_10d5m_v,
                    PSV_OOV_NI_0d5m_v, PSV_OOV_I_0d5m_v, PSV_OOV_NI_10d5m_v, PSV_OOV_I_10d5m_v])   

    sm[95] = mean([PSVL_NF_NI_10d10m_v, PSVL_NF_I_10d10m_v,
                    PSVL_F_NI_10d10m_v, PSVL_F_I_10d10m_v,
                    PSVL_FV_NI_10d10m_v, PSVL_FV_I_10d10m_v,
                    PSVL_OO_NI_10d10m_v, PSVL_OO_I_10d10m_v,
                    PSVL_OOV_NI_10d10m_v, PSVL_OOV_I_10d10m_v])/ mean([PSV_NF_NI_10d10m_v,
                    PSV_NF_I_10d10m_v,
                    PSV_F_NI_10d10m_v, PSV_F_I_10d10m_v,
                    PSV_FV_NI_10d10m_v, PSV_FV_I_10d10m_v,
                    PSV_OO_NI_10d10m_v, PSV_OO_I_10d10m_v,
                    PSV_OOV_NI_10d10m_v, PSV_OOV_I_10d10m_v])
                
    sm[96] = mean([PSVL_NF_NI_10d10m_8m_v, PSVL_NF_I_10d10m_8m_v,
                    PSVL_F_NI_10d10m_8m_v, PSVL_F_I_10d10m_8m_v,
                    PSVL_FV_NI_10d10m_8m_v, PSVL_FV_I_10d10m_8m_v,
                    PSVL_OO_NI_10d10m_8m_v, PSVL_OO_I_10d10m_8m_v,
                    PSVL_OOV_NI_10d10m_8m_v, PSVL_OOV_I_10d10m_8m_v])/ mean([PSV_NF_NI_10d10m_v,
                    PSV_NF_I_10d10m_v,
                    PSV_F_NI_10d10m_v, PSV_F_I_10d10m_v, PSV_FV_NI_10d10m_v, PSV_FV_I_10d10m_v,
                    PSV_OO_NI_10d10m_v, PSV_OO_I_10d10m_v,
                    PSV_OOV_NI_10d10m_v, PSV_OOV_I_10d10m_v]) 
                
    sm[97] = mean([PSVL_NF_NI_0d5m_nv, PSVL_NF_I_0d5m_nv, PSVL_NF_NI_10d5m_nv, PSVL_NF_I_10d5m_nv,
                    PSVL_F_NI_0d5m_nv, PSVL_F_I_0d5m_nv, PSVL_F_NI_10d5m_nv, PSVL_F_I_10d5m_nv,
                    PSVL_FV_NI_0d5m_nv, PSVL_FV_I_0d5m_nv, PSVL_FV_NI_10d5m_nv, PSVL_FV_I_10d5m_nv, 
                    PSVL_OO_NI_0d5m_nv, PSVL_OO_I_0d5m_nv, PSVL_OO_NI_10d5m_nv, PSVL_OO_I_10d5m_nv,
                    PSVL_OOV_NI_0d5m_nv, PSVL_OOV_I_0d5m_nv, PSVL_OOV_NI_10d5m_nv, PSVL_OOV_I_10d5m_nv])/ mean([PSV_NF_NI_0d5m_nv,
                    PSV_NF_I_0d5m_nv, PSV_NF_NI_10d5m_nv, PSV_NF_I_10d5m_nv,
                    PSV_F_NI_0d5m_nv, PSV_F_I_0d5m_nv, PSV_F_NI_10d5m_nv, PSV_F_I_10d5m_nv,
                    PSV_FV_NI_0d5m_nv, PSV_FV_I_0d5m_nv, PSV_FV_NI_10d5m_nv, PSV_FV_I_10d5m_nv,
                    PSV_OO_NI_0d5m_nv, PSV_OO_I_0d5m_nv,  PSV_OO_NI_10d5m_nv, PSV_OO_I_10d5m_nv,
                    PSV_OOV_NI_0d5m_nv, PSV_OOV_I_0d5m_nv, PSV_OOV_NI_10d5m_nv, PSV_OOV_I_10d5m_nv]) 
 
    sm[98] = mean([PSVL_NF_NI_0d5m_5d1m_nv, PSVL_NF_I_0d5m_5d1m_nv, PSVL_NF_NI_10d5m_5d1m_nv, PSVL_NF_I_10d5m_5d1m_nv, PSVL_F_NI_0d5m_5d1m_nv, PSVL_F_I_0d5m_5d1m_nv, PSVL_F_NI_10d5m_5d1m_nv, PSVL_F_I_10d5m_5d1m_nv, PSVL_FV_NI_0d5m_5d1m_nv, PSVL_FV_I_0d5m_5d1m_nv,PSVL_FV_NI_10d5m_5d1m_nv, PSVL_FV_I_10d5m_5d1m_nv,PSVL_OO_NI_0d5m_5d1m_nv, PSVL_OO_I_0d5m_5d1m_nv, PSVL_OO_NI_10d5m_5d1m_nv, PSVL_OO_I_10d5m_5d1m_nv,PSVL_OOV_NI_0d5m_5d1m_nv, PSVL_OOV_I_0d5m_5d1m_nv, PSVL_OOV_NI_10d5m_5d1m_nv, PSVL_OOV_I_10d5m_5d1m_nv ])/mean([PSV_NF_NI_0d5m_nv, PSV_NF_I_0d5m_nv, PSV_NF_NI_10d5m_nv, PSV_NF_I_10d5m_nv,
PSV_F_NI_0d5m_nv, PSV_F_I_0d5m_nv, PSV_F_NI_10d5m_nv, PSV_F_I_10d5m_nv,
PSV_FV_NI_0d5m_nv, PSV_FV_I_0d5m_nv, PSV_FV_NI_10d5m_nv, PSV_FV_I_10d5m_nv,
PSV_OO_NI_0d5m_nv, PSV_OO_I_0d5m_nv,  PSV_OO_NI_10d5m_nv, PSV_OO_I_10d5m_nv,
PSV_OOV_NI_0d5m_nv, PSV_OOV_I_0d5m_nv, PSV_OOV_NI_10d5m_nv, PSV_OOV_I_10d5m_nv])

    sm[99] =	mean([PSVL_NF_NI_10d10m_nv, PSVL_NF_I_10d10m_nv,
                    PSVL_F_NI_10d10m_nv, PSVL_F_I_10d10m_nv,
                    PSVL_FV_NI_10d10m_nv, PSVL_FV_I_10d10m_nv,
                    PSVL_OO_NI_10d10m_nv, PSVL_OO_I_10d10m_nv,
                    PSVL_OOV_NI_10d10m_nv, PSVL_OOV_I_10d10m_nv])/mean([PSV_NF_NI_10d10m_nv,
                    PSV_NF_I_10d10m_nv,
                    PSV_F_NI_10d10m_nv, PSV_F_I_10d10m_nv,
                    PSV_FV_NI_10d10m_nv, PSV_FV_I_10d10m_nv,
                    PSV_OO_NI_10d10m_nv, PSV_OO_I_10d10m_nv,
                    PSV_OOV_NI_10d10m_nv, PSV_OOV_I_10d10m_nv])
                
    sm[100] =	mean([PSVL_NF_NI_10d10m_8m_nv, PSVL_NF_I_10d10m_8m_nv,
                    PSVL_F_NI_10d10m_8m_nv, PSVL_F_I_10d10m_8m_nv,
                    PSVL_FV_NI_10d10m_8m_nv, PSVL_FV_I_10d10m_8m_nv,
                    PSVL_OO_NI_10d10m_8m_nv, PSVL_OO_I_10d10m_8m_nv,
                    PSVL_OOV_NI_10d10m_8m_nv, PSVL_OOV_I_10d10m_8m_nv])/mean([PSV_NF_NI_10d10m_nv,
                    PSV_NF_I_10d10m_nv,
                    PSV_F_NI_10d10m_nv, PSV_F_I_10d10m_nv,
                    PSV_FV_NI_10d10m_nv, PSV_FV_I_10d10m_nv,
                    PSV_OO_NI_10d10m_nv, PSV_OO_I_10d10m_nv,
                    PSV_OOV_NI_10d10m_nv, PSV_OOV_I_10d10m_nv])
	
	sm[101] = Turnout_control     # baseline turnout

    return sm
end

# ╔═╡ c6468ee0-3588-11ec-0d41-cbfce972d21c
begin

# 3. Estimation

# Define the number N of simulated individuals and the draws from a standard normal distribution to compute the simulated s, s_v, s_nv, eps in the function

sim_voters = 750000
rand1 = rand(Normal(0.0, 1.0), sim_voters)
rand2 = rand(Normal(0.0, 1.0), sim_voters)
rand3 = rand(Normal(0.0, 1.0), sim_voters)
rand4 = rand(Normal(0.0, 1.0), sim_voters)
rand_vec = [rand1,rand2,rand3,rand4]
	
end;

# ╔═╡ c616cc50-3588-11ec-25bc-1b52b9c8694d
# This is the function we want to minimize: the weighted sum of the squared differences between empirical and simulated moments

function criterion(parameters, rand_set)
	
    simMoments = voteSimEndogenousVoting_vary(parameters,rand_set)
	
    m = emp_moments .- simMoments
    y = m' * W * m
	
    return y 
end

# ╔═╡ 1bb22460-37c4-11ec-3897-732e028fc413
begin

# Define the quasi-random starting guesses for the parameters. These are found in part B of the appendix

h0_v_in = rand(Uniform(0.2, 0.4), 1)
h0_nv_in = rand(Uniform(0.2, 0.4), 1)
r_v_in = rand(Uniform(0.2, 0.4), 1)
r_nv_in = rand(Uniform(0.2, 0.4), 1)
eta_v_in = rand(Uniform(0.0, 0.5), 1)
eta_nv_in = rand(Uniform(0.0, 0.5), 1)
mu_s_v_in = rand(Uniform(-50.0, 0.0), 1)
mu_s_nv_in = rand(Uniform(-50.0, 0.0), 1)
sigma_s_v_in = rand(Uniform(0.0, 50.0), 1)
sigma_s_nv_in = rand(Uniform(0.0, 50.0), 1)
S_svy_v_in = rand(Uniform(0.0, 10.0), 1)
S_svy_nv_in = rand(Uniform(0.0, 10.0), 1)
timeval_v_in = rand(Uniform(0.0, 100.0), 1)
timeval_nv_in = rand(Uniform(0.0, 100.0), 1)
mu_sv_in = rand(Uniform(-20.0, 20.0), 1)
mu_sn_in = rand(Uniform(-30.0, 10.0), 1)
sigma_svn_in = rand(Uniform(0.0, 30.0), 1)
L_in = rand(Uniform(0.0, 20.0), 1)
mu_eps_in = rand(Uniform(-30.0, 100.0), 1)
sigma_eps_in = rand(Uniform(50.0, 200.0), 1)
	
params_init=[h0_v_in,h0_nv_in,r_v_in,r_nv_in,eta_v_in,eta_nv_in,mu_s_v_in,mu_s_nv_in,sigma_s_v_in,
sigma_s_nv_in,S_svy_v_in,S_svy_nv_in,timeval_v_in,timeval_nv_in,mu_sv_in,mu_sn_in,sigma_svn_in,L_in,mu_eps_in,sigma_eps_in,]
params_init = [item for sublist in params_init for item in sublist] # flatten array
	
end;

# ╔═╡ e33132a0-37c5-11ec-1079-e3c9bba57502
begin

# !!! Read Before Running !!!

# This cell computes the estimates. On our machines it takes around 10-13 minutes for each 500 iterations. The number of iterations needed for convergence depends a lot on the starting guesses, nonetheless 5000 iterations are usually sufficient to get a wsse (weighted sum of squared errors) of 160 which is the one found also by the authors. In their Matlab code they run 720 optimization routines starting from different initial parameters to avoid being stuck in local minima and take as the best estimates those that provide the min wsse. Given time constraints we run the algorithm a dozen times and took the estimates that provided the lower weighted sum of squared errors

# !!! We reduced the number of iterations since Pluto would start this cell automatically and does not support stopping cells on windows. To find the estimates change the number of iterations from 10 to 5000 !!!

sol = optimize(x->criterion(x,rand_vec),params_init,NelderMead(),Optim.Options(
		       iterations = 10))
res = Optim.minimizer(sol)   # vector containing the solutions
	
end;

# ╔═╡ 9b7b0c70-38c8-11ec-1484-1b1f791d99dd
begin

# Since the estimation procedure can take a lot of time and depends on the rand_vec and the starting guesses, here we provide the data to replicate our best estimates. The file estimates2 contains the estimates found by the authors, our best estimates and the initial parameters we used to found them. random_vector_used contains the rand_vec we used to compute our best estimates

data, header = readdlm(raw"../input/estimates2.csv",',',header=true)
best_estimates = DataFrame(data, vec(header))

be = best_estimates[:,2]                # our best estimates

data, header = readdlm(raw"../input/random_vector_used.csv",',',header=true)
best_random_vec = DataFrame(data, vec(header))
rand_vec_used = [best_random_vec[:,1],best_random_vec[:,2],best_random_vec[:,3],best_random_vec[:,4]]
params_init_used = best_estimates[:,3]  # The random draws used to compute the best                                             estimates
	
end;

# ╔═╡ 5b42d010-38c9-11ec-1474-e103040023f8
begin

# This cell computes our best estimates given the rand_vec and the initial parameters used

# !!! change the number of iterations to 10.000 to find our best estimates !!!
	
sol = optimize(x->criterion(x, rand_vec_used), params_init_used, NelderMead(),
Optim.Options(iterations = 10)) # We set 10.000 iterations even if the algorithm will                                   converge before 10.000 iterations
be = Optim.minimizer(sol)
	
end

# ╔═╡ c8ce3502-38c6-11ec-3b68-f1513e5f03b5
begin

# Compare our best estimates with the ones obtained by the authors

# These are the authors' estimates
parameters_authors = [0.38,0.36,0.38,0.30,0.14,0.16,-22.6,-27.7,26.9,24.7,1.6,1.2,42.7,23.9,-3.9,-11.3,9.5,7.6,64.1,318.7] 
	
# This difference is positive if the wsse obtained using our estimates is lower

dif = criterion(parameters_authors, rand_vec) - criterion(be, rand_vec)
	
end

# ╔═╡ 3f8c7950-38cb-11ec-1cf1-e79ce14eb776
# Below we compute the Jacobian using ForwardDiff but we do not use this jacobian to compute the standard errors since we get two columns of zeroes for the gradient with respect to mu_eps and sigma_eps and so the matrix is not invertible. This problem is either due to (1) mistakes in the function we wrote to compute the simulated moments or (2) the way in which julia computes the jacobian. Regarding the former possibility: using the same parameters and the same random vectors as the authors, we do see very small differences between ours and the authors' simulated moments (the largest difference being only 0.002); after inspecting it carefully, we cannot attribute these small discrepancies to mistakes in the function we wrote but it is possible that they are due to different default options between julia and matlab's operators (e.g., in rounding). Regarding the latter possibility: in julia, we computed the jacobian using the ForwardDiff Pkg while the authors used the Jacobianest function in Matlab; from a comment in the documentation of this Matlab function,  we learnt that "the error term on finite differences jacobian's estimates has a second order component, but also some 4th and 6th order terms in it" and the Jacobianest function uses Romberg extrapolation to improve the estimate, but ForwardDiff uses automatic differentiation so it should return a very precise jacobian.

jac_julia = ForwardDiff.jacobian(x->voteSimEndogenousVoting_vary(x, rand_vec_used), be);

# ╔═╡ 1f237b10-3e1f-11ec-2dde-092bc72714fb
# We load the jacobian evaluated in our best estimates computed in Matlab using the "Jacobianest" function (the same used by the authors).

jac_matlab = readdlm(raw"../input/jac_mat.csv",',', header=false);

# ╔═╡ 7a5da8fe-3b00-11ec-1f2c-7fc389373302
begin

# Compute the standard errors.
	
#For more informations on the formula for the standard errors we refer the reader to part B of the appendix or to the python notebook.
	
# We use the same notation as in their matlab code. sim_adjust is a scalar, DFDY_CSD_jacest is our jacobian (101x20 matrix) evaluated in the minimum, W is the weighting matrix (101x101 diagonal matrix), VCcontrol is the variance-covariance matrix of the empirical moments (a 101x101 matrix).

Jm_Js = 13197 / sim_voters # nr. empirical obs / nr. simulated obs
sim_adjust = 1 + Jm_Js 
DFDY_CSD_jacest = jac_matlab 
A = DFDY_CSD_jacest'*W*DFDY_CSD_jacest
B = DFDY_CSD_jacest'*W*(sim_adjust.*VCcontrol[1:101,1:101])*W*DFDY_CSD_jacest
VC = A\B/A
standard_errors = diag(VC).^0.5 # se are the square root of the main diagonal
	
end;

# ╔═╡ 427cce6e-37fa-11ec-2a6f-edc5c54a32aa
begin

# Create and save a csv file containing our estimates and their standard errors

# parameters
params = ["h0_v","h0_nv","r_v","r_nv","eta_v","eta_nv","mu_s_v","mu_s_nv","sigma_s_v","sigma_s_nv","S_svy_v","S_svy_nv","timeval_v","timeval_nv","mu_sv","mu_sn","sigma_si","L_in","mu_eps","sigma_eps"]
	
# authors' standard errors
se_authors= [0.0089,0.0092,0.0204,0.0183,0.1232,0.2987,2.9580,5.5938,5.5176,6.6687,1.2084,1.6149,9.4438,15.1539,1.4858,1.6422, 2.9750,1.5247,61.9195,248.6525]
	
table3=DataFrame()
table3.parameters = params
table3.our_point_est=round.(be, digits=2)
table3.our_se=round.(standard_errors, digits=4)
table3.authors_point_est=parameters_authors
table3.authors_se = se_authors

CSV.write(raw"../output/estimates_julia.csv", table3)
table3
	
end

# ╔═╡ Cell order:
# ╠═9c0d53a0-38c3-11ec-0d4d-1d549264f0f0
# ╠═92cfa5a2-5ce9-11ec-0e3e-61fe154ce557
# ╠═b34d84a0-5ce9-11ec-1316-39f597588bef
# ╠═73767a00-5db9-11ec-1bef-b9e6f6e711cc
# ╠═15e16bb0-3560-11ec-238e-e3018fb17cf9
# ╠═c6a68930-3588-11ec-08e1-7d0922b9816b
# ╠═c6731d20-3588-11ec-3127-51520038c40c
# ╠═c6468ee0-3588-11ec-0d41-cbfce972d21c
# ╠═c616cc50-3588-11ec-25bc-1b52b9c8694d
# ╠═1bb22460-37c4-11ec-3897-732e028fc413
# ╠═e33132a0-37c5-11ec-1079-e3c9bba57502
# ╠═9b7b0c70-38c8-11ec-1484-1b1f791d99dd
# ╠═5b42d010-38c9-11ec-1474-e103040023f8
# ╠═c8ce3502-38c6-11ec-3b68-f1513e5f03b5
# ╠═3f8c7950-38cb-11ec-1cf1-e79ce14eb776
# ╠═1f237b10-3e1f-11ec-2dde-092bc72714fb
# ╠═7a5da8fe-3b00-11ec-1f2c-7fc389373302
# ╠═427cce6e-37fa-11ec-2a6f-edc5c54a32aa
