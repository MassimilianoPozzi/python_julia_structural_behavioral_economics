
To reproduce the tables and figures in DellaVigna, List, Malmendier, and Rao 2016, use the code in this folder. 

The code in the “stata” folder produces the reduced form results and the moments used as inputs for the estimation code in the “benchmark” and “same” folder. These are Moments.mat for the benchmark moments and Moments_R_nodrop.mat for the robustness estimation using empirical moments calculated without dropping households that solicitors were not able to reach (e.g., had no solicitor sign or house was for sale).

Moments.mat created from samplemain_ctrl_moments1, samplemain_ctrl_moments2, and samplemain_ctrl_moments3. 

Moments_R_nodrop.mat created from samplenodrop_ctrl_moments1, samplenodrop_ctrl_moments2, and samplenodrop_ctrl_moments3. 

Baseline turnout and GOTV moments are added to the empirical moments in the Matlab code using the coefficient on the constant in column 1 of Table 6 and the GOTV effect (coefficient on indicator for receiving a flyer with an announcement about voting) in column 2 of Table 6. 

See below for an ordered list of the moments. 


STATA FILES
Input data files includes: 
ElectionData2011VotingHistoryCleanDec15 - voting record for households in the 2011 sample
exp2011Dec15 - 2011 data from the main field experiment
GOTVData2010Dec15 - 2010 GOTV data
GOTVData2012Dec15 - 2012 GOTV data

Stata do-file “electionMain” runs all reduced form analyses and generates the empirical moments for model estimation.  



MATLAB - model estimation in “benchmark” and “same” folders

The Matlab files take as inputs the empirical moments and variance covariance matrix from the Stata output (saved as Moments.mat and Moments_R_nodrop.mat).

The files in the “benchmark” folder produce estimations for the structural model allowing auxiliary parameters to differ across voters and non-voters. 
The files in the “same” folder do the corresponding estimations when auxiliary parameters are the same for voters and non-voters. 
The “same” folder also includes code for estimating the model with exogenous assignment into voters and non-voters. 

Files starting with “est” estimate the parameters that best fit the empirical moments. This includes the benchmark model “est_benchmark.m” as well as all alternative models considered in the robustness tables.

Robustness analyses may differ from the benchmark model by using a combination of alternative options set in the estimation code, alternative “voteSim” files, alternative “calc_implications” files, alternative “getSearchInits” files, and alternative empirical moments.
Robustness analyses include:
“hetL” - assume an exponential distribution of L
“UTalk" - allow for a utility of talking about politics
“meas_v10” - assume 10% of voters were mismeasured as non-voters
“halfN” - assume people are asked about voting half as much as benchmark
“fix_sige” - fixing sigma epsilon to a specific value
“meas_v20”- assume 20% of voters were mismeasured as non-voters
“meas_vnv10” - assume 10% of voters and non-voters were mismeasured
“2N” - assume people are asked about voting twice as often as benchmark
“nodrop” - use empirical moments calculated without dropping households that were not able to be reached
“GOTV” - estimate the model using the GOTV turnout moment
“noINI” - estimate the model without the I/NI moments
“nolie_fix_sigsvsn” - estimate the model without the lying incentive moments, fixing sigma sv/sn
“non_endog_v” - estimate the exogenous assignment model for voters only
“non_endog_nv” - estimate the exogenous assignment model for non-voters only
"best" - estimates the model using as starting points the best estimates from another estimation (used for the estimates in Column 2 of Table 3 and Columns 1 and 8 of Online Appendix Table 6). 


All other files are supporting functions, called when running the estimation files:

“voteSim” files take parameter values as inputs and outputs the simulated moments at those inputs. See “voteSimEndogenousVoting_vary.m” for detailed comments.

"minSearch" files perform the actual numerical minimization of the objective function to identify the minimum distance estimate.

"getSearchInits" files choose starting points for the minimization routine. This is usually done randomly within the starting range. 
The "getSearchInits_best" file chooses search inits using best estimates from other estimations. 
There are two corresponding .mat files in the "same" folder used to supply the start points: 
"fixsige490.mat" includes the best estimates from the estimation fixing sigma-epsilon to 490.6 
"bench_best.mat" includes the best estimates from the estimation of the benchmark model (with variable sigma-epsilon) using as start points the best estimates fixing sigma-epsilon to 490.6

"calc_implications" files calculate the following implications at the supplied parameter values (in this order). See “calc_implications.m” in the benchmark folder for more detail. 
	Mean Epsilon V then NV
	Mean Sv V then NV
	Mean Sn V then NV
	Mean value of voting to tell others V then NV
	Baseline Turnout
	Implied Change in Turnout if Never Asked About Voting
	Implied Change in Asked About Voting Twice as Often
	Mean Utility being Asked V then NV
	Implied GOTV Effect (Assuming N+1)
	Number people need to ask to get 1 more vote (N+1)
	Utility cost to get 1 more vote (N+1)
	Turnout if N=0...12 (Fig 9)
	Mean lying cost all sample, V, then NV (if distribution of L)

The implications reported for the exogenous voter assignment model are slightly different. There are calculated in the file ”calc_implications_vnv_sep" and include the following:
	Mean Value of saying voted
	Mean Value of saying didn't vote
	Implied value of voting to tell others (L=0,2,5,10)  
	Utility from being asked once

The Matlab files also include two supporting functions, written by John D'Errico, available at the Matlab File Exchange online. 
“fminsearchbnd” allows us to run fminsearch with bounds by transforming the parameter space.
“jacobianest” is used to numerically calculate derivatives for standard error calculations.  


Order of moments:
P(Answer Door), Voters, No Flyer, $0, 5minP(Answer Door), Voters, No Flyer, $10, 10 minP(Answer Door), Voters, No Flyer, $10, 5minP(Answer Door), Voters, Survey Flyer, $0, 5minP(Answer Door), Voters, Survey Flyer, $10, 10 minP(Answer Door), Voters, Survey Flyer, $10, 5minP(Answer Door), Voters, Election Flyer, $0, 5minP(Answer Door), Voters, Election Flyer, $10, 10 minP(Answer Door), Voters, Election Flyer, $10, 5minP(Answer Door), Voters, Survey Opt Out, $0, 5minP(Answer Door), Voters, Survey Opt Out, $10, 10 minP(Answer Door), Voters, Survey Opt Out, $10, 5minP(Answer Door), Voters, Election Opt Out, $0, 5minP(Answer Door), Voters, Election Opt Out, $10, 10 minP(Answer Door), Voters, Election Opt Out, $10, 5minP(Answer Door), Non-Voters, No Flyer, $0, 5minP(Answer Door), Non-Voters, No Flyer, $10, 10 minP(Answer Door), Non-Voters, No Flyer, $10, 5minP(Answer Door), Non-Voters, Survey Flyer, $0, 5minP(Answer Door), Non-Voters, Survey Flyer, $10, 10 minP(Answer Door), Non-Voters, Survey Flyer, $10, 5minP(Answer Door), Non-Voters, Election Flyer, $0, 5minP(Answer Door), Non-Voters, Election Flyer, $10, 10 minP(Answer Door), Non-Voters, Election Flyer, $10, 5minP(Answer Door), Non-Voters, Survey Opt Out, $0, 5minP(Answer Door), Non-Voters, Survey Opt Out, $10, 10 minP(Answer Door), Non-Voters, Survey Opt Out, $10, 5minP(Answer Door), Non-Voters, Election Opt Out, $0, 5minP(Answer Door), Non-Voters, Election Opt Out, $10, 10 minP(Answer Door), Non-Voters, Election Opt Out, $10, 5minP(Do Survey), Voters, No Flyer, $0, 5minP(Do Survey), Voters, No Flyer, $10, 10 minP(Do Survey), Voters, No Flyer, $10, 5minP(Do Survey), Voters, Survey Flyer, $0, 5minP(Do Survey), Voters, Survey Flyer, $10, 10 minP(Do Survey), Voters, Survey Flyer, $10, 5minP(Do Survey), Voters, Election Flyer, $0, 5minP(Do Survey), Voters, Election Flyer, $10, 10 minP(Do Survey), Voters, Election Flyer, $10, 5minP(Do Survey), Voters, Survey Opt Out, $0, 5minP(Do Survey), Voters, Survey Opt Out, $10, 10 minP(Do Survey), Voters, Survey Opt Out, $10, 5minP(Do Survey), Voters, Election Opt Out, $0, 5minP(Do Survey), Voters, Election Opt Out, $10, 10 minP(Do Survey), Voters, Election Opt Out, $10, 5minP(Do Survey), Non-Voters, No Flyer, $0, 5minP(Do Survey), Non-Voters, No Flyer, $10, 10 minP(Do Survey), Non-Voters, No Flyer, $10, 5minP(Do Survey), Non-Voters, Survey Flyer, $0, 5minP(Do Survey), Non-Voters, Survey Flyer, $10, 10 minP(Do Survey), Non-Voters, Survey Flyer, $10, 5minP(Do Survey), Non-Voters, Election Flyer, $0, 5minP(Do Survey), Non-Voters, Election Flyer, $10, 10 minP(Do Survey), Non-Voters, Election Flyer, $10, 5minP(Do Survey), Non-Voters, Survey Opt Out, $0, 5minP(Do Survey), Non-Voters, Survey Opt Out, $10, 10 minP(Do Survey), Non-Voters, Survey Opt Out, $10, 5minP(Do Survey), Non-Voters, Election Opt Out, $0, 5minP(Do Survey), Non-Voters, Election Opt Out, $10, 10 minP(Do Survey), Non-Voters, Election Opt Out, $10, 5minP(Opt Out), Voters, Opt Out, $0, 5minP(Opt Out), Voters, Opt Out, $10, 10 minP(Opt Out), Voters, Opt Out, $10, 5minP(Opt Out), Voters, Election Opt Out, $0, 5minP(Opt Out), Voters, Election Opt Out, $10, 10 minP(Opt Out), Voters, Election Opt Out, $10, 5minP(Opt Out), Non-Voters, Opt Out, $0, 5minP(Opt Out), Non-Voters, Opt Out, $10, 10 minP(Opt Out), Non-Voters, Opt Out, $10, 5minP(Opt Out), Non-Voters, Election Opt Out, $0, 5minP(Opt Out), Non-Voters, Election Opt Out, $10, 10 minP(Opt Out), Non-Voters, Election Opt Out, $10, 5minP(Do Survey), Voters, Not Informed, No FlyerP(Do Survey), Voters, Informed, No FlyerP(Do Survey), Voters, Not Informed, Survey FlyerP(Do Survey), Voters, Informed, Survey FlyerP(Do Survey), Voters, Not Informed, Election FlyerP(Do Survey), Voters, Informed, Election FlyerP(Do Survey), Voters, Not Informed, Opt OutP(Do Survey), Voters, Informed, Opt OutP(Do Survey), Voters, Not Informed, Election Opt OutP(Do Survey), Voters, Informed, Election Opt OutP(Do Survey), Non-Voters, Not Informed, No FlyerP(Do Survey), Non-Voters, Informed, No FlyerP(Do Survey), Non-Voters, Not Informed, Survey FlyerP(Do Survey), Non-Voters, Informed, Survey FlyerP(Do Survey), Non-Voters, Not Informed, Election FlyerP(Do Survey), Non-Voters, Informed, Election FlyerP(Do Survey), Non-Voters, Not Informed, Opt OutP(Do Survey), Non-Voters, Informed, Opt OutP(Do Survey), Non-Voters, Not Informed, Election Opt OutP(Do Survey), Non-Voters, Informed, Election Opt OutP(Lie|Survey), Voters, 5 min survey, controlP(Lie|Survey), Voters, 5 min survey, $5 incentiveP(Lie|Survey), Voters, 10 min survey, controlP(Lie|Survey), Voters, 10 min survey, 8 min incentiveP(Lie|Survey), Non-Voters, 5 min survey, controlP(Lie|Survey), Non-Voters, 5 min survey, $5 incentiveP(Lie|Survey), Non-Voters, 10 min survey, controlP(Lie|Survey), Non-Voters, 10 min survey, 8 min incentiveTurnout, baselineTurnout, GOTV intervention


MATLAB - "eratgneezy" folder includes all of the code necessary to generate the moments and estimates on Erat and Gneezy (2012) for Online Appendix Table 8.