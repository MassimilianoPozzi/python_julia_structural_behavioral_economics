
//change this!
cd "PUT YOUR FOLDER HERE"
set more off
cap log close 
log using Temp/log.txt, replace
version 11
do 01IndividualNonPara.do
************
** Constructs individual non-parametric measures of present bias and sophistication
** 
** INPUT:  
** - Data/decisions_data
** 
** OUTPUT:		
** - Data/decisions_data_w_non_para	
************

do 02individualMLE.do
************
** Find individual estimates which allows the identification of the main participant sample for the paper  
** This has to happen before all analysis because we limit to sample that can converge throughout main analysis 
** For all of the estimation strategies, this produces estimates for beta, beta-hat, etc.
** Note: This takes a while.
**
** OUTPUT:		
** - Data/ind_results
** - Note: I have the end result: ind_results file in package so that you can skip this if you want!
************

do 03MergeIndMLEAndConstructMainSample.do
**********
** Looks at individual estimations and finds issues/outliers, etc.
** Constructs main sample 
**
** INPUT:  
** - Data/decisions_data_w_non_para	
** - Data/ind_results
** 
** OUTPUT:		
** - Data/decisions_data_w_ind	
************

//now, all programs use Data/decisions_data_w_ind as the input
do 04stat_and_graphs.do
**********
** Constructs referenced statistics in the paper
** Constructs figures
************

do 05agggregate.do
**********
** Does a ton of aggregate ML estimations- some used in paper, some in appendix, some unused
** Creates tables from stored estimates
** Takes a while
************


do 06individualTable.do
**********
** Creates individual estimate table from estimates from 02individualMLE.do
************


do 07indBetaBetaHatTable.do
**********
** Creates individual beta/betahat correlation table from estimates from 02individualMLE.do
************


do 10app_learning.do
**********
** For Appendix: Stability of parameters across the experiment
************


do 11app_params_changing.do
**********
** For Appendix: Stability of parameters across the experiment
************


do 12app_reactiontobonus.do
**********
** For Appendix: Internal-Consistency Rewards
************


do 13app_plotsOfDroppedParticipants.do
**********
** For Appendix: Participants Not Included in the Main Sample - scatterplot
************


do 14app_effectofotherwages.do
**********
** Showing (null) effect of different wage options in a decision set: not discussed in detail in paper
************


