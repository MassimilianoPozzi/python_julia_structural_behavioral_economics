do setupAndDefinePrograms.do


**********************
******** Non-Monotonicities
**********************

use "Data/decisions_data_w_ind", clear

cap drop d_w_p p_non_mono p_opportunities

egen d_w_p=group(decisiondate workdate wid type)
bys d_w_p (wage): gen t_non_mono=(jobschosen<jobschosen[_n-1])  
bys d_w_p (wage): replace t_non_mono=. if _n==1
bys wid: egen p_non_mono=total(t_non_mono) 
gen p_more_than_5=(p_non_mono>5)

//QUOTED IN PAPER: On an individual level, 66 of the 72 participants have less than 5 total non-monotonicities...
tab p_more_than_5 if first==1 & $wid_sample2
//QUOTED IN PAPER: (given a total of 104 violation opportunities in adjacent decisions for the full experiment) = 4*26=104





**************************
******* Aggregate Non Parameteric Main Graphs comparing now vs. later and later vs. prediction
**************************


***********
** Figure 2: Present vs. future decisions (left) and predictions vs. future decisions (right)
***********

foreach s of numlist 2 6 { //s=2 is the main sample with MLE-problem-subjects removed, s=6 is everyone (in Appendix)
foreach wageGroups in "wage_group" "wage"  { //we use wageGroup=group wages by 3 in paper. Raw wage graph show in the Appendix for transparency

	local suffix=""
	local suffix2=""
	if ("`wageGroups'"=="wage") {
		local suffix="_all_wages"
	}
	if ("`s'"=="6") {
		local suffix2="_6"
	}
	
	use "Data/decisions_data_w_ind", clear
	cap drop t_*
	drop if wage>.315
	egen group = group(`wageGroups' type)
	xi: reg jobschosen i.group if ($data_sample1) & (${wid_sample`s'}), clu(wid)
	
	//could also do two-way clustering - all the same essentially. If anything, reduces SE
	//xi: ivreg2 jobschosen i.group if ($data_sample1) & (${wid_sample`s'}), clu(wid workdate)
	
	capture drop mean_jc
	predict mean_jc
	capture drop se_jc
	predict se_jc, stdp
	capture drop jc_low
	gen jc_low = mean_jc - 1.96*se_jc
	capture drop jc_high
	gen jc_high = mean_jc + 1.96*se_jc

	twoway (rcap jc_low jc_high `wageGroups' if type=="present", color(gs10)) ///
	(connected mean_jc `wageGroups' if type=="present", sort color(gs8))  ///
	(rcap jc_low jc_high `wageGroups' if type=="future", color(gs10) ) ///
	(connected mean_jc `wageGroups' if type=="future", lpattern(solid) lcolor(gs1) sort) ///
	if (${wid_sample`s'}) ///
	, legend(order(2 "Present Tasks" 4 "Future Tasks" 1 "95% CI") row(1) ) ///
	 xtitle("Wage") ytitle("Jobs Chosen") ylab(0(10)80) xlab(0(.05).31)  scheme(s1mono)
	graph save Output/t1, replace
	graph export "$savedir/nonpara_agg_pres_2`suffix'.png", replace width(2000)
	
	twoway (rcap jc_low jc_high `wageGroups' if type=="prediction", color(gs10)) ///
	(connected mean_jc `wageGroups' if type=="prediction", sort color(gs8))  ///
	(rcap jc_low jc_high `wageGroups' if type=="future", color(gs10) ) ///
	(connected mean_jc `wageGroups' if type=="future", lpattern(solid) lcolor(gs1) sort) ///
	if (${wid_sample`s'}) ///
	, legend(order(2 "Predicted Tasks" 4 "Future Tasks" 1 "95% CI") row(1) ) ///
	 xtitle("Wage") ytitle("Jobs Chosen") ylab(0(10)80) xlab(0(.05).31)  scheme(s1mono)
	graph save Output/t2, replace
	graph export "$savedir/nonpara_agg_pred_2`suffix'`suffix2'.png", replace width(2000)
	
	graph combine Output/t1.gph Output/t2.gph, xcommon ycommon ysize(2) altshrink 
	graph export "$savedir/nonpara_agg_2`suffix'`suffix2'.png", replace width(2000)
	
}
}






*****Density functions: In appendix.

use "Data/decisions_data_w_ind", clear
keep if ($data_sample1) & (${wid_sample2})
xi: reg jobschosen i.wage i.wid, 
predict res, res

cumul res if type=="present", gen(cum_pres)
cumul res if type=="future", gen(cum_future)
cumul res if type=="prediction", gen(cum_pred)

twoway (line cum_pred res, sort lcolor(gs10)) (line cum_future res, sort lcolor(gs1)) if res>-60 & res<40, ///
legend(order(1 "Predicted Tasks" 2 "Future Tasks") row(1) ) ///
xtitle("Residual Jobs Chosen") ytitle("Cumulative Probability") xlab(-60(20)40)  scheme(s1mono)
graph save Output/t1, replace 

twoway (line cum_pres res, sort lcolor(gs10)) (line cum_future res, sort lcolor(gs1)) if res>-60 & res<40, ///
legend(order(1 "Present Tasks" 2 "Future Tasks") row(1) ) ///
xtitle("Residual Jobs Chosen") ytitle("Cumulative Probability") xlab(-60(20)40)  scheme(s1mono)
graph save Output/t2, replace
	 
graph combine Output/t2.gph Output/t1.gph, xcommon ycommon ysize(2) altshrink 
graph export "$savedir/cumulative.png", replace width(2000)






*************
** Statistical tests for differences
*************

use "Data/decisions_data_w_ind", clear	

//NOW VS LATER
//USED -> basic
xi: reg jobschosen future if (type=="present" | type=="future")    & ($data_sample1) & (${wid_sample2}), clu(wid)
//USED -> decisiondate>2
xi: reg jobschosen future if (type=="present" | type=="future")  & decisiondatenum>2  & ($data_sample1) & (${wid_sample2}), clu(wid) 
//USED -> controlling for wage
xi: areg jobschosen future if (type=="present" | type=="future")    & ($data_sample1) & (${wid_sample2}), clu(wid) absorb(wage)

//PRED VS LATER	
//USED -> basic
xi: reg jobschosen future if (type=="prediction" | type=="future") & ($data_sample1) & (${wid_sample2}), clu(wid)
//USED -> decisiondate>2
xi: reg jobschosen future if (type=="prediction" | type=="future")  & decisiondatenum>2 & ($data_sample1) & (${wid_sample2}), clu(wid)
//USED -> controlling for wage
xi: areg jobschosen future if (type=="prediction" | type=="future") & ($data_sample1) & (${wid_sample2}), clu(wid) absorb(wage)

	
//PRED VS NOW
//USED -> basic
xi: reg jobschosen prediction if (type=="prediction" | type=="present") & ($data_sample1) & (${wid_sample2}), clu(wid)
//USED -> decisiondate>2
xi: reg jobschosen prediction if (type=="prediction" | type=="present")  & decisiondatenum>2 & ($data_sample1) & (${wid_sample2}), clu(wid)
//USED -> controlling for wage
xi: areg jobschosen prediction if (type=="prediction" | type=="present") & ($data_sample1) & (${wid_sample2}), clu(wid) absorb(wage)










**********************
******** Exploring Quasi-hyperbolic
**********************

***********
** Figure 3: Work decisions given different delays to the time of work
***********

foreach s of numlist 2 6 { //s=2 is the main sample with MLE-problem-subjects removed, s=6 is everyone (in Appendix)
forvalues a=0(1)1 { //this is whether we use days (4 days) or participation dates (1 date)

	use "Data/decisions_data_w_ind", clear
	
	drop if (type=="prediction")
	
	if (`a'==1) {
			gen days_to_work=workdate-decisiondate
			replace days_to_work=21 if days_to_work>21
			//note: one decision set violated rules (4 days) - occured because someone had limited time.
			drop if days_to_work<4 & days_to_work>0 
		} 
	else {
		gen days_to_work=workdatenum-decisiondatenum
	}
	
	gen days_to_work2=days_to_work
	local barwidth=3.5
	forvalues i=4(1)21 {
		replace days_to_work2=1.5+ceil((`i'-3)/4)*4 if days_to_work==`i'
	} 
	tab days_to_work2, gen(days_dummies_)
	
	reg jobschosen days_dummies_*  if ($data_sample1) & (${wid_sample`s'}), clu(wid) noconstant
	
	capture drop mean_jc
	predict mean_jc
	qui sum mean_jc if days_to_work>0
	local meanline=r(mean)
	capture drop se_jc
	predict se_jc, stdp
	capture drop jc_low
	gen jc_low = mean_jc - 1.96*se_jc
	capture drop jc_high
	gen jc_high = mean_jc + 1.96*se_jc

	if (`a'==1) {
		twoway  ///
		(bar mean_jc days_to_work2, sort color(gs8) barwidth(`barwidth'))  ///
		(rcap jc_low jc_high days_to_work2, color(gs1)) ///
		if (${wid_sample`s'}) ///
		, legend(off) ///
		xtitle("Days Until Work Occurs") ytitle("Jobs Chosen") ///
		xlab(0 "0"  5.5 "4-7"  9.5 "8-11" 13.5 "12-15" 17.5 "16-19" 21.5 "20-30") ///
		ylab(40(5)55) yline(`meanline', lwidth(vthin)) scheme(s1mono)
	}
	else {
		twoway  ///
		(bar mean_jc days_to_work2, sort color(gs8) barwidth(.8))  ///
		(rcap jc_low jc_high days_to_work2, color(gs1)) ///
		if (${wid_sample`s'}) ///
		, legend(off) ///
		xtitle("Participation Dates Until Work Occurs") ytitle("Jobs Chosen") xlab(0(1)3) ylab(40(5)55) yline(`meanline', lwidth(vthin)) scheme(s1mono)
	}
	cap mkdir "$direct/Output/hyperbolic"
	graph save "Output/hyperbolic/Hyperbolic`a'_`s'.gph", replace
	graph export "Output/hyperbolic/Hyperbolic`a'_`s'.png", replace width(2000)

}
}

//this is the main sample
graph combine "Output/hyperbolic/Hyperbolic0_2" "Output/hyperbolic/Hyperbolic1_2", xsize(10) ycommon altshrink 
graph export "$savedir/nonpara_hyperbolic_all.png", replace width(2000)

//this is everyone
graph combine "Output/hyperbolic/Hyperbolic0_6" "Output/hyperbolic/Hyperbolic1_6", xsize(10) ycommon altshrink 
graph export "$savedir/nonpara_hyperbolic_expanded.png", replace width(2000)





*************
** Statistical tests for quasi-hyperbolic.
*************

//for testing participation date differences to get Z scores
//used in paper after graphs above are shown.

use "Data/decisions_data_w_ind", clear
drop if (type=="prediction")
	
gen days_to_work=workdatenum-decisiondatenum
gen days_to_work2=days_to_work
tab days_to_work2, gen(days_dummies_)
	
//USED IN PAPER to get quasi-hyperbolic differences:
reg jobschosen days_dummies_2 days_dummies_3 days_dummies_4  if ($data_sample1) & (${wid_sample2}), clu(wid) 
//USED IN PAPER this is just used to see the highest piecewise signficance between other days:
reg jobschosen days_dummies_1 days_dummies_3 days_dummies_4  if ($data_sample1) & (${wid_sample2}), clu(wid) 
reg jobschosen days_dummies_1 days_dummies_2 days_dummies_4  if ($data_sample1) & (${wid_sample2}), clu(wid) 
reg jobschosen days_dummies_1 days_dummies_2 days_dummies_3  if ($data_sample1) & (${wid_sample2}), clu(wid) 

//USED IN PAPER (footnote) controlling for wage
areg jobschosen days_dummies_2 days_dummies_3 days_dummies_4  if ($data_sample1) & (${wid_sample2}), clu(wid) absorb(wage)

//USED IN PAPER (footnote) after decision date 2
reg jobschosen days_dummies_2 days_dummies_3 days_dummies_4  if ($data_sample1) & (${wid_sample2}) & decisiondatenum>2, clu(wid)




//for testing days-away-from-present differences to get Z scores
//this is referenced but not used explicitely in paper
use "Data/decisions_data_w_ind", clear
drop if (type=="prediction")

gen days_to_work=workdate-decisiondate
//one decision set violated rules (4 days) - occured because someone had limited time.
replace days_to_work=21 if days_to_work>21
drop if days_to_work<4 & days_to_work>0 
gen days_to_work2=days_to_work

forvalues i=4(1)21 {
	replace days_to_work2=1.5+ceil((`i'-3)/4)*4 if days_to_work==`i'
}
tab days_to_work2, gen(days_dummies_)
	
reg jobschosen days_dummies_2 days_dummies_3 days_dummies_4 days_dummies_5 days_dummies_6 if ($data_sample1) & (${wid_sample2}), clu(wid) 
reg jobschosen days_dummies_1 days_dummies_3 days_dummies_4 days_dummies_5 days_dummies_6  if ($data_sample1) & (${wid_sample2}), clu(wid) 
reg jobschosen days_dummies_1 days_dummies_2 days_dummies_4 days_dummies_5 days_dummies_6 if ($data_sample1) & (${wid_sample2}), clu(wid) 
reg jobschosen days_dummies_1 days_dummies_2 days_dummies_3 days_dummies_5 days_dummies_6 if ($data_sample1) & (${wid_sample2}), clu(wid) 
	




	
	
	
**************************
******** Individual Non Parameteric Reallocation - this is used below and combined with scatter to create Figure 4
**************************

foreach s of numlist 2 6 { //s=2 is the main sample with MLE-problem-subjects removed, s=6 is everyone (in Appendix)

	use "Data/decisions_data_w_ind", clear

	**** Histograms of participant-level differences between different types of decisions

	hist ind_diff_pres if first==1 & ${wid_sample`s'} & ind_diff_pres>-40 & ind_diff_pres<20, width(2.2) ///
	xtitle("Present Bias") xline(0, lpattern(dash) lwidth(1)) xsc(r(-30 10)) xlabel(-30(10)10)  ysc(r(0 .15)) ylabel(0(.05).15)  scheme(s2mono) 
	graph save Output/histnonpara1_`s', replace
	graph export "$savedir/nonpara_indhist_pres_`s'.png", replace width(2000)

	hist ind_diff_pred if first==1 & ${wid_sample`s'}, width(2.2) ///
	xtitle("Present Bias Perception") xline(0, lpattern(dash) lwidth(1)) xsc(r(-30 10)) xlabel(-30(10)10)  ysc(r(0 .15)) ylabel(0(.05).15) scheme(s2mono)
	graph save Output/histnonpara2_`s', replace
	graph export "$savedir/nonpara_indhist_pred_`s'.png", replace width(2000)

	graph combine Output/histnonpara1_`s'.gph Output/histnonpara2_`s'.gph, ysize(1) iscale(1.5)
	graph export "$savedir/nonpara_indhist_`s'.png", replace width(2000)
}
		




	


****** Statistics
use "Data/decisions_data_w_ind", clear
local cond="(first==1 & ml_prob==0)"
gen ind_diff_pres_l0=(ind_diff_pres<0)
gen ind_diff_pred_l0=(ind_diff_pred<0)

*** USED: number less than 0:
sum ind_diff_pres ind_diff_pred if `cond', detail
sum ind_diff_pres_l0 ind_diff_pred_l0 if `cond', detail

*** USED: SD of estimates
sdtest ind_diff_pres== ind_diff_pred if `cond'

*** USED: corr: 
pwcorr ind_diff_pres ind_diff_pred if `cond', sig

*** Not reported: spearman
spearman ind_diff_pres ind_diff_pred if `cond', stats(rho p)

*** Corr of bhi for comparison 
pwcorr bi1 bhi1 if `cond', sig

*** correlations
corr ind_diff_pres ind_diff_pred if  `cond'

***calculating the analagous "lambga" relative sophistication parameter
//gets significance for different robust regression techniques
reg ind_diff_pred ind_diff_pres  if  `cond'
robreg m   ind_diff_pred ind_diff_pres if  `cond'
robreg mm  ind_diff_pred ind_diff_pres if  `cond'
robreg s   ind_diff_pred ind_diff_pres if  `cond'
reg ind_diff_pred ind_diff_pres  [aweight=ind_diff_pred_var_inv ] if `cond', r nocons
eivreg ind_diff_pred ind_diff_pres if `cond' , r(ind_diff_pres .5)
//let's get reliability
sum ind_diff_pres_var if `cond'
local noiseVar=r(mean)
sum ind_diff_pres if `cond'
local totalVar=r(sd)^(2)
local reliability=(`totalVar'-`noiseVar')/(`totalVar')
eivreg ind_diff_pred ind_diff_pres if `cond' , r(ind_diff_pres `reliability')

*** In paper: are non parametric related to parametric 
pwcorr ind_diff_pres bi1 if first==1 & ml_prob1==0, sig
pwcorr ind_diff_pred bhi1 if first==1 & ml_prob1==0, sig







*************************************************
********Prediction-as-Commitment
*************************************************

**********
** The variation of predictions given different bonus amounts: Figure 5 in the paper
**********
foreach s of numlist 2 6 { //s=2 is the main sample with MLE-problem-subjects removed, s=6 is everyone (in Appendix)
	use "Data/decisions_data_w_ind", clear
	keep if (type=="prediction" | type=="future") 
	drop if wage>.315

	cap drop guessbonushigh
	gen guessbonushigh=0 if type=="future"
	replace guessbonushigh=1 if guessbonusamount<=1.01 & type=="prediction"
	replace guessbonushigh=2 if guessbonusamount>1.01 & guessbonusamount<=4 & type=="prediction"
	replace guessbonushigh=3 if guessbonusamount>4 & type=="prediction"
	egen group = group(wage_group guessbonushigh)
	xi: reg jobschosen i.group if ($data_sample1) & (${wid_sample`s'}), clu(wid)

	capture drop mean_jc
	predict mean_jc
	capture drop se_jc
	predict se_jc, stdp
	capture drop jc_low
	gen jc_low = mean_jc - 1.96*se_jc
	capture drop jc_high
	gen jc_high = mean_jc + 1.96*se_jc

	twoway (rcap jc_low jc_high wage_group if guessbonushigh==1, color(gs10) ) ///
	(connected mean_jc wage_group if guessbonushigh==1, lpattern(solid) lcolor(gs6)  mcolor(gs6) sort) ///
	(rcap jc_low jc_high wage_group if guessbonushigh==3, color(gs10)) ///
	(connected mean_jc wage_group if guessbonushigh==3, sort color(gs10)  mcolor(gs10))  ///
	(rcap jc_low jc_high wage_group if guessbonushigh==0, color(gs10) ) ///
	(connected mean_jc wage_group if guessbonushigh==0, lpattern(dashed) lcolor(gs1)  mcolor(gs1) msymbol(T) sort) ///
	if ($data_sample1) & (${wid_sample`s'})   ///
	, legend(order(2 "Predictions" "High Bonus" 4 "Predictions" "Low Bonus" 6 "Work Decision" "Future" 1 "95% CI") row(1) ) ///
	 xtitle("Wage") ytitle("Jobs Chosen") ylab(0(10)80) ysc(titlegap(3)) xlab(0(.05).31)  
	graph export "$savedir/nonpara_agg_bonus_future_`s'.png", replace width(2000)	
} 

*****************
****** Statistics
*****************

//USED: basic
xi: reg jobschosen i.guessbonushigh if ($data_sample1) & (${wid_sample2}), clu(wid)
xi: reg jobschosen i.guessbonushigh if future!=1 & guessbonushigh!=2 & ($data_sample1) & (${wid_sample2}), clu(wid)
//USED: controlling for wage
xi: areg jobschosen i.guessbonushigh if ($data_sample1) & (${wid_sample2}), clu(wid) absorb(wage)
xi: areg jobschosen i.guessbonushigh if future!=1 & guessbonushigh!=2 & ($data_sample1) & (${wid_sample2}), clu(wid) absorb(wage)
//USED: after decisiondate2
xi: reg jobschosen i.guessbonushigh if decisiondatenum>2 & ($data_sample1) & (${wid_sample2}), clu(wid)
xi: reg jobschosen i.guessbonushigh if decisiondatenum>2 & future!=1 & guessbonushigh!=2 & ($data_sample1) & (${wid_sample2}), clu(wid)







****************************************************************
********Projection Bias
****************************************************************

**********
** Figure 7: Comparison of decisions made before and after 10 mandatory tasks.
**********
foreach s of numlist 2 6 { //s=2 is the main sample with MLE-problem-subjects removed, s=6 is everyone (in Appendix)
	use "Data/decisions_data_w_ind", clear
	//keep if type=="present"
	cap drop t_* 
	drop if wage>.315
	gen projbias=(workdone1>0)
	egen group = group(wage_group projbias)
	xi: reg jobschosen i.group if ($data_sample1) & (${wid_sample`s'}), clu(wid)

	capture drop mean_jc
	predict mean_jc
	capture drop se_jc
	predict se_jc, stdp
	capture drop jc_low
	gen jc_low = mean_jc - 1.96*se_jc
	capture drop jc_high
	gen jc_high = mean_jc + 1.96*se_jc

	twoway (rcap jc_low jc_high wage_group if projbias==1, color(gs10)) ///
	(connected mean_jc wage_group if projbias==1, sort lpattern(dash) mcolor(gs8) lcolor(gs8))  ///
	(rcap jc_low jc_high wage_group if projbias==0, color(gs10) ) ///
	(connected mean_jc wage_group if projbias==0, msymbol(triangle) lpattern(solid) lcolor(gs1) mcolor(gs1) sort) ///
	if ($data_sample1) & (${wid_sample`s'}) ///
	, legend(order(2 "Decisions After Work" 4 "Decision Before Work" 1 "95% CI") row(1) ) ///
	 xtitle("Wage") ytitle("Jobs Chosen") ysc(titlegap(2)) ylab(0(10)80) xlab(0(.05).31)  
	graph export "$savedir/nonpara_agg_proj_`s'.png", replace width(2000)	
}

*****************
****** Statistics
*****************

use "Data/decisions_data_w_ind", clear
gen projbias=(workdone1>0)
cap drop group 
cap drop type_
encode type, gen(type_)
egen group = group(wage_group)

//USED: general differences
xi: reg jobschosen projbias i.type_ if  ($data_sample1) & (${wid_sample2}), clu(wid)
xi: areg jobschosen i.type_ projbias if  ($data_sample1) & (${wid_sample2}), clu(wid) absorb(wage)
xi: reg jobschosen i.type_ projbias if  decisiondatenum>2 & ($data_sample1) & (${wid_sample2}), clu(wid)

//USED: with tobit
xi: tobit jobschosen projbias i.type_ if  ($data_sample1) & (${wid_sample2}), clu(wid) ul(100) ll(0)
xi: tobit jobschosen i.wage projbias i.type_ if  ($data_sample1) & (${wid_sample2}), clu(wid) ul(100) ll(0)
xi: tobit jobschosen projbias i.type_ if  decisiondatenum>2 & ($data_sample1) & (${wid_sample2}), clu(wid) ul(100) ll(0)

//USED: corners:
gen corner0=(jobschosen==0)
gen corner100=(jobschosen==100)
ttest corner0, by(projbias) 
ttest corner100, by(projbias) 

//USED in footnote: separate
xi: areg jobschosen projbias if (type=="future") & ($data_sample1) & (${wid_sample2}), clu(wid) absorb(wage)
xi: areg jobschosen projbias if (type=="present") & ($data_sample1) & (${wid_sample2}), clu(wid) absorb(wage)
xi: areg jobschosen projbias if (type=="prediction") & ($data_sample1) & (${wid_sample2}), clu(wid) absorb(wage)




****************************************************************
********Individual histograms: beta beta-hat - will be used in Figure 6
****************************************************************
use "Data/decisions_data_w_ind", clear

hist bi1 if first==1 & ml_prob==0, width(.1) ///
xtitle("Individual {&beta} Estimates") xline(1, lpattern(dash) lwidth(1)) xsc(r(0 1.5)) xlabel(0(.5)1.5)  ysc(r(0 3)) ylabel(0(1)3) scheme(s2mono)
graph save Output/histbhb1, replace

hist bhi1 if first==1 & ml_prob==0, width(.1) ///
xtitle("Individual {&beta}{sub:h} Estimates") xline(1, lpattern(dash) lwidth(1)) xsc(r(0 1.5)) xlabel(0(.5)1.5)  ysc(r(0 3)) ylabel(0(1)3) scheme(s2mono)
graph save Output/histbhb2, replace

graph combine Output/histbhb1.gph Output/histbhb2.gph, ysize(1) iscale(1.5)
graph export "$savedir/para_indhist.png", replace width(2000)


	
	
**************************
******** Scatter plots: - used in Figures 4 and 6 - both for (prediction vs. present choices) and (bhi vs bi) 
**************************	
	

foreach s of numlist 2 6 { //s=2 is the main sample with MLE-problem-subjects removed, s=6 is everyone (in Appendix)
forvalues i=1(1)4 {
	local x="ind_diff_pred"
	local y="ind_diff_pres"
	local xText="Present Bias Perception"
	local yText="Present Bias"
	local cond="first==1"
	
	if (`i'>2) 	{	
		local x="bhi1"
		local y="bi1"
		local xText="Individual {&beta}{sub:h} Estimate"
		local yText="Individual {&beta} Estimate"
	}
	
	use "Data/decisions_data_w_ind", clear
	keep if ${wid_sample`s'}
	gen outlier=0
	local cond2="1==1"
	
	sum `x' if `cond'
	local xhigh=r(mean)+2.58*r(sd)
	local xlow=r(mean)-2.58*r(sd)
	sum `y'  if `cond'
	local yhigh=r(mean)+2.58*r(sd)
	local ylow=r(mean)-2.58*r(sd)	
	if (mod(`i',2)==0) {
		replace outlier=1 if `x'>`xhigh'
		replace outlier=1 if `x'<`xlow'
		replace outlier=1 if `y'>`yhigh'
		replace outlier=1 if `y'<`ylow'
		local cond2="outlier==0"
	}
	
	local addedlineCommand=""
	if (`i'==1) {
		local addedlineCommand="yline(`ylow' `yhigh', lcolor(gs13)) xline(`xlow' `xhigh', lcolor(gs13))" 
	}
	
	reg `y' `x' if `cond', r
	rreg `y' `x' if `cond'
	predict `y'_pred
	
	//(lfit `y' `x', lcolor(gs9) lpattern(dash)) 
	graph twoway (scatter `y' `x',  msize(small)) (line `y'_pred `x', lcolor(gs5) `addedlineCommand') ///
	if `cond' & `cond2', ///
	xtitle("`xText'") ytitle("`yText'") ///
	legend(off)
	graph save Output/scatter`i'_`s', replace
}
}

//this is main sample: non parametric
graph combine Output/histnonpara1_2.gph Output/histnonpara2_2.gph Output/scatter2_2.gph, ysize(1.7) iscale(1.5) col(3) imargin(small)
graph export "$savedir/nonpara_indhist_w_scatter.png", replace width(2000)

//this is everyone: non parametric
graph combine Output/histnonpara1_6.gph Output/histnonpara2_6.gph Output/scatter2_6.gph, ysize(1.7) iscale(1.5) col(3) imargin(small)
graph export "$savedir/nonpara_indhist_w_scatter_expanded.png", replace width(2000)

//this is main sample: parametric (ie bi vs bhi)
graph combine Output/histbhb1.gph Output/histbhb2.gph Output/scatter3_2.gph, ysize(1.7) iscale(1.5) col(3) imargin(small)
graph export "$savedir/para_indhist_w_scatter.png", replace width(2000)

//this is for the Appendix - shows how the scatter changes given removal of a few outliers
graph combine Output/scatter1_2.gph Output/scatter2_2.gph, col(2) ysize(2) iscale(1.5) 
graph export "$savedir/scatters_pred_pres_2.png", replace width(2000)


	
	
**************************
******** Coorelation and significance: B vs BH
**************************	
use "Data/decisions_data_w_ind", clear


foreach i of numlist 1/3 {
	local cond="(first==1 & ml_prob`i'==0)"
	
	//correlations
	corr bi`i' bhi`i' if `cond'
	//gets significance for different robust regression techniques
	reg bi`i' bhi`i'  if `cond'
	robreg m   bi`i' bhi`i' if `cond'
	robreg mm  bi`i' bhi`i' if `cond'
	robreg s   bi`i' bhi`i' if `cond'
}





****************************************************************
********Miscellaneous
***************************************************************



**************************
***Average time of task completion over time
**************************

use "Data/jobs", clear
//let's find the average for participation dates for people that completed
replace unixtime=unixtime-4*60*60*1000
gen day=dofc(unixtime)
egen group=group(wid day)
bys group: egen worktime=mean(timespent)
bys group: keep if _n==1
bys wid (day): gen pd=_n
bys wid (day): gen pdN=_N
keep if pdN==7
bys pd: sum worktime





********Understanding Individual Estimate Issues
use "Data/decisions_data_w_ind", clear

//four guys have no variance
bys wid: egen sd=sd(jobschosen)
tab wid if first==1 & ml_prob==1 & sd==0
//five more have very little variance
sum sd if first==1 & ml_prob==0
tab sd if first==1 & ml_prob==1
//these guys left early:
tab leftdatenum if first==1 & ml_prob==1
bys wid: gen N=_N
bys ml_prob: tab N if first==1

//these guys don't even get estimates
tab ml_prob if (bi1==. | bhi1==.) & first==1
//these guys have problems converging
tab ml_prob if iters1>=200 & (bi1==.) & first==1

//six participants had non-monotonicities in more than a third of their decision sets
cap drop d_w_p p_non_mono p_opportunities
egen d_w_p=group(decisiondate workdate wid type)
bys d_w_p (wage): gen t_non_mono=(jobschosen<jobschosen[_n-1])  
bys d_w_p (wage): replace t_non_mono=. if _n==1
bys d_w_p (wage): egen num_non_mono_in_d_w_p=total(t_non_mono)
bys d_w_p (wage): gen some_non_mono_in_d_w_p=num_non_mono_in_d_w_p>0
preserve
duplicates drop d_w_p, force
bys wid: gen num_d_w_p=_N
bys wid: egen num_non_mono_d_w_p=total(some_non_mono_in_d_w_p)
gen perc_d_w_p_w_non_mono=num_non_mono_d_w_p/num_d_w_p
duplicates drop wid, force
tab perc_d_w_p_w_non_mono if ml_prob==1
restore





******Payments
use "Data/decisions_data_w_ind", clear
sum payment if first==1 & left==0




******Summarizes Attrition 
use "Data/decisions_data_w_ind", clear
tab leftdatenum if first==1 & left==1



******Summarizes Date Modification
//some people "modified" during experimental day because button was on. don't count those.
use "Data/payments", clear
keep if type=="Modify"
keep if unixtime>clock("20oct2012 12:00:00", "DMY hms")
unique wid
duplicates tag wid, gen(dup)
tab dup
// 10 people changed date once. 2 people changed date twice.



********Average decision delay - for footnote on use of daily discount rate
use "Data/decisions_data_w_ind", clear
gen diff=workdate-decisiondate
sum diff if $wid_sample2 & $data_sample1

	
	





********Does getting more money affect future decisions?
//how to do this? See if person got random bonus in previous round.
//does it effect decisions now?
use "Data/decisions_data_w_ind", clear
//do decisions change with bonuses, condition on round
//so: how much money has person made on bonuses throughout?
//lets just get decision sets, get bonus level, and sum over time
bys group_decision_set: egen bonus_in_set=max(bonusamountifchosen)
duplicates drop group_decision_set, force
bys wid (unixtime): gen bonusesreceived=sum(bonus_in_set)
bys wid: gen final=(_n==_N)
//bonuses received at the end of the experiment: what is the variance?
sum bonusesreceived if final==1  & left==0, detail
keep group_decision_set bonusesreceived
save Data/bonusesreceived, replace
//how man 



use "Data/decisions_data_w_ind", clear
merge group_decision_set using Data/bonusesreceived, sort uniqusing
drop _merge
xi: areg jobschosen bonusesreceived i.decisiondatenum i.wid if  ($data_sample1) & (${wid_sample2}), clu(wid) absorb(wage)
xi: areg jobschosen bonusesreceived i.decisiondatenum i.wid if  ($data_sample1) & (${wid_sample2}), clu(wid) absorb(wage)











	
	
