do setupAndDefinePrograms.do


**********************************************
***How do people react to bonus decisions?
**********************************************



***** Note: problem with censored obs
use "Data/decisions_data_w_ind", clear
keep if diff_decision_bonus!=. & ml_prob==0
//clear that censored observations are an issue
//" This problem is clear in the data: for predictions of either 0 or 100, nearly 95% of the resultant work decisions perfectly match the prediction. As this e§ect is presumably due to censoring, we focus on predictions strictly between 0 and 100, where this confound does not exist. "
tab diff_decision_bonus if (ml_prob==0) & jobspredicted==0 | jobspredicted==100



***** summary statistics
use "Data/decisions_data_w_ind", clear
keep if diff_decision_bonus!=. & ml_prob==0
drop if jobspredicted==0 | jobspredicted==100

tab diff_decision_bonus
sum diff_decision_bonus
//is diff significantly diff from 0?
ttest diff_decision_bonus=0
//but this is not clustered! To cluster:
preserve
gen obs_id=_n
expand 2
bys obs_id: gen obs_n=_n
gen jobs=jobspredicted
replace jobs=jobschosen if obs_n==2
gen choice=(obs_n==2)
//replicate above
areg jobs choice, absorb(obs_id)
//adding in clusters
areg jobs choice, absorb(obs_id) clu(wid)
//diff across wages?
areg jobs choice if wage<.1, absorb(obs_id) clu(wid)
areg jobs choice if wage>.2, absorb(obs_id) clu(wid)
//check if increasing in wage?
reg diff_decision_bonus wage  			if ml_prob==0 & !(jobspredicted==0 | jobspredicted==100) & diff_decision_bonus>-6, clu(wid)
restore


***** histogram: in paper
use "Data/decisions_data_w_ind", clear
keep if diff_decision_bonus!=. & ml_prob==0
drop if jobspredicted==0 | jobspredicted==100
replace diff_decision_bonus=-6 if diff_decision_bonus<-5
gen diff_predicted=_n-7 if _n<14
gen perc_true=.
forvalues d=-6(1)6 {
	gen dummy=(diff_decision_bonus==`d')
	sum dummy
	replace perc_true=r(mean) if diff_predicted==`d'
	drop dummy
}
keep diff_predicted perc_true
keep if _n<14

gen freq1=round(perc_true*1000000)
twoway (hist diff_predicted [fweight=freq1], frac  barwidth(1) )

twoway (hist diff_predicted [fweight=freq1], frac  barwidth(.5) ), ///
xscale(r( -6 0 +6)) xlabel(-6(1)6) xline(-1.53, lpattern(solid ) lcolor(gs7) style(unextended)) xline(5.5 -5.5, lpattern(dash) lcolor(gs7) style(unextended))   ///
xlabel(-6 "<-5" -5 "-5" -4 "-4" -3 "-3" -2 "-2" -1 "-1" 0 "0" 6 ">5" 5 "5" 4 "4" 3 "3" 2 "2" 1 "1") ///
ytitle("Fraction of Data") xtitle("Work Decision minus Prediction") scheme(s2mono)
graph save "Output/t_hist1", replace
graph export "$savedir/pred_chosen_diff.png", replace



***consistency preferences are surprisingly "local."
use "Data/decisions_data_w_ind", clear
egen decision_set=group(wid unixtime type)
gen within1=0
gen within2=0
gen within3=0
gen within4=0
forvalues i=1(1)5 {
	bys decision_set (wage): gen t_p=wage if type=="present" & (bonusoffered==1) & _n==`i'
	bys decision_set (wage): egen p`i'=max(t_p)
	drop t_p
	order p`i'
	replace within1=1 if (wage<=(p`i'+.01) & wage>=(p`i'-.01)) 
	replace within2=1 if (wage<=(p`i'+.02) & wage>=(p`i'-.02)) 
	replace within3=1 if (wage<=(p`i'+.03) & wage>=(p`i'-.03)) 
	replace within4=1 if (wage<=(p`i'+.04) & wage>=(p`i'-.04)) 
}

forvalues i=1(1)4 {
	tab bonusoffered if within`i'==1 & (type=="present") & ($wid_sample2)
	xi: areg jobschosen bonusoffered if within`i'==1 & (type=="present") & ($wid_sample2), clu(wid) absorb(wage)
}


*******Rounding
use "Data/decisions_data_w_ind", clear
keep if diff_decision_bonus!=. & ml_prob==0
drop if jobspredicted==0 | jobspredicted==100
//see if rounded to 5 or 10
gen rounded5=(round(jobspredicted/5)-(jobspredicted/5)==0)
gen rounded10=(round(jobspredicted/10)-(jobspredicted/10)==0)
//lots of rounded predictions
sum rounded5
sum rounded10
//the predictions of a round number of more problematic
bys rounded5: tab diff_decision_bonus 
bys rounded10: tab diff_decision_bonus 



******SIMULATIONS: rough but gives an OK idea
//ok - so, let's look at predictions and think about the fact what happens when the person arrives given our parameters
use "Data/decisions_data_w_ind", clear
keep if ml_prob==0
drop if diff_decision_bonus==.

//let's think about the predictions we see: what do we expect people to do when they arrive?
gen jobs_when_arrive=((jobspredicted+10)*.83)-10
replace jobs_when_arrive=100 if jobspredicted==100
replace jobs_when_arrive=0 if jobspredicted==0
//this gives us very reasonable looking predictions
sum jobs_when_arrive jobspredicted

//but we dont want to use corners as discussed above
drop if jobspredicted==0 | jobspredicted==100

//so, we draw gammas. What gamma do we need to justify the jobs_when_arrive?
gen avg_gamma_when_arrive=(ln(10 + jobs_when_arrive) + ln(wage*.83*500))/ln(10 + jobs_when_arrive)
sum avg_gamma_when_arrive
//this is good -> has the right stats

//now, go through every prediction and simulate what happens when we arrive and draw gammas
//save a file to append to later
preserve
clear
gen perc=.
save "Output/sim", replace
restore

local N=_N
forvalues i=1(1)`N' {
disp `i'
qui {

	//store some vars
	local avg_gamma_when_arrive=avg_gamma_when_arrive[`i']
	local bonusamountifchosen=bonusamountifchosen[`i']
	local jobspredicted=jobspredicted[`i']
	local wage=wage[`i']
	
	preserve
	clear
	set obs 1000000
	gen jobspredicted=`jobspredicted'
	//this is the gamma we might get when we arrive -> its random and has SD (from ML simulation with error on gamma)
	gen gamma_when_arrive=`avg_gamma_when_arrive'+rnormal()*.10
	//given this gamma, what is optimal effort?
	gen optimal_effort=(500*.83*`wage')^(1/(gamma_when_arrive - 1)) - 10
	//what's the utility from this?
	gen utility_optimal=optimal_effort*500*`wage'*.83-(1/gamma_when_arrive)*(optimal_effort+10)^(gamma_when_arrive)
	//what's the utility from following through?
	gen utility_pred=`jobspredicted'*500*`wage'*.83-(1/gamma_when_arrive)*(`jobspredicted'+10)^(gamma_when_arrive)
	//how much will bonus add to utility?
	gen utility_pred_w_bonus=utility_pred+500*.83*`bonusamountifchosen'*.83
	//will the person follow through?
	gen follow_through=(utility_pred_w_bonus>utility_optimal)
	
	//now, think about what we would observe in terms of diff between pred and outcome
	//have -6 and +6 represent when don't follow through up or down
	replace optimal_effort=0 if optimal_effort<0
	replace optimal_effort=100 if optimal_effort>100
	gen diff=round(optimal_effort-jobspredicted)
	
	sum diff
	gen average_diff=r(mean)
	
	replace diff=-5 if follow_through==1 & diff<-5
	replace diff=5 if follow_through==1 & diff>5
	replace diff=-6 if follow_through==0 & diff<-5
	replace diff=6 if follow_through==0 & diff>5
	keep diff 
	
	//now store tabbed results only. Much easier
	gen diff_predicted=_n-7 if _n<14
	gen perc=.
	forvalues d=-6(1)6 {
		gen dummy=(diff==`d')
		sum dummy
		replace perc=r(mean) if diff_predicted==`d'
		drop dummy
	}
	keep diff_predicted perc 
	keep if _n<14
	append using "Output/sim"
	save "Output/sim", replace
	
	restore
}
}

//let's add in empirical data on 0/100 to see what happens
clear
set obs 13
gen diff_predicted=_n-7 if _n<14
gen perc=0
replace perc=.95 if diff_predicted==0
replace perc=.05 if diff_predicted==-5
gen dup=131
expand dup
drop dup
save "Output/sim_zeros", replace

use "Output/sim", clear
append using "Output/sim_zeros"
save "Output/sim_all", replace

//sim histogram -> without zeros
use "Output/sim", clear
collapse (mean) perc, by(diff_predicted)
gen freq1=round(perc*1000000)
twoway (hist diff_predicted [fweight=freq1], frac  barwidth(.5) ), ///
xscale(r( -6 0 +6)) xlabel(-6(1)6) xline(-3.4, lpattern(solid ) lcolor(gs7) style(unextended)) xline(5.5 -5.5, lpattern(dash) lcolor(gs7) style(unextended))   ///
xlabel(-6 "<-5" -5 "-5" -4 "-4" -3 "-3" -2 "-2" -1 "-1" 0 "0" 6 ">5" 5 "5" 4 "4" 3 "3" 2 "2" 1 "1") ///
ytitle("Fraction of Data") xtitle("Work Decision minus Prediction (simulated)")
graph save "Output/t_hist2", replace
graph export "$savedir/pred_chosen_diff_sim.png", replace

graph combine "Output/t_hist1" "Output/t_hist2",  ycommon xcommon  ysize(2) iscale(1.3)
graph export "$savedir/pred_chosen_diff_compare.png", replace







