***************************************************************************************************
* 
* Program: NLS_main.do
* Authors: Thomas Graeber, Avner Shlain
* Last Modified: 3/29/2017
* Purpose:
*       Run all NLS specifications - Table 5, 6 and Appendix Table 4
*		
* Files Used:
*     1. MTurkCleanedDataShort.dta 
*
***************************************************************************************************

* Preliminaries
global user ="ASS"

if("${user}"=="TG"){
	global main = "/Users/thomasgraeber/Dropbox/forecastExp/Data Analysis/Thomas Analysis"
}

if("${user}"=="ASS"){
	global main = "/Users/Avner/Dropbox/forecastExp/Data Analysis/Analysis/behavioral"
}

cd "$main"
set more off
use "$main/dtafiles/MTurkCleanedDataShort.dta", clear


* Additional variable indicating payoff per 100 button presses (p)
gen payoff_per_100 = 0
replace payoff_per_100 = 0.01 if treatment == "1.1"
replace payoff_per_100 = 0.1 if treatment == "1.2"
replace payoff_per_100 = 0.0 if treatment == "1.3"
replace payoff_per_100 = 0.001 if treatment == "2"
replace payoff_per_100 = 0.04 if treatment == "1.4"
replace payoff_per_100 = 0.01 if treatment == "4.1"
replace payoff_per_100 = 0.01 if treatment == "4.2"
replace payoff_per_100 = 0.02 if treatment == "6.2"
replace payoff_per_100 = 1 if treatment == "6.1"


* (***ALPHA***) Payoff per 100 to charity
gen payoff_charity_per_100 = 0
replace payoff_charity_per_100 = 0.01 if treatment == "3.1"
replace payoff_charity_per_100 = 0.1 if treatment == "3.2"

gen dummy_charity = 0
replace dummy_charity = 1 if treatment == "3.1"
replace dummy_charity = 1 if treatment == "3.2"

* (***BETA/DELTA***) Payoff per 100 delayed by 2 wks
gen delay_wks = 0
replace delay_wks = 2 if treatment == "4.1"
replace delay_wks = 4 if treatment == "4.2"

gen delay_dummy = 0
replace delay_dummy = 1 if treatment == "4.1"
replace delay_dummy = 1 if treatment == "4.2"

* Probability weights to back out curvature
gen prob = 1
replace prob = 0.5 if treatment == "6.2"
replace prob = 0.01 if treatment == "6.1"

gen weight_dummy = 0
replace weight_dummy = 1 if treatment == "6.1"

* Dummy for gift exchange
gen gift_dummy = 0
replace gift_dummy = 1 if treatment == "10"

capt log close
log using "$main/output/analyses", append

************************************************************
* Run NLS for both power and exponential cost function


** Models and treatments included in estimation
global models POWER EXPON
//global specifications Bench Bench_fix088 All_wo_prob_weighting No_pweight_fix088 All_fix_088 All_estim
global precision = 1e-15 // Default is 1e-5

// Table 5:
// A col 2 - power 3 (k,s,gamma)
// A col 4 - exp 3 + 3 specifications for the low pay
// B col 3 - power - k,s,gamma,a,alpha,s_ge,beta,delta
// B col 6 - exp   - k,s,gamma,a,alpha,s_ge,beta,delta


global specificationsT5 Bench All_nopweight
	preserve

	gen buttonpresses_nearest100 = round(buttonpresses, 100)
	replace buttonpresses_nearest100 = 25 if buttonpresses_nearest100 == 0
	label variable buttonpresses_nearest100 "Effort"
	* Log of rounded button presses
	gen logbuttonpresses_nearest100 = ln(buttonpresses_nearest100)	
	label variable logbuttonpresses_nearest100 "Log effort"
	
	foreach model of global models {
		
	* Cost function specifications		
		if ("`model'" == "POWER") {
			local depvar logbuttonpresses_`rounding'
			global st_values `"k 1.66306e-10 gamma 19.8117987 s 7.74996"'
			local k_scaler = 1e+57
			local s_scaler = 1e+6
		}
		
		if ("`model'" == "EXPON") {
			local depvar buttonpresses_`rounding'
			global st_values `"k 1.69443 gamma 0.015645717 s 3.69198"'
			local k_scaler = 1e+16
			local s_scaler = 1e+6
		}

	eststo clear

	foreach spec in $specificationsT5 {
		if "`spec'" == "Bench" { 
		// col 1
			local sample1 `""1.1", "1.2""'
			local sample2 `""1.3""'
			local sample Benchmark
			local alpha = 1
			local a = 1
			local beta = 1
			local delta = 1
			local p_weight = 1
			local curv = 1
			local gift = 1
			local stval_spec
		}				
	if "`spec'" == "All_nopweight" { 
		local sample1 `""1.1", "1.2", "1.3", "3.1", "3.2""'
		local sample2 `""4.1", "4.2", "10""'
		local sample Bench_Social_Time_Gift
		local alpha = `"{alpha}"'
		local a = `"{a}"'
		local beta = `"{beta}"'
		local delta = `"{delta}"'
		local p_weight = 1
		local curv = 1
		local gift `"{gift}"'
		local stval_spec `"alpha 0.003 a 0.13 beta 1.16 delta 0.75 gift 5e-6"'
	}	

	* Run model
	** Use nl default for standard errors (Gauss-Newton)

	tempvar sample_dummy
	g `sample_dummy' = 1 if inlist(treatment, `sample1')
	replace `sample_dummy' = 1 if inlist(treatment, `sample2')

	nl (`depvar' = -1/{gamma} * ln({k}/`k_scaler') + 1/ {gamma}*ln({s}/`s_scaler'+`gift'*0.4*gift_dummy + `beta'^delay_dummy*`delta'^delay_wks*`p_weight'^weight_dummy*prob*payoff_per_100^`curv'+`alpha'*payoff_charity_per_100+`a'*0.01*dummy_charity)) ///
	if `sample_dummy' == 1, initial($st_values `stval_spec') noc robust eps($precision)

			predict effort_hat_`spec'_`model'
							
			* Store results
			capt eststo `spec'_`model'
			capt estadd local sam `"`sample'"': `spec'_`model'
			capt estadd local mod `model': `spec'_`model'
			capt estadd local spec `spec': `spec'_`model'
				capt estadd scalar scaling_k = `k_scaler' : `spec'_`model'
				capt estadd scalar scaling_s = `s_scaler' : `spec'_`model'
			local prediction = 1
			foreach treatment in "1.1" "1.2" "1.3" "1.4" "2" "3.1" "3.2" "4.1" "4.2" "6.1" "6.2" {
				sum effort_hat_`spec'_`model' if treatment == "`treatment'", mean
				capt estadd scalar pred_`prediction' = r(mean): `spec'_`model'
				local ++prediction
			}
			//drop effort_hat
			
	}
* Report results

// levelsof treatment, l(treatments)
local tablenote1
local tablenote2
foreach t in "1.1" "1.2" "1.3" "1.4" "2" "3.1" "3.2" "4.1" "4.2" "6.1" "6.2" {
	sum buttonpresses_nearest100 if treatment == "`t'", mean
	local avg = int(`r(mean)')
	local tablenote1 `"`tablenote1' ''`t''': `avg', "'
}
foreach t in "1.1" "1.2" "1.3" "1.4" "2" "3.1" "3.2" "4.1" "4.2" "6.1" "6.2" {
	sum logbuttonpresses_nearest100 if treatment == "`t'", mean
	local avg = round(`r(mean)',0.001)
	local tablenote2 `"`tablenote2' ''`t''': `avg', "'

}

capt esttab _all using "${main}/output/NLS_results_Table5_`model'", nogaps replace ///
csv label b(a8) se(a8) r2(4) ar2(4) nostar ///
scalars("sam Sample" "mod Cost function" "spec Specification" "scaling_k Coef. k scaled up by" "scaling_s Coef. s scaled up by" "rss RSS" "rmse Root MSE" "msr Residual mean square" ///
"dev Residual deviance" "ic \# iterations" "pred_1 Predicted effort 1.1" /// 
 "pred_2 Predicted effort 1.2" "pred_3 Predicted effort 1.3" "pred_4 Predicted effort 1.4" "pred_5 Predicted effort 2" ///
 "pred_6 Predicted effort 3.1" "pred_7 Predicted effort 3.2" "pred_8 Predicted effort 4.1" ///
 "pred_9 Predicted effort 4.2" "pred_10 Predicted effort 6.1" "pred_11 Predicted effort 6.2") ///
 nonotes addnotes("Sample average effort: `tablenote1'" "Sample average log effort: `tablenote2'")

eststo clear
}

restore



// Table 6:
// col 1 - power - s,k,gamma, prob weighting, theta=1
// col 2 - power - s,k,gamma, prob weighting, theta=0.88
// col 3 - power - s,k,gamma, prob weighting, theta estimated
// col 4 - expon - s,k,gamma, prob weighting, theta=1
// col 2 - expon - s,k,gamma, prob weighting, theta=0.88
// col 3 - expon - s,k,gamma, prob weighting, theta estimated

global specificationsT6 Curv_fix1 Curv_fix088 Curv_estim



preserve
	gen buttonpresses_nearest100 = round(buttonpresses, 100)
	replace buttonpresses_nearest100 = 25 if buttonpresses_nearest100 == 0
	label variable buttonpresses_nearest100 "Effort"
	* Log of rounded button presses
	gen logbuttonpresses_nearest100 = ln(buttonpresses_nearest100)	
	label variable logbuttonpresses_nearest100 "Log effort"
	foreach model of global models {
		
* Cost function specifications
			
		if ("`model'" == "POWER") {
		local depvar logbuttonpresses_nearest100
		
			global st_values `"k 1.66306e-10 gamma 19.8117987 s 7.74996"'
			local k_scaler = 1e+57
			local s_scaler = 1e+6
		}
		
		if ("`model'" == "EXPON") {
		local depvar buttonpresses_nearest100
			global st_values `"k 1.69443 gamma 0.015645717 s 3.69198"'
			local k_scaler = 1e+16
			local s_scaler = 1e+6
		}

eststo clear

	foreach spec in $specificationsT6 {
		if "`spec'" == "Curv_fix088" { 
			local sample1 `""1.1", "1.2", "1.3""'
			local sample2 `""6.1", "6.2""'
			local sample Bench_Prob
			local alpha = 1
			local a = 1
			local beta = 1
			local delta = 1
			local p_weight = `"{p_weight}"'	
			local curv = 0.88
			local gift = 1
			local stval_spec `"p_weight 0.2"'
		}
		if "`spec'" == "Curv_fix1" { 
			local sample1 `""1.1", "1.2", "1.3""' 
			local sample2 `""6.1", "6.2""'
			local sample Bench_Prob
			local alpha = 1
			local a = 1
			local beta = 1
			local delta = 1
			local p_weight = `"{p_weight}"'	
			local curv = 1
			local gift = 1
			local stval_spec `"p_weight 0.2"'
		}
			if "`spec'" == "Curv_estim" { 
			local sample1 `""1.1", "1.2", "1.3""'
			local sample2 `""6.1", "6.2""'
			local sample Bench_Prob
			local alpha = 1
			local a = 1
			local beta = 1
			local delta = 1
			local p_weight = `"{p_weight}"'
			local curv `"{curv}"'
			local gift = 1
			local stval_spec `"p_weight 0.2 curv 0.5"'
		}	
* Run model
** Use nl default for standard errors (Gauss-Newton)

tempvar sample_dummy
g `sample_dummy' = 1 if inlist(treatment, `sample1')
replace `sample_dummy' = 1 if inlist(treatment, `sample2')

			nl (`depvar' = -1/{gamma} * ln({k}/`k_scaler') + 1/ {gamma}*ln({s}/`s_scaler'+`gift'*0.4*gift_dummy + `beta'^delay_dummy*`delta'^delay_wks*`p_weight'^weight_dummy*prob*payoff_per_100^`curv'+`alpha'*payoff_charity_per_100+`a'*0.01*dummy_charity)) ///
			if `sample_dummy' == 1, initial($st_values `stval_spec') noc robust eps($precision)

			predict effort_hat
							
			* Store results
			capt eststo `spec'_`model'
			capt estadd local sam `"`sample'"': `spec'_`model'
			capt estadd local mod `model': `spec'_`model'
			capt estadd local spec `spec': `spec'_`model'
				capt estadd scalar scaling_k = `k_scaler' : `spec'_`model'
				capt estadd scalar scaling_s = `s_scaler' : `spec'_`model'
			local prediction = 1
			foreach treatment in "1.1" "1.2" "1.3" "1.4" "2" "3.1" "3.2" "4.1" "4.2" "6.1" "6.2" {
				sum effort_hat if treatment == "`treatment'", mean
				capt estadd scalar pred_`prediction' = r(mean): `spec'_`model'
				local ++prediction
			}
			drop effort_hat
			
	}
	* Report results

// levelsof treatment, l(treatments)
local tablenote1
local tablenote2
foreach t in "1.1" "1.2" "1.3" "1.4" "2" "3.1" "3.2" "4.1" "4.2" "6.1" "6.2" {
	sum buttonpresses_`rounding' if treatment == "`t'", mean
	local avg = int(`r(mean)')
	local tablenote1 `"`tablenote1' ''`t''': `avg', "'
}
foreach t in "1.1" "1.2" "1.3" "1.4" "2" "3.1" "3.2" "4.1" "4.2" "6.1" "6.2" {
	sum logbuttonpresses_`rounding' if treatment == "`t'", mean
	local avg = round(`r(mean)',0.001)
	local tablenote2 `"`tablenote2' ''`t''': `avg', "'

}

capt esttab _all using "${main}/output/estimation_results_T6_`model'", nogaps replace ///
csv label b(a8) se(a8) r2(4) ar2(4) nostar ///
scalars("sam Sample" "mod Cost function" "spec Specification" "scaling_k Coef. k scaled up by" "scaling_s Coef. s scaled up by" "rss RSS" "rmse Root MSE" "msr Residual mean square" ///
"dev Residual deviance" "ic \# iterations" "pred_1 Predicted effort 1.1" /// 
 "pred_2 Predicted effort 1.2" "pred_3 Predicted effort 1.3" "pred_4 Predicted effort 1.4" "pred_5 Predicted effort 2" ///
 "pred_6 Predicted effort 3.1" "pred_7 Predicted effort 3.2" "pred_8 Predicted effort 4.1" ///
 "pred_9 Predicted effort 4.2" "pred_10 Predicted effort 6.1" "pred_11 Predicted effort 6.2") ///
 nonotes addnotes("Sample average effort: `tablenote1'" "Sample average log effort: `tablenote2'")

eststo clear
}

restore










// Appendix Table 4:
// col 1 - expon - s,k, gamma = 0.01
// col 2 - expon - s,k, gamma = 0.02
// col 3 - expon - s,k,gamma, theta = 0.88
// col 4 - expon - s,k,gamma, continuous points (no rounding)
local m=1
global models POWER

global specificationsAT4 All_g01 All_g02 All_fix088 All_noround
local Bench_Social_Time_Gift `""1.1", "1.2", "1.3", "3.1", "3.2", "4.1", "4.2", "10""'
local Bench `""1.1", "1.2", "1.3""'

eststo clear


* Round to nearest 100
		gen buttonpresses_nearest100 = round(buttonpresses, 100)
		replace buttonpresses_nearest100 = 25 if buttonpresses_nearest100 == 0
		label variable buttonpresses_nearest100 "Effort"

		gen buttonpresses_no_rounding = buttonpresses
		replace buttonpresses_no_rounding = 1 if buttonpresses_no_rounding == 0
		label variable buttonpresses_no_rounding "Effort"

	local k_scaler = 1e+16
	local s_scaler = 1e+6
	local sample Bench_Social_Time_Gift
	
	foreach spec in $specificationsAT4 {
		if "`spec'" == "All_g01" { 
			global st_values `"k 1.69443 s 3.69198"'
			local gamma = 0.01
			local alpha = `"{alpha}"'
			local a = `"{a}"'
			local beta = `"{beta}"'
			local delta = `"{delta}"'
			local p_weight = 1
			local curv = 1
			local gift `"{gift}"'
			local stval_spec `"alpha 0.003 a 0.13 beta 1.16 delta 0.75 gift 5e-6"'
			local rou = 1
		}				
		if "`spec'" == "All_g02" { 
			global st_values `"k 1.69443 s 3.69198"'
			local gamma = 0.02
			local alpha = `"{alpha}"'
			local a = `"{a}"'
			local beta = `"{beta}"'
			local delta = `"{delta}"'
			local p_weight = 1
			local curv = 1
			local gift `"{gift}"'
			local stval_spec `"alpha 0.003 a 0.13 beta 1.16 delta 0.75 gift 5e-6"'
			local rou = 1
		}				
		if "`spec'" == "All_fix088" { 
			global st_values `"k 1.69443 gamma 0.015 s 3.69198"'
			local gamma = `"{gamma}"'
			local alpha = `"{alpha}"'
			local a = `"{a}"'
			local beta = `"{beta}"'
			local delta = `"{delta}"'
			local p_weight = 1
			local curv = 0.88
			local gift `"{gift}"'
			local stval_spec `"alpha 0.003 a 0.13 beta 1.16 delta 0.75 gift 5e-6"'
			local rou = 1
		}

		if "`spec'" == "All_noround" { 
			global st_values `"k 1.69443 gamma 0.015 s 3.69198"'
			local gamma = `"{gamma}"'
			local alpha = `"{alpha}"'
			local a = `"{a}"'
			local beta = `"{beta}"'
			local delta = `"{delta}"'
			local p_weight = 1
			local curv = 1
			local gift `"{gift}"'
			local stval_spec `"alpha 0.003 a 0.13 beta 1.16 delta 0.75 gift 5e-6"'
			local rou = 0
		}


		if `rou' == 1 {
		 local depvar buttonpresses_nearest100
		 }
		 if `rou' == 0 {
		 local depvar buttonpresses_no_rounding
		 }
		
	nl (`depvar' = -1/`gamma' * ln({k}/`k_scaler') + 1/ `gamma'*ln({s}/`s_scaler'+`gift'*0.4*gift_dummy + `beta'^delay_dummy*`delta'^delay_wks*`p_weight'^weight_dummy*prob*payoff_per_100^`curv'+`alpha'*payoff_charity_per_100+`a'*0.01*dummy_charity)) ///
	if inlist(treatment, ``sample''), initial($st_values `stval_spec') noc robust eps($precision)

						predict effort_hat

			capt eststo m_`m'
			capt estadd local sam `"`sample'"': m_`m'
			capt estadd local spec `spec': m_`m'
				capt estadd scalar scaling_k = `k_scaler' : m_`m'
				capt estadd scalar scaling_s = `s_scaler' : m_`m'
			local prediction = 1
			foreach treatment in "1.1" "1.2" "1.3" "1.4" "2" "3.1" "3.2" "4.1" "4.2" "6.1" "6.2" {
				sum effort_hat if treatment == "`treatment'", mean
				capt estadd scalar pred_`prediction' = r(mean): m_`m'
				local ++prediction
			}
			drop effort_hat
			local ++m
	}
* Report results

// levelsof treatment, l(treatments)
local tablenote1
local tablenote2
foreach t in "1.1" "1.2" "1.3" "1.4" "2" "3.1" "3.2" "4.1" "4.2" "6.1" "6.2" {
	sum buttonpresses_no_rounding if treatment == "`t'", mean
	local avg = int(`r(mean)')
	local tablenote1 `"`tablenote1' ''`t''': `avg', "'
}



capt esttab _all using "${main}/output/estimation_results_AT4_`rounding'", nogaps replace ///
csv label b(a8) se(a8) r2(4) ar2(4) nostar ///
scalars("sam Sample" "mod Cost function" "spec Specification" "scaling_k Coef. k scaled up by" "scaling_s Coef. s scaled up by" "rss RSS" "rmse Root MSE" "msr Residual mean square" ///
"dev Residual deviance" "ic \# iterations" "pred_1 Predicted effort 1.1" /// 
 "pred_2 Predicted effort 1.2" "pred_3 Predicted effort 1.3" "pred_4 Predicted effort 1.4" "pred_5 Predicted effort 2" ///
 "pred_6 Predicted effort 3.1" "pred_7 Predicted effort 3.2" "pred_8 Predicted effort 4.1" ///
 "pred_9 Predicted effort 4.2" "pred_10 Predicted effort 6.1" "pred_11 Predicted effort 6.2") ///
 nonotes addnotes("Sample average effort: `tablenote1'" "Sample average log effort: `tablenote2'")

eststo clear

log close

