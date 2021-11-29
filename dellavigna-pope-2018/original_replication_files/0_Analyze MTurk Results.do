***************************************************************************************************
* 
* Program: Analyze MTurk Results.do
* Authors: Devin Pope
* Last Modified: 3/28/2017
* Purpose:
*       Upload and analyze the 10,000-person MTurk experiment data and related Forecast data
*		Runs analyses to output figures and tables for REStud paper 
*
* Files Used:
*     1. input/MTurkDataDemo
*     2. dtafiles/MTurkCleanedData
*     3. dtafiles/ExpertForecastCleanWide_all
*     4. dtafiles/ExpertForecastCleanWide
*     5. dtafiles/ExpertForecastCleanLong
*
***************************************************************************************************
cd your_foldername_here

********************************************
set more off
capture mkdir "output"

********************************************
*********** Load Packages ******************
********************************************
*ssc install estout
*ssc install confirmdir
*ssc install labutil
*ssc install egenmore

/* VERSION OF DISTPLOT FOR THIS CODE:
package name:  gr41_4.pkg
        from:  http://www.stata-journal.com/software/sj10-1/
NOTE: ssc install distplot will not work
*/

*************************************
*Define the treatment order
local tx_labels "No Payment" "1c PieceRate" "10c PieceRate" "4c PieceRate" "Very Low Pay" "1c RedCross" "10c RedCross" "Gift Exchange" "1c 2Wks" "1c 4Wks" "Gain 40c" "Loss 40c" "Gain 80c" "Prob.01 $1" "Prob.5 2c" "Social Comp" "Ranking" "Task Signif"

*************************************
* Analyze MTurk Effort - main analyses
*************************************


*************************************
*Summary Stats - Online Appendix Table 1

*button presses and duration
use dtafiles/MTurkCleanedData, replace
eststo clear
eststo, title("Mean"): estpost sum buttonpresses duration_time_min
esttab using output/TableSumStat_Mturk_noDemo.csv, mtitles csv main(mean) b(2) nostar nogaps unstack replace

*Demographic variables
use input/MTurkDataDemo, clear
keep if finished==1 & approved==1 & drop_flag==0
eststo clear
eststo, title("Mean"): estpost sum us_ip india_ip female ed_high ed_some_c ed_bach age_18 age_25 age_31 age_41 age_51 age_65 
esttab using output/TableSumStat_Mturk_Demo.csv, mtitles csv main(mean) b(2) nostar nogaps unstack replace


*************************************
*Summarize Effort

use dtafiles/MTurkCleanedData, clear	

**Let's see the average number of button presses by treatment group and also 
** create a effective_order which is a number from 1 to 18 indicating how effective each treatment 
** group was with 1 being the least effective and 18 being the most effective
bys treatmentname: egen avg_buttonpresses = mean(buttonpresses)
bys treatmentname: egen sd_buttonpresses = sd(buttonpresses)
bys treatmentname: egen obs_buttonpresses = count(_N)
list treatmentname avg_buttonpresses if treatmentname ~= treatmentname[_n-1]
sort avg_buttonpresses
list treatmentname avg_buttonpresses sd_buttonpresses obs_buttonpresses if treatmentname ~= treatmentname[_n-1]

egen presses_order_tx = group(avg_buttonpresses)
labmask presses_order_tx, value(treatmentname)
tabstat buttonpresses, by(presses_order_tx) stat(n mean semean)

* Online Appendix Figure 4
**Let's now look at the histogram of buttonpresses for all survey takers
histogram buttonpresses, width(25) start(0) xtitle("Button Presses") ytitle("Frequency") ///
	freq title("Distribution of Button Presses - All Treatments") graphregion(color(white)) lcolor(black) fcolor(black) ///
	xlabel(0(500)4000) xmtick(0(100)4000,grid) gap(75) ylabel(,nogrid)
graph export "output/Dist Button Presses All Treatments.png", replace


*Table 3 - Findings by treatment, Columns 3, 4
*Also used for Figure 3
tabstat buttonpresses, by(tx_order) stat(n mean semean)

collapse (count) n = buttonpresses (mean) mean_effort = buttonpresses ///
    (semean) se_mean = buttonpresses (sd) sd_effort = buttonpresses, by(tx_order treatmentname)
label values tx_order
export excel using output/TableSumStatByTreatment_c34.xlsx, firstrow(variables) replace

*************************************
*Bayesian shrinkage - online appendix figure 6

sum mean_effort
*generate a theta to equal the variance of the 18 treatment variables
gen theta = `r(Var)'
*generate a grand_avg variable that equals the average of the 18 treatment means
gen grand_avg = `r(mean)'
*gen se = sd_buttonpresses/sqrt(obs_buttonpresses) 

gen mean_bayesian_shrink = (theta/(theta+se_mean*se_mean))*mean_effort + (1-(theta/(theta+se_mean*se_mean)))*grand_avg
sort mean_effort
export excel using output/BayesianShrinkage.xlsx, firstrow(variables) replace

*************************************
*CDFs of effort for different treatments - Figure 4
use dtafiles/MTurkCleanedData, clear	
replace buttonpresses = 3000 if buttonpresses>3000 & buttonpresses!=.

*relabel treatments for these graphs
replace treatmentname="1 Cent for 100" if treatmentname=="1c PieceRate"
replace treatmentname="10 Cents for 100" if treatmentname=="10c PieceRate"
replace treatmentname="1 Cent for 1000" if treatmentname=="Very Low Pay"
replace treatmentname="Social Comparison" if treatmentname=="Social Comp"
replace treatmentname="Task Significance" if treatmentname=="Task Signif"
replace treatmentname="40 Cents Gain" if treatmentname=="Gain 40c"
replace treatmentname="80 Cents Gain" if treatmentname=="Gain 80c"
replace treatmentname="40 Cents Loss" if treatmentname=="Loss 40c"

labmask tx_order, value(treatmentname)

distplot buttonpresses if inlist(tx_order,1,2,3,5), ///
	legend(order(1 4 2 3)) over(tx_order) graphregion(color(white)) ///
	lcolor(navy forest_green dkorange maroon) ///
	lpattern(solid dash_dot shortdash longdash) ///
	ytitle("Cumulative Fraction") xtitle("Points in Task")
graph export "output/Fig4a.png", replace
graph export "output/Fig4a.eps", replace

distplot buttonpresses if inlist(tx_order,1,8,16,17,18), ///
	legend(order(1 2 3 4 5)) over(tx_order) graphregion(color(white)) ///
	lcolor(navy maroon forest_green dkorange edkblue) ///
	lpattern(solid longdash dash_dot shortdash shortdash_dot) ///
	ytitle("Cumulative Fraction") xtitle("Points in Task")
graph export "output/Fig4b.png", replace
graph export "output/Fig4b.eps", replace

distplot buttonpresses if inlist(tx_order,11,12,13), ///
	legend(order(1 2 3)) over(tx_order) graphregion(color(white)) ///
	lcolor(navy maroon forest_green) ///
	lpattern(solid longdash dash_dot) ///
	ytitle("Cumulative Fraction") xtitle("Points in Task")
graph export "output/Fig4c.png", replace
graph export "output/Fig4c.eps", replace


*************************************
* Analyze Forecasts
*************************************

**** Summary Statistics -- Table 2
	use dtafiles/ExpertForecastCleanWide_all, clear
	*median duration for Experts Completed All 15 Tx
	sum duration if sampleshort=="experts" & completed=="Yes", d
	gen median_duration = `r(p50)' if sampleshort=="experts" & completed=="Yes"

	eststo clear
	*beh_econ appl theory lab_econ psych
	eststo, title("All Experts Contacted"): estpost sum beh_econ_only beh_fin appl theory lab_econ psych_dm psych_only assistant associate professor other if sampleshort=="experts", listwise
	eststo, title("Experts Completed Surveys"): estpost sum beh_econ_only beh_fin appl theory lab_econ psych_dm psych_only assistant associate professor other if sampleshort=="experts" & inlist(completed,"Yes","Partial"), listwise
	eststo, title("Experts Completed All 15 Tx"): estpost sum beh_econ_only beh_fin appl theory lab_econ psych_dm psych_only assistant associate professor other median_duration click_practi click_instr heard used_turk if sampleshort=="experts" & completed=="Yes", listwise
	esttab using output/TableSumStat_Experts.csv, mtitles csv main(mean) b(2) nostar nogaps unstack replace	

*Figure 6 - Output raw data for expert CDF graphs
use dtafiles/ExpertForecastCleanWide, clear

outsheet treatment_t4_actual treatment_t5_actual treatment_t6_actual treatment_t7_actual ///
		treatment_t8_actual treatment_t9_actual treatment_t13_actual treatment_t14_actual ///
		treatment_t10_actual treatment_t11_actual treatment_t12_actual ///
		treatment_t15_actual treatment_t16_actual treatment_t17_actual treatment_t18_actual if _n==1 ///
		using output/expert_cdf_data_actual.xls, replace

keep if sampleshort=="experts"
gen N = _n
keep N treatment_t1_actual treatment_t2_actual treatment_t3_actual

tempfile cdf_output
save `cdf_output'

forval i=4/18 {
use dtafiles/ExpertForecastCleanWide, clear
keep if sampleshort=="experts"
sort treatment_t`i'
gen N = _n
keep N treatment_t`i'
merge 1:1 N using `cdf_output'
drop _merge
save `cdf_output', replace
}

order N treatment_t4 treatment_t5 treatment_t6 treatment_t7 ///
		treatment_t8 treatment_t9 treatment_t13 treatment_t14 ///
		treatment_t10 treatment_t11 treatment_t12 ///
		treatment_t15 treatment_t16 treatment_t17 treatment_t18 ///
		treatment_t1_actual treatment_t2_actual treatment_t3_actual

outsheet using output/expert_cdf_data.xls, replace

****** Summarize Forecasts by Treatment ********

	*Table 3- summary of findings by treatment, Col 5-8
	use dtafiles/ExpertForecastCleanLong, clear

	keep if sampleshort=="experts"
	collapse (mean) mean_forecast=forecast error_mean_forecast=WoC_abs_error mean_ind_error= abs_error outperform_WoC ///
	(sd) sd_ind_error= abs_error sd_ind_forecast= forecast, by(tx_order treatmentname)
	outsheet tx_order treatmentname mean_forecast sd_ind_forecast using output/TableSumStatByTreatment_c56.xls, replace

	
	*Figures 5 and 7 - forecasts by treatment and by specialty
	use dtafiles/ExpertForecastCleanLong, clear

	keep if sampleshort=="experts"
	collapse (mean) mean_=forecast, by(tx_order treatmentname actual WoC_forecast std_econ lab_econ beh_econ psych)
	gen group = "std_econ" if std_econ==1
	replace group = "lab_econ" if lab_econ==1
	replace group = "beh_econ" if beh_econ==1
	replace group = "psych" if psych==1
	drop std-psych
	reshape wide mean_, i(tx_order treat actual WoC_forecast) j(group) string
	rename WoC_forecast all_experts
	outsheet tx_order treatmentname actual all_experts mean_std_econ mean_lab_econ mean_beh_econ mean_psych using output/MeanForecastsByTx.xls, replace
	
***************************************
***** Graph of SD MTurk Effort X SD Forecast
 

*save MTurk SD for each treatment 
	use dtafiles/MTurkCleanedData, clear
	collapse (sd) sd_treat = buttonpresses, by(tx_order treatmentname)
	tempfile MTurk_sds
	save `MTurk_sds'

use dtafiles/ExpertForecastCleanLong, clear
	
collapse (sd) sd_ex_forecast=forecast, by(tx_order treatno)  
	
merge 1:1 tx_order using `MTurk_sds'
 
 *reformat for graph
 foreach var in sd_ex_forecast sd_treat  {
 forvalues i = 4/18 {
 gen `var'`i' = .
 replace `var'`i' = `var' if treatno == `i'
 }
 }
 
 
*Online Appendix Figure 8
 foreach y in sd_ex_forecast {
 foreach x in sd_treat {
 local titlesd_ex_forecast "Standard Deviation of Expert Forecast"
 local titlesd_treat "Standard Deviation of MTurker Effort"
 local xtitle `title`x''
 local ytitle `title`y''
twoway (lfit `y' `x') ///
(scatter `y'4 `x'4, mcolor(forest_green) mlabel("treatment") mlabcolor(black) mlabposition(10)) ///
(scatter `y'5 `x'5, mcolor(forest_green) mlabel("treatment") mlabcolor(black) mlabposition(9)) ///
(scatter `y'6 `x'6, mcolor(brown) mlabel("treatment") mlabcolor(black) mlabposition(9)) ///
(scatter `y'7 `x'7, mcolor(brown) mlabel("treatment") mlabcolor(black) mlabposition(7)) ///
(scatter `y'8 `x'8, mcolor(lavender) mlabel("treatment") mlabcolor(black)) ///
(scatter `y'9 `x'9, mcolor(lavender) mlabel("treatment") mlabcolor(black) mlabposition(2)) ///
(scatter `y'10 `x'10, mcolor(maroon) mlabel("treatment") mlabcolor(black)) ///
(scatter `y'11 `x'11, mcolor(maroon) mlabel("treatment") mlabcolor(black)) ///
(scatter `y'12 `x'12, mcolor(maroon) mlabel("treatment") mlabcolor(black)) ///
(scatter `y'13 `x'13, mcolor(dkorange) mlabel("treatment") mlabcolor(black)) ///
(scatter `y'14 `x'14, mcolor(dkorange) mlabel("treatment") mlabcolor(black)) ///
(scatter `y'15 `x'15, mcolor(navy) mlabel("treatment") mlabcolor(black) mlabposition(6)) ///
(scatter `y'16 `x'16, mcolor(navy) mlabel("treatment") mlabcolor(black) mlabposition(2)) ///
(scatter `y'17 `x'17, mcolor(navy) mlabel("treatment") mlabcolor(black)) ///
(scatter `y'18 `x'18, mcolor(navy) mlabel("treatment") mlabcolor(black) mlabposition(6)), ///
graphregion(color(white)) ///
legend(off) ///
xtitle("`xtitle'") ytitle("`ytitle'")
graph export "output/`xtitle' vs `ytitle'.png", replace 
}
}


*************************************
* Analyze MTurk Effort - Minute by Minute
*************************************

*Minute by minute - recorded times are the remaining seconds left

* Online Appendix Figure 5

use dtafiles/MTurkCleanedData, clear	

keep tx_order times button
tostring times, replace
*find the first (doesn't have space in the beginning)
gen first = substr(times,1,strpos(times, ",")-1)
destring first, replace

*split out number of buttonpresses per seconds left
forval s=0/600 {
	egen sleft`s'=noccur(times), string(" `s',")
	quietly: replace sleft`s'= sleft`s'+1 if first==`s'
}

*One subject missing 10 button presses
*string of characters between 574 and 573
*assign half to 574 and half to 573 (characters in between)
egen total_check = rowtotal(sleft*)
gen check = abs( total_check - buttonpresses)
tab check
replace sleft574 = sleft574 + 5 if check==10
replace sleft573 = sleft573 + 5 if check==10
drop total_check check

*sum to minute level 
local min_minute = 0
forval m=1/10 {
	local max_minute = `m'*60
	local actual_min = 11 - `m'
	egen min`actual_min'= rowtotal(sleft`min_minute' - sleft`max_minute')
	local min_minute = `max_minute'+1
}

*save so we can look at this later
save dtafiles/MTurk_sec_data, replace	


*graph by minute
collapse (mean) min*, by(tx_order)
reshape long min, i(tx_order) j(minute)
rename min bp

reshape wide bp, i(minute) j(tx_order)

local txnum = 1
foreach t_label in "`tx_labels'" {
label var bp`txnum' "`t_label'"
local txnum = `txnum' + 1
}


label var minute "Minute in the Task"

graph twoway line bp1 bp5 bp2 bp4 bp3 bp8 bp18 bp17 bp16 minute, ///
	xlabel(1(1)10) ytitle("Button Presses") title("Minute-by-minute Effort by Treatment") ///
	legend(position(3) cols(1)) xsize(6.5) graphregion(color(white)) ///
	lpattern(dash dash dash dash dash solid solid solid solid) ///
	lcolor(navy*0.2 navy*0.4 navy*0.6 navy*0.8 navy dkorange*0.4 dkorange*0.6 dkorange*0.8 dkorange)
graph export "output/minute1.png", replace

graph twoway line bp1 bp6 bp7 bp14 bp15 bp9 bp10 bp11 bp12 bp13 minute, ///
	xlabel(1(1)10) ytitle("Button Presses") title("Minute-by-minute Effort by Treatment") ///
	legend(position(3) cols(1)) xsize(6.5) graphregion(color(white)) ///
	lpattern(dash dash_dot dash_dot solid solid dash dash solid solid solid) ///
	lcolor(navy*0.2 purple*0.5 purple dkgreen*0.5 dkgreen navy*0.5 navy dkorange*0.4 dkorange*0.7 dkorange)
graph export "output/minute2.png", replace
