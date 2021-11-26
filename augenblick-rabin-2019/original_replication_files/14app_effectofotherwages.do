do setupAndDefinePrograms.do

*** Asked by refereee: does wage mix (other wages in the decision set) effect decisions? 
*** The results are just mentioned in a footnote!
** Multiple possible measures of other wages
** Below basically does the same thing four times with different measures

**1) Measure 1: avg of other decisions
use "Data/decisions_data_w_ind", clear
cap drop t_*
bys group_decision_set (wage): egen avg_of_wages_in_ds=mean(wage)
gen avg_of_other_wages_in_ds=(avg_of_wages_in_ds*5-wage)/4

egen wageGroup=cut(avg_of_other_wages_in_ds),group(10)
bys wageGroup: egen wageGroupAvg=mean(avg_of_other_wages_in_ds)

//Problem: wages are not completely random - for $.02 wage, then need higher wage in set so throw out sets with only lower wages 
//=> wages and avg_wage_other negatively correlated
corr avg_of_other_wages_in_ds wage
reg avg_of_other_wages_in_ds wage

//so, need to control for wage by looking at the residuals after controlling for wage
xi: reg jobschosen i.wage if ($data_sample1) & (${wid_sample2}), clu(wid)
predict jobschosen_wageRes, residuals

//just linear reg: no effect
xi: reg jobschosen_wageRes avg_of_other_wages_in_ds if ($data_sample1) & (${wid_sample2}), clu(wid)

//non parametrically: look at graph
xi: reg jobschosen_wageRes i.wageGroup if ($data_sample1) & (${wid_sample2}), clu(wid)

capture drop mean_jc
predict mean_jc
capture drop se_jc
predict se_jc, stdp
capture drop jc_low
gen jc_low = mean_jc - 1*se_jc
capture drop jc_high
gen jc_high = mean_jc + 1*se_jc

twoway (rcap jc_low jc_high wageGroupAvg, color(gs10) ) ///
(connected mean_jc wageGroupAvg, sort color(gs8) lpattern(solid))  ///
(hist avg_of_other_wages_in_ds, bin(20) yaxis(2) )  ///
if (${wid_sample2}) ///
, legend(off) ///
 xtitle("Average of Other Wages in Set") yscale(range(-10 5) axis(1)) yscale(range(0 50) axis(2))  ///
 ytitle("Jobs Chosen (Residual cont. for wages)")  scheme(s1mono)
graph save Output/t1, replace
graph export "$savedir/effect_other_wages1.png", replace width(2000)

//does controlling for avg effect results - given above, probably not:
xi: reg jobschosen future if (type=="present" | type=="future")    & ($data_sample1) & (${wid_sample2}), clu(wid)
xi: reg jobschosen future i.wageGroup if (type=="present" | type=="future")    & ($data_sample1) & (${wid_sample2}), clu(wid)





**2) Measure 2: what is the highest of the other wages?
use "Data/decisions_data_w_ind", clear
cap drop t_*
bys group_decision_set (wage): gen t_1=wage[5]
bys group_decision_set (wage): gen t_2=wage[4]
order  group_decision_set t_1 t_2 wage
gen highest_other_wages=t_1
//if this wage is highest, choose 2nd highest
replace highest_other_wages=t_2 if wage==t_1

//so, need to control for wage by looking at the residuals
xi: reg jobschosen i.wage if ($data_sample1) & (${wid_sample2}), clu(wid)
predict jobschosen_wageRes, residuals

//just linear reg
xi: reg jobschosen_wageRes highest_other_wages if ($data_sample1) & (${wid_sample2}), clu(wid)

//non parametrically
xi: reg jobschosen_wageRes i.highest_other_wages if ($data_sample1) & (${wid_sample2}), clu(wid)

capture drop mean_jc
predict mean_jc
capture drop se_jc
predict se_jc, stdp
capture drop jc_low
gen jc_low = mean_jc - 1*se_jc
capture drop jc_high
gen jc_high = mean_jc + 1*se_jc

twoway (rcap jc_low jc_high highest_other_wages, color(gs10) ) ///
(connected mean_jc highest_other_wages, sort color(gs8) lpattern(solid))  ///
(hist highest_other_wages, discrete yaxis(2) )  ///
if (${wid_sample2}) & highest_other_wages>.10 & highest_other_wages<.31 ///
, legend(off) ///
 xtitle("Highest of other wages in set") yscale(range(-10 5) axis(1)) yscale(range(0 60) axis(2))  ///
 ytitle("Jobs Chosen (Residual cont. for wages)")  scheme(s1mono)
graph save Output/t2, replace
graph export "$savedir/effect_other_wages2.png", replace  width(2000)


//NOW VS LATER
//USED -> basic
xi: reg jobschosen future if (type=="present" | type=="future")    & ($data_sample1) & (${wid_sample2}), clu(wid)
xi: reg jobschosen future i.highest_other_wages if (type=="present" | type=="future")    & ($data_sample1) & (${wid_sample2}), clu(wid)





**3) what is the lowest of the other wages?
use "Data/decisions_data_w_ind", clear
cap drop t_*
bys group_decision_set (wage): gen t_1=wage[1]
bys group_decision_set (wage): gen t_2=wage[2]
order  group_decision_set t_1 t_2 wage
gen lowest_other_wages=t_1
//if this wage is highest, choose 2nd highest
replace lowest_other_wages=t_2 if wage==t_1

//need to control for wage by looking at the residuals
xi: reg jobschosen i.wage if ($data_sample1) & (${wid_sample2}), clu(wid)
predict jobschosen_wageRes, residuals

//just linear reg
xi: reg jobschosen_wageRes lowest_other_wages if ($data_sample1) & (${wid_sample2}), clu(wid)

//non parametrically
xi: reg jobschosen_wageRes i.lowest_other_wages if ($data_sample1) & (${wid_sample2}), clu(wid)

capture drop mean_jc
predict mean_jc
capture drop se_jc
predict se_jc, stdp
capture drop jc_low
gen jc_low = mean_jc - 1*se_jc
capture drop jc_high
gen jc_high = mean_jc + 1*se_jc

twoway (rcap jc_low jc_high lowest_other_wages, color(gs10) ) ///
(connected mean_jc lowest_other_wages, sort color(gs8) lpattern(solid))  ///
(hist lowest_other_wages, discrete yaxis(2) )  ///
if (${wid_sample2}) & lowest_other_wages<.21 ///
, legend(off) ///
 xtitle("Lowest of other wages in set") yscale(range(-10 5) axis(1)) yscale(range(0 60) axis(2))  ///
 ytitle("Jobs Chosen (Residual cont. for wages)")  scheme(s1mono)
graph save Output/t3, replace
graph export "$savedir/effect_other_wages3.png", replace  width(2000)



//NOW VS LATER
//USED -> basic
xi: reg jobschosen future if (type=="present" | type=="future")    & ($data_sample1) & (${wid_sample2}), clu(wid)
xi: reg jobschosen future i.lowest_other_wages if (type=="present" | type=="future")    & ($data_sample1) & (${wid_sample2}), clu(wid)








**4) where does the wage rank?
//so, if we have the same wage, but sometimes it is relatively large or small, does that effect outcome?
use "Data/decisions_data_w_ind", clear
cap drop t_*


bys group_decision_set (wage): gen rank_of_wage_in_ds=_n
egen g_wage_rank=group(wage rank)
bys g_wage_rank: gen N=_N
//scatter rank wage [w=N], msymbol(circle_hollow)

//so, need to control for wage by looking at the residuals
xi: reg jobschosen i.wage if ($data_sample1) & (${wid_sample2}), clu(wid)
predict jobschosen_wageRes, residuals

//just linear reg
xi: reg jobschosen_wageRes rank_of_wage_in_ds if ($data_sample1) & (${wid_sample2}), clu(wid)

//non parametrically
xi: reg jobschosen_wageRes i.rank_of_wage_in_ds if ($data_sample1) & (${wid_sample2}), clu(wid)

capture drop mean_jc
predict mean_jc
capture drop se_jc
predict se_jc, stdp
capture drop jc_low
gen jc_low = mean_jc - 1*se_jc
capture drop jc_high
gen jc_high = mean_jc + 1*se_jc

twoway (rcap jc_low jc_high rank_of_wage_in_ds, color(gs10) ) ///
(connected mean_jc rank_of_wage_in_ds, sort color(gs8) lpattern(solid))  ///
(hist rank_of_wage_in_ds, bin(20) yaxis(2) )  ///
if (${wid_sample2}) ///
, legend(off) ///
 xtitle("Rank in Decision Set") yscale(range(-10 5) axis(1)) yscale(range(0 5) axis(2))  ///
 ytitle("Jobs Chosen (Residual cont. for wages)")  scheme(s1mono)
graph save Output/t4, replace
graph export "$savedir/effect_other_wages4.png", replace  width(2000)


//NOW VS LATER
//USED -> basic
xi: reg jobschosen future if (type=="present" | type=="future")    & ($data_sample1) & (${wid_sample2}), clu(wid)
xi: reg jobschosen future i.rank_of_wage_in_ds if (type=="present" | type=="future")    & ($data_sample1) & (${wid_sample2}), clu(wid)















