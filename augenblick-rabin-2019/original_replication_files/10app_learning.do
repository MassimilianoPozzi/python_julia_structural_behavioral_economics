do setupAndDefinePrograms.do

****************************
***Cost curves over time: these are just hand taken from the specification in column (3) with date FE
****************************

clear

set obs 7
gen ddn=_n
gen gamma=2.117939
replace gamma=gamma+.1354211 if _n==1
replace gamma=gamma+.0166715 if _n==2
replace gamma=gamma-.0273651 if _n==3
replace gamma=gamma-.0435105 if _n==4
replace gamma=gamma-.0257543 if _n==5
replace gamma=gamma-.0167649 if _n==6
replace gamma=gamma-.0386978 if _n==7
gen phi=686.88
replace phi=phi+533.3638 if _n==1
replace phi=phi+80.39228 if _n==2
replace phi=phi-72.98571 if _n==3
replace phi=phi-143.7833 if _n==4
replace phi=phi-121.9545 if _n==5
replace phi=phi-114.7756 if _n==6
replace phi=phi-160.257 if _n==7
expand 100
bys ddn: gen tasks=_n

gen totalcost=(1/phi)*(1/gamma)*tasks^(gamma)
gen marginalcost=(1/phi)*tasks^(gamma-1)

local graph_total="twoway "
local graph_marginal="twoway "
forvalues i=1(1)7 {
	local colornum=((8-`i')*2)+0
	local graph_total="`graph_total' (line totalcost tasks if ddn==`i', sort lc(gs`colornum'))"
	local graph_marginal="`graph_marginal' (line marginalcost tasks if ddn==`i', sort lc(gs`colornum'))"
}
`graph_total', legend(off) xtitle("Tasks") ytitle("Total Cost (Dollars)") ///
text( 15.3 10 "Estimated aggregate cost curves across" "participation dates (darker=>later)" , fcolor(white) size(big) place(se) box just(center) margin(l+4 t+1 b+1))  
graph save t_total, replace
graph export "$savedir/total_cost.png", replace





**************************
********Gets the average # of tasks in different situations -> raw data -> creates pics
*************************
//drop attrition people -> otherwise creates issue
forvalues s=2(1)2 {

//forvalues s=1(1)3 {	
	use "Data/decisions_data_w_ind", clear
	//drop if decisiondatenum==1
	egen group = group(decisiondatenum type)
	qui tab wage, gen(wageFE_)
	xi: reg jobschosen i.group wageFE_* if ($data_sample1) & (${wid_sample`s'}), clu(wid) nocons

	//this doesn't use wage in predictions
	foreach var of varlist wageFE_* {
		qui sum `var'
		replace `var'=r(mean) 
	}
	capture drop mean_jc
	predict mean_jc
	capture drop se_jc
	predict se_jc, stdp
	capture drop jc_low
	gen jc_low = mean_jc - 1*se_jc
	capture drop jc_high
	gen jc_high = mean_jc + 1*se_jc

	//, legend() ///
	twoway (connected mean_jc decisiondatenum if type=="present", sort lpattern(dash) mcolor(gs12) color(gs12))  ///
	(connected mean_jc decisiondatenum if type=="future", lpattern(dash_dot) lcolor(gs9) mcolor(gs9) sort) ///
	(connected mean_jc decisiondatenum if type=="prediction", sort lpattern(solid) mcolor(gs1) color(gs1))  ///
	if (${wid_sample`s'}), ///
	legend(order(1 "Present" 2 "Future" 3 "Prediction") row(1) ring(0) position(1)) ///
	 xtitle("Decision Date Num") ytitle("Jobs Chosen")  xlab(1(1)7)   ysc(r(40(5)55)) ylabel(40(5)55) 

	 graph save t1, replace	
	

}

graph combine t1.gph t_total.gph, xsize(13)  iscale(1.3)    

graph export "$savedir/raw_est_learning.png", replace  width(2500)
