do setupAndDefinePrograms.do

//just create file I will append to
gen x=.
save Data/changingParams, replace

foreach i of numlist 4014 4015 4016 {

	use "Data/decisions_data_w_ind", clear
	genVarsNeededForMl
	setupMlGlobals

	global vars="(betahat:) (delta:) (phi:) (gamma:) (beta:) (alpha:)"
	global setupvars="sigma betahat delta phi gamma beta alpha"
	global wid_sample="${wid_sample2}"
	global data_sample="${data_sample1}"
	global init="ml init beta:_cons=1 betahat:_cons=1 gamma:_cons=2 delta:_cons=1 phi:_cons=200 sigma:_cons=50"
	
	//fixed effects
	cap drop ddn_*
	cap drop wdn_*
	cap drop wid_*
	qui tab decisiondatenum, gen(ddn_)
	qui tab workdatenum, gen(wdn_)
	qui tab wid if $wid_sample2, gen(wid_) 
	
	// Getting individual estimates for beta for each day - present to current future decisions. 
	// Can't do week 1 or week 7 because need both future and past decisions
	if (`i'==4014) {
		local est="beta"
		local startcount=90
		global data_sample="$data_sample1 & decisiondatenum>1 & decisiondatenum<7"
		global vars="(betahat:) (delta:) (phi: wid_* ddn_*) (gamma: ddn_*, nocons ) (beta: ddn_*, nocons) (alpha:)"
		global setupvars="sigma betahat delta phi gamma beta alpha"
		global init="ml init betahat:_cons=1 beta:ddn_2=1 beta:ddn_3=1 beta:ddn_4=1 beta:ddn_5=1 beta:ddn_6=1 delta:_cons=1 gamma:ddn_2=2 gamma:ddn_3=2 gamma:ddn_4=2 gamma:ddn_5=2 gamma:ddn_6=2 phi:_cons=500 sigma:_cons=50"
	}
	
	// Getting individual estimates for alpha for each day - present to current future decisions. 
	// Can't do week 1 because need decisions before and after mandatory tasks
	if (`i'==4015) {
		local est="alpha"
		local startcount=91
		global data_sample="$data_sample1 & decisiondatenum>1"
		global vars="(betahat:) (delta:) (phi: wid_* ddn_*) (gamma: ddn_*, nocons ) (beta:) (alpha: ddn_*, nocons)"
		global setupvars="sigma betahat delta phi gamma beta alpha"
		global init="ml init betahat:_cons=1 beta:_cons=1 delta:_cons=1 gamma:ddn_2=2 gamma:ddn_3=2 gamma:ddn_4=2 gamma:ddn_5=2 gamma:ddn_6=2 gamma:ddn_7=2 phi:_cons=500 sigma:_cons=50"
	}

	// Getting individual estimates for beta for each day - present to current future decisions. 
	// Can't do week 1 or 7 
	if (`i'==4016) {
		local est="betahat"
		local startcount=1
		global data_sample="$data_sample1 & decisiondatenum>1 & decisiondatenum<7"
		global vars="(betahat: ddn_*, nocons) (delta:) (phi: wid_* ddn_*) (gamma: ddn_*, nocons ) (beta:) (alpha:)"
		global setupvars="sigma betahat delta phi gamma beta alpha"
		global init="ml init betahat:ddn_2=1 betahat:ddn_3=1 betahat:ddn_4=1 betahat:ddn_5=1 betahat:ddn_6=1 beta:_cons=1 delta:_cons=1 gamma:ddn_2=2 gamma:ddn_3=2 gamma:ddn_4=2 gamma:ddn_5=2 gamma:ddn_6=2 gamma:ddn_7=2 phi:_cons=500 sigma:_cons=50"
	}

	ml model lf NBetaHatStructural1_agg (sigma: $netdistance effort wage prediction guessbonusamount today pb group_for_test sophfixedweight distanceToPayment distanceToWork=) $vars ///
	if $wid_sample & $data_sample, clu(wid) technique(nr) 
	
	$init
	ml max

	clear
	set obs 8
	gen decisiondatenum_=_n 
	gen `est'=.
	gen `est'h=.
	gen `est'l=.
	forvalues w=2(1)7 {
		matrix V=e(V)
		local se=V[(`startcount'+`w'),(`startcount'+`w')]^(1/2)
		disp `se'
		replace `est'=[`est']ddn_`w' if _n==`w'
		replace `est'h=[`est']ddn_`w'+1.96*`se' if _n==`w'
		replace `est'l=[`est']ddn_`w'-1.96*`se' if _n==`w'
	}	
	cap append using Data/changingParams
	duplicates drop
	save Data/changingParams, replace
}



use Data/changingParams, clear
local est="beta"
twoway (rcap `est'l `est'h decisiondatenum_, color(gs10)) ///
(connected `est' decisiondatenum_, sort lpattern(solid) color(gs1))  ///
if decisiondatenum_<7 & decisiondatenum_>1 ///
, legend(off ) ///
 xtitle("Decision Date Num") ytitle("Estimate of Beta") xlab(2(1)6) ylab(.5(.1)1.1)  yline(1, lpattern(-) lcolor(gs5))
graph save Temp/t4, replace	

local est="betahat"
twoway (rcap `est'l `est'h decisiondatenum_, color(gs10)) ///
(connected `est' decisiondatenum_, sort lpattern(solid) color(gs1))  ///
if decisiondatenum_<7 & decisiondatenum_>1 ///
, legend(off ) ///
 xtitle("Decision Date Num") ytitle("Estimate of Beta-hat") xlab(2(1)6) ylab(.5(.1)1.1)  yline(1, lpattern(-) lcolor(gs5)) 
graph save Temp/t5, replace	

local est="alpha"
twoway (rcap `est'l `est'h decisiondatenum_, color(gs10)) ///
(connected `est' decisiondatenum_, sort lpattern(solid) color(gs1))  ///
if decisiondatenum_<8 & decisiondatenum_>1 ///
, legend(off ) ///
 xtitle("Decision Date Num") ytitle("Estimate of Alpha") xlab(2(1)6) yline(0, lpattern(-) lcolor(gs5))
graph save Temp/t6, replace	
	
	
graph combine Temp/t4.gph Temp/t5.gph Temp/t6.gph, ///
rows(1) ysize(4) graphregion(margin(zero)) xsize(13) iscale(1.3) 
graph export "$savedir/stability_all.png", replace  width(2500)



