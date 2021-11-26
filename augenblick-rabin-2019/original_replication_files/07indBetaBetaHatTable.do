do setupAndDefinePrograms.do
use "Data/decisions_data_w_ind", clear

*******************************************
********Get relationship between bi and bhi
*******************************************

foreach i of numlist 1 2 3 4 5 {
	cap drop biSE_inv`i' bhiSE_inv`i'
	gen biSE_inv`i'=1/(biSE`i'^2)
	gen bhiSE_inv`i'=1/(bhiSE`i'^2)
}

//these care clearly wrong because SE abnormally low or high due to stopping convergence early
//unfortunately, the way-too-big ones are exactly the worse!
replace bhiSE_inv5=10 if bhiSE_inv5>100000


//get everything
eststo clear
local estadded=""
foreach i of numlist 1 2 3 4 5 {

	local cond="(first==1 & ml_prob`i'==0)"

	cap drop bi`i'dev bhi`i'dev
	gen bi`i'dev=1-bi`i'
	gen bhi`i'dev=1-bhi`i'
	
	eststo est`i': estpost sum bi`i' if `cond'
	unique wid if `cond'
	estadd local N_clust="`r(N)'", replace
	
	est store base
	reg bhi`i'dev bi`i'dev  if `cond', r nocons
	getEstimates 1 `i'

	est store base
	reg bhi`i'dev bi`i'dev [aweight=bhiSE_inv`i'] if `cond', r nocons
	getEstimates 2 `i'

	est store base
	robreg mm bhi`i'dev bi`i'dev  if `cond', 
	getEstimates 3 `i'
	
	est store base
	robreg m bhi`i'dev bi`i'dev  if `cond', 
	getEstimates 4 `i'
	
	est store base
	robreg s bhi`i'dev bi`i'dev  if `cond', 
	getEstimates 5 `i'
	
	est store base
	eivreg bhi`i'dev bi`i'dev if `cond', r(bi`i'dev .5)
	getEstimates 6 `i'
	
	//let's get reliability
	gen biVAR`i'=biSE`i'^2
	qui sum biVAR`i' if `cond'
	local noiseVar=r(mean)
	sum bi`i' if `cond'
	local totalVar=r(sd)^(2)
	local reliability=(`totalVar'-`noiseVar')/(`totalVar')
	
	est store base
	eivreg bhi`i'dev bi`i'dev if `cond' , r(bi`i'dev `reliability')
	getEstimates 7 `i'
	local rel: display %9.2f `reliability'
	local rel_formatted="\emph{"+"`rel'"+"}"
	estadd local rel `rel_formatted'
	

}

//make table

esttab `estadded' /// 
using "$savedir/betabetahatcorr_new.tex", ///
stats(beta1 betaSD1 beta2 betaSD2 beta3 betaSD3 beta4 betaSD4 beta5 betaSD5 beta6 betaSD6 beta7 betaSD7 rel N_clust, star labels( ///
"\emph{OLS}"  " "  ///
"\emph{GLS}"  " "  ///
"\emph{MM-estimation}"  " "  ///
"\emph{M-estimation}"  " "  ///
"\emph{S-estimation}"  " "  ///
"\emph{EIV \$\rho_{xx'}\$:=0.5}"  " "  ///
"\emph{EIV w/ est \$\rho_{xx'}\$}"  " "  ///
"\emph{(est \$\rho_{xx'}\$)}"  ///
"\midrule Observations"  ///
)) ///
mtitle( "\specialcell{Primary \\ Estimation}" ///
"\specialcell{Early \\ Decisions}" ///
"\specialcell{Later \\ Decisions}" ///
"\specialcell{Projection\\ Bias}" ///
"\specialcell{With \\ Sophistication}" ///
) ///
cells( b(fmt(%9.3f)) se(par fmt(%9.3f))) ///
width(1.2\hsize) replace booktabs gaps f collabels(none) eqlabels(none) ///
style(tex) legend label compress substitute(\_ _) 

