do setupAndDefinePrograms.do
use "Data/decisions_data_w_ind", clear

**************************
********Individual Estimates Analysis
**************************


foreach i of numlist 1/5 {
	sum bhi`i' bi`i' gi`i' ai`i' pi`i' if ml_prob`i'==0 & first==1
	//sum bhi`i' bi`i' gi`i' ai`i' if iters`i'!=200 & first==1
	//sum bhi`i' bi`i' gi`i' ai`i' if iters`i'!=200 & first==1 &ml_prob`i'==1
	//pwcorr bi`i' bhi`i' if ml_prob`i'==0 & first==1, sig obs
	//pwcorr bi`i' bhi`i' if ml_prob`i'==0 & first==1 & bhi`i'<1, sig obs
}

********Weights based on SD of estimates

eststo clear
local estadded=""
foreach i of numlist 1/5 {
	local estadded="`estadded' est`i'" 
	
	local cond="(first==1 & ml_prob`i'==0)"

	eststo est`i': estpost sum bi`i' 

	local stats "bi bhi di gi"
	if (`i'==4) local stats "bi bhi di gi ai"
     
	foreach x of local stats {
	      
		sum `x'`i' if `cond', detail

		local temp: display %9.3f r(mean)		
		estadd local `x'mean="`temp'"
		local temp: display %9.3f r(p50)		
		estadd local `x'med="`temp'"
		local temp: display %5.3f r(sd)		
		estadd local `x'sd="("+"`temp'"+")"
      	
	}	

	estadd local N="`r(N)'", replace

	//correlations
	corr bi`i' bhi`i' if first==1 & ml_prob`i'==0
	local temp: display %5.3f r(rho)
	estadd local corr_bi_bhi="`temp'"
	//pwcorr doesn't store?
	local p=min(2*ttail(r(N)-2,abs(r(rho))*sqrt(r(N)-2)/sqrt(1-r(rho)^2)),1)
	local temp: display %5.3f `p'
	estadd local corr_sig="`temp'"
	
	//percentage less than zero
	gen bi`i'_l1=(bi`i'<1)
	gen bhi`i'_l1=(bhi`i'<1)
	sum bi`i'_l1  if first==1 & ml_prob`i'==0
	local temp: display %9.2f r(mean)		
	estadd local bi_l0="`temp'"
	sum bhi`i'_l1  if first==1 & ml_prob`i'==0
	local temp: display %9.2f r(mean)		
	estadd local bhi_l0="`temp'"
								
}


local stats "bi bhi di gi ai"
local statsadded=""
foreach x of local stats {
	if ("`x'"=="bi" | "`x'"=="bhi") {
		local statsadded="`statsadded' `x'mean `x'med `x'sd" 
	} 
	else {
		local statsadded="`statsadded' `x'mean `x'med `x'sd" 
	}
}
local statsadded="`statsadded' bi_l0 bhi_l0 corr_bi_bhi corr_sig" 
	

esttab `estadded' /// 
using "$savedir/individual_new.tex", ///
 ///
stats(`statsadded' N, labels( ///
"mean(\$\hat{\beta_{i}}\$)" "median(\$\hat{\beta_{i}}\$)"  "sd(\$\hat{\beta_{i}}\$)"  ///
"\\ mean(\$\hat{\beta_{h,i}}\$)" "median(\$\hat{\beta_{h,i}}\$)"  "sd(\$\hat{\beta_{h,i}}\$)"  ///
"\\ mean(\$\hat{\delta_{i}}\$)" "median(\$\hat{\delta_{i}}\$)"  "sd(\$\hat{\delta_{i}}\$)"  ///
"\\ mean(\$\hat{\gamma_{i}}\$)" "median(\$\hat{\gamma_{i}}\$)"  "sd(\$\hat{\gamma_{i}}\$)"  ///
"\\ mean(\$\hat{\alpha_{i}}\$)" "median(\$\hat{\alpha_{i}}\$)"  "sd(\$\hat{\alpha_{i}}\$)"  ///
"\midrule P[\$\hat{\beta_{i}}\$]<1" "P[\$\hat{\beta_{h,i}}\$]<1" "r(\$\hat{\beta_{i}}\$, \$\hat{\beta_{h,i}}\$)" ///
"p-value r(\$\hat{\beta_{i}}\$, \$\hat{\beta_{h,i}}\$)"  ///
"\midrule Observations" ///
fmt(2))) ///
mtitle( "\specialcell{Primary \\ Estimation}" ///
"\specialcell{Early \\ Decisions}" ///
"\specialcell{Later \\ Decisions}" ///
"\specialcell{Proj. \\ Bias}" ///
"\specialcell{Pred. \\ Soph.}" ///
) ///
cells( b(fmt(%9.3f)) se(par fmt(%9.3f))) ///
width(1.2\hsize) replace booktabs gaps f collabels(none) eqlabels(none) ///
style(tex) legend label compress substitute(\_ _) 



