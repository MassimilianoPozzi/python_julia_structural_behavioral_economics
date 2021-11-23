do setupAndDefinePrograms.do

**************************
********Identifying Observations w/ ML Issues
**************************
use "Data/ind_results", clear


foreach i of numlist 1/5 {
	cap drop ml_prob`i' 
	gen ml_prob`i'=0
	label variable ml_prob`i' "Is there an issue with this person's ml estimates for routine `i'?"
	
	//so, we have a problem if we don't get estimates
	replace ml_prob`i'=1 if bi`i'==. | iters`i'>=200 | biSE`i'==0 | bhiSE`i'==0 
	
	//then, see if estimates are very far outliers (grubbs test)
	local stats "bi bhi gi ai"
	disp "------------> i=`i'"
	tab ml_prob`i'
	sum bhi`i' if ml_prob`i'==0
    foreach stat of local stats {
		label variable `stat'`i' "individual estimate : `stat' running routine `i'"
		//used to do: grubbs1 `stat'`i' if ml_prob`i'==0 & first==1, level(99.99)
		//but fails in stata 13 for some reason
		//so just write own code:
		qui gen grubbs_`stat'`i'=0
		label variable grubbs_`stat'`i' "Is ind stat rejected by grubbs test?"
		local conf=(100-99.99)/100
		forvalues  iters=1(1)6 {
			qui sum `stat'`i' if ml_prob`i'==0 
			qui gen t_centered = (abs(`stat'`i' -r(mean)))/r(sd) 
			qui local cutoff = (r(N)-1)*sqrt(invttail(r(N)-2,`conf'/(2*r(N)))^2/(r(N)*(r(N)-2+invttail(r(N)-2,`conf'/(2*r(N)))^2)))
			qui replace grubbs_`stat'`i'=1 if t_centered>`cutoff'
			qui drop t_centered
			cap replace ml_prob`i'=1 if grubbs_`stat'`i'==1
		}		
	}
	tab ml_prob`i'
	sum bhi`i' if ml_prob`i'==0

}

//main one we use
gen ml_prob=ml_prob1
label variable ml_prob "Is there an issue with this person's ml estimates?"
tab ml_prob

//merge results
save "Data/temp", replace
use "Data/decisions_data_w_non_para", replace
merge wid using "Data/temp", sort uniqusing 
drop _merge
save "Data/decisions_data_w_ind", replace 
