do setupAndDefinePrograms.do


use "Data/decisions_data", clear
cap drop ind_diff_pres ind_diff_pred
bys wid: egen ind_var=sd(jobschosen)
replace ind_var=ind_var^2
gen ind_diff_pres=.
gen ind_diff_pred=.
gen ind_diff_pred_var=.
gen ind_diff_pres_var=.
gen ind_diff_pres_var_inv=.
gen ind_diff_pred_var_inv=.
egen group = group(wage_group)
	
**** Calculations of non-parametric measure for each individual

forvalues wid=1(1)100 {
	
	xi: qui reg jobschosen today i.group if (type=="present" | type=="future") & wid==`wid' & $data_sample1
	matrix define coeffs=e(b)
	local present=coeffs[1,1]
	matrix define V=e(V)
	local var=V[1,1]
	replace ind_diff_pres=`present' if wid==`wid'
	replace ind_diff_pres_var=`var' if wid==`wid'
	replace ind_diff_pres_var_inv=1/`var' if wid==`wid'
}

forvalues wid=1(1)100 {
	xi: qui reg jobschosen prediction i.group if (type=="prediction" | type=="future") & wid==`wid' & $data_sample1
	matrix define coeffs=e(b)
	local prediction=coeffs[1,1]
	matrix define V=e(V)
	local var=V[1,1]
	replace ind_diff_pred=`prediction' if wid==`wid'
	replace ind_diff_pred_var=`var' if wid==`wid'
	replace ind_diff_pred_var_inv=1/`var' if wid==`wid'
}

label variable ind_var "Individual variance of jobsChosen"
label variable ind_diff_pres "Individual avg diff b/t present and future (controlling for wagegroup)"
label variable ind_diff_pres_var "var: ind_diff_pres"
label variable ind_diff_pres_var_inv "inverse var: ind_diff_pres"
label variable ind_diff_pred "Individual avg diff b/t predict and future (controlling for wagegroup)"
label variable ind_diff_pred_var "var: ind_diff_pred"
label variable ind_diff_pred_var_inv "inverse var: ind_diff_pred"
cap drop group
cap drop _I*
save "Data/decisions_data_w_non_para", replace
