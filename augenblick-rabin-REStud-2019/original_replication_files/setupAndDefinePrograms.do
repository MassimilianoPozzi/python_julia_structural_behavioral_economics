************
** NOTE: Change these directories to match your own 
************
clear
clear programs 
macro drop _all 
set more off
set matsize 11000

//change these directories!
global direct="PUT YOUR FOLDER HERE"
global savedir="PUT YOUR SAVE FOLDER HERE"
cd "$direct"

** First, define some subsamples we will use over and over
global wid_sample1 "(1==1)"								//everyone
global wid_sample2 "(ml_prob==0)"						//main sample: people without ML estimation issues
global wid_sample3 "(ml_prob==0) & (left==0)"			//main sample minus attritors
global wid_sample4 "(left==0)"							//everyone minus attritors
global wid_sample5 "(ml_prob==1)"						//people with ML estimation issues
global wid_sample6 "(1==1)"								//everyone


global data_sample1 "(bonusoffered==0)"											//main decision sample: can't use work decision data if person offered bonus for that decision
global data_sample2 "(bonusoffered==0) & (bonusoffered_in_set==0)"				//main sample minus any decision set in which a bonus was offered for any wage
global data_sample3 "(bonusoffered==0) & (first_p==1 | first_f==1 | today==1)"	//main sample but only includes first prediction set and first future set
global data_sample4 "(bonusoffered==0) & (first_pf==1 | today==1)"				//main sample but only includes first set that is either prediction or future


program define setupMlGlobals
	global loud1=0
	global loud2=0
	global tobit=1
	global ln=0
	global soph=0
	global sophfixed=0
	global pb=1
	global structuralpb1=0
	global structuralpb2=0
	global left=0
	global uncertainty=0
	global separatedelta=0
	global thirddegpoly=0
	global seconddegpoly=0
	global riskaversionalpha=0
	global fixedbeta1=0
	global fixeddelta1=0
	global error_gamma=0
	global error_phi=0
	global log_error=0
	global double_error=0
	global nonpara=0
	global nonrtol=0
	global later=0
	global testdifference=0
	global monetarydelta=0
	global widfe=1
	global dayfe=0
	global netdistance="netdistance"
end


program define genVarsNeededForMl
	gen pb=workdone1
	gen group_for_test=1
	gen sophfixedweight=1
	gen distanceToPayment=paymentdate-decisiondate
	gen distanceToWork=workdate-decisiondate
end


program define setIndInitUsingBaseResults
	local i=`1'
	local wid=`2'
	//if we are not in baseline, let's use past estimates as initial to speed up convergence and not hit bad numerical derivatives
	if (`i'!=1) {
		preserve
		use "Data/ind_results", replace
		
		//if we didn't find an answer last time, let's use default
		sum iters1 if wid==`wid'
		if (r(N)==0 | r(mean)==200) {
			//hold at default
		} 
		else {
			//use past results
			local storedVars "bi1 bhi1 gi1 di1 pi1 si1"
			local initVars "beta betahat gamma delta phi sigma"
			global init="ml init "
			forvalues n=1(1)6 {
				local storedVar `: word `n' of `storedVars''
				local initVar `:word `n' of `initVars''
				qui sum `storedVar' if wid==`wid'
				local mean=r(mean)
				global init="$init `initVar':_cons=`mean'"
			}
		}
		restore
	}
end

program define createIndResultsFileIfMissing
	capture confirm file "Data/ind_results.dta"
	if _rc!=0 {
		preserve
		clear
		set obs 100
		gen wid=_n
		sort wid
		save "Data/ind_results", replace
		restore
	}
end

program define getAndStoreIndEst
	local i=`1'
	local wid=`2'
	
	matrix define b=e(b)
	matrix list b
	local si=b[1,1]
	local bhi=b[1,2]
	local di=b[1,3]
	local gi=b[1,4]
	local bi=b[1,5]
	local pi=b[1,6]
	cap local ai=b[1,7]
	
	matrix define V=e(V)
	local siSE=V[1,1]^(1/2)
	local bhiSE=V[2,2]^(1/2)
	local diSE=V[3,3]^(1/2)
	local giSE=V[4,4]^(1/2)
	local biSE=V[5,5]^(1/2)
	local piSE=V[6,6]^(1/2)
	cap local aiSE=V[7,7]^(1/2)

	local iters=e(ic)
	
	local vars "gi di bi bhi pi ai si giSE diSE biSE bhiSE piSE aiSE siSE iters" 
	foreach var of local vars {
		cap gen `var'`i'=``var'' if wid==`wid'	
	}
 
	//store estimates
	keep wid gi`i'-iters`i'
	duplicates drop wid, force
	save "Data/temp", replace
	cap use "Data/ind_results", clear
	//sometimes have a weird problem where computer is running to fast and cant open this file. Dropbox?
	if _rc!=0 {
		disp "******************Computer running too fast?"
		sleep 3000
		use "Data/ind_results", clear
	}
	merge wid using "Data/temp", sort uniqusing update replace
	drop _merge
	cap save "Data/ind_results", replace
	//sometimes have a weird problem where computer is running to fast and cant open this file. Dropbox?
	if _rc!=0 {
		disp "******************Computer running too fast?"
		sleep 3000
		save "Data/ind_results", replace
	}
end


program define generateFE
	local s=`1'
	cap drop ddn_*
	cap drop wdn_*
	qui tab decisiondatenum, gen(ddn_)
	qui tab workdatenum, gen(wdn_)
	//can't use individual FE for people not in main sample because their individual FE won't converge-> these guys have to be combined
	qui tab wid if $wid_sample2, gen(wid_) 
	foreach var of varlist wid_* {
		qui replace `var' = 0 if `var'==.
    }
	gen wid_all=(!($wid_sample2))
	if (`s'==2) {
		//needed because we use collinear to get average FE level below
		drop wid_all
	}
end

program define generateConstraintsSoConsIsAvg
	****This inelegantly forces constant to be the average FE - not sure if any other way
	//note need "collinear" for this to work in the ml command
	foreach type in "ddn" "wdn" {
		foreach param in "gamma" "phi" "psi1" "psi2" "psi3" {
			local constraint_`type'_`param'="0"
			foreach var of varlist `type'_* {
				local constraint_`type'_`param'="`constraint_`type'_`param''+[`param']`var'"
			}
			local constraint_`type'_`param'="`constraint_`type'_`param''=0"
		}
	}
	foreach param in "gamma" "phi"  "psi1" "psi2" "psi3" {
		local constraint_wid_`param'="0"
		foreach var of varlist wid_* {
			local constraint_wid_`param'="`constraint_wid_`param''+[`param']`var'"
		}
		local constraint_wid_`param'="`constraint_wid_`param''=0"
	}
	
	constraint drop _all
	constraint define 1 `constraint_ddn_gamma'
	constraint define 2 `constraint_ddn_phi'
	constraint define 3 `constraint_wid_phi'
	//only used for robustness checks
	constraint define 4 `constraint_wid_gamma'
	//only used in appendix
	constraint define 5 `constraint_wdn_gamma'
	constraint define 6 `constraint_wdn_phi'	
	constraint define 7 `constraint_delta_m_e'
	
	constraint define 8 `constraint_ddn_psi1'
	constraint define 9 `constraint_ddn_psi2'
	constraint define 10 `constraint_ddn_psi3'
	constraint define 11 `constraint_wid_psi2'
	constraint define 12 `constraint_wid_psi1'
end


program define storeEstVarsForTables
	local i=`1'

	if ($tobit==1) {
		estadd local tobit "X", replace
	}
	if ($widfe==1) {
		estadd local widfe "X", replace
	}
	if ($dayfe==1) {
		estadd local dayfe "X", replace
	}
	if ($soph==1) {
		estadd local soph "X", replace
	}
	if ($later==1) {
		estadd local later "X", replace
	}
	
	
	//just for FE regressions
	if (`i'==11001) { 
		estadd local gamma_wid "X", replace
	}
	if (`i'==11002) { 
		estadd local phi_wid "X", replace
		estadd local phi_day "X", replace
	}
	if (`i'==11003) { 
		estadd local phi_wid "X", replace
		estadd local gamma_day "X", replace
	}
	if (`i'==11004) { 
		estadd local gamma_wid "X", replace
		estadd local phi_day "X", replace
	}
	if (`i'==11005) { 
		estadd local gamma_wid "X", replace
		estadd local gamma_day "X", replace
	}	
	
	if (`i'==11006) { 
		estadd local gamma_wid "X", replace
		estadd local gamma_day "X", replace
		estadd local phi_day "X", replace	
	}	
	
	
	cap {
		test [beta]_cons==1
		if `r(p)'>=.001 {
			local temp: display %9.3f r(p)
			local temp "p=`temp'"
		} 
		else {
			local temp "p<0.001"
		}
		estadd local betatest `temp', replace

		test [betahat]_cons==1
		local temp: display %9.2f r(p)
		local temp "p=`temp'"
		estadd local betahattest `temp', replace

		test [delta]_cons==1
		local temp: display %9.2f r(p)
		local temp "p=`temp'"
		estadd local deltatest `temp', replace
	}
	
	cap {
		test [alpha]_cons==0
		if `r(p)'>=.001 {
			local temp: display %9.3f r(p)
			local temp "p=`temp'"
		} 
		else {
			local temp "p<0.001"
		}
		estadd local alphatest `temp', replace
	}

	
	if ($monetarydelta) {	
		
		cap {
			test [delta_e]_cons==1
			if `r(p)'>=.001 {
				local temp: display %9.3f r(p)
				local temp "p=`temp'"
			} 
			else {
				local temp "p<0.001"
			}
			estadd local deltaetest `temp', replace
		}
			
		if ($includeDeltaM1TestInTable) {
			cap {
				test [delta_m]_cons==1
				if `r(p)'>=.001 {
					local temp: display %9.3f r(p)
					local temp "p=`temp'"
				} 
				else {
					local temp "p<0.001"
				}
				estadd local deltamtest `temp', replace
			}
		}

		cap {
			if ($includeDeltaMDeltaETestInTable) {
				test [delta_m]_cons==[delta_e]_cons
				if `r(p)'>=.001 {
					local temp: display %9.3f r(p)
					local temp "p=`temp'"
				} 
				else {
					local temp "p<0.001"
				}
				estadd local deltametest `temp', replace
			}
		}
	}
	
	cap {
		test [theta]_cons==0
	}
	
	if ($separatedelta) {
		test [delta1]_cons==1
		test [delta2]_cons==1
		test [delta3]_cons==1
		test [delta4]_cons==1
		
		test [delta1]_cons==[delta2]_cons
		test [delta1]_cons==[delta3]_cons
		test [delta1]_cons==[delta4]_cons
		test [delta2]_cons==[delta3]_cons
		test [delta2]_cons==[delta4]_cons
		test [delta3]_cons==[delta4]_cons
	}
end
	
	
	


program define formatData, rclass
	local b=b[1,`1']
	local V=V[`1',`1']^(1/2)
	local stars ""
	if (abs(`b'/`V')>=2.58) local stars "***"
	if (abs(`b'/`V')<2.58 & abs(`b'/`V')>=1.96) local stars "**"
	if (abs(`b'/`V')<1.96 & abs(`b'/`V')>=1.65) local stars "*"
	local b_: display %9.3f `b'
	local V_: display %4.3f `V'
	return local b_formatted="`b_'"+"`stars'"
	return local V_formatted="("+"`V_'"+")"
end


program define getEstimates
	matrix define b=e(b)
	matrix define V=e(V)

	est restore base
	eststo est`2'

	formatData 1
	estadd local beta`1' `r(b_formatted)'
	estadd local betaSD`1' `r(V_formatted)'
	
	//esttab `estadded', stats($statsadded)
end

