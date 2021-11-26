do setupAndDefinePrograms.do

createIndResultsFileIfMissing

foreach i of numlist 1/5 {		//different estimations / columns in table
foreach wid of numlist  1(1)100 {			//different individuals

	use "Data/decisions_data_w_non_para", clear
	disp "`i'-----wid=`wid'"
	keep if wid==`wid'
	
	**Maybe Needed vars - these are used by the aggregate and sometimes not needed
	genVarsNeededForMl
	
	**Setup vars for ML - lots of these mean nothing for the individual estimation. They mean something for aggregate estimation. If there is an error, likely due to not having an aggregate global defined.
	//default values
	setupMlGlobals
	global pb=0  						//this is on for aggregate by default but not for individual
	global condition="(1==1)"			//this is used for date sample
	//setup vars used in default
	global vars="(betahat:) (delta:) (gamma:) (beta:) (phi:)"
	global setupvars="sigma betahat delta gamma beta phi"
	global init="ml init beta:_cons=1 betahat:_cons=1 gamma:_cons=2 delta:_cons=1 phi:_cons=250 sigma:_cons=50"

	//how do estimations change with each i?
	if (`i'==2) {
		global init="ml init beta:_cons=1 betahat:_cons=1 gamma:_cons=2 delta:_cons=1 phi:_cons=200 sigma:_cons=50"
		global condition="(decisiondatenum<4)"
	}
	if (`i'==3) {
		global init="ml init beta:_cons=1 betahat:_cons=1 gamma:_cons=2 delta:_cons=1 phi:_cons=200 sigma:_cons=50"  //this gets more people
		global condition="(decisiondatenum>=4)"
	}
	if (`i'==4) {
		global init="ml init beta:_cons=1 betahat:_cons=1 gamma:_cons=2 delta:_cons=1 phi:_cons=200 sigma:_cons=50"  //this gets more people
		global pb=1
		global vars="(betahat:) (delta:) (gamma:) (beta:) (phi:) (alpha:)"
		global setupvars="sigma betahat delta gamma beta phi alpha"
	}
	if (`i'==5) {
		global init="ml init beta:_cons=1 betahat:_cons=1 gamma:_cons=2 delta:_cons=1 phi:_cons=250 sigma:_cons=50"  //this gets more people
		global soph=1
		global sophfixed=1
		replace sophfixedweight=1
		global vars="(betahat:) (delta:) (gamma:) (beta:) (phi:)"
		global setupvars="sigma betahat delta gamma beta phi"
		//global init="ml init beta:_cons=.375 betahat:_cons=1.179 gamma:_cons=1.64 delta:_cons=.98 phi:_cons=807 sigma:_cons=130"
	}	
	
	
	cap noisily {
		
		ml model lf NBetaHatStructural1_agg (sigma: netdistance effort wage prediction guessbonusamount today pb group_for_test sophfixedweight distanceToPayment distanceToWork=) $vars if wid==`wid' & ($data_sample1) & $condition 
		$init
	
		local ml_iters=200
		if ($soph==1) local ml_iters=100 			//gets stuck on soph - takes soooo long and doesn't change
		ml max, iter(`ml_iters') 
		
		**Get and save needed vars
		getAndStoreIndEst `i' `wid'
	
	}
	disp "bh:`bhi' g:`gi' d:`di' p:`pi' b:`bi' a: `ai' -> iters=`iters'"

}
}

