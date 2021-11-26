do setupAndDefinePrograms.do

****************
******Main Results
****************

foreach s of numlist 2 1 4 {		// this changes the subject sample
foreach ds of numlist 1 2 3 4 {		// this changes the data sample
foreach i of numlist  1/5 101/103 1001/1005  2001/2004 4001/4005 5001/5002 6001/6004 7001/7004 8001/8004 9001/9005 10001/10004 11001/11011 12001/12004 13001/13004 14001/14006 15001/15004  16001/16005  17001/17005 18001/18005    {  //all types of estimations



	
	
	disp "****************************** s=`s' ds=`ds' i=`i' "
	
	//this skips combinations of data/subjects we don't need
	if ( (`ds'!=1 & `i'>10) | (`s'!=2 & `i'>10) ) 	continue    //only do the alternative participant samples and data samples for main specification
	if ((`ds'!=1 & `s'!=2))  						continue	//combinations of alternative p samples and data samples are fine, but not estimated
	
	use "Data/decisions_data_w_ind", clear
	
	genVarsNeededForMl

	setupMlGlobals
	global vars="(betahat:) (delta:) (phi:) (gamma:) (beta:) (alpha:)"
	global setupvars="sigma betahat delta phi gamma beta alpha"
	global wid_sample="${wid_sample`s'}"
	global data_sample="${data_sample`ds'}"
	global init="ml init beta:_cons=1 betahat:_cons=1 gamma:_cons=2 delta:_cons=1 phi:_cons=200 sigma:_cons=50"
	if (`s'==1) {
		global init="ml init beta:_cons=.8 betahat:_cons=1 gamma:_cons=2 delta:_cons=1 phi:_cons=500 sigma:_cons=40"
	}
	
	************Main Table	
	// baseline estimation
	if (`i'==1) {
		global widfe=0
	}
	
	// baseline with wid fixed effects
	if (`i'==2) {
		global vars="(betahat:) (delta:) (phi: wid_*) (gamma:) (beta:) (alpha:)"
	}
		
	// baseline with fixed effects for decisiondatenum
	if (`i'==3) { 
		global dayfe=1
		global vars="(betahat:) (delta:) (phi: wid_* ddn_*) (gamma: ddn_*) (beta:) (alpha:)"	
	}

	//sophistication forced to be 1
	if (`i'==4) {
		global soph=1	
		global sophfixed=1	
		global dayfe=1			
		global vars="(betahat:) (delta:) (phi: wid_* ddn_*) (gamma: ddn_*) (beta:) (alpha:)"
		//can start from anywhere, but takes a while. Hard w/ all FE
		global init="ml init beta:_cons=.86 betahat:_cons=1 gamma:_cons=2.07 delta:_cons=1.002 phi:_cons=480 sigma:_cons=43"
		/*
		if (`s'==4) {
			global init="ml init beta:_cons=.81 betahat:_cons=1 gamma:_cons=2.07 delta:_cons=1 phi:_cons=550 sigma:_cons=45"
		}
		if (`s'==1) {
			global init="ml init beta:_cons=.74 betahat:_cons=1 gamma:_cons=2.07 delta:_cons=1 phi:_cons=550 sigma:_cons=45"
		}
		*/
		if (`ds'==4) {
			global init="ml init beta:_cons=.74 betahat:_cons=1 gamma:_cons=2.07 delta:_cons=1 phi:_cons=550 sigma:_cons=45"
		}
		
	}
	//later date
	if (`i'==5) {
		global later=1
		global dayfe=1
		global vars="(betahat:) (delta:) (phi: wid_* ddn_*) (gamma: ddn_*) (beta:) (alpha:)" 
	}

		

	
		
		
	*************To get quasi-hyperbolic -> allow delta1, delta2, etc. on distance from decisiondatenum to workdatenum
	//Following the discussion of the appropriateness of the quasi-hyperbolic discounting model in the reduced-form analysis, we separately estimate the relative weight placed on the disutility from tasks that must be completed on future dates in comparison to tasks that must be completely immediately. We find that the weights for tasks that must be completed one, two, and three dates away are 17%, 20%, and 17% higher than the weight placed on immediate tasks. All of these estimates are significantly different from 0 (?²(1)=13.66, p<0.001; ?²(1)=13.69, p<0.001; ?²(1)=9.85, p= 0.002).
	if (`i'==101) {
		cap drop days_to_payment
		gen netdistanceNum=workdatenum-decisiondatenum
		global separatedelta=1
		global fixedbeta1=1
		global fixeddelta1=1
		global netdistance="netdistanceNum"
		global setupvars="sigma betahat delta1 delta2 delta3 delta4 phi gamma alpha"
		global vars="(betahat:) (delta1: ) (delta2: ) (delta3: ) (delta4: ) (phi:) (gamma:) (alpha:)"	
		global init="ml init betahat:_cons=1 gamma:_cons=2 delta1:_cons=1 delta2:_cons=1 delta3:_cons=1 delta4:_cons=1 phi:_cons=500 sigma:_cons=40"
		
	}	
	if (`i'==102) {
		cap drop days_to_payment
		gen netdistanceNum=workdatenum-decisiondatenum
		global separatedelta=1
		global fixedbeta1=1
		global fixeddelta1=1
		global netdistance="netdistanceNum"
		global setupvars="sigma betahat delta1 delta2 delta3 delta4 phi gamma alpha"
		global vars="(betahat:) (delta1: ) (delta2: ) (delta3: ) (delta4: ) (phi: wid_*) (gamma:) (alpha:)"	
		global init="ml init betahat:_cons=1 gamma:_cons=2 delta1:_cons=1 delta2:_cons=1 delta3:_cons=1 delta4:_cons=1 phi:_cons=500 sigma:_cons=40"
		
	}	
	if (`i'==103) {
		global dayfe=1
		cap drop days_to_payment
		gen netdistanceNum=workdatenum-decisiondatenum
		global separatedelta=1
		global fixedbeta1=1
		global fixeddelta1=1
		global netdistance="netdistanceNum"
		global setupvars="sigma betahat delta1 delta2 delta3 delta4 phi gamma alpha"
		global vars="(betahat:) (delta1: ) (delta2: ) (delta3: ) (delta4: ) (phi: wid_* ddn_*) (gamma: ddn_*) (alpha:)"	
		global init="ml init betahat:_cons=1 gamma:_cons=2 delta1:_cons=1 delta2:_cons=1 delta3:_cons=1 delta4:_cons=1 phi:_cons=500 sigma:_cons=40"
		
	}	
	
	
		
		
		
	************ Non-tobit
	// baseline estimation
	if (`i'==1001) {
		global tobit=0
		global widfe=0
	}
	
	// baseline with wid fixed effects
	if (`i'==1002) {
		global tobit=0
		global vars="(betahat:) (delta:) (phi: wid_*) (gamma:) (beta:) (alpha:)"
	}
		
	// baseline with fixed effects for decisiondatenum
	if (`i'==1003) { 
		global tobit=0
		global dayfe=1
		global vars="(betahat:) (delta:) (phi: wid_* ddn_*) (gamma: ddn_*) (beta:) (alpha:)"	
		//oddly required given constraint to make average fixed effect work.
		//not needed without constraint (converges to same place though)
		global init="ml init beta:_cons=.9 betahat:_cons=1.003 gamma:_cons=2.13 delta:_cons=1.002 phi:_cons=680 sigma:_cons=43 alpha:_cons=6.8"

	}

	//sophistication forced to be 1
	if (`i'==1004) {
		global tobit=0
		global soph=1	
		global sophfixed=1	
		global dayfe=1			
		global vars="(betahat:) (delta:) (phi: wid_* ddn_*) (gamma: ddn_*) (beta:) (alpha:)"
		//can start from anywhere, but takes a while. Hard w/ all FE
		global init="ml init beta:_cons=.86 betahat:_cons=1 gamma:_cons=2.07 delta:_cons=1.002 phi:_cons=480 sigma:_cons=43"	
		global init="ml init beta:_cons=.9 betahat:_cons=1 gamma:_cons=2.7 delta:_cons=1.002 phi:_cons=5000 sigma:_cons=30"	

		//if (`s'==4) {
		//	global init="ml init beta:_cons=.81 betahat:_cons=1 gamma:_cons=2.07 delta:_cons=1 phi:_cons=550 sigma:_cons=45"
		//}
		//if (`s'==1) {
		//		global init="ml init beta:_cons=1 betahat:_cons=1 gamma:_cons=2 delta:_cons=1 phi:_cons=500 sigma:_cons=50"
		//}
	}
	//later
	if (`i'==1005) {
		global tobit=0
		global later=1
		global dayfe=1
		global vars="(betahat:) (delta:) (phi: wid_* ddn_*) (gamma: ddn_*) (beta:) (alpha:)"
	}
	

	

	************Forcing Sophistication
	// baseline estimation
	if (`i'==2001) {
	
		global soph=1	
		global sophfixed=1
		global widfe=0
		global init="ml init beta:_cons=.83 betahat:_cons=1.003 gamma:_cons=2.13 delta:_cons=1.002 phi:_cons=680 sigma:_cons=43 alpha:_cons=6.8"

	}
	// baseline with wid fixed effects
	if (`i'==2002) {
		global soph=1	
		global sophfixed=1
		global vars="(betahat:) (delta:) (phi: wid_*) (gamma:) (beta:) (alpha:)"
		global init="ml init beta:_cons=.86 betahat:_cons=1 gamma:_cons=2.07 delta:_cons=1.002 phi:_cons=480 sigma:_cons=43 alpha:_cons=4"
	}
	// baseline with fixed effects for decisiondatenum?? This is the one done above. Drop this?
	if (`i'==2003) { 
		global soph=1	
		global sophfixed=1
		global dayfe=1
		global vars="(betahat:) (delta:) (phi: wid_* ddn_*) (gamma: ddn_*) (beta:) (alpha:)"
		global init="ml init beta:_cons=.86 betahat:_cons=1 gamma:_cons=2.07 delta:_cons=1.002 phi:_cons=480 sigma:_cons=26"
	
	}

	//later
	if (`i'==2004) {
		
		global soph=1	
		global sophfixed=1
		global later=1
		global dayfe=1
		global vars="(betahat:) (delta:) (phi: wid_* ddn_*) (gamma: ddn_*) (beta:) (alpha:)"	
		global init="ml init beta:_cons=.86 betahat:_cons=1 gamma:_cons=2.2 delta:_cons=1.002 phi:_cons=800 sigma:_cons=30"
	}	
	//later - no FE (NOT USED - > just testing
	if (`i'==2005) {
		global soph=1	
		global sophfixed=1
		global later=1
		global dayfe=1
		global vars="(betahat:) (delta:) (phi: wid_*) (gamma:  ddn_*) (beta:) (alpha:)"	
		global init="ml init beta:_cons=.86 betahat:_cons=1 gamma:_cons=2.2 delta:_cons=1.002 phi:_cons=750 sigma:_cons=45"
	}

	
	
	**************** Delta Adjusted
	
	// baseline estimation
	if (`i'==4001) {
		global netdistance="netdistance_adj"
		global widfe=0
	}
	// baseline with wid fixed effects
	if (`i'==4002) {
		global netdistance="netdistance_adj"
		global vars="(betahat:) (delta:) (phi: wid_*) (gamma:) (beta:) (alpha:)"
	}
	// baseline with fixed effects for decisiondatenum
	if (`i'==4003) { 
		global netdistance="netdistance_adj"
		global dayfe=1
		global vars="(betahat:) (delta:) (phi: wid_* ddn_*) (gamma: ddn_*) (beta:) (alpha:)"	
	}
	//sophistication forced to be 1
	if (`i'==4004) {
		global netdistance="netdistance_adj"
		global soph=1	
		global sophfixed=1	
		global dayfe=1			
		global vars="(betahat:) (delta:) (phi: wid_* ddn_*) (gamma: ddn_*) (beta:) (alpha:)"
		//can start from anywhere, but takes a while.
		global init="ml init beta:_cons=.8589 betahat:_cons=1 gamma:_cons=2.066 delta:_cons=1.002 phi:_cons=482.39 sigma:_cons=42.65"
	}
	//later
	if (`i'==4005) {
		global netdistance="netdistance_adj"
		global later=1
		global dayfe=1
		global vars="(betahat:) (delta:) (phi: wid_* ddn_*) (gamma: ddn_*) (beta:) (alpha:)"	
	}



	
	**************** Uncertainty
	if (`i'==5001) {
		global widfe=0
		global double_error=1
		global error_gamma=1
		global vars="(betahat:) (delta:) (phi:) (gamma:) (beta:) (sigma2:) (alpha:)"
		global setupvars="sigma betahat delta phi gamma beta sigma2 alpha"
		global init "ml init beta:_cons=.8437 betahat:_cons=1.000 gamma:_cons=2.351 delta:_cons=1.005 phi:_cons=1041 sigma:_cons=.2156 sigma2:_cons=.1367"
	}
	if (`i'==5002) {
		global widfe=0
		do uncertainty_setup.do
		global uncertainty=1
		global pb=0
		global error_gamma=1
		global double_error=1
		global vars="(betahat:) (delta:) (phi:) (gamma:) (beta:) (sigma2:)"
		global setupvars="sigma betahat delta phi gamma beta sigma2"
		global init="ml init beta:_cons=.8 betahat:_cons=1 gamma:_cons=2 delta:_cons=1 phi:_cons=1000 sigma:_cons=.2 sigma2:_cons=.1"
		global init="ml init beta:_cons=.767 betahat:_cons=.906 gamma:_cons=2.33 delta:_cons=1.0055 phi:_cons=1066.656 sigma:_cons=.217 sigma2:_cons=.1228"
	}	
	
		
	**************** Error on gamma
	
	//error gamma
	if (`i'==6001) {
		global widfe=0
		global error_gamma=1
		global vars="(betahat:) (delta:) (phi:) (gamma:)  (beta:) (alpha:)"	
		global init="ml init beta:_cons=1 betahat:_cons=1 gamma:_cons=2 delta:_cons=1 phi:_cons=1000 sigma:_cons=.22"
	}
	if (`i'==6002) {
		global error_gamma=1
		global vars="(betahat:) (delta:) (phi: wid_*) (gamma:)  (beta:) (alpha:)"	
		global init="ml init beta:_cons=1 betahat:_cons=1 gamma:_cons=2 delta:_cons=1 phi:_cons=1000 sigma:_cons=.22"
	}
	if (`i'==6003) {
		global dayfe=1
		global error_gamma=1
		global vars="(betahat:) (delta:) (phi: wid_* ddn_*) (gamma: ddn_*) (beta:) (alpha:)"
		global init="ml init beta:_cons=1 betahat:_cons=1 gamma:_cons=2 delta:_cons=1 phi:_cons=1000 sigma:_cons=.22"
	}
	if (`i'==6004) {
		global dayfe=1
		global later=1
		global error_gamma=1
		global vars="(betahat:) (delta:) (phi: wid_* ddn_*) (gamma: ddn_*)  (beta:) (alpha:)"	
		global init="ml init beta:_cons=1 betahat:_cons=1 gamma:_cons=2 delta:_cons=1 phi:_cons=1000 sigma:_cons=.22"
	}

	
	**************** Error on ln(phi)
	
	if (`i'==7001) {
		global error_phi=1
		global ln=1
		global widfe=0
		global setupvars="sigma betahat delta phi gamma beta alpha"
		global vars="(betahat:) (delta:) (phi:) (gamma:) (beta:) (alpha:)"	
		global init="ml init beta:_cons=1 betahat:_cons=1 gamma:_cons=2.2 delta:_cons=1 phi:_cons=750 sigma:_cons=.7"
	}
	if (`i'==7002) {
		global error_phi=1
		global ln=1
		global setupvars="sigma betahat delta phi gamma beta alpha"
		global vars="(betahat:) (delta:) (phi: wid_*) (gamma: ) (beta:) (alpha:)"	
		global init="ml init beta:_cons=1 betahat:_cons=1 gamma:_cons=2.2 delta:_cons=1 phi:_cons=750 sigma:_cons=.7"

	}
	if (`i'==7003) {
		global error_phi=1
		global ln=1
		global dayfe=1
		global setupvars="sigma betahat delta phi gamma beta alpha"
		global vars="(betahat:) (delta:) (phi: wid_* ddn_*) (gamma: ddn_*) (beta:) (alpha:)"	
		global init="ml init beta:_cons=1 betahat:_cons=1 gamma:_cons=2.2 delta:_cons=1 phi:_cons=750 sigma:_cons=.7"

	}
	if (`i'==7004) {
		global later=1
		global error_phi=1
		global ln=1
		global dayfe=1
		global setupvars="sigma betahat delta phi gamma beta alpha"
		global vars="(betahat:) (delta:) (phi: wid_* ddn_*) (gamma: ddn_*) (beta:) (alpha:)"	
		global init="ml init beta:_cons=1 betahat:_cons=1 gamma:_cons=2.2 delta:_cons=1 phi:_cons=750 sigma:_cons=.7"

	}
	
	
	**************** Error on log of effort
	

	if (`i'==8001) {
		global ln=1
		global vars="(betahat:) (delta:) (phi:) (gamma:)  (beta:) (alpha:)"	
		global init="ml init beta:_cons=1 betahat:_cons=1 gamma:_cons=2 delta:_cons=1 phi:_cons=500 sigma:_cons=.7"
	}
	if (`i'==8002) {
		global ln=1
		global vars="(betahat:) (delta:) (phi: wid_*) (gamma:)  (beta:) (alpha:)"	
		global init="ml init beta:_cons=1 betahat:_cons=1 gamma:_cons=2 delta:_cons=1 phi:_cons=500 sigma:_cons=.7"
	}
	//takes a while - gets stuck on FE
	if (`i'==8003) {
		global dayfe=1
		global ln=1
		global vars="(betahat:) (delta:) (phi: wid_* ddn_*) (gamma: ddn_*) (beta:) (alpha:)"
		global init="ml init beta:_cons=1 betahat:_cons=1 gamma:_cons=2 delta:_cons=1 phi:_cons=500 sigma:_cons=.7"
	}
	if (`i'==8004) {
		global dayfe=1
		global later=1
		global ln=1
		global vars="(betahat:) (delta:) (phi: wid_* ddn_*) (gamma: ddn_*)  (beta:) (alpha:)"	
		global init="ml init beta:_cons=1 betahat:_cons=1 gamma:_cons=2 delta:_cons=1 phi:_cons=500 sigma:_cons=.7"
	}

	
	**************** Projection Bias 1

	// baseline estimation
	if (`i'==9001) {
		global pb=0
		global structuralpb1=1
		global widfe=0
	}
	
	// baseline with wid fixed effects
	if (`i'==9002) {
		global pb=0
		global structuralpb1=1
		global vars="(betahat:) (delta:) (phi: wid_*) (gamma:) (beta:) (alpha:)"
	}
		
	// baseline with fixed effects for decisiondatenum
	if (`i'==9003) { 
		global pb=0
		global structuralpb1=1
		global dayfe=1
		global vars="(betahat:) (delta:) (phi: wid_* ddn_*) (gamma: ddn_*) (beta:) (alpha:)"	
	}

	//sophistication forced to be 1
	if (`i'==9004) {
		global pb=0
		global structuralpb1=1
		global soph=1	
		global sophfixed=1	
		global dayfe=1			
		global vars="(betahat:) (delta:) (phi: wid_* ddn_*) (gamma: ddn_*) (beta:) (alpha:)"
		global init="ml init beta:_cons=.86 betahat:_cons=1 gamma:_cons=2.07 delta:_cons=1.002 phi:_cons=480 sigma:_cons=43"	
	}
	
	
	//later
	if (`i'==9005) {
		global pb=0
		global structuralpb1=1
		global later=1
		global dayfe=1
		global vars="(betahat:) (delta:) (phi: wid_* ddn_*) (gamma: ddn_*) (beta:) (alpha:)"	
	}

	
	**************** Projection Bias 2 - issues with convergence?

	// baseline estimation
	if (`i'==10001) {
		global nonrtol=1
		global pb=0
		global structuralpb2=1
		global widfe=0
	}
	
	// baseline with wid fixed effects
	if (`i'==10002) {
		global nonrtol=1
		global pb=0
		global structuralpb2=1
		global vars="(betahat:) (delta:) (phi: wid_*) (gamma:) (beta:) (alpha:)"
	}
		
	// baseline with fixed effects for decisiondatenum
	if (`i'==10003) { 
		global nonrtol=1
		global pb=0
		global structuralpb2=1
		global dayfe=1
		global vars="(betahat:) (delta:) (phi: wid_* ddn_*) (gamma: ddn_*) (beta:) (alpha:)"	
	}
	if (`i'==10004) {
		global pb=0
		global structuralpb2=2
		global later=1
		global dayfe=1
		global vars="(betahat:) (delta:) (phi: wid_* ddn_*) (gamma: ddn_*) (beta:) (alpha:)"	
	}

	
	
	
	*******Different forms of FE
	// this is already done:
	if (`i'==11000) { 
		global widfe=1
		global dayfe=0
		global vars="(betahat:) (delta:) (phi: wid_*) (gamma:) (beta:) (alpha:)"	
	}
	if (`i'==11001) { 
		global widfe=1
		global dayfe=0
		global vars="(betahat:) (delta:) (phi:) (gamma: wid_*) (beta:) (alpha:)"	
	}
	if (`i'==11002) { 
		global widfe=1
		global dayfe=1
		global vars="(betahat:) (delta:) (phi: wid_* ddn_*) (gamma:) (beta:) (alpha:)"	
	}
	if (`i'==11003) { 
		global widfe=1
		global dayfe=1
		global vars="(betahat:) (delta:) (phi: wid_*) (gamma: ddn_*) (beta:) (alpha:)"	
	}
	if (`i'==11004) { 
		global widfe=1
		global dayfe=1
		global vars="(betahat:) (delta:) (phi: ddn_*) (gamma: wid_*) (beta:) (alpha:)"	
	}	
	if (`i'==11005) { 
		global widfe=1
		global dayfe=1
		global vars="(betahat:) (delta:) (phi:) (gamma: wid_* ddn_*) (beta:) (alpha:)"	
	}
	if (`i'==11006) { 
		global widfe=1
		global dayfe=1
		global vars="(betahat:) (delta:) (phi: ddn_*) (gamma: ddn_* wid_* ) (beta:) (alpha:)"	
	}
	if (`i'==11007) { 
		global widfe=1
		global dayfe=1
		global vars="(betahat:) (delta:) (phi: ddn_*) (gamma: ddn_*) (beta:) (alpha:)"	
	}
	

	*******Second Degree Polynomial
	if (`i'==12001) {
		global seconddegpoly=1
		global widfe=0
		global setupvars="sigma betahat delta psi1 psi2 beta alpha"
		global vars="(betahat:) (delta:) (psi1:) (psi2:)  (beta:) (alpha:)"
		global init="ml init beta:_cons=.8 betahat:_cons=1 delta:_cons=1 psi1:_cons=0 psi2:_cons=.0028 sigma:_cons=26"
	}
	if (`i'==12002) {
		global seconddegpoly=1
		global setupvars="sigma betahat delta psi1 psi2 beta alpha"
		global vars="(betahat:) (delta:) (psi1:) (psi2: wid_*)  (beta:) (alpha:)"
		global init="ml init beta:_cons=.8 betahat:_cons=1 delta:_cons=1 psi1:_cons=0 psi2:_cons=.0028 sigma:_cons=26"
	}
	if (`i'==12003) { 
		global seconddegpoly=1
		global dayfe=1
		global setupvars="sigma betahat delta psi1 psi2 beta alpha"
		global vars="(betahat:) (delta:) (psi1: ddn_*) (psi2: wid_* ddn_*)  (beta:) (alpha:)"
		global init="ml init beta:_cons=.8 betahat:_cons=1 delta:_cons=1 psi1:_cons=0 psi2:_cons=.0028 sigma:_cons=26"
	}
	if (`i'==12004) {
		global seconddegpoly=1
		global later=1
		global dayfe=1
		global setupvars="sigma betahat delta psi1 psi2 beta alpha"
		global vars="(betahat:) (delta:) (psi1: ddn_*) (psi2: wid_* ddn_*)  (beta:) (alpha:)"
		global init="ml init beta:_cons=.8 betahat:_cons=1 delta:_cons=1 psi1:_cons=0 psi2:_cons=.0028 sigma:_cons=26"
	}

		


	*******Third Degree Polynomial
	if (`i'==13001) {
		drop if wid==2 
		global thirddegpoly=1
		global widfe=0
		global setupvars="sigma betahat delta psi1 psi2 psi3 beta alpha"
		global vars="(betahat:) (delta:) (psi1:) (psi2:)  (psi3:) (beta:) (alpha:)"
		global init="ml init beta:_cons=.8 betahat:_cons=1 delta:_cons=1 psi1:_cons=0 psi2:_cons=.004 psi3:_cons=.0000004 sigma:_cons=25"
	}
	if (`i'==13002) {
		//this is a pain to get to converge given all of the params and FE
		drop if wid==2 
		global thirddegpoly=1
		global setupvars="sigma betahat delta psi1 psi2 psi3 beta alpha"
		global vars="(betahat:) (delta:) (psi1: ) (psi2: wid_*)  (psi3: ) (beta:) (alpha:)"
		global init="ml init beta:_cons=.8 betahat:_cons=1 delta:_cons=1 psi1:_cons=0 psi2:_cons=.0045 psi3:_cons=.00000025 sigma:_cons=25"
	}
	if (`i'==13003) { 
		drop if wid==2 
		global thirddegpoly=1
		global dayfe=1
		global setupvars="sigma betahat delta psi1 psi2 psi3 beta alpha"
		global vars="(betahat:) (delta:) (psi1: ddn_*) (psi2:  wid_* ddn_*)  (psi3:) (beta:) (alpha:)"
		global init="ml init beta:_cons=.8 betahat:_cons=1 delta:_cons=1 psi1:_cons=0 psi2:_cons=.004 psi3:_cons=.0000004 sigma:_cons=25"
	}
	if (`i'==13004) {
		drop if wid==2 
		global thirddegpoly=1
		global later=1
		global dayfe=1
		global setupvars="sigma betahat delta psi1 psi2 psi3 beta alpha"
		global vars="(betahat:) (delta:) (psi1: ddn_*) (psi2:  wid_* ddn_*)  (psi3:) (beta:) (alpha:)"
		global init="ml init beta:_cons=.8 betahat:_cons=1 delta:_cons=1 psi1:_cons=0 psi2:_cons=.004 psi3:_cons=.0000004 sigma:_cons=25"
	}

	
	*******
	// Baseline: not used - just use with FE
	local constraint_delta_m_e=""
	if (`i'==14001) {
		global widfe=0
		local constraint_delta_m_e="[delta_e]_cons=[delta_m]_cons"
		global monetarydelta=1
		global includeDeltaM1TestInTable=1
		global includeDeltaMDeltaETestInTable=0
		global setupvars="sigma betahat delta_m delta_e phi gamma beta alpha"
		global vars="(betahat:) (delta_e:) (delta_m:) (phi:) (gamma:) (beta:) (alpha:)"
		global init="ml init beta:_cons=1 betahat:_cons=1 gamma:_cons=2 delta_m:_cons=1 delta_e:_cons=1 phi:_cons=200 sigma:_cons=50"
	}
	if (`i'==14002) {
		global widfe=0
		global monetarydelta=1
		global includeDeltaM1TestInTable=1
		global includeDeltaMDeltaETestInTable=1
		global setupvars="sigma betahat delta_m delta_e phi gamma beta alpha"
		global vars="(betahat:) (delta_e:) (delta_m:) (phi:) (gamma:) (beta:) (alpha:)"
		global init="ml init beta:_cons=1 betahat:_cons=1 gamma:_cons=2 delta_m:_cons=1 delta_e:_cons=1 phi:_cons=200 sigma:_cons=50"
	}
	if (`i'==14003) {
		global widfe=0
		local constraint_delta_m_e="[delta_m]_cons=1"
		global monetarydelta=1
		global includeDeltaM1TestInTable=0
		global includeDeltaMDeltaETestInTable=0
		global setupvars="sigma betahat delta_m delta_e phi gamma beta alpha"
		global vars="(betahat:) (delta_e:) (delta_m:) (phi:) (gamma:) (beta:) (alpha:)"
		global init="ml init beta:_cons=1 betahat:_cons=1 gamma:_cons=2 delta_m:_cons=1 delta_e:_cons=1 phi:_cons=200 sigma:_cons=50"
	}
	if (`i'==14004) {
		global dayfe=1
		local constraint_delta_m_e="[delta_e]_cons=[delta_m]_cons"
		global monetarydelta=1
		global includeDeltaM1TestInTable=1
		global includeDeltaMDeltaETestInTable=0
		global setupvars="sigma betahat delta_m delta_e phi gamma beta alpha"
		global vars="(betahat:) (delta_e:) (delta_m:) (phi: wid_* ddn_*) (gamma: ddn_*) (beta:) (alpha:)"
		global init="ml init beta:_cons=1 betahat:_cons=1 gamma:_cons=2 delta_m:_cons=1 delta_e:_cons=1 phi:_cons=200 sigma:_cons=50"
	}
	if (`i'==14005) {
		global dayfe=1
		global monetarydelta=1
		global includeDeltaM1TestInTable=1
		global includeDeltaMDeltaETestInTable=1
		global setupvars="sigma betahat delta_m delta_e phi gamma beta alpha"
		global vars="(betahat:) (delta_e:) (delta_m:) (phi: wid_* ddn_*) (gamma: ddn_*) (beta:) (alpha:)"
		global init="ml init beta:_cons=1 betahat:_cons=1 gamma:_cons=2 delta_m:_cons=1 delta_e:_cons=1 phi:_cons=200 sigma:_cons=50"
	}
	if (`i'==14006) {
		global dayfe=1
		local constraint_delta_m_e="[delta_m]_cons=1"
		global monetarydelta=1
		global includeDeltaM1TestInTable=0
		global includeDeltaMDeltaETestInTable=0
		global setupvars="sigma betahat delta_m delta_e phi gamma beta alpha"
		global vars="(betahat:) (delta_e:) (delta_m:) (phi: wid_* ddn_*) (gamma: ddn_*) (beta:) (alpha:)"
		global init="ml init beta:_cons=1 betahat:_cons=1 gamma:_cons=2 delta_m:_cons=1 delta_e:_cons=1 phi:_cons=200 sigma:_cons=50"
	}
	
	
	*******Risk aversion
	if (`i'==15001) {
		global riskaversionalpha=0.01
		global widfe=0
		global vars="(betahat:) (delta:) (phi:) (gamma:) (beta:) (alpha:)"
	}
	if (`i'==15002) {
		//this is a pain to get to converge given all of the params and FE
		global riskaversionalpha=0.01
		global vars="(betahat:) (delta:) (phi: wid_*) (gamma:) (beta:) (alpha:)"
	}
	if (`i'==15003) { 
		global riskaversionalpha=0.01
		global dayfe=1
		global vars="(betahat:) (delta:) (phi: wid_* ddn_*) (gamma: ddn_*) (beta:) (alpha:)"	
	}
	if (`i'==15004) {
		global riskaversionalpha=0.01
		global later=1
		global dayfe=1
		global vars="(betahat:) (delta:) (phi: wid_* ddn_*) (gamma: ddn_*) (beta:) (alpha:)"	
	}

	
	
	*******Projection Bias leading to bias if ignored. Note: dont estimate alpha then
	if (`i'>=16000 & `i'<16010) {
		drop if type=="future" & workdone1==10
		drop if type=="prediction" & workdone1==10
		drop if type=="present" & workdone1==0
		global pb=0
		global setupvars="sigma betahat delta phi gamma beta"
	}
	// baseline estimation
	if (`i'==16001) {
		global widfe=0
		global vars="(betahat:) (delta:) (phi:) (gamma:) (beta:)"
	}
	
	// baseline with wid fixed effects
	if (`i'==16002) {
		global vars="(betahat:) (delta:) (phi: wid_*) (gamma:) (beta:)"
	}
		
	// baseline with fixed effects for decisiondatenum
	if (`i'==16003) { 
		global dayfe=1
		global vars="(betahat:) (delta:) (phi: wid_* ddn_*) (gamma: ddn_*) (beta:)"	
	}

	//sophistication forced to be 1
	if (`i'==16004) {
		global soph=1	
		global sophfixed=1	
		global dayfe=1			
		global vars="(betahat:) (delta:) (phi: wid_* ddn_*) (gamma: ddn_*) (beta:)"
		//can start from anywhere, but takes a while. Hard w/ all FE
		global init="ml init beta:_cons=.86 betahat:_cons=1 gamma:_cons=2.07 delta:_cons=1.002 phi:_cons=480 sigma:_cons=43"
		/*
		if (`s'==4) {
			global init="ml init beta:_cons=.81 betahat:_cons=1 gamma:_cons=2.07 delta:_cons=1 phi:_cons=550 sigma:_cons=45"
		}
		if (`s'==1) {
			global init="ml init beta:_cons=.74 betahat:_cons=1 gamma:_cons=2.07 delta:_cons=1 phi:_cons=550 sigma:_cons=45"
		}
		*/
		if (`ds'==4) {
			global init="ml init beta:_cons=.74 betahat:_cons=1 gamma:_cons=2.07 delta:_cons=1 phi:_cons=550 sigma:_cons=45"
		}
		
	}
	//later date
	if (`i'==16005) {
		global later=1
		global dayfe=1
		global vars="(betahat:) (delta:) (phi: wid_* ddn_*) (gamma: ddn_*) (beta:)" 
	}

	
	
	*******With no projection bias parameter
	if (`i'>=17000 & `i'<17010) {
		global pb=0
		global setupvars="sigma betahat delta phi gamma beta"
	}
	// baseline estimation
	if (`i'==17001) {
		global widfe=0
		global vars="(betahat:) (delta:) (phi:) (gamma:) (beta:)"
	}
	
	// baseline with wid fixed effects
	if (`i'==17002) {
		global vars="(betahat:) (delta:) (phi: wid_*) (gamma:) (beta:)"
	}
		
	// baseline with fixed effects for decisiondatenum
	if (`i'==17003) { 
		global dayfe=1
		global vars="(betahat:) (delta:) (phi: wid_* ddn_*) (gamma: ddn_*) (beta:)"	
	}

	//sophistication forced to be 1
	if (`i'==17004) {
		global soph=1	
		global sophfixed=1	
		global dayfe=1			
		global vars="(betahat:) (delta:) (phi: wid_* ddn_*) (gamma: ddn_*) (beta:)"
		global init="ml init beta:_cons=.86 betahat:_cons=1 gamma:_cons=2.07 delta:_cons=1.002 phi:_cons=480 sigma:_cons=43"
		if (`ds'==4) {
			global init="ml init beta:_cons=.74 betahat:_cons=1 gamma:_cons=2.07 delta:_cons=1 phi:_cons=550 sigma:_cons=45"
		}
		
	}
	//later date
	if (`i'==17005) {
		global later=1
		global dayfe=1
		global vars="(betahat:) (delta:) (phi: wid_* ddn_*) (gamma: ddn_*) (beta:)" 
	}

		
	********work date FE instead of decision date FE
	if (`i'==18001) {
		global widfe=0
	}
	
	// baseline with wid fixed effects
	if (`i'==18002) {
		global vars="(betahat:) (delta:) (phi: wid_*) (gamma:) (beta:) (alpha:)"
	}
		
	// baseline with fixed effects for decisiondatenum
	if (`i'==18003) { 
		global dayfe=1
		global vars="(betahat:) (delta:) (phi: wid_* wdn_*) (gamma: wdn_*) (beta:) (alpha:)"	
	}

	//sophistication forced to be 1
	if (`i'==18004) {
		global soph=1	
		global sophfixed=1	
		global dayfe=1			
		global vars="(betahat:) (delta:) (phi: wid_* wdn_*) (gamma: wdn_*) (beta:) (alpha:)"
		//can start from anywhere, but takes a while. Hard w/ all FE
		global init="ml init beta:_cons=.86 betahat:_cons=1 gamma:_cons=2.07 delta:_cons=1.002 phi:_cons=480 sigma:_cons=43"
		/*
		if (`s'==4) {
			global init="ml init beta:_cons=.81 betahat:_cons=1 gamma:_cons=2.07 delta:_cons=1 phi:_cons=550 sigma:_cons=45"
		}
		if (`s'==1) {
			global init="ml init beta:_cons=.74 betahat:_cons=1 gamma:_cons=2.07 delta:_cons=1 phi:_cons=550 sigma:_cons=45"
		}
		*/
		if (`ds'==4) {
			global init="ml init beta:_cons=.74 betahat:_cons=1 gamma:_cons=2.07 delta:_cons=1 phi:_cons=550 sigma:_cons=45"
		}
		
	}
	
	//later date
	if (`i'==18005) {
		global later=1
		global dayfe=1
		global vars="(betahat:) (delta:) (phi: wid_* wdn_*) (gamma: wdn_*) (beta:) (alpha:)" 
	}

	
	
	
	****Get sample we are using
	if ($later) {
		global data_sample="$data_sample & decisiondatenum>2"
	}
	keep if $wid_sample & $data_sample
	
	
	
	****Generate Fixed Effects and constraints so that constant=avg of FE
	generateFE `s'
	generateConstraintsSoConsIsAvg
	constraint define 13 `constraint_delta_m_e' 	//this is the monetary delta constraint

	***** Run model
	ml model lf NBetaHatStructural1_agg (sigma: $netdistance effort wage prediction guessbonusamount today pb group_for_test sophfixedweight distanceToPayment distanceToWork=) $vars ///
	if $wid_sample & $data_sample, clu(wid) technique(nr) constraints(1 2 3 4 5 6 7 8 9 10 11 12 13) collinear
	$init

	if ($soph & (`s'!=2)) | ($nonrtol==1) {
		ml max, nonrtolerance  iter(200)
	} 	
	else {
		ml max, iter(200)
	}

	eststo est`i'_`s'_`ds'
	
	storeEstVarsForTables `i'

}
}
}






*********** These are tables of main results for a variety of samples
** 1 = full sample (in appendix)
** 2 = main sample (in paper)
** 3 = main sample without attrition (not used)
** 4 = without attrition (in appendix)

foreach s of numlist 2 1 4 {

esttab est1_`s'_1 est2_`s'_1 est3_`s'_1 est5_`s'_1 est4_`s'_1 /// 
using "c:/graphics/workovertime2Final/aggresults_`s'.tex", ///
nolines ///
stats(widfe dayfe soph later N N_clust ll betatest betahattest alphatest deltatest, ///
labels("\\[-6pt] Participant FE" "\\[-6pt] Day FE" "\\[-6pt] Prediction Soph." "\\[-6pt] Later Decisions" "\midrule Observations" "Participants" "Log Likelihood" "\midrule \$H_{0}\$(\$\hat{\beta}\$ =1)" "\$H_{0}\$(\$\widehat{\beta_{h}}\$ =1)" "\$H_{0}\$(\$\widehat{\alpha}\$ =0)" "\$H_{0}\$(\$\widehat{\delta}\$ =1)") fmt(0)) ///
order(beta:_cons betahat:_cons delta:_cons gamma:_cons phi:_cons alpha:_cons) ///
cells( b(fmt(%9.3f %9.3f %9.3f %9.3f %9.0f %9.3f)) se(par fmt(%9.3f %9.3f %9.3f %9.3f %9.0f %9.3f ))) ///
width(0.8\hsize) drop(sigma:_cons wid_* ddn*)  replace booktabs gaps f collabels(none) eqlabels(none) ///
varlabels(beta:_cons "Present Bias \$\beta\$" ///
betahat:_cons "Naive Pres. Bias \$\beta_{h}\$"  ///
delta:_cons "Discount Factor \$\delta\$"  ///
gamma:_cons "Cost Curvature \$\gamma\$"  ///
phi:_cons "Cost Slope \$\varphi\$"  ///
alpha:_cons "Proj Task Reduction \$\tilde{\alpha}\$") ///
style(tex) legend label compress substitute(\_ _) ///
mtitle("\specialcell{Initial \\ Estimation}"  ///
"\specialcell{Participant \\ FE}"  ///
"\specialcell{Decision \\ Day FE}"  ///
"\specialcell{Later \\ Decisions}"  ///
"\specialcell{Pred. \\ Soph.}"  ///
)

}



*********** Other Samples 
* Later Choices
* No Bonus in decision set
* Only first in decision set
* Only first decision set


foreach ds of numlist 2 3 4 {

esttab est1_2_`ds' est2_2_`ds' est3_2_`ds' est5_2_`ds' est4_2_`ds' /// 
using "c:/graphics/workovertime2Final/datasample_`ds'.tex", ///
nolines ///
stats(widfe dayfe soph later N N_clust ll betatest betahattest alphatest deltatest, ///
labels("\\[-6pt] Participant FE" "\\[-6pt] Day FE" "\\[-6pt] Prediction Soph." "\\[-6pt] Later Decisions" "\midrule Observations" "Participants" "Log Likelihood" "\midrule \$H_{0}\$(\$\hat{\beta}\$ =1)" "\$H_{0}\$(\$\widehat{\beta_{h}}\$ =1)" "\$H_{0}\$(\$\widehat{\alpha}\$ =0)" "\$H_{0}\$(\$\widehat{\delta}\$ =1)") fmt(0)) ///
order(beta:_cons betahat:_cons delta:_cons gamma:_cons phi:_cons alpha:_cons) ///
cells( b(fmt(%9.3f %9.3f %9.3f %9.3f %9.0f %9.3f)) se(par fmt(%9.3f %9.3f %9.3f %9.3f %9.0f %9.3f ))) ///
width(0.8\hsize) drop(sigma:_cons wid_* ddn*)  replace booktabs gaps f collabels(none) eqlabels(none) ///
varlabels(beta:_cons "Present Bias \$\beta\$" ///
betahat:_cons "Naive Pres. Bias \$\beta_{h}\$"  ///
delta:_cons "Discount Factor \$\delta\$"  ///
gamma:_cons "Cost Curvature \$\gamma\$"  ///
phi:_cons "Cost Slope \$\varphi\$"  ///
alpha:_cons "Proj Task Reduction \$\tilde{\alpha}\$") ///
style(tex) legend label compress substitute(\_ _) ///
mtitle("\specialcell{Initial \\ Estimation}"  ///
"\specialcell{Participant \\ FE}"  ///
"\specialcell{Decision \\ Day FE}"  ///
"\specialcell{Later \\ Decisions}"  ///
"\specialcell{Pred. \\ Soph.}"  ///
)

}





***SAME AS First BUT NON TOBIT
esttab est1001_2_1 est1002_2_1 est1003_2_1 est1005_2_1 est1004_2_1 /// 
using "c:/graphics/workovertime2Final/non_tobit_2.tex", ///
nolines ///
stats(widfe dayfe soph later N N_clust ll betatest betahattest alphatest deltatest, ///
labels("\\[-6pt] Participant FE" "\\[-6pt] Day FE" "\\[-6pt] Prediction Soph." "\\[-6pt] Later Decisions" "\midrule Observations" "Participants" "Log Likelihood" "\midrule \$H_{0}\$(\$\hat{\beta}\$ =1)" "\$H_{0}\$(\$\widehat{\beta_{h}}\$ =1)" "\$H_{0}\$(\$\widehat{\alpha}\$ =0)" "\$H_{0}\$(\$\widehat{\delta}\$ =1)") fmt(0)) ///
order(beta:_cons betahat:_cons delta:_cons gamma:_cons phi:_cons alpha:_cons) ///
cells( b(fmt(%9.3f %9.3f %9.3f %9.3f %9.0f %9.3f)) se(par fmt(%9.3f %9.3f %9.3f %9.3f %9.0f %9.3f ))) ///
width(0.8\hsize) drop(sigma:_cons wid_* ddn*)  replace booktabs gaps f collabels(none) eqlabels(none) ///
varlabels(beta:_cons "Present Bias \$\beta\$" ///
betahat:_cons "Naive Pres. Bias \$\beta_{h}\$"  ///
delta:_cons "Discount Factor \$\delta\$"  ///
gamma:_cons "Cost Curvature \$\gamma\$"  ///
phi:_cons "Cost Slope \$\varphi\$"  ///
alpha:_cons "Proj Task Reduction \$\tilde{\alpha}\$") ///
style(tex) legend label compress substitute(\_ _) ///
mtitle("\specialcell{Initial \\ Estimation}"  ///
"\specialcell{Participant \\ FE}"  ///
"\specialcell{Decision \\ Day FE}"  ///
"\specialcell{Later \\ Decisions}"  ///
"\specialcell{Pred. \\ Soph.}"  ///
)




***SAME AS First BUT INCLUDES SOPHISTICATION
esttab est2001_2_1 est2002_2_1 est2003_2_1 est2004_2_1 /// 
using "c:/graphics/workovertime2Final/sophistication_2.tex", ///
nolines ///
stats(widfe dayfe soph later N N_clust ll betatest betahattest alphatest deltatest, ///
labels("\\[-6pt] Participant FE" "\\[-6pt] Day FE" "\\[-6pt] Prediction Soph." "\\[-6pt] Later Decisions" "\midrule Observations" "Participants" "Log Likelihood" "\midrule \$H_{0}\$(\$\hat{\beta}\$ =1)" "\$H_{0}\$(\$\widehat{\beta_{h}}\$ =1)" "\$H_{0}\$(\$\widehat{\alpha}\$ =0)" "\$H_{0}\$(\$\widehat{\delta}\$ =1)") fmt(0)) ///
order(beta:_cons betahat:_cons delta:_cons gamma:_cons phi:_cons alpha:_cons) ///
cells( b(fmt(%9.3f %9.3f %9.3f %9.3f %9.0f %9.3f)) se(par fmt(%9.3f %9.3f %9.3f %9.3f %9.0f %9.3f ))) ///
width(0.8\hsize) drop(sigma:_cons wid_* ddn*)  replace booktabs gaps f collabels(none) eqlabels(none) ///
varlabels(beta:_cons "Present Bias \$\beta\$" ///
betahat:_cons "Naive Pres. Bias \$\beta_{h}\$"  ///
delta:_cons "Discount Factor \$\delta\$"  ///
gamma:_cons "Cost Curvature \$\gamma\$"  ///
phi:_cons "Cost Slope \$\varphi\$"  ///
alpha:_cons "Proj Task Reduction \$\tilde{\alpha}\$") ///
style(tex) legend label compress substitute(\_ _) ///
mtitle("\specialcell{Initial \\ Estimation}"  ///
"\specialcell{Participant \\ FE}"  ///
"\specialcell{Decision \\ Day FE}"  ///
"\specialcell{Later \\ Decisions}"  ///
)



*****************Alt Discount
esttab est4001_2_1 est4002_2_1 est4003_2_1 est4005_2_1 est4004_2_1 /// 
using "c:/graphics/workovertime2Final/aggaltdelta_2.tex", ///
nolines ///
stats(widfe dayfe soph later N N_clust ll betatest betahattest alphatest deltatest, ///
labels("\\[-6pt] Participant FE" "\\[-6pt] Day FE" "\\[-6pt] Prediction Soph." "\\[-6pt] Later Decisions" "\midrule Observations" "Participants" "Log Likelihood" "\midrule \$H_{0}\$(\$\hat{\beta}\$ =1)" "\$H_{0}\$(\$\widehat{\beta_{h}}\$ =1)" "\$H_{0}\$(\$\widehat{\alpha}\$ =0)" "\$H_{0}\$(\$\widehat{\delta}\$ =1)") fmt(0)) ///
order(beta:_cons betahat:_cons delta:_cons gamma:_cons phi:_cons alpha:_cons) ///
cells( b(fmt(%9.3f %9.3f %9.3f %9.3f %9.0f %9.3f)) se(par fmt(%9.3f %9.3f %9.3f %9.3f %9.0f %9.3f ))) ///
width(0.8\hsize) drop(sigma:_cons wid_* ddn*)  replace booktabs gaps f collabels(none) eqlabels(none) ///
varlabels(beta:_cons "Present Bias \$\beta\$" ///
betahat:_cons "Naive Pres. Bias \$\beta_{h}\$"  ///
delta:_cons "Discount Factor \$\delta\$"  ///
gamma:_cons "Cost Curvature \$\gamma\$"  ///
phi:_cons "Cost Slope \$\varphi\$"  ///
alpha:_cons "Proj Bias \$\alpha\$") ///
style(tex) legend label compress substitute(\_ _) ///
mtitle("\specialcell{Initial \\ Estimation}"  ///
"\specialcell{Participant \\ FE}"  ///
"\specialcell{Decision \\ Day FE}"  ///
"\specialcell{Later \\ Decisions}"  ///
"\specialcell{Pred. \\ Soph.}"  ///
)




*********** This is table for uncertainty
esttab est5001_2_1 est5002_2_1 ///  
using "c:/graphics/workovertime2Final/agguncertainty_2.tex", ///
nolines ///
stats(N N_clust ll betatest deltatest, ///
labels("\midrule Observations" "Participants" "Log Likelihood" "\midrule \$H_{0}\$(\$\hat{\beta}\$ =1)" "\$H_{0}\$(\$\widehat{\delta}\$ =1)") fmt(0)) ///
order(beta:_cons delta:_cons gamma:_cons sigma:_cons sigma2:_cons) ///
cells( b(fmt(%9.3f %9.3f %9.3f %9.3f %9.3f %9.3f)) se(par fmt(%9.3f %9.3f %9.3f %9.3f %9.3f %9.3f ))) ///
width(0.8\hsize)  replace booktabs gaps f collabels(none) eqlabels(none) ///
drop(betahat:_cons phi:_cons alpha:_cons delta:_cons) ///
varlabels(beta:_cons "Present Bias \$\hat{\beta}\$" ///
gamma:_cons "Cost Curvature \$\hat{\gamma}\$"  ///
sigma:_cons "Decision Error \$\widehat{\sigma }(\varepsilon _{\gamma })\$" ///
sigma2:_cons "Preference Shock \$\widehat{\sigma }(\eta _{\gamma })\$" ///
) ///
style(tex) legend label compress substitute(\_ _) ///
mtitle("\specialcell{Shock on \$\gamma\$ \\ No Reaction }"  ///
"\specialcell{Shock on \$\gamma\$ \\ Optimal Reaction }"  ///
)



*********** Error on gamma
 
esttab est6001_2_1 est6002_2_1 est6003_2_1 est6004_2_1  /// 
using "c:/graphics/workovertime2Final/errorongamma_2.tex", ///
nolines ///
stats(widfe dayfe soph later N N_clust ll betatest betahattest alphatest deltatest, ///
labels("\\[-6pt] Participant FE" "\\[-6pt] Day FE" "\\[-6pt] Prediction Soph." "\\[-6pt] Later Decisions" "\midrule Observations" "Participants" "Log Likelihood" "\midrule \$H_{0}\$(\$\hat{\beta}\$ =1)" "\$H_{0}\$(\$\widehat{\beta_{h}}\$ =1)" "\$H_{0}\$(\$\widehat{\alpha}\$ =0)" "\$H_{0}\$(\$\widehat{\delta}\$ =1)") fmt(0)) ///
order(beta:_cons betahat:_cons delta:_cons gamma:_cons phi:_cons alpha:_cons) ///
cells( b(fmt(%9.3f %9.3f %9.3f %9.3f %9.0f %9.3f)) se(par fmt(%9.3f %9.3f %9.3f %9.3f %9.0f %9.3f ))) ///
width(0.8\hsize) drop(sigma:_cons wid_* ddn*)  replace booktabs gaps f collabels(none) eqlabels(none) ///
varlabels(beta:_cons "Present Bias \$\beta\$" ///
betahat:_cons "Naive Pres. Bias \$\beta_{h}\$"  ///
delta:_cons "Discount Factor \$\delta\$"  ///
gamma:_cons "Cost Curvature \$\gamma\$"  ///
phi:_cons "Cost Slope \$\varphi\$"  ///
alpha:_cons "Proj Task Reduction \$\tilde{\alpha}\$") ///
style(tex) legend label compress substitute(\_ _) ///
mtitle("\specialcell{Initial \\ Estimation}"  ///
"\specialcell{Participant \\ FE}"  ///
"\specialcell{Decision \\ Day FE}"  ///
"\specialcell{Later \\ Decisions}"  ///
)


*********** Error on ln(phi)
 
esttab est7001_2_1 est7002_2_1 est7003_2_1 est7004_2_1 /// 
using "c:/graphics/workovertime2Final/erroronphi_2.tex", ///
nolines ///
stats(widfe dayfe soph later N N_clust ll betatest betahattest alphatest deltatest, ///
labels("\\[-6pt] Participant FE" "\\[-6pt] Day FE" "\\[-6pt] Prediction Soph." "\\[-6pt] Later Decisions" "\midrule Observations" "Participants" "Log Likelihood" "\midrule \$H_{0}\$(\$\hat{\beta}\$ =1)" "\$H_{0}\$(\$\widehat{\beta_{h}}\$ =1)" "\$H_{0}\$(\$\widehat{\alpha}\$ =0)" "\$H_{0}\$(\$\widehat{\delta}\$ =1)") fmt(0)) ///
order(beta:_cons betahat:_cons delta:_cons gamma:_cons phi:_cons alpha:_cons) ///
cells( b(fmt(%9.3f %9.3f %9.3f %9.3f %9.0f %9.3f)) se(par fmt(%9.3f %9.3f %9.3f %9.3f %9.0f %9.3f ))) ///
width(0.8\hsize) drop(sigma:_cons wid_* ddn*)  replace booktabs gaps f collabels(none) eqlabels(none) ///
varlabels(beta:_cons "Present Bias \$\beta\$" ///
betahat:_cons "Naive Pres. Bias \$\beta_{h}\$"  ///
delta:_cons "Discount Factor \$\delta\$"  ///
gamma:_cons "Cost Curvature \$\gamma\$"  ///
phi:_cons "Cost Slope \$\varphi\$"  ///
alpha:_cons "Proj Task Reduction \$\tilde{\alpha}\$") ///
style(tex) legend label compress substitute(\_ _) ///
mtitle("\specialcell{Initial \\ Estimation}"  ///
"\specialcell{Participant \\ FE}"  ///
"\specialcell{Decision \\ Day FE}"  ///
"\specialcell{Later \\ Decisions}"  ///
)


*********** Error on ln(effort)
 
esttab est8001_2_1 est8002_2_1 est8003_2_1 est8004_2_1 /// 
using "c:/graphics/workovertime2Final/errorlogeffort_2.tex", ///
nolines ///
stats(widfe dayfe soph later N N_clust ll betatest betahattest alphatest deltatest, ///
labels("\\[-6pt] Participant FE" "\\[-6pt] Day FE" "\\[-6pt] Prediction Soph." "\\[-6pt] Later Decisions" "\midrule Observations" "Participants" "Log Likelihood" "\midrule \$H_{0}\$(\$\hat{\beta}\$ =1)" "\$H_{0}\$(\$\widehat{\beta_{h}}\$ =1)" "\$H_{0}\$(\$\widehat{\alpha}\$ =0)" "\$H_{0}\$(\$\widehat{\delta}\$ =1)") fmt(0)) ///
order(beta:_cons betahat:_cons delta:_cons gamma:_cons phi:_cons alpha:_cons) ///
cells( b(fmt(%9.3f %9.3f %9.3f %9.3f %9.0f %9.3f)) se(par fmt(%9.3f %9.3f %9.3f %9.3f %9.0f %9.3f ))) ///
width(0.8\hsize) drop(sigma:_cons wid_* ddn*)  replace booktabs gaps f collabels(none) eqlabels(none) ///
varlabels(beta:_cons "Present Bias \$\beta\$" ///
betahat:_cons "Naive Pres. Bias \$\beta_{h}\$"  ///
delta:_cons "Discount Factor \$\delta\$"  ///
gamma:_cons "Cost Curvature \$\gamma\$"  ///
phi:_cons "Cost Slope \$\varphi\$"  ///
alpha:_cons "Proj Task Reduction \$\tilde{\alpha}\$") ///
style(tex) legend label compress substitute(\_ _) ///
mtitle("\specialcell{Initial \\ Estimation}"  ///
"\specialcell{Participant \\ FE}"  ///
"\specialcell{Decision \\ Day FE}"  ///
"\specialcell{Later \\ Decisions}"  ///
)




*********** Projection Bias 1

esttab est9001_2_1 est9002_2_1 est9003_2_1 est9005_2_1 est9004_2_1 /// 
using "c:/graphics/workovertime2Final/structural_pb1_2.tex", ///
nolines ///
stats(widfe dayfe soph later N N_clust ll betatest betahattest alphatest deltatest, ///
labels("\\[-6pt] Participant FE" "\\[-6pt] Day FE" "\\[-6pt] Prediction Soph." "\\[-6pt] Later Decisions" "\midrule Observations" "Participants" "Log Likelihood" "\midrule \$H_{0}\$(\$\hat{\beta}\$ =1)" "\$H_{0}\$(\$\widehat{\beta_{h}}\$ =1)" "\$H_{0}\$(\$\widehat{\alpha}\$ =0)" "\$H_{0}\$(\$\widehat{\delta}\$ =1)") fmt(0)) ///
order(beta:_cons betahat:_cons delta:_cons gamma:_cons phi:_cons alpha:_cons) ///
cells( b(fmt(%9.3f %9.3f %9.3f %9.3f %9.0f %9.3f)) se(par fmt(%9.3f %9.3f %9.3f %9.3f %9.0f %9.3f ))) ///
width(0.8\hsize) drop(sigma:_cons wid_* ddn*)  replace booktabs gaps f collabels(none) eqlabels(none) ///
varlabels(beta:_cons "Present Bias \$\beta\$" ///
betahat:_cons "Naive Pres. Bias \$\beta_{h}\$"  ///
delta:_cons "Discount Factor \$\delta\$"  ///
gamma:_cons "Cost Curvature \$\gamma\$"  ///
phi:_cons "Cost Slope \$\varphi\$"  ///
alpha:_cons "Proj Bias \$\alpha\$") ///
style(tex) legend label compress substitute(\_ _) ///
mtitle("\specialcell{Initial \\ Estimation}"  ///
"\specialcell{Participant \\ FE}"  ///
"\specialcell{Decision \\ Day FE}"  ///
"\specialcell{Later \\ Decisions}"  ///
"\specialcell{Pred. \\ Soph.}"  ///
)


*********** Projection Bias 2

esttab est10001_2_1 est10002_2_1 est10003_2_1 est10004_2_1 /// 
using "c:/graphics/workovertime2Final/structural_pb2_2.tex", ///
nolines ///
stats(widfe dayfe soph later N N_clust ll betatest betahattest alphatest deltatest , ///
labels("\\[-6pt] Participant FE" "\\[-6pt] Day FE" "\\[-6pt] Prediction Soph." "\\[-6pt] Later Decisions" "\midrule Observations" "Participants" "Log Likelihood" "\midrule \$H_{0}\$(\$\hat{\beta}\$ =1)" "\$H_{0}\$(\$\widehat{\beta_{h}}\$ =1)" "\$H_{0}\$(\$\widehat{\alpha}\$ =0)" "\$H_{0}\$(\$\widehat{\delta}\$ =1)") fmt(0)) ///
order(beta:_cons betahat:_cons delta:_cons gamma:_cons phi:_cons alpha:_cons) ///
cells( b(fmt(%9.3f %9.3f %9.3f %9.3f %9.0f %9.3f)) se(par fmt(%9.3f %9.3f %9.3f %9.3f %9.0f %9.3f ))) ///
width(0.8\hsize) drop(sigma:_cons wid_* ddn*)  replace booktabs gaps f collabels(none) eqlabels(none) ///
varlabels(beta:_cons "Present Bias \$\beta\$" ///
betahat:_cons "Naive Pres. Bias \$\beta_{h}\$"  ///
delta:_cons "Discount Factor \$\delta\$"  ///
gamma:_cons "Cost Curvature \$\gamma\$"  ///
phi:_cons "Cost Slope \$\varphi\$"  ///
alpha:_cons "Proj Bias \$\alpha\$") ///
style(tex) legend label compress substitute(\_ _) ///
mtitle("\specialcell{Initial \\ Estimation}"  ///
"\specialcell{Participant \\ FE}"  ///
"\specialcell{Decision \\ Day FE}"  ///
"\specialcell{Later \\ Decisions}"  ///
)


//Differing locations of FE

esttab est11001_2_1 est11002_2_1 est11003_2_1 est11004_2_1 est11005_2_1 est11006_2_1 /// 
using "c:/graphics/workovertime2Final/structural_firstfe_2.tex", ///
nolines ///
stats(gamma_wid phi_wid gamma_day phi_day N N_clust ll betatest betahattest alphatest deltatest, ///
labels("\\[-6pt] \$\gamma\$ Participant FE" "\\[-6pt] \$\varphi\$ Participant FE" "\\[-6pt] \$\gamma\$ Day FE" "\\[-6pt] \$\varphi\$ Day FE" "\midrule Observations" "Participants" "Log Likelihood" "\midrule \$H_{0}\$(\$\hat{\beta}\$ =1)" "\$H_{0}\$(\$\widehat{\beta_{h}}\$ =1)" "\$H_{0}\$(\$\widehat{\alpha}\$ =0)" "\$H_{0}\$(\$\widehat{\delta}\$ =1)") fmt(0)) ///
order(beta:_cons betahat:_cons delta:_cons gamma:_cons phi:_cons alpha:_cons) ///
cells( b(fmt(%9.3f %9.3f %9.3f %9.3f %9.0f %9.3f)) se(par fmt(%9.3f %9.3f %9.3f %9.3f %9.0f %9.3f ))) ///
width(0.8\hsize) drop(sigma:_cons wid_* ddn*)  replace booktabs gaps f collabels(none) eqlabels(none) ///
varlabels(beta:_cons "Present Bias \$\beta\$" ///
betahat:_cons "Naive Pres. Bias \$\beta_{h}\$"  ///
delta:_cons "Discount Factor \$\delta\$"  ///
gamma:_cons "Cost Curvature \$\gamma\$"  ///
phi:_cons "Cost Slope \$\varphi\$"  ///
alpha:_cons "Proj Bias \$\tilde{\alpha}\$") ///
style(tex) legend label compress substitute(\_ _) ///
mtitle("" "" "" "" "" "" ///
)



*********** Different functional form: 2nd degree poly

esttab est12001_2_1 est12002_2_1 est12003_2_1 est12004_2_1 /// 
using "c:/graphics/workovertime2Final/secdegpoly_2.tex", ///
nolines ///
stats(widfe dayfe soph later N N_clust ll betatest betahattest alphatest deltatest, ///
labels("\\[-6pt] Participant FE" "\\[-6pt] Day FE" "\\[-6pt] Prediction Soph." "\\[-6pt] Later Decisions" "\midrule Observations" "Participants" "Log Likelihood" "\midrule \$H_{0}\$(\$\hat{\beta}\$ =1)" "\$H_{0}\$(\$\widehat{\beta_{h}}\$ =1)" "\$H_{0}\$(\$\widehat{\alpha}\$ =0)" "\$H_{0}\$(\$\widehat{\delta}\$ =1)") fmt(0)) ///
order(beta:_cons betahat:_cons delta:_cons psi1:_cons psi2:_cons alpha:_cons) ///
cells( b(fmt(%9.3f %9.3f %9.3f %9.3f %9.4f %9.3f)) se(par fmt(%9.3f %9.3f %9.3f %9.3f %9.4f %9.3f ))) ///
width(0.8\hsize) drop(sigma:_cons wid_* ddn*)  replace booktabs gaps f collabels(none) eqlabels(none) ///
varlabels(beta:_cons "Present Bias \$\beta\$" ///
betahat:_cons "Naive Pres. Bias \$\beta_{h}\$"  ///
delta:_cons "Discount Factor \$\delta\$"  ///
psi1:_cons "Poly Var 1 \$\psi_{1}\$"  ///
psi2:_cons "Poly Var 2 \$\psi_{2}\$"  ///
alpha:_cons "Proj Task Reduction \$\tilde{\alpha}\$") ///
style(tex) legend label compress substitute(\_ _) ///
mtitle("\specialcell{Initial \\ Estimation}"  ///
"\specialcell{Participant \\ FE}"  ///
"\specialcell{Decision \\ Day FE}"  ///
"\specialcell{Later \\ Decisions}"  ///
)



*********** Different functional form: 3rd degree poly

esttab est13001_2_1 est13002_2_1 est13003_2_1 est13004_2_1 /// 
using "c:/graphics/workovertime2Final/thirddegpoly_2.tex", ///
nolines ///
stats(widfe dayfe soph later N N_clust ll betatest betahattest alphatest deltatest, ///
labels("\\[-6pt] Participant FE" "\\[-6pt] Day FE" "\\[-6pt] Prediction Soph." "\\[-6pt] Later Decisions" "\midrule Observations" "Participants" "Log Likelihood" "\midrule \$H_{0}\$(\$\hat{\beta}\$ =1)" "\$H_{0}\$(\$\widehat{\beta_{h}}\$ =1)" "\$H_{0}\$(\$\widehat{\alpha}\$ =0)" "\$H_{0}\$(\$\widehat{\delta}\$ =1)") fmt(0)) ///
order(beta:_cons betahat:_cons delta:_cons psi1:_cons psi2:_cons  psi3:_cons alpha:_cons) ///
cells( b(fmt(%9.3f %9.3f %9.3f %9.3f %9.4f %9.6f %9.3f)) se(par fmt(%9.3f %9.3f %9.3f %9.3f %9.4f %9.6f %9.3f ))) ///
width(0.8\hsize) drop(sigma:_cons wid_* ddn*)  replace booktabs gaps f collabels(none) eqlabels(none) ///
varlabels(beta:_cons "Present Bias \$\beta\$" ///
betahat:_cons "Naive Pres. Bias \$\beta_{h}\$"  ///
delta:_cons "Discount Factor \$\delta\$"  ///
psi1:_cons "Poly Var 1 \$\psi_{1}\$"  ///
psi2:_cons "Poly Var 2 \$\psi_{2}\$"  ///
psi3:_cons "Poly Var 2 \$\psi_{3}\$"  ///
alpha:_cons "Proj Task Reduction \$\tilde{\alpha}\$") ///
style(tex) legend label compress substitute(\_ _) ///
mtitle("\specialcell{Initial \\ Estimation}"  ///
"\specialcell{Participant \\ FE}"  ///
"\specialcell{Decision \\ Day FE}"  ///
"\specialcell{Later \\ Decisions}"  ///
)


*****Separate Deltas 

esttab est14001_2_1 est14002_2_1 est14003_2_1 est14004_2_1 est14005_2_1 est14006_2_1 /// 
using "$savedir/aggresults_delta_m_e.tex", ///
nolines ///
stats(widfe dayfe soph later N N_clust ll betatest betahattest alphatest deltamtest deltaetest deltametest, ///
labels("\\[-6pt] Participant FE" "\\[-6pt] Work Day FE" "\\[-6pt] Prediction Soph." "\\[-6pt] Later Decisions" "\midrule Observations" "Participants" "Log Likelihood" "\midrule \$H_{0}\$(\$\hat{\beta}\$ =1)" "\$H_{0}\$(\$\widehat{\beta_{h}}\$ =1)" "\$H_{0}\$(\$\widehat{\alpha}\$ =0)" "\$H_{0}\$(\$\widehat{\delta_{m}}\$ =1)" "\$H_{0}\$(\$\widehat{\delta_{e}}\$ =1)" "\$H_{0}\$(\$\widehat{\delta_{m}}\$ =\$\widehat{\delta_{e}}\$)") fmt(0)) ///
order(beta:_cons betahat:_cons delta_m:_cons delta_e:_cons gamma:_cons phi:_cons alpha:_cons) ///
cells( b(fmt(%9.3f %9.3f %9.3f %9.3f %9.3f %9.0f %9.3f)) se(par fmt(%9.3f %9.3f %9.3f %9.3f %9.3f %9.0f %9.3f ))) ///
width(0.8\hsize) drop(sigma:_cons wid_* ddn* )  replace booktabs gaps f collabels(none) eqlabels(none) ///
varlabels(beta:_cons "Present Bias \$\beta\$" ///
betahat:_cons "Naive Pres. Bias \$\beta_{h}\$"  ///
delta_m:_cons "Money Discount \$\delta_{m}\$"  ///
delta_e:_cons "Effort Discount \$\delta_{e}\$"  ///
gamma:_cons "Cost Curvature \$\gamma\$"  ///
phi:_cons "Cost Slope \$\varphi\$"  ///
alpha:_cons "Proj Task Reduct. \$\tilde{\alpha}\$") ///
style(tex) legend label compress substitute(\_ _) ///
mtitle("\specialcell{Initial \\ Estimation \\ \$\delta_{m}\$=\$\delta_{e}\$}"  ///
"\specialcell{Initial \\ Estimation \\ \$\delta_{m}\$, \$\delta_{e}\$}"  ///
"\specialcell{Initial \\ Estimation \\ \$\delta_{m}\$=1, \$\delta_{e}\$}"  ///
"\specialcell{Decision \\ Day FE \\ \$\delta_{m}\$=\$\delta_{e}\$}"   ///
"\specialcell{Decision \\ Day FE \\  \$\delta_{m}\$, \$\delta_{e}\$}"   ///
"\specialcell{Decision \\ Day FE \\ \$\delta_{m}\$=1, \$\delta_{e}\$}"   ///
)






***Risk Aversion

esttab est15001_2_1 est15002_2_1 est15003_2_1 est15004_2_1 /// 
using "c:/graphics/workovertime2Final/risk_aversion_2.tex", ///
nolines ///
stats(widfe dayfe soph later N N_clust ll betatest betahattest alphatest deltatest, ///
labels("\\[-6pt] Participant FE" "\\[-6pt] Day FE" "\\[-6pt] Prediction Soph." "\\[-6pt] Later Decisions" "\midrule Observations" "Participants" "Log Likelihood" "\midrule \$H_{0}\$(\$\hat{\beta}\$ =1)" "\$H_{0}\$(\$\widehat{\beta_{h}}\$ =1)" "\$H_{0}\$(\$\widehat{\alpha}\$ =0)" "\$H_{0}\$(\$\widehat{\delta}\$ =1)") fmt(0)) ///
order(beta:_cons betahat:_cons delta:_cons gamma:_cons phi:_cons alpha:_cons) ///
cells( b(fmt(%9.3f %9.3f %9.3f %9.3f %9.0f %9.3f)) se(par fmt(%9.3f %9.3f %9.3f %9.3f %9.0f %9.3f ))) ///
width(0.8\hsize) drop(sigma:_cons wid_* ddn*)  replace booktabs gaps f collabels(none) eqlabels(none) ///
varlabels(beta:_cons "Present Bias \$\beta\$" ///
betahat:_cons "Naive Pres. Bias \$\beta_{h}\$"  ///
delta:_cons "Discount Factor \$\delta\$"  ///
gamma:_cons "Cost Curvature \$\gamma\$"  ///
phi:_cons "Cost Slope*Wealth \$\varphi \cdot exp(-a\cdot y)\$"  ///
alpha:_cons "Proj Task Reduction \$\tilde{\alpha}\$") ///
style(tex) legend label compress substitute(\_ _) ///
mtitle("\specialcell{Initial \\ Estimation}"  ///
"\specialcell{Participant \\ FE}"  ///
"\specialcell{Decision \\ Day FE}"  ///
"\specialcell{Later \\ Decisions}"  ///
"\specialcell{Pred. \\ Soph.}"  ///
)



***If we screwed up and confounded present bias and projection bias because of the timing

esttab est16001_2_1 est16002_2_1 est16003_2_1 est16004_2_1 est16005_2_1 /// 
using "c:/graphics/workovertime2Final/structural_confound_pb_2.tex", ///
nolines ///
stats(widfe dayfe soph later N N_clust ll betatest betahattest deltatest, ///
labels("\\[-6pt] Participant FE" "\\[-6pt] Day FE" "\\[-6pt] Prediction Soph." "\\[-6pt] Later Decisions" "\midrule Observations" "Participants" "Log Likelihood" "\midrule \$H_{0}\$(\$\hat{\beta}\$ =1)" "\$H_{0}\$(\$\widehat{\beta_{h}}\$ =1)" "\$H_{0}\$(\$\widehat{\delta}\$ =1)") fmt(0)) ///
order(beta:_cons betahat:_cons delta:_cons gamma:_cons phi:_cons) ///
cells( b(fmt(%9.3f %9.3f %9.3f %9.3f %9.0f)) se(par fmt(%9.3f %9.3f %9.3f %9.3f %9.0f))) ///
width(0.8\hsize) drop(sigma:_cons wid_* ddn*)  replace booktabs gaps f collabels(none) eqlabels(none) ///
varlabels(beta:_cons "Present Bias \$\beta\$" ///
betahat:_cons "Naive Pres. Bias \$\beta_{h}\$"  ///
delta:_cons "Discount Factor \$\delta\$"  ///
gamma:_cons "Cost Curvature \$\gamma\$"  ///
phi:_cons "Cost Slope \$\varphi\$") ///
style(tex) legend label compress substitute(\_ _) ///
mtitle("\specialcell{Initial \\ Estimation}"  ///
"\specialcell{Participant \\ FE}"  ///
"\specialcell{Decision \\ Day FE}"  ///
"\specialcell{Later \\ Decisions}"  ///
"\specialcell{Pred. \\ Soph.}"  ///
)





***No projection bias

esttab est17001_2_1 est17002_2_1 est17003_2_1 est17004_2_1 est17005_2_1 /// 
using "c:/graphics/workovertime2Final/structural_remove_pb_2.tex", ///
nolines ///
stats(widfe dayfe soph later N N_clust ll betatest betahattest deltatest, ///
labels("\\[-6pt] Participant FE" "\\[-6pt] Day FE" "\\[-6pt] Prediction Soph." "\\[-6pt] Later Decisions" "\midrule Observations" "Participants" "Log Likelihood" "\midrule \$H_{0}\$(\$\hat{\beta}\$ =1)" "\$H_{0}\$(\$\widehat{\beta_{h}}\$ =1)" "\$H_{0}\$(\$\widehat{\delta}\$ =1)") fmt(0)) ///
order(beta:_cons betahat:_cons delta:_cons gamma:_cons phi:_cons) ///
cells( b(fmt(%9.3f %9.3f %9.3f %9.3f %9.0f)) se(par fmt(%9.3f %9.3f %9.3f %9.3f %9.0f))) ///
width(0.8\hsize) drop(sigma:_cons wid_* ddn*)  replace booktabs gaps f collabels(none) eqlabels(none) ///
varlabels(beta:_cons "Present Bias \$\beta\$" ///
betahat:_cons "Naive Pres. Bias \$\beta_{h}\$"  ///
delta:_cons "Discount Factor \$\delta\$"  ///
gamma:_cons "Cost Curvature \$\gamma\$"  ///
phi:_cons "Cost Slope \$\varphi\$") ///
style(tex) legend label compress substitute(\_ _) ///
mtitle("\specialcell{Initial \\ Estimation}"  ///
"\specialcell{Participant \\ FE}"  ///
"\specialcell{Decision \\ Day FE}"  ///
"\specialcell{Later \\ Decisions}"  ///
"\specialcell{Pred. \\ Soph.}"  ///
)


*****Work date num

esttab est18001_2_1 est18002_2_1 est18003_2_1 est18005_2_1 est18004_2_1 /// 
using "$savedir/aggresults_wdn.tex", ///
nolines ///
stats(widfe dayfe soph later N N_clust ll betatest betahattest alphatest deltatest, ///
labels("\\[-6pt] Participant FE" "\\[-6pt] Work Day FE" "\\[-6pt] Prediction Soph." "\\[-6pt] Later Decisions" "\midrule Observations" "Participants" "Log Likelihood" "\midrule \$H_{0}\$(\$\hat{\beta}\$ =1)" "\$H_{0}\$(\$\widehat{\beta_{h}}\$ =1)" "\$H_{0}\$(\$\widehat{\alpha}\$ =0)" "\$H_{0}\$(\$\widehat{\delta}\$ =1)") fmt(0)) ///
order(beta:_cons betahat:_cons delta:_cons gamma:_cons phi:_cons alpha:_cons) ///
cells( b(fmt(%9.3f %9.3f %9.3f %9.3f %9.0f %9.3f)) se(par fmt(%9.3f %9.3f %9.3f %9.3f %9.0f %9.3f ))) ///
width(0.8\hsize) drop(sigma:_cons wid_* wdn*)  replace booktabs gaps f collabels(none) eqlabels(none) ///
varlabels(beta:_cons "Present Bias \$\beta\$" ///
betahat:_cons "Naive Pres. Bias \$\beta_{h}\$"  ///
delta:_cons "Discount Factor \$\delta\$"  ///
gamma:_cons "Cost Curvature \$\gamma\$"  ///
phi:_cons "Cost Slope \$\varphi\$"  ///
alpha:_cons "Proj Task Reduction \$\tilde{\alpha}\$") ///
style(tex) legend label compress substitute(\_ _) ///
mtitle("\specialcell{Initial \\ Estimation}"  ///
"\specialcell{Participant \\ FE}"  ///
"\specialcell{Decision \\ Day FE}"  ///
"\specialcell{Later \\ Decisions}"  ///
"\specialcell{Pred. \\ Soph.}"  ///
)

