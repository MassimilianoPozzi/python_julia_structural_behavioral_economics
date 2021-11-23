program NBetaHatStructural1_agg
	version 9.1
	args lnf $setupvars
	tempvar netdistance effort wage prediction bonus today pb sigma_total group_for_test sophfixedweight ///
	        predgamma truegamma predphi truephi added e_future e_present p_future p_utility cur_effort ///
		  diff bound1 bound2 cur_utility cur_utility_diff best_utility_diff best_effort past_utility ///
		  derivative weight predchoice truechoice prob bound1 e_future_uncertainty ///
		  tasks mc expected_mc param_realized mc_realized add diff lastdiff param_with_uncertainty ///
		  candidate_base_param mb bound1 diff lastdiff candidate_pred gamma_l gamma_h prob_h prob_l derv_h ///
		  lbound ubound cur_effort discountFactor ///			
		  derv_l predpred_uncertainty alpha_ ///
		  mb_nonpara mc0 mc1 mc2 mc3 mc4 mc5 mc6 mc7 mc8 mc9 mc10 mc11 mc12 mc13 mc14 mc15 mc16 mc17 mc18 mc19 mc20 ///
		  delta_m_var distanceToPayment distanceToWork ///
		  t_x1 t_x2 t_x3 t_x t_w t_v t_f t_e t_stat productLog
		  


	local quietcommand ""
	if ($loud1!=1) {
		local quietcommand "qui"
	}
			
		
	`quietcommand' {
	

		***********
		**** Setup
		***********

		gen double `netdistance'=$ML_y1
		gen double `effort'=$ML_y2
		gen double `wage'=$ML_y3
		gen double `prediction'=$ML_y4
		gen double `bonus'=$ML_y5
		gen double `today'=$ML_y6
		gen double `pb'=$ML_y7
		gen double `group_for_test'=$ML_y8
		gen double `sophfixedweight'=$ML_y9
		gen double `distanceToPayment'=$ML_y10
		gen double `distanceToWork'=$ML_y11



		*********fixes effect
		if ($fixedbeta1==1) {
			tempvar beta
			gen double `beta'=1
		}
		if ($fixeddelta1==1) {
			tempvar delta
			gen double `delta'=1
		}
		
		
		***********
		**** Base predictions depending on location of error
		***********

		**** We have to determine the predicted effort choices given the data (if error is on effort)
		**** Note: effort *includes* mandatory work as "effort" is +10
		**** Present, future, and predictions
		
		if ($monetarydelta) {   			//appendix: separate monetary and effort delta
			gen double `e_present'=((`phi'*((`delta_m'^(`distanceToPayment'))/((`delta_e'^(`distanceToWork'))))*(`beta')*`wage')^(1/(`gamma'-1)))	
			gen double `e_future'=((`phi'*((`delta_m'^(`distanceToPayment'))/((`delta_e'^(`distanceToWork'))))*`wage')^(1/(`gamma'-1)))	
			gen double `p_future'=((`phi'*((`delta_m'^(`distanceToPayment'))/((`delta_e'^(`distanceToWork'))))*`betahat'*`wage')^(1/(`gamma'-1)))
		} 
		else if ($seconddegpoly) {			//appendix: different cost specification
			gen double `e_present'=(-(1/(`beta'*(`delta'^(`netdistance'))))*`psi1'+`wage')/((1/(`beta'*(`delta'^(`netdistance'))))*`psi2')						
			gen double `e_future'=(-(1/(`delta'^(`netdistance')))*`psi1'+`wage')/((1/(`delta'^(`netdistance')))*`psi2')						
			gen double `p_future'=(-(1/(`betahat'*(`delta'^(`netdistance'))))*`psi1'+`wage')/((1/(`betahat'*(`delta'^(`netdistance'))))*`psi2')					
		}
		else if ($thirddegpoly) {			//appendix: different cost specification
			//qui gen double  `predchoice'=(-`discountFactor'*(`psi2')+(((`discountFactor'^2)*((`psi2'^2)-4*`psi1'*`psi3')+`discountFactor'*(4*`wage'*`psi3'))^(1/2)))/(2*`discountFactor'*`psi3')			
			gen double `e_present'=(-(1/(`beta'*(`delta'^(`netdistance'))))*(`psi2')+((((1/(`beta'*(`delta'^(`netdistance'))))^2)*((`psi2'^2)-4*`psi1'*`psi3')+(1/(`beta'*(`delta'^(`netdistance'))))*(4*`wage'*`psi3'))^(1/2)))/(2*(1/(`beta'*(`delta'^(`netdistance'))))*`psi3')			
			gen double `e_future'=(-(1/(1*(`delta'^(`netdistance'))))*(`psi2')+((((1/(1*(`delta'^(`netdistance'))))^2)*((`psi2'^2)-4*`psi1'*`psi3')+(1/(1*(`delta'^(`netdistance'))))*(4*`wage'*`psi3'))^(1/2)))/(2*(1/(1*(`delta'^(`netdistance'))))*`psi3')			
			gen double `p_future'=(-(1/(`betahat'*(`delta'^(`netdistance'))))*(`psi2')+((((1/(`betahat'*(`delta'^(`netdistance'))))^2)*((`psi2'^2)-4*`psi1'*`psi3')+(1/(`betahat'*(`delta'^(`netdistance'))))*(4*`wage'*`psi3'))^(1/2)))/(2*(1/(`betahat'*(`delta'^(`netdistance'))))*`psi3')			

		}
		else { 								//this is the default: used almost universally
			gen double `e_present'=((`phi'*(`delta'^(`netdistance'))*(`beta')*`wage')^(1/(`gamma'-1)))	
			gen double `e_future'=((`phi'*(`delta'^(`netdistance'))*`wage')^(1/(`gamma'-1)))			
			gen double `p_future'=((`phi'*(`delta'^(`netdistance'))*`betahat'*`wage')^(1/(`gamma'-1)))
		}
		
		
		gen double `predchoice'=`e_future'		
		qui replace `predchoice'=`e_present' if `today'==1
		qui replace `predchoice'=`p_future' if `prediction'==1 	
		
		
		
		//this is when we have monetary curvature
		if ($riskaversionalpha!=0) {
			gen double `t_x1'=exp((10*$riskaversionalpha*`wage')/(`gamma'-1))*(($riskaversionalpha*`wage')/(`gamma'-1))*((`phi'*(1/(1/(`beta'*(`delta'^(`netdistance')))))*`wage')^(1/(`gamma'-1)))
			gen double `t_x2'=exp((10*$riskaversionalpha*`wage')/(`gamma'-1))*(($riskaversionalpha*`wage')/(`gamma'-1))*((`phi'*(1/(1/(1*(`delta'^(`netdistance')))))*`wage')^(1/(`gamma'-1)))
			gen double `t_x3'=exp((10*$riskaversionalpha*`wage')/(`gamma'-1))*(($riskaversionalpha*`wage')/(`gamma'-1))*((`phi'*(1/(1/(`betahat'*(`delta'^(`netdistance')))))*`wage')^(1/(`gamma'-1)))
			
			gen double `t_x'=	`t_x2'		
			qui replace `t_x'=	`t_x1' if `today'==1
			qui replace `t_x'=	`t_x3' if `prediction'==1 	
			
			//issue: need to use Haley's method for Lambert's w
			//starting guess
			gen double `t_w' = 1 
			gen double `t_v'=.
			gen double `t_f'=.
			gen double `t_e'=.
			gen double `t_stat'=.
	
			local farthest=1
			local count=1 		//this just stops loops that aren't going anywhere. normally converges in like 5 steps
			while (`farthest' > 1e-8 & `count'<100) {
			   local count=`count'+1
			   qui replace `t_v' = `t_w'
			   qui replace `t_e' = exp(`t_w')
			   qui replace `t_f' = `t_w'*`t_e' - `t_x'  // Iterate to make this quantity zero
			   qui replace `t_w' = `t_w' - `t_f'/((`t_e'*(`t_w'+1) - (`t_w'+2)*`t_f'/(2*`t_w'+2)))
			   qui replace `t_stat'=abs(`t_w' - `t_v')/abs(`t_w')
			   qui sum `t_stat'
			   local farthest=r(max)
			   disp "`farthest'"
			}
			qui rename `t_w' `productLog'
			qui replace `predchoice'=((`gamma'-1)/($riskaversionalpha*`wage'))*`productLog'		
		}
		
		
	


		//get the separate delta for each day
		if ($separatedelta) {
			forvalues num=1(1)4 {
				//future
				replace `predchoice'=((`phi'*(`delta`num'')*`wage')^(1/(`gamma'-1))) if `netdistance'==`num' 
				//predictions
				replace `predchoice'=((`phi'*(`delta`num'')*`betahat'*`wage')^(1/(`gamma'-1))) if `netdistance'==`num' & `prediction'==1
			}
		}

		
	
	

		
		***Projection Bias type 2
		if ($structuralpb2) {
		
			//need to force alpha or will have real issues around alpha=1
			//without constraint, clearly heading positive
			//gen double `alpha_'=1/(1+exp(-`alpha'))
		
			replace `e_present'=(((`phi'*(`delta'^(`netdistance'))*(`beta')*`wage')-(`alpha')*(`pb')^(`gamma'-1))/(1-`alpha'))^(1/(`gamma'-1))
			replace `e_future'=(((`phi'*(`delta'^(`netdistance'))*`wage')-(`alpha')*(`pb')^(`gamma'-1))/(1-`alpha'))^(1/(`gamma'-1))			
			replace `p_future'=(((`phi'*(`delta'^(`netdistance'))*`betahat'*`wage')-(`alpha')*(`pb')^(`gamma'-1))/(1-`alpha'))^(1/(`gamma'-1))
		
			replace `predchoice'=`e_future'		
			replace `predchoice'=`e_present' if `today'==1
			replace `predchoice'=`p_future' if `prediction'==1
			
			//prob is that negatives to a power are imaginary -> missing
			replace `predchoice'=0 if `predchoice'==.
			
		}
		
		gen double `truechoice'=(`effort')


		**** We have to determine the predicted gammma "choices" given the data (if error is on gamma)
		**** Present, future, and predictions
		if ($error_gamma) {			
			qui gen double `predgamma'=1+(ln(`phi'*(`delta'^(`netdistance'))*`wage')/ln(`effort'))
			qui replace `predgamma'=1+(ln(`phi'*(`delta'^(`netdistance'))*(`beta')*`wage')/ln(`effort')) if `today'==1
			qui replace `predgamma'=1+(ln(`phi'*(`delta'^(`netdistance'))*(`betahat')*`wage')/ln(`effort')) if `prediction'==1
			
			if ($pb) {
				qui replace `predgamma'=1+(ln(`phi'*(`delta'^(`netdistance'))*`wage')/ln(`effort'+`alpha')) if `pb'>0
				qui replace `predgamma'=1+(ln(`phi'*(`delta'^(`netdistance'))*(`beta')*`wage')/ln(`effort'+`alpha')) if `today'==1 & `pb'>0
				qui replace `predgamma'=1+(ln(`phi'*(`delta'^(`netdistance'))*(`betahat')*`wage')/ln(`effort'+`alpha')) if `prediction'==1 & `pb'>0
			}
			
			qui gen double `truegamma'=`gamma'
			
		}


		**** We have to determine the predicted phi "choices" given the data (if error is on phi)
		**** Present, future, and predictions
		if ($error_phi) {			
			qui gen double `predphi'=((`effort')^(`gamma'-1))/((`delta'^(`netdistance'))*`wage')
			qui replace `predphi'=((`effort')^(`gamma'-1))/((`delta'^(`netdistance'))*(`beta')*`wage') if `today'==1
			qui replace `predphi'=((`effort')^(`gamma'-1))/((`delta'^(`netdistance'))*(`betahat')*`wage') if `prediction'==1
			
			if ($pb) {
				qui replace `predphi'=((`effort'-`alpha')^(`gamma'-1))/((`delta'^(`netdistance'))*`wage') if `pb'>0
				qui replace `predphi'=((`effort'-`alpha')^(`gamma'-1))/((`delta'^(`netdistance'))*(`beta')*`wage') if `today'==1 & `pb'>0
				qui replace `predphi'=((`effort'-`alpha')^(`gamma'-1))/((`delta'^(`netdistance'))*(`betahat')*`wage') if `prediction'==1 & `pb'>0
			}
			
			gen double `truephi'=`phi'
			sum `predphi' `phi'
		}



		***********
		**** Sophistication 
		***********
		if ($soph) {
			**guy's percieved future self's utility from predicted effort [this is what guy will do with no further incentives]
			gen double `p_utility'=(`phi'*(`delta'^(`netdistance'))*`betahat'*`wage'*`p_future')-(1/`gamma')*((`p_future')^(`gamma'))
			**need to find the effort level, that when combined to bonus, gives closest to the p_utility
			gen double `cur_utility'=0
			gen double `cur_utility_diff'=0
			gen double `best_effort'=0
			gen double `best_utility_diff'=100000000
			gen double `past_utility'=0
			gen double `derivative'=10
			
			qui gen double `lbound'=1
			//need to span a large range.
			qui gen double `ubound'=801
			qui gen double `cur_effort'=.
			
			
			foreach interval of numlist 80 40 20 10 5 1 .1 .01 .001 {
				local reps=((`ubound'[1]-`lbound'[1])/`interval') 
				qui replace `cur_utility'=0
				qui replace `cur_utility_diff'=0
				qui replace `best_effort'=0
				qui replace `best_utility_diff'=100000000
				qui replace `past_utility'=0
				qui replace `derivative'=10
				forvalues repnum=0(1)`reps' {
					
					qui replace `cur_effort'=`lbound'+`interval'*`repnum'
					qui replace `cur_utility'=((`phi'*(`delta'^(`netdistance'))*`betahat'*(`wage'*`cur_effort'+`bonus'))-(1/`gamma')*((`cur_effort')^(`gamma')))
					qui replace `cur_utility_diff'=`p_utility'-`cur_utility'
					qui replace `derivative'=`cur_utility'-`past_utility'
					qui replace `best_effort'=`cur_effort'-`interval' if abs(`cur_utility_diff')<abs(`best_utility_diff')  & `derivative'<0
					qui replace `best_utility_diff'=`cur_utility_diff' if abs(`cur_utility_diff')<abs(`best_utility_diff') 	  & `derivative'<0
					qui replace `past_utility'=`cur_utility'
					//qui replace `best_effort'=`e_future' if `best_effort'<1 & `e_future'!=.
					
				}
				qui replace `lbound'=`best_effort'-2*`interval'
				qui replace `ubound'=`best_effort'+2*`interval'	
			}
			qui replace `best_effort'=`e_future' if `e_future'<`best_effort'
			//qui replace `best_effort'=`e_future' if `e_future'>`best_effort' 
			
			//correction
			gen double `added'=(1-`betahat')*100
			replace `added'=5 if `added'>5
			replace `added'=-5 if `added'<-5	
			//smooth out added amount - its +5 if b>1 and -5 if b<1
			replace `added'=5-10*(1/(1+2.71828^-(100*(`betahat'-1))))
			replace `best_effort'=`best_effort'+`added' 
			
			//now, I have the predicted choice given future decision and prediction.
			if ($sophfixed) {
				gen `weight'=`sophfixedweight'
			} 
			else { 				//this is for trying to estimate theta: never used
				gen double `weight'=(1/(1+exp(-`theta'))) 				
				//gen double `weight'=`theta'
			}
			replace `predchoice'=`weight'*`best_effort'+(1-`weight')*`p_future' if `prediction'==1 
		}

		
		***********
		**** Uncertainty
		***********
		if ($uncertainty) {
			
			***********
			**** Future Decisions
			***********
			**** we are trying to find the new level of effort with uncertainty that gives the same MU as the old level without uncertainty
			local cond1=" if (`today'==0) & (`prediction'==0) "
			local cond2=" if (`prediction'==1)"
			local bigiters=8
			local smalliters=8
			
			qui gen double `bound2'=`predchoice' `cond1'
			qui gen double `param_with_uncertainty'=`predchoice' `cond1'	 
			qui gen double `mb'=(`phi'*(`delta'^(`netdistance'))*`wage') `cond1'	 
			local vars "candidate_base_param bound1 expected_mc param_realized mc_realized add diff lastdiff "
		      foreach var of local vars {
		      	qui gen double ``var''=0 `cond1'
		    	}
			if ($error_gamma) {			
				replace `bound2'=`predgamma' `cond1'
				replace `param_with_uncertainty'=`predgamma' `cond1'
				replace `bound1'=1 `cond1'
			}

			**** Double cycle through potential parameter guesses
			forvalues bigiter=1(1)`bigiters' {
				qui replace `lastdiff'=-100 `cond1'
				forvalues smalliter=0(1)`smalliters' {
					qui replace `candidate_base_param'=`bound1'+`smalliter'*((`bound2'-`bound1')/`smalliters') `cond1'

					**** The "parameter" is either gamma or effort, depending on the error structure
					**** Given subject chooses param=candidate_base, they will have error - we need the expected MC = MB

					qui replace `expected_mc'=0 `cond1'
					forvalues i=0(1)$num1 {
						qui replace `param_realized'=`candidate_base_param'+gridmid[(`i'+1),1]*`sigma2' `cond1'
						if ($error_gamma) {
							qui replace `mc_realized'=`effort'^(`param_realized'-1) `cond1'
						} 
						else {
							qui replace `mc_realized'=`param_realized'^(`gamma'-1) `cond1'
						}
						qui replace `add'=PDF[(`i'+1),1]*`mc_realized' `cond1'
						qui replace `expected_mc'=`expected_mc'+`add' `cond1'			
					} 
					qui replace `diff'=`mb'-`expected_mc' `cond1'
					qui replace `param_with_uncertainty'=`candidate_base_param' if (`diff'<=0  & `lastdiff'>0)
					qui replace `lastdiff'=`diff' `cond1'
				}
				qui replace `bound1'=`param_with_uncertainty'-(`bound2'-`bound1')/`smalliters' `cond1'
				qui replace `bound2'=`param_with_uncertainty' `cond1'
			}
			if ($error_gamma) {
				qui replace `predgamma'=`param_with_uncertainty' `cond1'
			} 
			else {
				qui replace `predchoice'=`param_with_uncertainty' `cond1'
			}

		
		}

		*****
		***Things that will change predicted choices - different methods
		*****

		***Projection Bias
		if ($pb) {
			replace `predchoice'=`predchoice'-`alpha' if `pb'>0
		}
		if ($structuralpb1) {
			replace `predchoice'=`predchoice'-`alpha'*10 if `pb'>0
		}
		
		***Testing for some group to have a different B
		if ($testdifference) {
			replace `predchoice'=((`phi'*(`delta'^(`netdistance'))*(`beta'+`betaleft')*`wage')^(1/(`gamma'-1))) if `group_for_test'==1 & `today'==1
		}
		
		***Logged decisions
		if ($ln) {
			replace `predchoice'=ln(`predchoice'+10)
			replace `truechoice'=ln(`truechoice')
			if ($error_phi) {
				replace `predphi'=ln(`predphi')
				replace `truephi'=ln(`truephi')
			}
		}
		
		***Error
		gen double `sigma_total'=`sigma'
		if ($double_error) {
			replace `sigma_total'=(`sigma'^2+`sigma2'^2)^(1/2) if `today'==1 
		}

		***Log Error
		if (($double_error) & ($log_error)) {
			replace `sigma_total'=(`sigma'/ln(`effort'))
			replace `sigma_total'=((`sigma'/ln(`effort'))^2+`sigma2'^2)^(1/2) if `today'==1 
			bys `today': sum `sigma_total'
		}
		

		*****
		***Converting into probabilities
		*****

		gen double `prob'=normalden(`truechoice',`predchoice',`sigma_total')
		***If Tobit
		if ($tobit) {
			replace `prob'=1-normal((`predchoice'-(`truechoice'))/`sigma_total') if `effort'==10
			replace `prob'=normal((`predchoice'-(`truechoice'))/`sigma_total') if `effort'==110
		}
		
		***When Error is on Gamma
		if ($error_gamma) {
			replace `prob'=normalden(`truegamma',`predgamma',`sigma_total')
			if ($tobit) {
				//**************IS THIS WRONG??? I THINK I MIGHT HAVE TO REVERSE THESE 1-normal and then normal 
				replace `prob'=normal((`predgamma'-(`truegamma'))/`sigma_total') if `effort'==10
				replace `prob'=1-normal((`predgamma'-(`truegamma'))/`sigma_total') if `effort'==110
			}
		}

		***When Error is on phi
		if ($error_phi) {
			replace `prob'=normalden(`truephi',`predphi',`sigma_total')
			if ($tobit) {
				replace `prob'=1-normal((`predphi'-(`truephi'))/`sigma_total') if `effort'==10
				replace `prob'=normal((`predphi'-(`truephi'))/`sigma_total') if `effort'==110
			}
		}


	}

	quietly replace `lnf' = ln(`prob')

	
	//RETURN RESULTS?
	if ($loud2==1) {
		local parameters "lnf prob beta betahat gamma phi delta delta_m delta_e alpha theta sigma sigma2"
        	foreach param of local parameters {
			capture confirm variable ``param''
			if (!_rc) {
                		qui sum ``param''
				local `param'1=round(r(mean),.001) 
				local N=r(N) 
			}
        		}
		local command ""
		foreach param of local parameters {
			capture confirm variable ``param''
			if (!_rc) {
				local command "`command' `param'=``param'1'" 
			}
        		}
		disp "`command' N=`N'"	
		qui sum `lnf'
		disp r(N)*r(mean)
	}


end





