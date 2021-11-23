do setupAndDefinePrograms.do


use "Data/decisions_data_w_ind", clear
tab leftdatenum if first==1


//graphs all the rest of the guys
use "Data/decisions_data_w_ind", clear
//standard relationship
reg jobschosen wage if ml_prob==1
reg jobschosen wage if ml_prob==0

keep if ml_prob==1
levelsof wid, local(levels)
mat coeffs_wids = J(28,2,1)
local count=0
foreach wid of local levels {
	local count=`count'+1
	cap {
	local conditions "(1==1)" 
	foreach condition of local conditions {
        reg jobschosen wage if wid==`wid' & `condition'
		cap matrix define coeff=e(b) 
		
		matrix coeffs_wids[`count',1]=coeff[1,1]
		matrix coeffs_wids[`count',2]=`wid'
	}
	}
}

mata : st_matrix("coeffs_wids", sort(st_matrix("coeffs_wids"), 1))
noisily matrix list coeffs_wids

//here are the scattters of all of the ml_prob=0 people in order
forvalues i=1(1)28 {
	local wid=coeffs_wids[`i',2]
	sum ml_prob if wid==`wid'
	local removed=""
	if r(mean)>0 {
		local removed="" 
	}
	sum left if wid==`wid'
	local left=""
	if r(mean)>0 {
		local left="[Left Early]" 
	}
	
	//have to do it this way to force between 0 and 100
	cap drop pred
	qui reg jobschosen wage if wid==`wid'
	predict pred
	replace pred=. if pred<0
	replace pred=. if pred>100
	replace wage=.3 if wage>.3 & wid==`wid' 
	
	twoway (scatter jobschosen wage if wid==`wid' & type=="present", m(circle_hollow) mc(gs5)) (scatter jobschosen wage if wid==`wid' & type=="future", m(diamond_hollow) mc(gs5)) (scatter jobschosen wage if wid==`wid' & type=="prediction", m(triangle_hollow) mc(gs5)) (line pred wage if wid==`wid', sort lc(gs5)), legend(off) ytitle("") xtitle("") title("ID: `wid' `removed' `left'", size(medium))
	graph save Temp/t`wid'.gph, replace
}


local command="graph combine "
forvalues i=1(1)28 {
	local wid=coeffs_wids[`i',2]
	local command="`command' Temp/t`wid'.gph"
}

`command', ycommon xcommon cols(4) ysize(7) scale(1.2) imargin(0)
graph export "$savedir/removed_scatters.png", replace width(3000)


