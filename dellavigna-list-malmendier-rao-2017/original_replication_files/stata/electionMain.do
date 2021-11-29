clear
set more off

*****************************************************************************
** "Voting To Tell Others"
** by Stefano DellaVigna, John A. List, Ulrike Malmendier, and Gautam Rao
**
** This do file loads the data from the field experiments,
** It generates the reduced-form results, as well as the moments used in the
** structural estimation
*****************************************************************************

* Additional packages needed to run this code: 
*net describe dm79 //(includes svmat2)  STB-56 dm79
*ssc install distplot
*ssc install outreg2
*net from http://www.stata-press.com/data/gps/ml3/   // mysureg from maximum likelihood package
*net install ml3_prog.pkg

*cd into folder
*cd "your-directory-here"

capture mkdir cleaneddata
capture mkdir output

******************************************************************
*** STEP 0. Merge Experimental Results To Voting History ***
******************************************************************

** Load file with voting behavior for 2011
use inputdata/ElectionData2011VotingHistoryCleanDec15, clear
su

* Generate a proxy of party affiliation using primary voting, prioratizing the most recent affiliation
gen reg="N"
foreach x in mar2004 mar2006 feb2008 feb2010 {
	replace reg="R" if `x'=="R"
	replace reg="D" if `x'=="D"
}
* Generate turnout indicator for all the elections in the record, including primaries
foreach x in mar2004 nov2004 feb2005 mar2006 nov2006 apr2007 feb2008 nov2008 apr2009 feb2010 nov2010 {
rename `x' `x'old
	gen `x'=.
	replace `x'=1 if `x'old~="N" & `x'old~="" 
	replace `x'=0 if `x'old=="N"
}
* There are missing values for older voting records, since they come from different data set from the County
* Affects variables mar2004 nov2004 feb2005 nov2006 
* These variables have the same missing observations, replace to 0, and add indicator for missing value
gen nov2006dmiss=(nov2006old=="" & nov2010old~="")
foreach x in mar2004 nov2004 feb2005 nov2006 {
replace `x'=0 if nov2006dmiss==1
}
drop *old
su


* Collapse to household level
keep row household_id reg
bys household_id: gen no=_n
sum no
local max_no = `r(max)'
reshape wide reg, i(household_id) j(no)
gen dhhRep=0
gen dhhDem=0
forvalues x = 1/`max_no' {
replace dhhRep=1 if reg`x'=="R"
replace dhhDem=1 if reg`x'=="D"
}
label var dhhRep "Indicator for whether a household includes a voter who voted in a Rep. primary"
label var dhhDem "Indicator for whether a household includes a voter who voted in a Dem. primary"

keep row dhhRep dhhDem

*Merge in the main survey data file
merge 1:1 row using inputdata/exp2011Dec15
ta _m
drop _m
su


* Order variables
order row treat* location year date hour solicitor charity nosol toel /*
*/ voter* answer saidyes dnd voted* *ask*
su
compress
desc,f
save cleaneddata/exp2011WithVoting, replace




******************************************************************
*** STEP 1. Reduced-Form Data Analysis of Experimental Results ***
******************************************************************

use cleaneddata/exp2011WithVoting, replace

** Initial Sample
count
ta date

*For main analysis, drop two groups of observations which the solicitors did not contact
* I. Households with no solicitor sign
* II. Households where solicitor could not knock on door (big dog barking, house for sale, etc.)
* Important to eliminate these observations because these households are disproportionally
* in No-Warning treatment (since some are excluded after flyering), 
* and have no giving/no survey completion
* nosolsign==1 indicates that these households have a sign that they do not accept solicitors
gen samplemain = 1
ta treatment nosol, row
replace samplemain=0 if nosol==1

ta treatment toel,row
replace samplemain=0 if toeliminate==1

drop nosol toeliminate
* Final Sample
count if samplemain==1

* Identify full sample (without dropping)
gen samplenodrop=1


*** Create treatment flags 
gen treatmindoll=""
foreach min in "5" "10" {
foreach doll in "0" "5" "10" {
	replace treatmindoll="`doll'd`min'm" if treatmmin==`min' & treatmdollar==`doll'
}
}

** Generate fixed effects
egen sodate=concat(solicitor date)
egen byte grsol=group(solicitor)
egen byte grdatloc=group(date location)
egen byte grdatlocsol=group(date location solicitor)
egen byte grhour=group(hour)
* drop solicitor area* hour 

* Define variables which apply both to presidential and congressional elections
gen voter=voterin2010
gen voted=votedin2010
label var voted "Self-reported voting"
gen lie=(voted~=voter) if voter~=. & voted~=.
gen period="Summer" if date<mdy(9,1,2011)
replace period="Fall" if date>=mdy(9,1,2011)



* Party affiliation in the primary
gen primary="Other"
replace primary="Rep" if dhhRep==1 & dhhDem==0
replace primary="Dem" if dhhDem==1 & dhhRep==0

* Save final file
save cleaneddata/exp2011WithVotingFinalSample, replace


*************
* Figures and Tables for Main Sample
keep if samplemain==1

*** Document sample sizes in different treatments for main sample
* Figure 3a
ta treatflyer,m
ta treatmindoll,m
ta treatmindoll treatv,m
* Figure 3b
ta treatmmin treatinfo,m
ta treatmmin treatinfo if saidyes==1,m


*** Document voter and non-voter sampling
* Figures 4, 5, 7
ta voterin2010,m
* Figure 6 (only included if either dhhRep=1 or dhhDem=1, not both)
ta dhhRep dhhDem if voterin2010==1,m
* Figure 8
ta voterin2010 if saidyes==1 & votedin2010~=.,m


*** Coherence checks
ta answer saidyes,m
ta answer votedin2010,m
* there are 5 cases of surveys answered with no recorded response to the voting 
*     question (but responses to the other questions) 
ta votedin2010 saidyes,m


** Generate treatment variables for figures and tables
foreach x in "Nw" "W" "We" "Oo" "Ooe" {
	gen dflyer`x'=treatflyer=="`x'"
	}
foreach z in "0d5m" "10d10m" "10d5m" {
	gen dmd`z'=treatmindoll=="`z'"
	}
foreach z in "Nv" "V" {
	gen dnv`z'=treatv=="`z'"
	}
foreach z in "No" "V" {
	gen dlt`z'=treatinfo=="`z'"
	}
gen dlie=treatinfo=="V"
* Generate alternate treatment variable for flyer treatments
gen flyerelect=(treatflyer=="We"|treatflyer=="Ooe")
gen flyer=(treatflyer=="W"|treatflyer=="We")
gen optout=(treatflyer=="Oo"|treatflyer=="Ooe")


*** FIGURES

* Figure 4 -- by survey duration and payment
xi: reg answer if voter==1, robust cluster(sodate)
outreg2 using output/Figure5, excel side label ctitle(placeholder,delete,"") se noaster bdec(4) rdec(4) replace
foreach y in "answer" "saidyes" {
foreach voter in "1" "0" {
xi: reg `y' dmd* if voter==`voter', nocons robust cluster(sodate)
outreg2 dmd* using output/Figure4, excel side label ctitle(`y',Voter`voter',noCtrl) se noaster bdec(4) rdec(4) append
}
}

* Figure 5 -- by flyer type
xi: reg answer if voter==1, robust cluster(sodate)
outreg2 using output/Figure6, excel side label ctitle(placeholder,delete,"") se noaster bdec(4) rdec(4) replace
foreach y in "answer" "saidyes" dnd {
foreach voter in "1" "0" {
xi: reg `y' dflyer* if voter==`voter', nocons robust cluster(sodate)
outreg2 dflyer* using output/Figure5, excel side label ctitle(`y',Voter`voter',noCtrl) se noaster bdec(4) rdec(4) append
}
}

* Figure 6 -- by voter type
xi: reg answer if voter==1, robust cluster(sodate)
outreg2 using output/Figure7, excel side label ctitle(placeholder,delete,"") se noaster bdec(4) rdec(4) replace
foreach y in "answer" "saidyes" dnd {
foreach voter in "1" {
foreach party in "Rep" "Dem" "Other" {
xi: reg `y' dflyer* if voter==`voter' & primary=="`party'", nocons robust cluster(sodate)
outreg2 dflyer* using output/Figure6, excel side label ctitle(`y',`party',Voter`voter',noCtrl) se noaster bdec(4) rdec(4) append
}
}
}

* Figure 7 -- by announcement at the door
xi: reg answer if voter==1, robust cluster(sodate)
outreg2 using output/Figure8, excel side label ctitle(placeholder,delete,"") se noaster bdec(4) rdec(4) replace
foreach y in "saidyes" {
foreach voter in "1" "0" {
xi: reg `y' dnv* if voter==`voter', nocons robust cluster(sodate)
outreg2 dnv* using output/Figure7, excel side label ctitle(`y',Voter`voter',noCtrl) se noaster bdec(4) rdec(4) append
}
}

* Figure 8 -- by lying incentives
xi: reg lie if voter==1, robust cluster(sodate)
outreg2 using output/Figure9, excel side label ctitle(placeholder,delete,"") se noaster bdec(4) rdec(4) replace
foreach y in "lie" {
foreach voter in "1" "0" {
foreach min in 10 5 {
xi: reg `y' dlt* if voter==`voter' & treatmmin==`min', nocons robust cluster(sodate)
outreg2 dlt* using output/Figure8, excel side label ctitle(`y',Voter`voter',min`min',noCtrl) se noaster bdec(4) rdec(4) append
}
}
}





*** TABLES
* drop omitted group variables
drop *dflyerNw *dmd0d5m


** TABLE 1
** Generate Survey Table on answering the door and survey completion
xi: reg answer if voter==1, robust cluster(sodate)
outreg2 using output/Table1, excel label ctitle(Place-Holder-Delete) se alpha(.01,.05,.10) bdec(4) rdec(4) replace
* Answering the door
foreach voter in "1" "0" {
foreach ctr in "" "i.grsol i.grdatloc i.grhour" {
xi: reg answer dmd* flyer optout flyerelect `ctr' if voter==`voter', robust cluster(sodate)
local contr="YES"
if "`ctr'"=="" {
local contr="NO"
}
outreg2 dmd* flyer optout flyerelect using output/Table1, excel label ctitle(Voter`voter',noCtrl) addtext(Fixed Effects, `contr') se alpha(.01,.05,.10) bdec(4) rdec(4) append
}
}
* Completing the survey
foreach voter in "1" "0" {
foreach ctr in "" "i.grsol i.grdatloc i.grhour" {
xi: reg saidyes dmd* flyer optout flyerelect dnvV `ctr' if voter==`voter', robust cluster(sodate)
local contr="YES"
if "`ctr'"=="" {
local contr="NO"
}
outreg2 dmd* flyer optout flyerelect dnvV using output/Table1, excel label ctitle(Voter`voter',noCtrl) addtext(Fixed Effects, `contr') se alpha(.01,.05,.10) bdec(4) rdec(4) append
}
}


** TABLE 2
** Generate Survey Table -- lying in voting
* Panel A -- Pooled across different tduration-payment treatments for the survey
xi: reg lie dlie if voter==1, robust cluster(sodate)
outreg2 dlie using output/Table2, excel label ctitle(Placeholder-Delete) se alpha(.01,.05,.10) bdec(4) rdec(4) replace
foreach voter in "1" "0" {
xi: reg lie dlie if voter==`voter', robust cluster(sodate)
outreg2 dlie using output/Table2, excel label ctitle(Voter`voter',noCtrl) addtext(Panel, All) se alpha(.01,.05,.10) bdec(4) rdec(4) append
xi: reg lie dlie i.grdatloc if voter==`voter', robust cluster(sodate)
outreg2 dlie using output/Table2, excel label ctitle(Voter`voter',Ctrl) addtext(Panel, All) se alpha(.01,.05,.10) bdec(4) rdec(4) append
}
* Panels B-C-D -- By tduration-payment treatments for the survey
foreach x in "0d5m" "10d5m" "10d10m" {
foreach voter in "1" "0" {
xi: reg lie dlie if voter==`voter' & treatmindoll=="`x'", robust cluster(sodate)
outreg2 dlie using output/Table2, excel label ctitle(Voter`voter',noCtrl) addtext(Panel, `x') se alpha(.01,.05,.10) bdec(4) rdec(4) append
xi: reg lie dlie i.grdatloc if voter==`voter' & treatmindoll=="`x'", robust cluster(sodate)
outreg2 dlie using output/Table2, excel label ctitle(Voter`voter',Ctrl) addtext(Panel, `x') se alpha(.01,.05,.10) bdec(4) rdec(4) append
}
}

** Online Appendix Table 1 -- Robustness Checks to Table 1
xi: reg answer if voter==1, robust cluster(sodate)
outreg2 dmd* dflyer* using output/OnlAppTable1, excel label ctitle(Place-Holder-Delete) se alpha(.01,.05,.10) bdec(4) rdec(4) replace
* Answering the door
foreach voter in "1" "0" {
foreach ctr in "i.grsol i.grdatloc i.grhour" "i.grdatlocsol i.grhour" {
xi: reg answer dmd* flyer optout flyerelect `ctr' if voter==`voter', robust cluster(sodate)
local contr="YES-Benchmark"
if "`ctr'"=="i.grdatlocsol i.grhour" {
local contr="YES-Loc*Date*Sol"
}
outreg2 dmd* flyer optout flyerelect using output/OnlAppTable1, excel label ctitle(answer,Voter`voter') addtext(Fixed Effects, `contr') se alpha(.01,.05,.10) bdec(4) rdec(4) append
}
}
* Completing the survey
foreach voter in "1" "0" {
foreach ctr in "i.grsol i.grdatloc i.grhour" "i.grdatlocsol i.grhour" {
xi: reg saidyes dmd* flyer optout flyerelect dnvV `ctr' if voter==`voter', robust cluster(sodate)
local contr="YES-Benchmark"
if "`ctr'"=="i.grdatlocsol i.grhour" {
local contr="YES-Loc*Date*Sol"
}
outreg2 dmd* flyer optout flyerelect dnvV using output/OnlAppTable1, excel label ctitle(saidyes,Voter`voter') addtext(Fixed Effects, `contr') se alpha(.01,.05,.10) bdec(4) rdec(4) append
}
}


** Online Appendix Table 2 -- Table 1 and 2, split by time period
** Generate Survey Table on answering the door and survey completion
xi: reg answer if voter==1, robust cluster(sodate)
outreg2 using output/OnlAppTable2, excel label ctitle(Place-Holder-Delete) se alpha(.01,.05,.10) bdec(4) rdec(4) replace
* Answering the door
foreach voter in "1" "0" {
foreach ctr in "i.grsol i.grdatloc i.grhour" {
foreach time in "Summer" "Fall" {
xi: reg answer dmd* flyer optout flyerelect `ctr' if voter==`voter' & period=="`time'", robust cluster(sodate)
outreg2 dmd* dmd* flyer optout flyerelect using output/OnlAppTable2, excel label ctitle(answer,Voter`voter',Ctrl) addtext(Period, `time') se alpha(.01,.05,.10) bdec(4) rdec(4) append
}
}
}
* Completing the survey
foreach voter in "1" "0" {
foreach ctr in "i.grsol i.grdatloc i.grhour" {
foreach time in "Summer" "Fall" {
xi: reg saidyes dmd* flyer optout flyerelect dnvV `ctr' if voter==`voter' & period=="`time'", robust cluster(sodate)
outreg2 dmd* flyer optout flyerelect dnvV using output/OnlAppTable2, excel label ctitle(saiyes,Voter`voter',Ctrl) addtext(Period, `time') se alpha(.01,.05,.10) bdec(4) rdec(4) append
}
}
}
* Lying
foreach voter in "1" "0" {
foreach time in "Summer" "Fall" {
xi: reg lie dlie i.grdatloc if voter==`voter' & period=="`time'", robust cluster(sodate)
outreg2 dlie using output/OnlAppTable2, excel label ctitle(lie,Voter`voter',Ctrl)  addtext(Period, `time') se alpha(.01,.05,.10) bdec(4) rdec(4) append
}
}


** Online Appendix Table 3 -- Table 1 and 2, split by Rep. and Dem.
** Generate Survey Table on answering the door and survey completion
xi: reg answer if voter==1, robust cluster(sodate)
outreg2 dmd* dflyer* using output/OnlAppTable3, excel label ctitle(Place-Holder-Delete) se alpha(.01,.05,.10) bdec(4) rdec(4) replace
* Answering the door
foreach voter in "1" "0" {
foreach ctr in "i.grsol i.grdatloc i.grhour" {
foreach party in  "Rep" "Dem" "Other"  {
xi: reg answer dmd* flyer optout flyerelect `ctr' if voter==`voter' & primary=="`party'", robust cluster(sodate)
outreg2 dmd* flyer optout flyerelect using output/OnlAppTable3, excel label ctitle(answer,Voter`voter',Ctrl) addtext(Registered,`party') se alpha(.01,.05,.10) bdec(4) rdec(4) append
}
}
* Completing the survey
foreach ctr in "i.grsol i.grdatloc i.grhour" {
foreach party in  "Rep" "Dem" "Other"  {
xi: reg saidyes dmd* flyer optout flyerelect dnvV `ctr' if voter==`voter' & primary=="`party'", robust cluster(sodate)
outreg2 dmd* flyer optout flyerelect dnvV using output/OnlAppTable3, excel label ctitle(saidyes,Voter`voter',Ctrl) addtext(Registered,`party') se alpha(.01,.05,.10) bdec(4) rdec(4) append
}
}
* Lying
foreach party in  "Rep" "Dem" "Other"  {
xi: reg lie dlie i.grdatloc if voter==`voter' & primary=="`party'", robust cluster(sodate)
outreg2 dlie using output/OnlAppTable3, excel label ctitle(lie,Voter`voter',Ctrl)  addtext(Registered,`party') se alpha(.01,.05,.10) bdec(4) rdec(4) append
}
}

** Online Appendix Table 4 -- Robustness Checks to Table 2
egen grdat = group(date)
egen grloc = group(location)

* Panel A -- Pooled across different tduration-payment treatments for the survey
*placeholder column
xi: reg lie dlie if voter==1, robust cluster(sodate)
outreg2 dlie using output/OnlAppTable4, excel label ctitle(Placeholder-Delete) se alpha(.01,.05,.10) bdec(4) rdec(4) replace
foreach voter in "1" "0" {
* Date-Location FE
xi: reg lie dlie i.grdatloc if voter==`voter', robust cluster(sodate)
outreg2 dlie using output/OnlAppTable4, excel label ctitle(Voter`voter',Dat-Loc FE) addtext(Panel, All) se alpha(.01,.05,.10) bdec(4) rdec(4) append
* Sol + Date-Location + Hour FE
xi: reg lie dlie i.grdatloc i.grhour i.grsol if voter==`voter', robust cluster(sodate)
outreg2 dlie using output/OnlAppTable4, excel label ctitle(Voter`voter',Dat-Loc + Hour + Sol FE) addtext(Panel, All) se alpha(.01,.05,.10) bdec(4) rdec(4) append
}






*****************************************
*** STEP 2. Evidence on number of times voters were asked about voting
*****************************************

** 2011 survey data
use cleaneddata/exp2011WithVotingFinalSample, replace
keep if samplemain==1

* Cap each category at 20 and add
foreach x in friends relatives coworkers people {
	foreach y in 2008 2010 {
		cap replace `x'ask`y'=20 if `x'ask`y'>20 & `x'ask`y'~=.
		}
	}
su *ask*
egen totask2010=rsum (friendsask2010 relativesask2010 coworkersask2010 peopleask2010) if ~(friendsask2010==. & relativesask2010==. & coworkersask2010==. & peopleask2010==.)
egen totask2008=rsum (friendsask2008 relativesask2008) if ~(friendsask2008==. & relativesask2008==.)
drop friends* relatives* coworkers* people* 
* Summary stats on how often people are asked whether they voted
su totask*
su totask*,d
* Small differences depending on voting status
bys voterin2010: su totask*
bys voterin2010: su totask*,d
* Small differences depending on primary voter affiliation
bys primary: su totask*,d

*asked at least once in 2010
gen asked2010least1 = totask2010>= 1 if totask2010!=.
tab asked2010least1

*asked more than 10 times in 2010
gen asked2010more10 = totask2010> 10 if totask2010!=.
tab asked2010more10

*asked more than 10 times in 2008
gen asked2008more10 = totask2008> 10 if totask2008!=.
tab asked2008more10

*Note: number of times asked reported in the text and Table 4, 5, Online Appendix 5-7


* Online Appendix Figure 1 
distplot line totask2010 totask2008, clp(l -) title("CDFs of number of times asked about turnout") xtitle("Number of times asked") xline(5) ///
	legend(label(1 "2010 Congr. Elections") label(2 "2008 Presid. Elections")) saving(output/TimesAsked, replace)
graph export output/TimesAsked.png, replace
	



*****************************************
*** STEP 3. Generation of moments and variance-covariance matrix for structural estimation
* Estimation is done in Matlab
*****************************************

* Loop through the main sample and sample without dropping (for robustness)
foreach sample in main nodrop {

use cleaneddata/exp2011WithVotingFinalSample, replace

* Run only on relevant sample
keep if sample`sample'==1

* Rename treatment variable to remove the dash
foreach x in "0" 5 10 {
	foreach y in 5 10 {
		foreach z in "Nw" "W" "We" "Oo" "Ooe" {
			replace treatment="`z'`x'd`y'm" if treatment=="`z'-`x'd`y'm"
			}
		}
	}
ta treatment,m

levelsof treatment
foreach v in `r(levels)' {
	gen tr`v' = treatment=="`v'"
}

//Generate lying treatments
gen lyingTreatment = string(treatmmin) + treatinfo
levelsof lyingTreatment
foreach v in `r(levels)' {
	gen vtrLying`v' = lyingTreatment=="`v'"
} 
drop lyingTreatment 

//Gen labels for {Nw W We Oo Ooe}*treatV (whether informed at door)
gen treatCat=""
foreach x in "0" 5 10 {
	foreach y in 5 10 {
		foreach z in "Nw" "W" "We" "Oo" "Ooe" {
				replace treatCat="`z'" if treatment=="`z'`x'd`y'm"
			}
		}
	}
gen SVContent = treatCat+"_"+treatv
levelsof SVContent
foreach v in `r(levels)' {
	gen vtrSVContent`v' = SVContent=="`v'"
}
drop SVContent

* Number of unique values
quietly: tab sol 
local sol_num = `r(r)'
quietly: tab grdatloc
local datloc_num = `r(r)'
quietly: tab hour
local hour_num = `r(r)'

*** Demeaning the control variables to get the moments evaluated at the mean value of the parameters
foreach x in sol datloc hour {
	local j=`sol_num'*("`x'"=="sol")+`datloc_num'*("`x'"=="datloc")+`hour_num'*("`x'"=="hour")
	forvalues i=1/`j'{
			gen byte gr`x'`i' = (gr`x'==`i')
			egen mgr`x'`i' = mean(gr`x'`i')
			gen dgr`x'`i' = gr`x'`i'-mgr`x'`i'
			replace dgr`x'`i'=0 if dgr`x'`i'==.
			** if this isn't a real category (all zeroes), just drop the new control variable
         summ gr`x'`i'
         if r(max)==0 {
            drop dgr`x'`i'
		}
	}
}



//order moments
aorder tr*
order trW* trNw* trOo*, alpha
aorder vtr*
drop mgr* gr*
aorder dgr*

compress
save cleaneddata/NewVoting2011`sample'.dta, replace


// Set of moments on answering the door and completing the survey
use cleaneddata/NewVoting2011`sample'.dta, replace

foreach vars of varlist trNw0d5m-trWe10d5m {
	gen V`vars'= `vars'*voterin2010
	gen NV`vars'= `vars'*(1-voterin2010)
}
order VtrNw* VtrW* VtrWe* VtrOo* VtrOoe* NVtrNw* NVtrW* NVtrWe* NVtrOo* NVtrOoe*

* sureg to get VC matrix
#delimit ;

mysureg 
(answer VtrNw0d5m-NVtrOoe10d5m dgrdatloc* dgrhour* dgrsol*, noconst)
(saidyes VtrNw0d5m-NVtrOoe10d5m dgrdatloc* dgrhour* dgrsol*, noconst)
(dnd VtrOo0d5m-VtrOoe10d5m NVtrOo0d5m-NVtrOoe10d5m dgrdatloc* dgrhour* dgrsol*, noconst)
, cluster(sodate);
 
#delimit cr

outreg2 Vtr* NVtr* using "output/sample`sample'_ctrl_moments1", excel ctitle(Mom,Ctrl) nose bdec(4) rdec(4) auto(6) replace noaster sideway nopa
matrix list e(V)
matrix VCE=e(V)
drop *
svmat2 VCE, rnames(rownames) full
drop VCE*
svmat VCE, names(eqcol)
drop *sigma*
drop if regexm(rownames,".*dgr.*") | regexm(rownames,".*sigma.*")
drop *dgr*
ds, has(type numeric)
format `r(varlist)' %15.0g
outsheet using "output/sample`sample'_ctrl_moments1_VarCov.csv", comma replace


//Separate regression for survey info treatments
use cleaneddata/NewVoting2011`sample'.dta, replace

foreach vars of varlist vtrSVContentNw_Nv-vtrSVContentWe_V{
	gen V`vars'= `vars'*voterin2010
	gen NV`vars'= `vars'*(1-voterin2010)
}
order VvtrSVContentNw_* VvtrSVContentW_* VvtrSVContentWe_* VvtrSVContentOo_* VvtrSVContentOoe_* NVvtrSVContentNw_* NVvtrSVContentW_* NVvtrSVContentWe_* NVvtrSVContentOo_* NVvtrSVContentOoe_* 


reg saidyes VvtrSVContent* NVvtrSVContent* dgrdatloc* dgrhour* dgrsol*, noconst vce(cluster sodate)
outreg2 tr* Vvtr* NVvtr* using "output/sample`sample'_ctrl_moments2", excel nose bdec(4) rdec(4) auto(6) replace noaster sideway nopa
matrix list e(V)
matrix VCE=e(V)
drop *
svmat2 VCE, rnames(rownames) full 
drop VCE*
svmat VCE, names(string)

*Note: only want to keep the 20 relevant treatment variables
local max_n = _N
drop if regexm(rownames,".*dgr.*") | regexm(rownames,".*sigma.*")
drop string21-string`max_n'

ds, has(type numeric)
format `r(varlist)' %15.0g
outsheet using "output/sample`sample'_ctrl_moments2_VarCov.csv", comma replace


//2011 lying about vote moments with date-location FE 
use cleaneddata/NewVoting2011`sample'.dta, replace

foreach vars of varlist vtrLying5No-vtrLying10V {
	gen V`vars'= `vars'*voterin2010
	gen NV`vars'= `vars'*(1-voterin2010)
}

* Demeaned date-loc variables for the whole sample
reg lie VvtrLying* NVvtrLying* dgrdatloc*, noconst vce(cluster sodate)
outreg2 Vvtr* NVvtr* using "output/sample`sample'_ctrl_moments3", excel nose bdec(4) rdec(4) auto(6) replace noaster sideway nopa

matrix list e(V)
matrix VCE=e(V)
drop *
svmat2 VCE, rnames(rownames) full
drop VCE*
svmat VCE, names(string)

*Note: only want to keep the 8 relevant treatment variables
local max_n = _N
drop if regexm(rownames,".*dgr.*") | regexm(rownames,".*sigma.*")
drop string9-string`max_n'

ds, has(type numeric)
format `r(varlist)' %15.0g
outsheet using "output/sample`sample'_ctrl_moments3_VarCov.csv", comma replace

}





*****************************************
*** STEP 4. Get-out-the-vote experiment
*****************************************

*** 2010 GOTV
use inputdata/GOTVData2010Dec15, clear
su
* Generate turnout indicator for all the elections in the record, including primaries
foreach x in mar2004 nov2004 feb2005 mar2006 nov2006 apr2007 feb2008 nov2008 apr2009 feb2010 nov2010 {
rename `x' `x'old
	gen `x'=.
	replace `x'=1 if `x'old~="N" & `x'old~="" 
	replace `x'=0 if `x'old=="N"
}

*summarize past + present voting by treatment
tabstat mar2004 nov2004 feb2005 mar2006 nov2006 apr2007 feb2008 nov2008 apr2009 feb2010 nov2010, stat(mean n) by(treatment)

* There are missing values for older voting records, since they come from different data set from the County
* These variables have the same missing observations, replace to 0, and add indicator for missing value
gen nov2006dmiss=(nov2006old=="" & nov2010old~="")
foreach x in mar2004 nov2004 feb2005 nov2006 {
replace `x'=0 if nov2006dmiss==1
}
su
drop *old
su

* Level of randomization is route*city
egen route_id=group(route city)

** Treatments
gen tf = treatment=="TF"
gen cf = treatment=="CF"
tab treatment


*** Table 6
* Note: coeff and SE from these regressions used for turnout and GOTV moments
xi: reg nov2010 cf tf, robust cluster(route_id)
test cf-tf=0
outreg2 cf tf using output/Table6,excel label ctitle(2010,noctr) se alpha(.01,.05,.10) dec(4) replace
xi: reg nov2010 mar2004-feb2010 nov2006dmiss cf tf, robust cluster(route_id)
test cf-tf=0
outreg2 cf tf using output/Table6, excel label ctitle(2010,Control) se alpha(.01,.05,.10) dec(4) append


*** 2012 GOTV
use inputdata/GOTVData2012Dec15, clear
su
* Generate turnout indicator for all the elections in the record, including primaries
foreach x in mar2004 nov2004 feb2005 mar2006 nov2006 apr2007 feb2008 nov2008 apr2009 feb2010 nov2010 apr2011 mar2012 nov2012 {
rename `x' `x'old
	gen `x'=.
	replace `x'=1 if `x'old~="N" & `x'old~="" 
	replace `x'=0 if `x'old=="N"
}

*summarize past + present voting by treatment
tabstat mar2004 nov2004 feb2005 mar2006 nov2006 apr2007 feb2008 nov2008 apr2009 feb2010 nov2010 apr2011 mar2012 nov2012, stat(mean n) by(treatment)

* There are missing values for older voting records, since they come from different data set from the County
* These variables have the same missing observations, replace to 0, and add indicator for missing value
gen nov2006dmiss=(nov2006old=="" & nov2010old~="")
foreach x in mar2004 nov2004 feb2005 mar2006 nov2006 apr2007 {
replace `x'=0 if nov2006dmiss==1
}
su
drop *old
su

* Level of randomization is route*city
egen route_id=group(route city)

** Treatments
gen tf = treatment=="TF"
gen cf = treatment=="CF"
tab treatment

*** Table 6
xi: reg nov2012 cf tf, robust cluster(route_id)
test cf-tf=0
outreg2 cf tf using output/Table6, excel label ctitle(2012,noctr) se alpha(.01,.05,.10) dec(4) append
xi: reg nov2012 mar2004-mar2012 nov2006dmiss cf tf, robust cluster(route_id)
test cf-tf=0
outreg2 cf tf using output/Table6, excel label ctitle(2012,Control) se alpha(.01,.05,.10) dec(4) append

