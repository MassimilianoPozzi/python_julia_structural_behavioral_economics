***************************************************************************************************
* 
* Program: Clean Forecast Data.do
* Authors: Devin Pope & Stefano DellaVigna
* Purpose:
*       Clean the Expert Forecasts
*
* Files Used:
*     1. input/ExpertForecast

* Last Modified: 3/31/2017
***************************************************************************************************


********************************************
************ Setup Analysis ****************
********************************************
cd your_foldername_here

********************************************
set more off
set matsize 800

capture mkdir "dtafiles"

*************************************
*Define the treatment order
local tx_labels "No Payment" "1c PieceRate" "10c PieceRate" "4c PieceRate" "Very Low Pay" "1c RedCross" "10c RedCross" "Gift Exchange" "1c 2Wks" "1c 4Wks" "Gain 40c" "Loss 40c" "Gain 80c" "Prob.01 $1" "Prob.5 2c" "Social Comp" "Ranking" "Task Signif"

********************************************
******WORK BEFORE RESHAPING DATA************
********************************************

use input/ExpertForecast, clear

*****************
*Clean the Demographic Vars, Duration, and Confidence

	*Rank variables for experts
	gen assistant = (Rank == "Assistant Professor") & sample=="experts"
	gen associate = (Rank == "Associate Professor") & sample=="experts"
	gen professor = (Rank == "Professor") & sample=="experts"
	gen other = (assistant == 0 & associate == 0 & professor == 0) & sample=="experts"

	*Field variables
	gen std_econ = (Field1 == "Applied Microeconomics" | Field1 == "Economic Theory") & sample=="experts"
	gen appl = (Field1 == "Applied Microeconomics") & sample=="experts"
	gen theory = (Field1 == "Economic Theory") & sample=="experts"
	gen lab_econ = (Field1 == "Lab Experiments") & sample=="experts"
	gen beh_econ = (Field1 == "Behavioral Economics" | Field1 == "Behavioral Finance") & sample=="experts"
	gen psych = (Field1 == "Decision-Making" | Field1 == "Psychology") & sample=="experts"

	gen psych_dm = Field1 == "Decision-Making" & sample=="experts"
	gen psych_only = Field1 == "Psychology" & sample=="experts"
	gen beh_econ_only = Field1 == "Behavioral Economics" & sample=="experts"
	gen beh_fin = Field1 == "Behavioral Finance" & sample=="experts"

	*Identify Subsamples
	gen expert=sampleshort=="experts"

	*top code the duration variable to the 90th percentile
	replace duration=min(max(0,duration),50) if duration~=.
	rename duration duration


*****************
*Clean the Treatment Vars	

	* Create variables that contain the actual effort levels for each treatment
	gen treatment_t1_actual = 1521
	gen treatment_t2_actual = 2029
	gen treatment_t3_actual = 2175
	gen treatment_t4_actual = 2132 
	gen treatment_t5_actual = 1883
	gen treatment_t6_actual = 1907
	gen treatment_t7_actual = 1918
	gen treatment_t8_actual = 2004
	gen treatment_t9_actual = 1970
	gen treatment_t10_actual = 2136
	gen treatment_t11_actual = 2155
	gen treatment_t12_actual = 2188
	gen treatment_t13_actual = 1896
	gen treatment_t14_actual = 1977
	gen treatment_t15_actual = 1848
	gen treatment_t16_actual = 1761
	gen treatment_t17_actual = 1740
	gen treatment_t18_actual = 1602

	forvalues i=1/3 {
		gen actual`i'=treatment_t`i'_actual 
	}
	forvalues i=4/18 {
		gen forecast`i'=treatment_t`i' 
		gen actual`i'=treatment_t`i'_actual 
	}

* wide data with partial responses - only for experts paper
saveold dtafiles/ExpertForecastCleanWide_all, replace version(12)

* drop the non-matching experts
drop if completed!="Yes"
drop completed

*** Save Wide data set
compress
saveold dtafiles/ExpertForecastCleanWide, replace version(12)

********************************************
******WORK AFTER RESHAPING DATA************
********************************************

	reshape long actual forecast, i(id expert sampleshort) j(treatno)

	gen abs_error=abs(actual-forecast)

	*treatment order for table
	gen tx_order=treatno
	replace tx_order=tx_order+1 if tx_order>7
	replace tx_order = 8 if tx_order==19
	
	*define treatment names
	gen treatmentname=""
	local txnum = 1
	foreach t_label in "`tx_labels'" {
	replace treatmentname = "`t_label'" if tx_order==`txnum'
	local txnum = `txnum' + 1
	}

	* Give labels to treatment groups
	gen treatname_group=""
	replace treatname_group="4c" if treatno==4
	replace treatname_group="Co" if treatno==5
	replace treatname_group="Char" if treatno==6|treatno==7
	replace treatname_group="Time" if treatno==8|treatno==9
	replace treatname_group="RD" if treatno==10|treatno==11|treatno==12
	replace treatname_group="PWt" if treatno==13|treatno==14
	replace treatname_group="Psych" if treatno==15|treatno==16|treatno==17
	replace treatname_group="Ge" if treatno==18

	* Wisdom of Crowds for each subsample 
	*   calculate error metrics using average forecast
	bysort sampleshort treatno: egen WoC_forecast = mean(forecast)
	gen WoC_abs_error = abs(WoC_forecast - actual)
	gen WoC_squared_error = (WoC_forecast - actual)^2
	gen outperform_WoC = (abs_error<WoC_abs_error) if sampleshort=="experts"
	

*** Save Long data set
compress
saveold dtafiles/ExpertForecastCleanLong, replace version(12)
