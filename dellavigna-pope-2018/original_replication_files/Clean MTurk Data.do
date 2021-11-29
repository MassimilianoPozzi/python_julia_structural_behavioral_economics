
***************************************************************************************************
* 
* Program: Clean MTurk Data.do
* Authors: Devin Pope
* Last Modified: 3/7/2017
* Purpose:
*       Upload and clean the 10,000-person MTurk experiment data
*
* Files Used:
*     1. input/MTurkData.dta
*
***************************************************************************************************

cd your_foldername_here

*ssc install labutil

********************************************
set more off
capture mkdir "dtafiles"

*************************************
*Define the treatment order
local tx_labels "No Payment" "1c PieceRate" "10c PieceRate" "4c PieceRate" "Very Low Pay" "1c RedCross" "10c RedCross" "Gift Exchange" "1c 2Wks" "1c 4Wks" "Gain 40c" "Loss 40c" "Gain 80c" "Prob.01 $1" "Prob.5 2c" "Social Comp" "Ranking" "Task Signif"

*************************************
*************************************

**************Loading and cleaning the data********************

use "input/MTurkData", clear

*generate a new variable that are nice labels for the treatments
gen stringtreatment = string(treatment)
drop treatment
rename stringtreatment treatment

*define order of treatments for tables - same for forecasters and MTurk
gen tx_order = 1 if treatment == "1.3"
replace tx_order = 2 if treatment == "1.1" 
replace tx_order = 3 if treatment == "1.2" 
replace tx_order = 4 if treatment == "1.4" 
replace tx_order = 5 if treatment == "2" 
replace tx_order = 6 if treatment == "3.1" 
replace tx_order = 7 if treatment == "3.2" 
replace tx_order = 8 if treatment == "10" 
replace tx_order = 9 if treatment == "4.1" 
replace tx_order = 10 if treatment == "4.2" 
replace tx_order = 11 if treatment == "5.1" 
replace tx_order = 12 if treatment == "5.2" 
replace tx_order = 13 if treatment == "5.3" 
replace tx_order = 14 if treatment == "6.1" 
replace tx_order = 15 if treatment == "6.2" 
replace tx_order = 16 if treatment == "7" 
replace tx_order = 17 if treatment == "8" 
replace tx_order = 18 if treatment == "9" 

*define treatment names
gen treatmentname=""
local txnum = 1
foreach t_label in "`tx_labels'" {
replace treatmentname = "`t_label'" if tx_order==`txnum'
local txnum = `txnum' + 1
}

labmask tx_order, value(treatmentname)

***********Correcting for Qualtrics Error**************************
*While the survey was in progress, Qualtrics encountered an error which affected
*a large number of survey takers. We will filter out these results

*Time variables
gen double starttime = clock(startdate, "YMDhms")
format starttime %tc

*how many are we dropping?
count if starttime >= tc(18may2015 15:34:00) & starttime < tc(19may2015 10:52:00)

*drop them
drop if starttime >= tc(18may2015 15:34:00) & starttime < tc(19may2015 10:52:00)

gen double endtime = clock(enddate, "YMDhms")
format endtime %tc

gen duration_time_ms = endtime - starttime
gen duration_time_minutes = duration_time_ms/60000


*real quick, I want to create a counter for number of observations per treatment (this will be used further below for tests of attrition)
bys treatmentname: egen initial_obs_count = count(_N) 


****************Dealing with quitters and cheaters*******************

*We want to drop anybody who cheated and did more than 4,000 button pushes
drop if buttonpresses > 4000 & buttonpresses < 9000

*we drop them if they dropped out before even getting a treatment
drop if treatment == "."

*We drop if mturkid_code == . because this means they didn't enter in their id number and we don't know if they came 
*in and out multiple times to get a good treatment
drop if mturkid_code == .

*We want to drop anybody who didn't finish the survey and then started it again - and deleted their cookies in the 
*meantime (so that qualtrics didn't just start them again from where they left off). So, we are going to delete 
*anybody with more than 2 observations for the same mturkid
bys mturkid_code: egen counter = count(_n)
drop if counter >= 2
drop counter

***Now, let's restrict the sample to people who completed the survey.
keep if finished == 1

*We want to drop anyone with an NA button push count (this is equivalent to 0 button pushes)
drop if buttonpresses == .

**restrict further to people who had their HIT accepted
**HITs might be rejected if they entered an invalid survey code
**or did not submit in time
keep if approved == 1


*********We now have a final sample of 9,861 people*************
**************Let's analyze the results*************************
compress
saveold dtafiles/MTurkCleanedData, replace

keep buttonpresses treatment*
saveold dtafiles/MTurkCleanedDataShort, replace version(12)
