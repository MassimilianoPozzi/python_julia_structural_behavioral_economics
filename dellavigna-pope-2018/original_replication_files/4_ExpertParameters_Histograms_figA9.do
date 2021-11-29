* histograms.do
* This file plots histograms for implied parameters of the experts' forecasts
clear all
cd "C:\Users\Avner\Dropbox\forecastExp\Data Analysis\Analysis\behavioral"
* Import Data
insheet using "dtafiles\HistogramData.csv", comma names

gen beta = 1.2
replace beta = expertbeta if expertbeta <= 1.2
qui sum beta, d
local median = r(p50)
hist beta, fraction width(0.05) xlabel(0(0.1)1.1 "") text(-0.0075 1.25 ">1.2", place(s)) /*xline(1.1301)
*/ addplot(pci 0 `median' 0.25 `median', lwidth(medthick) lcolor(green) lpattern(solid) || pci 0 1.16911162737392 0.25 1.16911162737392, lwidth(medthick) lcolor(red) lpattern(dash)) xline(1.2, lpattern(dash) lcolor(gs12) style(foreground)) xtitle("Implied Beta") title("Distribution of Experts' Parameters - Beta Discounting") note("Red dashed line indicates the parameter value implied by the actual mTurker mean effort." "Green solid line indicates the median parameter value of the experts." "Observations are truncated at 1.2.") graphregion(color(white)) bgcolor(white) fcolor(navy) lcolor(gs16) legend(off)
graph export "output\expertParametersPowerBetaDiscounting.png", width(2000) height(1500) replace

gen delta = 1.2
replace delta = expertdelta if expertdelta <= 1.2
qui sum delta, d
local median = r(p50)
hist delta, fraction width(0.05) xlabel(0(0.1)1.1 "") text(-0.007 1.23 ">1.2", place(s)) /*xline(0.7614)
*/ addplot(pci 0 `median' 0.25 `median', lwidth(medthick) lcolor(green) lpattern(solid) || pci 0 0.753103713210399 0.25 0.753103713210399, lwidth(medthick) lcolor(red) lpattern(dash)) xline(1.2, lpattern(dash) lcolor(gs12) style(foreground)) xtitle("Implied Delta") title("Distribution of Experts' Parameters - Delta Discounting") note("Red dashed line indicates the parameter value implied by the actual mTurker mean effort." "Green solid line indicates the median parameter value of the experts." "Observations are truncated at 1.2") graphregion(color(white)) bgcolor(white) fcolor(navy) lcolor(gs16) legend(off)
graph export "output\expertParametersPowerDeltaDiscounting.png", width(2000) height(1500) replace

gen lambda = -0.2
replace lambda = expertlambda if expertlambda >-0.2 & expertlambda <= 5
replace lambda = 5.05 if expertlambda > 5
replace lambda = . if expertlambda == .
* drop observations with negative lambda (12 observations)
replace lambda = . if expertlambda < 0
qui sum lambda, d
local median = r(p50)
hist lambda, fraction width(0.2) xlabel(0(1)5 "") start(0) text(-0.006 5.2 ">5", place(s)) /*xline(1.7307)
*/ addplot(pci 0 `median' 0.2 `median', lwidth(medthick) lcolor(green) lpattern(solid) || pci 0 1.7307 0.2 1.7307, lwidth(medthick) lcolor(red) lpattern(dash)) xline(5, lpattern(dash) lcolor(gs12) style(foreground)) xtitle("Implied Lambda") title("Distribution of Experts' Parameters - Loss Aversion") note("Red dashed line indicates the parameter value implied by the actual mTurker mean effort." "Green solid line indicates the median parameter value of the experts." "Observations are truncated at 5. Observations with a difference of less than 10 in the two gain-" "treatments and negative parameter values are dropped (i.e., 58 and 7 observations)") graphregion(color(white)) bgcolor(white) fcolor(navy) lcolor(gs16) legend(off)
graph export "output\expertParametersPowerLossAversion.png", width(2000) height(1500) replace

* warm glow and altruism jointly backed-out
gen warmglow = .
replace warmglow = expertwarmglow if expertwarmglow>-0.01 & expertwarmglow <= 0.5
replace warmglow = 0.52 if expertwarmglow > 0.5
replace warmglow = -0.02 if expertwarmglow < 0
* multiply warmglow times 100 to faciliate interpretation as "The warm-glow corresponds to paying a piece rate of x cents per hundred effort units
*replace warmglow = warmglow*100
qui sum warmglow, d
local median = r(p50)
hist warmglow, fraction width(0.025) xlabel(0(0.1)0.45) start(-0.025) text(-0.0115 0.52 ">0.5", place(s)) text(-0.0115 -0.022 "<0", place(s)) /*xline(0.1408, 
style(foreground))*/ addplot(pci 0 `median' 0.3 `median', lwidth(medthick) lcolor(green) lpattern(solid) || pci 0 0.1250594499496 0.3 0.1250594499496, lwidth(medthick) lcolor(red) lpattern(dash)) xline(0.5, lpattern(dash) lcolor(gs12) style(foreground)) xline(0, lpattern(dash) lcolor(gs12) style(foreground)) xtitle("Implied Warm Glow*100") title("Distribution of Experts' Parameters - Warm Glow") note("Red dashed line indicates the parameter value implied by the actual mTurker mean effort." "Green solid line indicates the median parameter value of the experts." "Parameter values can be interpreted as an intrinsic value" "corresponding to paying an additional piece rate of x per 100 effort units." "Observations are truncated at 0 and 0.5.") graphregion(color(white)) bgcolor(white) fcolor(navy) lcolor(gs16) legend(off)
graph export "output\expertParametersPowerWarmglow2.png", width(2000) height(1500) replace

gen altruism = .
replace altruism = expertaltruism if expertaltruism >= 0 & expertaltruism <= 1.2
replace altruism = -0.02 if expertaltruism < 0
replace altruism = 0.52 if expertaltruism > 0.5
qui sum altruism, d
local median = r(p50)
hist altruism, fraction width(0.025) start(-0.025) xlabel(0(0.1)0.45 "") text(-0.0135 0.52 ">0.5", place(s)) text(-0.0135 -0.022 "<0", place(s)) /*
xline(0.1408, style(foreground))*/ addplot(pci 0 `median' 0.42 `median', lwidth(medthick) lcolor(green) lpattern(solid) || pci 0 0.00298894400368287 0.42 0.00298894400368287, lwidth(medthick) lcolor(red) lpattern(dash)) xline(0.5, lpattern(dash) lcolor(gs12) style(foreground)) xline(0, lpattern(dash) lcolor(gs12) style(foreground)) xtitle("Implied Altruism") title("Distribution of Experts' Parameters - Altruism") note("Red dashed line indicates the parameter value implied by the actual mTurker mean effort." "Green solid line indicates the median parameter value of the experts." "Observations are truncated at 0 and 0.5.") graphregion(color(white)) bgcolor(white) fcolor(navy) lcolor(gs16) legend(off)
graph export "output\expertParametersPowerAltruism2.png", width(2000) height(1500) replace


gen curv1a = curv1
replace curv1a = 0 if curv1<0
replace curv1a = 5 if curv1>5

gen curv88a = curv88
replace curv88a = 0 if curv88<0
replace curv88a = 5 if curv88>5

gen curvEsta = curv44
replace curvEsta = 0 if curv44<0
replace curvEsta = 15 if curv44>15

qui sum curv1, d
local median1 = r(p50)
qui sum curv88, d
local median88 = r(p50)
qui sum curv44, d
local median44 = r(p50)
hist curv1a, fraction width(.2) start(0) xlabel(0(1)5, labsize(vsmall)) ylabel(0(0.1)0.4, labsize(vsmall)) addplot(pci 0 `median1' 0.4 `median1', lwidth(thick) lcolor(green) || pci 0 0.23897314 0.4 0.23897314, lwidth(thick) lcolor(red) lpattern(dash)) xtitle("Implied 1% Probability Weighting (in %) - Curvature of 1.00 (assumed)",size(small)) ytitle("") graphregion(color(white)) bgcolor(white) fcolor(navy) lcolor(gs16) legend(off) name(curv1hist)
hist curv88a, fraction width(.2) start(0) xlabel(0(1)5, labsize(vsmall)) ylabel(0(0.1)0.4, labsize(vsmall)) addplot(pci 0 `median88' 0.4 `median88', lwidth(thick) lcolor(green) || pci 0 0.46599932 0.4 0.46599932, lwidth(thick) lcolor(red) lpattern(dash)) xtitle("Implied 1% Probability Weighting (in %) - Curvature of 0.88 (assumed)",size(small)) ytitle("") graphregion(color(white)) bgcolor(white) fcolor(navy) lcolor(gs16) legend(off) name(curv88hist)
hist curvEsta, fraction width(.5) start(0) xlabel(0(1)15, labsize(vsmall)) ylabel(0(0.1)0.4, labsize(vsmall)) addplot(pci 0 `median44' 0.4 `median44', lwidth(thick) lcolor(green) || pci 0 4.29606134 0.4 4.29606134, lwidth(thick) lcolor(red) lpattern(dash)) xtitle("Implied 1% Probability Weighting (in %) - Curvature of 0.47 (estimated)",size(small)) ytitle("") graphregion(color(white)) bgcolor(white) fcolor(navy) lcolor(gs16) legend(off) name(curvEsthist)

graph combine curv1hist curv88hist curvEsthist, col(1) graphregion(color(white)) note("Implied probability weighting under three different specifications of the curvature" "over the piece rate (1, 0.88, 0.47). The red line indicates the actual probability" "weighting from the Mturkers' effort, the green line denotes the median forecast of" "the experts. Observations are truncated at 5, 5, and 15, respectively.", size(vsmall))
graph export "output\expNLLS_ImpliedProbWeightCombined.png", width(1600) height(2000) replace

graph combine curv1hist curv88hist curvEsthist, col(1) graphregion(color(white))
graph export "output\expNLLS_ImpliedProbWeightCombinedNoNotes.png", width(1600) height(2000) replace

