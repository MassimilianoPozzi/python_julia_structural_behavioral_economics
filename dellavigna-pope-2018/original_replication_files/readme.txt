This folder contains the files necessary to replicate the tables and figures in “What Motivates Effort? Evidence and Expert Forecasts” by Stefano DellaVigna and Devin Pope.

Raw Data:
The input folder contains the three raw data files needed for the analysis. 
ExpertForecast.dta contains forecasts from academic experts. 
MTurkData.dta includes effort of all MTurkers who participated in the task. 
MTurkDataDemo.dta is a separate file containing demographic information on the MTurkers who completed the task. We intentionally de-link MTurker demographics from their effort. 

Data Cleaning:
Clean Forecast Data.do and Clean MTurk Data.do create cleaned data files (saved in “dtafiles”) that we use in the analyses. These need to be run first to generate the data used in subsequent analysis files. 

Analysis of MTurk Data and Forecasts: 
0_Analyze MTurk Results.do creates the reduced form figures and tables. 

1_NLS_main.do, 2_Tables_5_6_A4_figs_2_a3_a9.R, 3_Table5_panelA_Col4_LowPaySimulations.m, and 4_ExpertParameters_Histograms_figA9.do generate parameter estimates by NLS and minimum distance and generate the corresponding tables and charts. 

The output from these analyses are saved in the output folder. In some cases, figures or table calculations are done in Excel, using corresponding data saved in the output folder. All Excel calculations are all saved in the file BehavioralTablesFigures.xlsx.

The MetaAnalysisCode folder contains code for figures 8 and online appendix figure 7 in the paper, as well as for some calculations in online appendix table 3. ExpertLitPlot_Fig8_AppFig7.R is an R script that reads in input data (ExpertsCohens.csv) based on extra columns of the posted version of online appendix table 2, and produces figure 8 and online appendix figure 7 of the paper. PNG and EPS versions of these plots are also included in the folder. OnlineAppTable3Calculations.R contains calculations underlying the implied effort values under the average probability weights for p=0.01 and p=0.5 at the bottom of online appendix table 3.
