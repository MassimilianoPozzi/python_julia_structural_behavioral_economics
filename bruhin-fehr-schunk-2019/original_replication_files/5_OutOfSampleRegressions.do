* If you use this code or parts of it, please cite the following paper
* Bruhin A., E. Fehr, D. Schunk (2018): "The Many Faces of Human Sociality: Uncovering the Distribution and Stability of Social Preferences", 
* Journal of the European Economic Association, forthcoming.
* ===========================================================================================================================================


clear
cls

* ========= Trust Games ========= 
insheet using "OutOfSample_trustgames.csv", comma
 
gen diff = tworthy_indivpred - tworthy_typepred

reg trustworthy                                     bf_* cogabil pe_d1_* pe_monthinc pe_age pe_female, cl(sid)
reg trustworthy                   tworthy_typepred  bf_* cogabil pe_d1_* pe_monthinc pe_age pe_female, cl(sid)
reg trustworthy tworthy_indivpred                   bf_* cogabil pe_d1_* pe_monthinc pe_age pe_female, cl(sid)
reg trustworthy tworthy_typepred  diff              bf_* cogabil pe_d1_* pe_monthinc pe_age pe_female, cl(sid)

clear
* ========= Reward Punishment Games ========= 
insheet using "OutOfSample_rpgames.csv", comma

gen diff = reward_punish_typepred - reward_punish_indivpred

reg reward_punish                                                 bf_* cogabil pe_d1_* pe_monthinc pe_age pe_female , cl(sid)
reg reward_punish                         reward_punish_typepred  bf_* cogabil pe_d1_* pe_monthinc pe_age pe_female , cl(sid)
reg reward_punish reward_punish_indivpred                         bf_* cogabil pe_d1_* pe_monthinc pe_age pe_female , cl(sid)
reg reward_punish reward_punish_typepred  diff                    bf_* cogabil pe_d1_* pe_monthinc pe_age pe_female , cl(sid)
