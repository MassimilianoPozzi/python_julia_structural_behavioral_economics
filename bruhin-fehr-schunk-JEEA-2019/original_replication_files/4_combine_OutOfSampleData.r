# If you use this code or parts of it, please cite the following paper
# Bruhin A., E. Fehr, D. Schunk (2018): "The Many Faces of Human Sociality: Uncovering the Distribution and Stability of Social Preferences",
# Journal of the European Economic Association, forthcoming.
# ===========================================================================================================================================

rm(list=ls()) # Clear the memory

# Function to prepare data sets for out of sample predictions
f.combine <- function(choices, subjects, ipe, finmixobj) {

	d1 <- merge(choices, subjects, by="sid", all.x=T, all.y=F)
	d2 <- merge(d1, ipe, by="sid", all.x=T, all.y=F)

	# ex-post probabilities of type membershipt => classification of subjects into types
	tau       <- finmixobj$tau
	typenam   <- paste("type_", 1:ncol(tau), sep="")
	class     <- matrix(0, nrow(tau), ncol(tau)+1, dimnames=list(c(), c("sid", typenam)))
	class[,1] <- as.numeric(rownames(tau))
	for (i in 1:nrow(tau)) {
		class[i, which.max(tau[i,])+1] <- 1
	}

	d3 <- merge(d2, class, by="sid", all.x=T, all.y=F)

	# add type specific parameter estimates
	d3$alpha_type      <- NA
	d3$beta_type       <- NA
	d3$gamma_type      <- NA
	d3$delta_type      <- NA
	d3$choicesens_type <- NA
	pars <- paste(c("alpha","beta","gamma","delta","choicesens"),"_type", sep="")


	for (i in 1:nrow(d3)) {

		typeind     <- which.max(d3[i, typenam])  # which type?
		d3[i, pars] <- finmixobj$pe[2:6, typeind] # add type-specific estimates

	}

	d3

}

# Function to compute the behavioral predictions for the reward punishment games
f.rppredict <- function(d, betaind, alphaind, gammaind, deltaind, sigmaind) {
	pay.r1 <- d[,"payoff_r"] - d[,"c1"]
	pay.r2 <- d[,"payoff_r"] - d[,"c2"]
	pay.r3 <- d[,"payoff_r"] - d[,"c3"]
	pay.r4 <- d[,"payoff_r"] - d[,"c4"]
	pay.r5 <- d[,"payoff_r"] - d[,"c5"]
	pay.r6 <- d[,"payoff_r"] - d[,"c6"]
	pay.r7 <- d[,"payoff_r"]

	pay.p1 <- d[,"payoff_p"] + d[,"r1"]
	pay.p2 <- d[,"payoff_p"] + d[,"r2"]
	pay.p3 <- d[,"payoff_p"] + d[,"r3"]
	pay.p4 <- d[,"payoff_p"] + d[,"r4"]
	pay.p5 <- d[,"payoff_p"] + d[,"r5"]
	pay.p6 <- d[,"payoff_p"] + d[,"r6"]
	pay.p7 <- d[,"payoff_p"]

	r1 <- as.numeric(pay.r1 > pay.p1)
	r2 <- as.numeric(pay.r2 > pay.p2)
	r3 <- as.numeric(pay.r3 > pay.p3)
	r4 <- as.numeric(pay.r4 > pay.p4)
	r5 <- as.numeric(pay.r5 > pay.p5)
	r6 <- as.numeric(pay.r6 > pay.p6)
	r7 <- as.numeric(pay.r7 > pay.p7)

	s1 <- as.numeric(pay.r1 < pay.p1)
	s2 <- as.numeric(pay.r2 < pay.p2)
	s3 <- as.numeric(pay.r3 < pay.p3)
	s4 <- as.numeric(pay.r4 < pay.p4)
	s5 <- as.numeric(pay.r5 < pay.p5)
	s6 <- as.numeric(pay.r6 < pay.p6)
	s7 <- as.numeric(pay.r7 < pay.p7)

	par1 <- d[,betaind]*r1 + d[,alphaind]*s1 + d[,gammaind]*d[,"pgoesleft_kind"] + d[,deltaind]*(1-d[,"pgoesleft_kind"])
	par2 <- d[,betaind]*r2 + d[,alphaind]*s2 + d[,gammaind]*d[,"pgoesleft_kind"] + d[,deltaind]*(1-d[,"pgoesleft_kind"])
	par3 <- d[,betaind]*r3 + d[,alphaind]*s3 + d[,gammaind]*d[,"pgoesleft_kind"] + d[,deltaind]*(1-d[,"pgoesleft_kind"])
	par4 <- d[,betaind]*r4 + d[,alphaind]*s4 + d[,gammaind]*d[,"pgoesleft_kind"] + d[,deltaind]*(1-d[,"pgoesleft_kind"])
	par5 <- d[,betaind]*r5 + d[,alphaind]*s5 + d[,gammaind]*d[,"pgoesleft_kind"] + d[,deltaind]*(1-d[,"pgoesleft_kind"])
	par6 <- d[,betaind]*r6 + d[,alphaind]*s6 + d[,gammaind]*d[,"pgoesleft_kind"] + d[,deltaind]*(1-d[,"pgoesleft_kind"])
	par7 <- d[,betaind]*r7 + d[,alphaind]*s7 + d[,gammaind]*d[,"pgoesleft_kind"] + d[,deltaind]*(1-d[,"pgoesleft_kind"])

	umat <- matrix(NA, nrow(d), 7)
	umat[,1] <- (1-par1)*pay.r1 + par1*pay.p1
	umat[,2] <- (1-par2)*pay.r2 + par2*pay.p2
	umat[,3] <- (1-par3)*pay.r3 + par3*pay.p3
	umat[,4] <- (1-par4)*pay.r4 + par4*pay.p4
	umat[,5] <- (1-par5)*pay.r5 + par5*pay.p5
	umat[,6] <- (1-par6)*pay.r6 + par6*pay.p6
	umat[,7] <- (1-par7)*pay.r7 + par7*pay.p7

	umaxind   <- max.col(umat)
	pred.hard <- d[,"r1"]*(umaxind == 1) + d[,"r2"]*(umaxind == 2) + d[,"r3"]*(umaxind == 3) + d[,"r4"]*(umaxind == 4) + d[,"r5"]*(umaxind == 5) + d[,"r6"]*(umaxind == 6)

	expu  <- exp(d[,sigmaind]*umat)
	probs <- expu / matrix(rowSums(expu), nrow(expu), ncol(expu))
	pred.prob <- d[,"r1"]*probs[,1] + d[,"r2"]*probs[,2] + d[,"r3"]*probs[,3] + d[,"r4"]*probs[,4] + d[,"r5"]*probs[,5] + d[,"r6"]*probs[,6]
	pred.prob[is.na(pred.prob)] <- pred.hard[is.na(pred.prob)]
	pred.prob
}

# Function to compute the behavioral predictions for the trust games
f.trustpredict <- function(d, betaind, alphaind, gammaind, sigmaind) {
	trustX <- 1200-d[,"cost"]  #trustee's payoff when trustworthy
	trustR <- as.numeric(trustX > 900) #indicator for advantageous inequality
	trustS <- as.numeric(trustX < 900) #indicator for disadvantageous inequality

	umat <- matrix(NA, nrow(d), 2)
	umat[,1] <- (1 - d[,betaind]*trustR - d[,alphaind]*trustS - d[,gammaind])*trustX + (d[,betaind]*trustR + d[,alphaind]*trustS + d[,gammaind])*900
	umat[,2] <- (1 - d[,betaind] - d[,gammaind])*1200

	umaxind   <- max.col(umat)
	pred.hard <- as.numeric(umaxind == 1)

	expu  <- exp(d[,sigmaind]*umat)
	probs <- expu / matrix(rowSums(expu), nrow(expu), ncol(expu))
	pred.prob <- probs[,1]
	pred.prob[is.na(pred.prob)] <- pred.hard[is.na(pred.prob)]
	pred.prob

}

# Load the data that needs to be combined
subjects   <- read.table("subjects.csv",            sep=",", head=T)
trustgames <- read.table("choices_trustgames.csv",  sep=",", head=T)
rpgames    <- read.table("choices_rpgames.csv",     sep=",", head=T)
ipe        <- read.table("individual_est_exp2.csv", sep=",", head=T)

load("finmix_est_3cl.RData") # load the finite mixture estimates

# Prepare the data sets
trg <- f.combine(trustgames, subjects, ipe, pe2)
rpg <- f.combine(rpgames, subjects, ipe, pe2)

# Add the predictions based on the individual- and type-specific estimates (see Section 4.5)
trg$tworthy_indivpred <- f.trustpredict(trg, "beta", "alpha", "gamma", "choicesens")
trg$tworthy_typepred  <- f.trustpredict(trg, "beta_type", "alpha_type", "gamma_type", "choicesens_type")

rpg$reward_punish_indivpred <- f.rppredict(rpg, "beta", "alpha", "gamma", "delta", "choicesens")
rpg$reward_punish_typepred  <- f.rppredict(rpg, "beta_type", "alpha_type", "gamma_type", "delta_type", "choicesens_type")

# Wirte data to files
write.table(trg, "OutOfSample_trustgames.csv", sep=",", col.names=T, row.names=F)
write.table(rpg, "OutOfSample_rpgames.csv",    sep=",", col.names=T, row.names=F)
