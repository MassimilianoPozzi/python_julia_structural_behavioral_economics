# If you use this code or parts of it, please cite the following paper
# Bruhin A., E. Fehr, D. Schunk (2018): "The Many Faces of Human Sociality: Uncovering the Distribution and Stability of Social Preferences",
# Journal of the European Economic Association, forthcoming.
# ===========================================================================================================================================

library(compiler)

# Log Likelihood
fdummy.ll <- function(v,y,self_x, other_x, self_y, other_y,indicators_x, indicators_y) {
	beta  <- v[1:(length(v)-1)]
	gamma <- exp(v[length(v)])

	lli    <- indicators_x%*%beta
	rli    <- indicators_y%*%beta
	uleft  <- gamma*((1-lli)*self_x+lli*other_x)
	uright <- gamma*((1-rli)*self_y+rli*other_y)

	probs  <- (exp(uleft)/(exp(uleft)+exp(uright)))^y * (exp(uright)/(exp(uleft)+exp(uright)))^(1-y)

	sum(log(probs))
}
f.ll <- cmpfun(fdummy.ll)

# Log Likelihood of a model without indicators
fdummy.ll.nopar <- function(v,y,self_x,self_y) {
	gamma <- exp(v)

	uleft  <- gamma*self_x
	uright <- gamma*self_y

	probs  <- (exp(uleft)/(exp(uleft)+exp(uright)))^y * (exp(uright)/(exp(uleft)+exp(uright)))^(1-y)

	sum(log(probs))
}
f.ll.nopar <- cmpfun(fdummy.ll.nopar)

# Gradient of Log Likelihood
fdummy.ll.gr <- function(v, y, self_x, other_x, self_y, other_y, indicators_x, indicators_y) {
	beta  <- v[1:(length(v)-1)]
	gamma <- exp(v[length(v)])

	lli <- indicators_x%*%beta
	rli <- indicators_y%*%beta
	utl <- gamma*((1-lli)*self_x+lli*other_x)
	utr <- gamma*((1-rli)*self_y+rli*other_y)
	probs  <- (exp(utl)/(exp(utl)+exp(utr)))^y * (exp(utr)/(exp(utl)+exp(utr)))^(1-y)
	probsm <- matrix(probs, nrow(indicators_x), ncol(indicators_x))

	u   <- matrix(exp(utl)            , nrow(indicators_x), ncol(indicators_x))
	w   <- matrix(exp(utl) + exp(utr) , nrow(indicators_x), ncol(indicators_x))
	up  <- gamma * indicators_x * matrix(other_x-self_x, nrow(indicators_x), ncol(indicators_x)) * u
	wp  <- up + gamma * indicators_y * matrix((other_y-self_y)*exp(utr), nrow(indicators_x), ncol(indicators_x))

	prebgrad <- (-1)^(1-y)*((up*w - u*wp)/w^2)/probsm
	bgrad    <- colSums(prebgrad)

	up2   <- utl * u[,1]
	wp2   <- up2 + utr * exp(utr)

	ggrad <- sum((-1)^(1-y) * (((up2*w[,1] - u[,1]*wp2)/w[,1]^2)/probs))

	c(bgrad, ggrad)
}
f.ll.gr <- cmpfun(fdummy.ll.gr)

# Gradient of Log Likelihood of a model without indicators
fdummy.ll.gr.nopar <- function(v, y, self_x,self_y) {
	gamma <- exp(v)

	utl <- gamma*self_x
	utr <- gamma*self_y
	probs  <- (exp(utl)/(exp(utl)+exp(utr)))^y * (exp(utr)/(exp(utl)+exp(utr)))^(1-y)

	u   <- exp(utl)
	w   <- exp(utl) + exp(utr)
	up2   <- utl * u
	wp2   <- up2 + utr * exp(utr)

	ggrad <- sum((-1)^(1-y) * (((up2*w - u*wp2)/w^2)/probs))
	ggrad
}
f.ll.gr.nopar <- cmpfun(fdummy.ll.gr.nopar)

# Function to compute cluster robust standard errors
f.clusterse <- function(estobj, dat, label_choice_x, label_self_x, label_other_x, label_self_y, label_other_y, label_indicators_x, label_indicators_y, clustervar) {
	y         <- as.vector(dat[,label_choice_x])
	self_x    <- as.matrix(dat[,label_self_x])
	other_x   <- as.matrix(dat[,label_other_x])
	self_y    <- as.matrix(dat[,label_self_y])
	other_y   <- as.matrix(dat[,label_other_y])
	indicators_x <- as.matrix(dat[,label_indicators_x])
	indicators_y <- as.matrix(dat[,label_indicators_y])
	clusters  <- as.vector(dat[,clustervar])

	ll1 <- f.ll(estobj$retobj$par, y, self_x, other_x, self_y, other_y, indicators_x, indicators_y)
	if (abs(ll1 - estobj$retobj$value)>0.001) stop(paste("Saved loglik", ll1, "and computed loglik", estobj$retobj$value, "differ"))

	sandwich <- f.sandwich(estobj$retobj$par, y, self_x, other_x, self_y, other_y,indicators_x, indicators_y, clusters) # sandwich part; see function below

	k <- length(estobj$retobj$par)

	vcm <- solve(-estobj$retobj$hessian)
	vcm <- sandwich$qc * (vcm %*% sandwich$sandwich %*% vcm)

	se <- sqrt(diag(vcm))
	se[k] <- sqrt( exp(estobj$retobj$par[k])^2 * se[k]^2 )

	pe <- estobj$retobj$par
	pe[k] <- exp(pe[k])
	out <- matrix(c(pe, se, pe/se, 2*(1-pnorm(abs(pe)/se))), k, 4, dimnames=list(c(paste(label_indicators_x,"label_idicators_y"),"choicesens"),c("point_estimates","se_cluster","z-stat_cluster","p-val_cluster")));

	cat("\nNumber of clusters", sandwich$j,"\n")
	print(round(out,3))

	list(out=out, vcm=vcm[1:(k-1), 1:(k-1)])
}

# Sandwitch part for calculating cluster robust standard errors
f.sandwich <- function(v, y, self_x, other_x, self_y, other_y,indicators_x, indicators_y, clusters) {

	#compute individual gradient contributions
	beta  <- v[1:(length(v)-1)]
	gamma <- exp(v[length(v)])

	lli <- indicators_x%*%beta;
	rli <- indicators_y%*%beta;
	utl <- gamma*((1-lli)*self_x+lli*other_x)
	utr <- gamma*((1-rli)*self_y+rli*other_y)
	probs  <- (exp(utl)/(exp(utl)+exp(utr)))^y * (exp(utr)/(exp(utl)+exp(utr)))^(1-y);
	probsm <- matrix(probs, nrow(indicators_x), ncol(indicators_x))

	u   <- matrix(exp(utl)            , nrow(indicators_x), ncol(indicators_x))
	w   <- matrix(exp(utl) + exp(utr) , nrow(indicators_x), ncol(indicators_x))
	up  <- gamma * indicators_x * matrix(other_x-self_x, nrow(indicators_x), ncol(indicators_x)) * u
	wp  <- up + gamma * indicators_y * matrix((other_y-self_y)*exp(utr), nrow(indicators_x), ncol(indicators_x))

	bgradi <- (-1)^(1-y)*((up*w - u*wp)/w^2)/probsm

	up2   <- utl * u[,1]
	wp2   <- up2 + utr * exp(utr)

	ggradi <- (-1)^(1-y) * (((up2*w[,1] - u[,1]*wp2)/w[,1]^2)/probs)

	gradi <- cbind(bgradi, ggradi)

	#Sandwich part of the sandwich VC estimator
	cl <- unique(clusters)
	j  <- length(cl)
	k  <- ncol(gradi)

	sandwich <- matrix(0, k, k)
	for (i in 1:j) {
		sel <- clusters == cl[i]
		gradsel <- colSums(gradi[sel, ])
		sandwich <- sandwich + (gradsel %o% gradsel)
	}

	qc <- (length(y)-1) / (length(y)-k) * j/(j-1) #df adjustment

	list(qc=qc, sandwich=sandwich, j=j)
}

# Estimation Function (to be called directly for aggregate estimations and by f.indivest for individual estimations)
fdummy.estim <- function(dat, label_choice_x, label_indicators_x, label_indicators_y, label_self_x, label_other_x, label_self_y, label_other_y, silent=F, comp.se=T, v0=NULL) {
	y         <- as.vector(dat[,label_choice_x])
	self_x    <- as.matrix(dat[,label_self_x])
	other_x   <- as.matrix(dat[,label_other_x])
	self_y    <- as.matrix(dat[,label_self_y])
	other_y   <- as.matrix(dat[,label_other_y])
	indicators_x <- as.matrix(dat[,label_indicators_x])
	indicators_y <- as.matrix(dat[,label_indicators_y])

	k <- length(label_indicators_x)+1

	if (is.null(v0)) {v0  <- c(rep(0,k-1), log(runif(1,.05,.8)/mean(c(mean(self_x),mean(other_x),mean(self_y),mean(other_y)))))} # Random start values if start values are not provided
	ret <- optim(v0, f.ll, f.ll.gr, method="BFGS", control=list(maxit=100000,trace=as.numeric(silent==F),REPORT=50, fnscale=-1), hessian=comp.se, y=y, self_x=self_x, other_x=other_x, self_y=self_y, other_y=other_y, indicators_x=indicators_x, indicators_y=indicators_y)

	rawpe <- ret$par
	pe    <- rawpe
	pe[k] <- exp(rawpe[k]) # Choice sensitivity > 0

	se <- rep(NA, k)
	overall.pval <- NA
	if(comp.se) {
		rawse <- sqrt(diag( solve(-ret$hessian) ))
		se    <- rawse
		se[k] <- sqrt( exp(rawpe[k])^2 * rawse[k]^2 ) #Delta method

		ret2 <- optim(log(0.5/1000), f.ll.nopar, f.ll.gr.nopar, method="BFGS", control=list(maxit=100000,trace=F,REPORT=50, fnscale=-1), hessian=F, y=y, self_x=self_x, self_y=self_y)
		overall.pval <- 1-pchisq(2*(ret$value-ret2$value),4)
	}

	out <- matrix(c(pe, se, pe/se, 2*(1-pnorm(abs(pe)/se))), k, 4, dimnames=list(c(paste(label_indicators_x,label_indicators_y),"choicesens"),c("point_estimates","se_noCluster","z-stat_noCluster","p-val_noCluster")))

	if (silent==F) {
		print(paste("LogLik:",ret$value), quote=F)
		print(round(out,4))
	}

	list(out=out, retobj=ret, overall.pval=overall.pval)
}
f.estim <- cmpfun(fdummy.estim)


# Function for running individual estimations
f.indivest <- function(dat){

	# Set up grid of start values (individual estimations are more fragile than pooled ones => test a lot of start values for each subject)
	v01   <- c(0,0,0,0)
	v02   <- c(.1,-.1,0,0)
	v03   <- c(.1,-.1,.01,-.01)
	v04   <- c(.1,-.1,-.01,.01)
	v05   <- c(-.1,.1,0,0)
	v06   <- c(-.1,.1,.01,-.01)
	v07   <- c(-.1,.1,-.01,.01)
	v00   <- matrix(c(v01,v02,v03,v04,v05,v06,v07),7,4,byrow=T)
	errv0 <- log(seq(0.0001, .001, length.out=20))
	irps  <- 7*20 #140 estimation per subject with varying start values

	ids <- sort(unique(dat$sid)) #Individual IDs
	j   <- length(ids)           #Number of individuals

	parset <- c("beta", "alpha", "gamma", "delta", "choicesens")
	parmatrix    <- matrix(NA, j,    length(parset)+2, dimnames=list(ids, c(parset, "loglik", "success"))) # Matrix for saving the final estimation results for every subject
	parsubmatrix <- matrix(NA, irps, length(parset)+2, dimnames=list(NULL,c(parset, "loglik", "success"))) # Matrix for saving the estimation results for every start value combination in the grid of start values

	# Run the individual Estimations including all irps subreplications
	for (i in 1:j) {
		idat <- dat[dat$sid == ids[i], ]
		cat("\nEstimating individual:", i,"/",j)

		# Cycle through all 140 start value combinations
		for (s in 1:irps) {
			verr <- errv0[((s-1) %% 20) + 1] #start values for the choice sensitivity
			v    <- v00[((s-1) %/% 20) +1,] #start values of the parameters

			#estimtion
			pe <- f.estim(idat, "choice_x", c("r_x","s_x","q","v"), c("r_y","s_y","q","v"), "self_x", "other_x", "self_y", "other_y", silent=T, comp.se=F, v0=c(v,verr))
			ipar <- pe$out[,1]
			ill  <- pe$retobj$value
			isuc <- as.numeric(pe$retobj$convergence == 0)
			parsubmatrix[s,] <- c(ipar, ill, isuc)
		}

		# Extract start values that led to the best estimation (i.e. the estimation with the highest log likelihood)
		v0 <- as.matrix(parsubmatrix[order(parsubmatrix[,"loglik"], decreasing=T),])[1,1:length(parset)]
		v0[length(v0)] <- log(v0[length(v0)]) # Choice sensitivty (last parameter in v0) is > 0

		# Redo best estimation
		pe <- f.estim(idat, "choice_x", c("r_x","s_x","q","v"), c("r_y","s_y","q","v"), "self_x", "other_x", "self_y", "other_y", silent=T, comp.se=F, v0=v0)

		parmatrix[i, 1:length(parset)] <- pe$out[,1] # save point estimates
		parmatrix[i, length(parset)+1] <- pe$retobj$value # save log likelihood
		parmatrix[i, length(parset)+2] <- as.numeric(pe$retobj$convergence == 0) # save success indicator
	}

	parmatrix
}

# Funtion that identifies subjects that behaved erraticly and exhibit individual parameter estimates outside the identifiable range
f.droperratic <- function(par1, par2) {

	range <- c(-3,1) # Identifiable parameter range

	drop <- matrix(NA, nrow(par1), 4, dimnames=list(rownames(par1),c()))
	for (i in 1:4) {
		drop[,i] <- (par1[,i] < range[1]) | (par1[,i] > range[2]) | (par2[,i] < range[1]) | (par2[,i] > range[2])
	}

	rowSums(drop)
}
