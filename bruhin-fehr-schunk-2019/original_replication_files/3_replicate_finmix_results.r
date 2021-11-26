# If you use this code or parts of it, please cite the following paper
# Bruhin A., E. Fehr, D. Schunk (2018): "The Many Faces of Human Sociality: Uncovering the Distribution and Stability of Social Preferences",
# Journal of the European Economic Association, forthcoming.
# ===========================================================================================================================================

rm(list=ls()) # Clear the memory
library(compiler)
set.seed(51)

# Data-Preparation-Function
f.dataprep <- function(dat, persid, label_choice_x, label_indicators_x, label_indicators_y, label_self_x, label_other_x, label_self_y, label_other_y, nafiller=1) {
	tid <- table(dat[,persid])    #Table of individuals
	pid <- as.numeric(names(tid)) #Vector of individual identifiers
	j   <- length(tid)            #Number of individuals
	n   <- max(tid)               #Maximum number of choices per individual
	lk  <- length(label_indicators_x)        #Number of "left" regressors
	rk  <- length(label_indicators_y)        #Number of "right" regressors
	{if (lk != rk) {stop("Length of label_indicators_x and label_indicators_y differ!")}}
	k   <- lk                     #Number of regressors

	indicators_x  <- matrix(1, n*j, k)  #"left" regressor matrix
	indicators_y <- matrix(1, n*j, k)  #"right" regressor matrix

	nlist      <- list(c(), pid)
	self_x     <- matrix(0, n, j, dimnames=nlist) # Payoff matrix, player A, "X"
	other_x    <- matrix(0, n, j, dimnames=nlist) # Payoff matrix, player B, "X"
	self_y     <- matrix(0, n, j, dimnames=nlist) # Payoff matrix, player A, "Y"
	other_y    <- matrix(0, n, j, dimnames=nlist) # Payoff matrix, player B, "Y"
	choice_x <- matrix(0, n, j, dimnames=nlist) # Matrix indicating the choice (1=left, 0=right)
	sel        <- matrix(0, n, j, dimnames=nlist) # Matrix indicating whether an observation is present (panel may be unbalanced!)
	rm(nlist)

	#Fill the data and selection matrices individual by individual (could handle an unbalanced panel)
	for (i in 1:j) {
		sel[,i]        <- c(rep(1, tid[i]), rep(0, n-tid[i]))

		idx            <- dat[,persid] == pid[i]
		self_x[,i]     <- c(dat[idx, label_self_x],   rep(nafiller, n-tid[i]))
		other_x[,i]    <- c(dat[idx, label_other_x],  rep(nafiller, n-tid[i]))
		self_y[,i]     <- c(dat[idx, label_self_y],   rep(nafiller, n-tid[i]))
		other_y[,i]    <- c(dat[idx, label_other_y],  rep(nafiller, n-tid[i]))
		choice_x[,i] <- c(dat[idx, label_choice_x], rep(nafiller, n-tid[i]))

		{if (k > 1) {
			m1 <- as.matrix(dat[idx, label_indicators_x])
			m2 <- matrix(0, n-nrow(m1), ncol(m1))
			m  <- rbind(m1, m2)
			indicators_x[(n*(i-1)+1):(i*n),] <- m
			rm (m1,m2)

			m1 <- as.matrix(dat[idx, label_indicators_y])
			m2 <- matrix(0, n-nrow(m1), ncol(m1))
			m  <- rbind(m1, m2)
			indicators_y[(n*(i-1)+1):(i*n),] <- m
			rm(m1,m2)
		}}
	}
	list(tid=tid, pid=pid, n=n, j=j, k=k, indicators_x=indicators_x, indicators_y=indicators_y, self_x=self_x, other_x=other_x, self_y=self_y, other_y=other_y, choice_x=choice_x, sel=sel)
}
f.dataprep <- cmpfun(f.dataprep)

# Component Density Function
f.compdens <- function(theta.v, choice_x, self_x, other_x, self_y, other_y, indicators_x, indicators_y, sel, nc ,k, n, j) {
	theta.m <- matrix(theta.v, k+1, nc)

	linindex_x <- indicators_x  %*% theta.m[1:k,] # Linear index of X-choice, a (n*j x nc)-matrix
	linindex_y <- indicators_y %*% theta.m[1:k,]  # Linear index of Y-choice, a (n*j x nc)-matrix

	ret     <- matrix(0, j, nc)

	# Loop over all nc preference types
	for (c in 1:nc) {
		li_x   <- matrix(linindex_x[,c], n, j)
		li_y   <- matrix(linindex_y[,c], n, j)
		choicesens <- theta.m[k+1,c]
		cs_li_x   <- choicesens * ((1-li_x)*self_x + li_x*other_x)  #Type-specific deterministic utility of chosing "X" times choice sensitivity
		cs_li_y   <- choicesens * ((1-li_y)*self_y + li_y*other_y)  #Type-specific deterministic utility of chosing "Y" times choice sensitivity

		probs <- ( (exp(cs_li_x)/(exp(cs_li_x)+exp(cs_li_y)))^choice_x * (exp(cs_li_y)/(exp(cs_li_x)+exp(cs_li_y)))^(1-choice_x) )^sel

		for (i in 1:j) {
			ret[i,c] <- prod(probs[,i])
		}
	}
	ret
}
f.compdens <- cmpfun(f.compdens)

# Log Likelihood of the Mixture Model and its Gradient
f.totloglik <- function(psi, choice_x, self_x, other_x, self_y, other_y, indicators_x, indicators_y, sel, nc, k, n, j) {
	lpsi    <- length(psi)
	theta.v <- psi[1:(lpsi-nc+1)]
	mix     <- psi[(lpsi-nc+2):lpsi]
	mix     <- c(mix, 1-sum(mix))
	m       <- matrix(mix, j, nc, byrow=T)

	sum(log(rowSums( m*f.compdens(theta.v, choice_x, self_x, other_x, self_y, other_y, indicators_x, indicators_y, sel, nc, k, n, j))))
}
f.totloglik <- cmpfun(f.totloglik)

f.totloglik.gradient <- function(psi, choice_x, self_x, other_x, self_y, other_y, indicators_x, indicators_y, sel, nc, k, n, j) {
	lpsi    <- length(psi)
	theta.v <- psi[1:(lpsi-nc+1)]
	theta.m <- matrix(theta.v, k+1, nc)
	mix     <- psi[(lpsi-nc+2):lpsi]
	mix     <- c(mix, 1-sum(mix))
	m       <- matrix(mix, j, nc, byrow=T)

	nj           <- n*j
	choice_x.k   <- matrix(as.vector(choice_x), nj, k)
	sel.k        <- matrix(as.vector(sel)     , nj, k)
	self_x.k     <- matrix(as.vector(self_x)  , nj, k)
	other_x.k    <- matrix(as.vector(other_x) , nj, k)
	self_y.k     <- matrix(as.vector(self_y)  , nj, k)
	other_y.k    <- matrix(as.vector(other_y) , nj, k)

	cdens.m <- f.compdens(theta.v, choice_x, self_x, other_x, self_y, other_y, indicators_x, indicators_y, sel, nc ,k, n, j)
	denom.outerderiv   <- 1/matrix(rowSums(m*cdens.m), n, j, byrow=T)
	denom.outerderiv.k <- matrix(as.vector(denom.outerderiv), nj, k)

	ret.theta <- rep(0, nc*(k+1))
	ret.mix   <- rep(0, nc)

	for (c in 1:nc) {
		choicesens <- theta.m[k+1, c]
		beta  <- theta.m[1:k, c]

		li_x   <- matrix(as.vector(indicators_x %*% beta), n, j) # linear index X-choice
		li_y   <- matrix(as.vector(indicators_y %*% beta), n, j) # linear index Y-choice

		ul    <- (1-li_x)*self_x + li_x*other_x # Deterministic utility from X-choice
		ur    <- (1-li_y)*self_y + li_y*other_y # Deterministic utility from Y-choice
		ul.k  <- matrix(as.vector(ul), nj, k)
		ur.k  <- matrix(as.vector(ur), nj, k)

		innerderiv   <- mix[c]*matrix(cdens.m[,c],n,j,byrow=T)/((exp(choicesens*ul)/(exp(choicesens*ul)+exp(choicesens*ur)))^choice_x * (exp(choicesens*ur)/(exp(choicesens*ul)+exp(choicesens*ur)))^(1-choice_x))
		innerderiv.k <- matrix(as.vector(innerderiv), nj, k)

		choicesens.grad <- sel * denom.outerderiv * innerderiv * ((1/(exp(choicesens*ul)+exp(choicesens*ur))^2) * ( (exp(choicesens*ul)/(exp(choicesens*ul)+exp(choicesens*ur)))^choice_x * exp(choicesens*(ul + 2*ur))*(-ur + ul) * (exp(choicesens*ur)/(exp(choicesens*ul)+exp(choicesens*ur)))^(-choice_x) * (choice_x * exp(-choicesens*ul) - exp(-choicesens*ur) + choice_x * exp(-choicesens*ur)) ))
		choicesens.grad <- sum(choicesens.grad)

		beta.grad <- sel.k * denom.outerderiv.k * innerderiv.k * ((-1/(exp(choicesens*ul.k)+exp(choicesens*ur.k))^2) * ( (exp(choicesens*ul.k)/(exp(choicesens*ul.k)+exp(choicesens*ur.k)))^choice_x.k * choicesens*exp(choicesens*(ul.k + 2*ur.k))*(indicators_x*(self_x.k-other_x.k) + indicators_y*(other_y.k-self_y.k)) * (exp(choicesens*ur.k)/(exp(choicesens*ul.k)+exp(choicesens*ur.k)))^(-choice_x.k) * (choice_x.k * exp(-choicesens*ul.k) - exp(-choicesens*ur.k) + choice_x.k * exp(-choicesens*ur.k)) ))
		beta.grad <- colSums(beta.grad)

		mix.grad  <- sum(cdens.m[,c] * denom.outerderiv[1,])

		ret.theta[((c-1)*(k+1)+1):(c*(k+1))] <- c(beta.grad, choicesens.grad)
		ret.mix[c]                           <- mix.grad
	}
	ret.mix <- ret.mix - ret.mix[nc]
	ret     <- c(ret.theta, ret.mix[1:(nc-1)])

	ret
}
f.totloglik.gradient <- cmpfun(f.totloglik.gradient)

# Log Likelihood Part in the M-Step and its Gradient
f.mloglik <- function(theta.v, choice_x, self_x, other_x, self_y, other_y, indicators_x, indicators_y, sel, nc, k, n, j, tau) {
	sum(tau * log(f.compdens(theta.v, choice_x, self_x, other_x, self_y, other_y, indicators_x, indicators_y, sel, nc, k, n, j)))
}
f.mloglik <- cmpfun(f.mloglik)

f.mloglik.pregrad <- function(theta.c, weight, choice_x, self_x, other_x, self_y, other_y, indicators_x, indicators_y, sel, nc, k, n, j) {
	choicesens <- theta.c[k+1]
	beta       <- theta.c[1:k]

	weights <- matrix(weight, n, j, byrow=T) #The individual weights (used for the taus)

	li_x   <- matrix(as.vector(indicators_x %*% beta), n, j) # linear index X-choice
	li_y   <- matrix(as.vector(indicators_y %*% beta), n, j) # linear index Y-choice

	ul    <- (1-li_x)*self_x + li_x*other_x # Deterministic utility from X-choice
	ur    <- (1-li_y)*self_y + li_y*other_y # Deterministic utility from Y-choice

	denom <- 1/((exp(choicesens*ul)/(exp(choicesens*ul)+exp(choicesens*ur)))^choice_x * (exp(choicesens*ur)/(exp(choicesens*ul)+exp(choicesens*ur)))^(1-choice_x))

	#The gradient of choicesens (scalar)
	choicesens.grad <- sel * weights * denom * ((1/(exp(choicesens*ul)+exp(choicesens*ur))^2) * ( (exp(choicesens*ul)/(exp(choicesens*ul)+exp(choicesens*ur)))^choice_x * exp(choicesens*(ul + 2*ur))*(-ur + ul) * (exp(choicesens*ur)/(exp(choicesens*ul)+exp(choicesens*ur)))^(-choice_x) * (choice_x * exp(-choicesens*ul) - exp(-choicesens*ur) + choice_x * exp(-choicesens*ur)) ))
	choicesens.grad <- sum(choicesens.grad)

	#The gradient (k x 1)-vector of beta
	nj           <- n*j
	choice_x.k <- matrix(as.vector(choice_x), nj, k)
	denom.k      <- matrix(as.vector(denom)     , nj, k)
	weights.k    <- matrix(as.vector(weights)   , nj, k)
	ul.k         <- matrix(as.vector(ul)        , nj, k)
	ur.k         <- matrix(as.vector(ur)        , nj, k)
	sel.k        <- matrix(as.vector(sel)       , nj, k)
	self_x.k   <- matrix(as.vector(self_x)  , nj, k)
	other_x.k   <- matrix(as.vector(other_x)  , nj, k)
	self_y.k  <- matrix(as.vector(self_y) , nj, k)
	other_y.k  <- matrix(as.vector(other_y) , nj, k)

	beta.grad <- sel.k * weights.k * denom.k * ((-1/(exp(choicesens*ul.k)+exp(choicesens*ur.k))^2) * ( (exp(choicesens*ul.k)/(exp(choicesens*ul.k)+exp(choicesens*ur.k)))^choice_x.k * choicesens*exp(choicesens*(ul.k + 2*ur.k))*(indicators_x*(self_x.k-other_x.k) + indicators_y*(other_y.k-self_y.k)) * (exp(choicesens*ur.k)/(exp(choicesens*ul.k)+exp(choicesens*ur.k)))^(-choice_x.k) * (choice_x.k * exp(-choicesens*ul.k) - exp(-choicesens*ur.k) + choice_x.k * exp(-choicesens*ur.k)) ))
	beta.grad <- colSums(beta.grad)

	c(beta.grad, choicesens.grad)
}
f.mloglik.pregrad <- cmpfun(f.mloglik.pregrad)

f.mloglik.grad <- function(theta.v, choice_x, self_x, other_x, self_y, other_y, indicators_x, indicators_y, sel, nc, k, n, j, tau) {
	ret <- rep(0, (k+1)*nc)

	for (c in 1:nc) {
		ret[((k+1)*(c-1)+1):((k+1)*c)] <- f.mloglik.pregrad(theta.v[((c-1)*(k+1)+1):(c*(k+1))], tau[,c], choice_x, self_x, other_x, self_y, other_y, indicators_x, indicators_y, sel, nc, k, n, j)
	}

	ret
}
f.mloglik.grad <- cmpfun(f.mloglik.grad)

# E-Step Function
f.estep <- function(theta.v, choice_x, self_x, other_x, self_y, other_y, indicators_x, indicators_y, sel, nc, k, n, j, mix) {
	m <- matrix(mix, j, nc, byrow=T)

	numerator   <- m * f.compdens(theta.v, choice_x, self_x, other_x, self_y, other_y, indicators_x, indicators_y, sel, nc, k, n, j)
	denominator <- matrix(rowSums(numerator), j, nc)

	ret <- matrix(numerator/denominator, j, nc, dimnames=list(colnames(sel), paste("comp",1:nc,sep="")))
	ret
}
f.estep <- cmpfun(f.estep)

# S-Step-Function
f.sstep <- function(em.tau,j,nc){
	sem.tau <- matrix(0,j,nc)
	for (i in 1:j){
		sem.tau[i,] <- t(rmultinom(1,1,em.tau[i,]))
	}
	sem.tau
}
f.sstep <- cmpfun(f.sstep)

# M-Step Function
f.mstep <- function(theta.v0, choice_x, self_x, other_x, self_y, other_y, indicators_x, indicators_y, sel, nc, k, n, j, tau) {
	ret <- optim(theta.v0, f.mloglik, f.mloglik.grad, method="BFGS", control=list(maxit=10000, trace=0, fnscale=-1), hessian=F, choice_x=choice_x, self_x=self_x, other_x=other_x, self_y=self_y, other_y=other_y, indicators_x=indicators_x, indicators_y=indicators_y, sel=sel, nc=nc, k=k, n=n, j=j, tau=tau)

	mix <- colMeans(tau)

	list(mix=mix, theta.v=ret$par)
}
f.mstep <- cmpfun(f.mstep)

# Annealing-Part-Function for the SAEM-Algorithm
f.annealing <- function(em.mix, em.vm, sem.mix, sem.vm, iter, b=20){
	{if (iter <= b)
		{q <- cos((iter/b)*acos(b/100))}
	 else
		{q <- (b/100)*sqrt(b/iter)}
	}

	mix <- q*sem.mix+(1-q)*em.mix
	vm  <- q*sem.vm+(1-q)*em.vm

	list(mix=mix, vm=vm)
}
f.annealing <- cmpfun(f.annealing)

# Sandwitch part of the clustered standard errors
f.sandwich <- function(psi, choice_x, self_x, other_x, self_y, other_y, indicators_x, indicators_y, sel, nc, k, n, j) {
	lpsi    <- length(psi);
	theta.v <- psi[1:(lpsi-nc+1)];
	theta.m <- matrix(theta.v, k+1, nc);
	mix     <- psi[(lpsi-nc+2):lpsi];
	mix     <- c(mix, 1-sum(mix));
	m       <- matrix(mix, j, nc, byrow=T);

	# Individual gradient
	nj           <- n*j;
	choice_x.k <- matrix(as.vector(choice_x), nj, k);
	sel.k        <- matrix(as.vector(sel)       , nj, k);
	self_x.k   <- matrix(as.vector(self_x)  , nj, k);
	other_x.k   <- matrix(as.vector(other_x)  , nj, k);
	self_y.k  <- matrix(as.vector(self_y) , nj, k);
	other_y.k  <- matrix(as.vector(other_y) , nj, k);

	cdens.m <- f.compdens(theta.v, choice_x, self_x, other_x, self_y, other_y, indicators_x, indicators_y, sel, nc ,k, n, j);
	denom.outerderiv   <- 1/matrix(rowSums(m*cdens.m), n, j, byrow=T);
	denom.outerderiv.k <- matrix(as.vector(denom.outerderiv), nj, k);

	ret.theta <- matrix(0, nj, nc*(k+1));
	ret.mix   <- matrix(0, nj, nc);

	for (c in 1:nc) {
		gamma <- theta.m[k+1, c];
		beta  <- theta.m[1:k, c];

		lli   <- matrix(as.vector(indicators_x  %*% beta), n, j); #Left linear index
		rli   <- matrix(as.vector(indicators_y %*% beta), n, j); #Right linear index

		ul    <- (1-lli)*self_x  + lli*other_x;  #Utility from "left" choice
		ur    <- (1-rli)*self_y + rli*other_y; #Utility from "right" choice
		ul.k  <- matrix(as.vector(ul), nj, k);
		ur.k  <- matrix(as.vector(ur), nj, k);

		innerderiv   <- mix[c]*matrix(cdens.m[,c],n,j,byrow=T)/((exp(gamma*ul)/(exp(gamma*ul)+exp(gamma*ur)))^choice_x * (exp(gamma*ur)/(exp(gamma*ul)+exp(gamma*ur)))^(1-choice_x));
		innerderiv.k <- matrix(as.vector(innerderiv), nj, k);

		gamma.grad <- sel * denom.outerderiv * innerderiv * ((1/(exp(gamma*ul)+exp(gamma*ur))^2) * ( (exp(gamma*ul)/(exp(gamma*ul)+exp(gamma*ur)))^choice_x * exp(gamma*(ul + 2*ur))*(-ur + ul) * (exp(gamma*ur)/(exp(gamma*ul)+exp(gamma*ur)))^(-choice_x) * (choice_x * exp(-gamma*ul) - exp(-gamma*ur) + choice_x * exp(-gamma*ur)) ));

		beta.grad <- sel.k * denom.outerderiv.k * innerderiv.k * ((-1/(exp(gamma*ul.k)+exp(gamma*ur.k))^2) * ( (exp(gamma*ul.k)/(exp(gamma*ul.k)+exp(gamma*ur.k)))^choice_x.k * gamma*exp(gamma*(ul.k + 2*ur.k))*(indicators_x*(self_x.k-other_x.k) + indicators_y*(other_y.k-self_y.k)) * (exp(gamma*ur.k)/(exp(gamma*ul.k)+exp(gamma*ur.k)))^(-choice_x.k) * (choice_x.k * exp(-gamma*ul.k) - exp(-gamma*ur.k) + choice_x.k * exp(-gamma*ur.k)) ));

		mix.grad  <- cdens.m[,c] * denom.outerderiv[1,]

		ret.theta[,((c-1)*(k+1)+1):(c*(k+1))] <- cbind(beta.grad, as.vector(gamma.grad));
		ret.mix[,c]                           <- as.vector(matrix(mix.grad/n,n,j,byrow=T));
	};

	ret.mix <- ret.mix - matrix(ret.mix[,nc], nj, nc)
	ret     <- cbind(ret.theta, ret.mix[,1:(nc-1)]);

	#Sandwich part
	sandwich <- matrix(0, lpsi, lpsi)
	for (i in 1:j) {
		gradi <- colSums(ret[((i-1)*n+1):(i*n),])
		sandwich <- sandwich + (gradi %o% gradi)
	}

	qc <- (nj-1) / (nj-lpsi) * j/(j-1) #df adjustment

	list(qc=qc, sandwich=sandwich)
}
f.sandwich <- cmpfun(f.sandwich)

# Function to compute clustered standard errors, relies on f.sandwich
f.clusterse <- function(estobj, dat, persid, label_choice_x, label_indicators_x, label_indicators_y, self_x, other_x, self_y, other_y, nc, sortrow=1) {
	prep         <- f.dataprep(dat, persid, label_choice_x, label_indicators_x, label_indicators_y, self_x, other_x, self_y, other_y, nafiller=1);
	n            <- prep$n;
	j            <- prep$j;
	k            <- prep$k;
	indicators_x <- prep$indicators_x;
	indicators_y <- prep$indicators_y;
	self_x       <- prep$self_x;
	other_x      <- prep$other_x;
	self_y       <- prep$self_y;
	other_y      <- prep$other_y;
	choice_x     <- prep$choice_x;
	sel          <- prep$sel;
	rm(prep);

	ll1 <- f.totloglik(estobj$ret$par, choice_x, self_x, other_x, self_y, other_y, indicators_x, indicators_y, sel, nc, k, n, j)
	if (abs(ll1 - estobj$ret$value)>0.001) stop(paste("Saved loglik", ll1, "and computed loglik", estobj$value, "differ"))


	sandwich <- f.sandwich(estobj$ret$par, choice_x, self_x, other_x, self_y, other_y, indicators_x, indicators_y, sel, nc, k, n, j)

	lpsi <- length(estobj$ret$par)
	vcm <- solve(-estobj$ret$hessian)
	vcm <- sandwich$qc * (vcm %*% sandwich$sandwich %*% vcm)

	se <- sqrt(diag(vcm))
	lastmix.se <- sqrt( rep(1,nc-1) %*% vcm[(lpsi-(nc-2)):lpsi, (lpsi-(nc-2)):lpsi] %*% rep(1,nc-1) )
	psi.se <- c(se,lastmix.se)
	psi    <- estobj$ret$par

	#Extract the point estimates
	theta.v    <- psi[1:(lpsi-nc+1)];
	mix        <- psi[(lpsi-nc+2):lpsi];
	mix        <- c(mix, 1-sum(mix));

	mix.m      <- matrix(mix,       1, nc, dimnames=list(c("mixture"), paste("comp",1:nc,sep="")));
	theta.m    <- matrix(theta.v, k+1, nc, dimnames=list(c(paste(label_indicators_x,label_indicators_y),"choicesens"), paste("comp",1:nc,sep="")));

	#Extract the standard errors
	theta.v.se <- psi.se[1:(lpsi-nc+1)];
	mix.se     <- psi.se[(lpsi-nc+2):(lpsi+1)];

	mix.m.se   <- matrix(mix.se,       1, nc, dimnames=list(c("mixture"), paste("comp",1:nc,sep="")));
	theta.m.se <- matrix(theta.v.se, k+1, nc, dimnames=list(c(paste(label_indicators_x,label_indicators_y),"choicesens"), paste("comp",1:nc,sep="")));

	pe.out <- rbind(mix.m,    theta.m);
	se.out <- rbind(mix.m.se, theta.m.se);

	#Sort the groups according to their size to avoid problems of label switching
	sortidx  <- order(pe.out[sortrow,]);
	pe.out   <- pe.out[, sortidx];
	se.out   <- se.out[, sortidx];

	p.out  <- 2*(1-pnorm(abs(pe.out)/se.out));

	cat("POINT ESTIMATES\n")
	print(round(pe.out,3));
	cat("\n------------------------------------------------\n");
	cat("ROBUST STANDARD ERRORS\n")
	print(round(se.out,3));
	cat("\n------------------------------------------------\n");
	cat("ROBUST P-VALUES\n")
	print(round(p.out,3));

	list(se=se.out, p=p.out, vcm=vcm)
}
f.clusterse <- cmpfun(f.clusterse)

# Main Estimation Function
f.estim <- function(dat, v0path, persid, nc, label_choice_x, label_indicators_x, label_indicators_y, label_self_x, label_other_x, label_self_y, label_other_y, saem=F, loglik1=NA, maxemiter=50, diffcrit=5e-3, sortrow=1) {

	# Prepares the data and stores it
	prep       <- f.dataprep(dat, persid, label_choice_x, label_indicators_x, label_indicators_y, label_self_x, label_other_x, label_self_y, label_other_y, nafiller=1)
	n          <- prep$n
	j          <- prep$j
	k          <- prep$k
	indicators_x      <- prep$indicators_x
	indicators_y     <- prep$indicators_y
	self_x   <- prep$self_x
	other_x   <- prep$other_x
	self_y  <- prep$self_y
	other_y  <- prep$other_y
	choice_x <- prep$choice_x
	sel        <- prep$sel
	rm(prep)

	# Prepares the start values
	v0      <- read.table(v0path, sep=",", header=F)[,1]
	theta.v <- runif(nc*(1+k),0.8,1.2) * rep(v0,nc) # Random variation of the start values from the pooled model
	aa      <- runif(nc, 0.5, 2)
	mix     <- aa/sum(aa) # Initial mixture
	rm(aa)

	# Initial EM control flow parameters
	iter  <- 1
	diff  <- 1
	llold <- 0
	llnew <- 1

	# EM/SAEM algorithm (see Dempster et al. (1977) and book by McLachlan (2000))
	cat("Starting EM type maximization...\n")
	while (iter < maxemiter+1 & abs(diff) > diffcrit) {
		cat("\n  EM-Type Iteration:", iter,"of max.",maxemiter,"\n")
		cat("\n     working on E-step.\n")
		tau.em <- f.estep(theta.v, choice_x, self_x, other_x, self_y, other_y, indicators_x, indicators_y, sel, nc, k, n, j, mix)
		cat("     working on M-step (EM algorithm).\n")
		mstep <- f.mstep(theta.v, choice_x, self_x, other_x, self_y, other_y, indicators_x, indicators_y, sel, nc, k, n, j, tau.em)
		theta.v.em <- mstep$theta.v
		mix.em     <- mstep$mix

		if (saem==T) {
			cat("\n     working on S-step.\n")
			tau.sem <- f.sstep(tau.em, j, nc)

			cat("     working on M-step (SEM algorithm).\n")
			mstep <- f.mstep(theta.v, choice_x, self_x, other_x, self_y, other_y, indicators_x, indicators_y, sel, nc, k, n, j, tau.sem)
			theta.v.sem <- mstep$theta.v
			mix.sem     <- mstep$mix

			ann <- f.annealing(mix.em, theta.v.em, mix.sem, theta.v.sem, iter, b=20)
			mix     <- ann$mix
			theta.v <- ann$v
		 } else {
		 	mix     <- mix.em
		 	theta.v <- theta.v.em
		 }

		# Update the control flow parameters
		iter  <- iter+1
		llold <- llnew
		psi   <- c(theta.v, mix[1:(length(mix)-1)])
		llnew <- f.totloglik(psi, choice_x, self_x, other_x, self_y, other_y, indicators_x, indicators_y, sel, nc, k, n, j)
		diff  <- llnew-llold
		cat("\n  Achieved Log Likelihood:", round(llnew,4),"\n\n  ----------\n")
	}
	cat("\n...done. Parameters after EM-type estimation:\n")
	mix.m   <- matrix(mix,       1, nc, dimnames=list(c("mixture"), paste("comp",1:nc,sep="")))
	theta.m <- matrix(theta.v, k+1, nc, dimnames=list(c(paste(label_indicators_x,label_indicators_y),"choicesens"), paste("comp",1:nc,sep="")))
	out     <- rbind(mix.m, theta.m)

	#Direct ML estimation
	ret  <- optim(psi,f.totloglik, f.totloglik.gradient, method="BFGS", control=list(maxit=100000,trace=1,REPORT=20, fnscale=-1), hessian=T, choice_x=choice_x, self_x=self_x, other_x=other_x, self_y=self_y, other_y=other_y, indicators_x=indicators_x, indicators_y=indicators_y, sel=sel, nc=nc, k=k, n=n, j=j)

	psi     <- ret$par #Extract the MLE
	lpsi    <- length(psi)
	psi.se <- sqrt(diag(solve(-ret$hessian))) #Standard errors

	#Extract the point estimates
	theta.v    <- psi[1:(lpsi-nc+1)]
	mix        <- psi[(lpsi-nc+2):lpsi]
	mix        <- c(mix, 1-sum(mix))

	mix.m      <- matrix(mix,       1, nc, dimnames=list(c("mixture"), paste("comp",1:nc,sep="")))
	theta.m    <- matrix(theta.v, k+1, nc, dimnames=list(c(paste(label_indicators_x,label_indicators_y),"choicesens"), paste("comp",1:nc,sep="")))

	#Extract the standard errors
	theta.v.se <- psi.se[1:(lpsi-nc+1)]
	mix.se   <- psi.se[(lpsi-nc+2):lpsi]
	mix.se   <- c(mix.se,NA)

	mix.m.se   <- matrix(mix.se,       1, nc, dimnames=list(c("mixture"), paste("comp",1:nc,sep="")))
	theta.m.se <- matrix(theta.v.se, k+1, nc, dimnames=list(c(paste(label_indicators_x,label_indicators_y),"choicesens"), paste("comp",1:nc,sep="")))

	tau     <- f.estep(theta.v, choice_x, self_x, other_x, self_y, other_y, indicators_x, indicators_y, sel, nc, k, n, j, mix)

	aic <- -2*ret$value + 2*lpsi
	bic <- -2*ret$value + log(sum(sel))*lpsi

	tauv   <- as.vector(tau)
	nec    <- sum(na.omit(-tauv*log(tauv)))/(ret$value-loglik1)

	i.out  <- matrix(c(ret$value, lpsi, sum(sel), sortrow, aic, bic, nec), , 1, dimnames=list(c("Loglik:", "# of parameters:", "# of observations:", "Sortrow (to avoid label switching):", "AIC:", "BIC:", "NEC:"), c()))
	pe.out <- rbind(mix.m,    theta.m)
	se.out <- rbind(mix.m.se, theta.m.se)

	#Sort the groups according to their size to avoid problems of label switching
	sortidx  <- order(pe.out[sortrow,])
	pe.out   <- pe.out[, sortidx]
	se.out   <- se.out[, sortidx]
	tau      <- tau[, sortidx]

	p.out  <- 2*(1-pnorm(abs(pe.out)/se.out))

	cat("\n------------------------------------------------\n")
	cat("GENERAL INFORMATION\n")
	print(round(i.out,3))
	cat("\n------------------------------------------------\n")
	cat("POINT ESTIMATES (NOT CLUSTERED)\n")
	print(round(pe.out,3))
	cat("\n------------------------------------------------\n")
	cat("STANDARD ERRORS (NOT CLUSTERED)\n")
	print(round(se.out,3))
	cat("\n------------------------------------------------\n")
	cat("P-VALUES\n")
	print(round(p.out,3))
	cat("\n------------------------------------------------\n")

	list(pe=pe.out, se=se.out, p=p.out, info=i.out, ret=ret, tau=tau, success=T)
}
f.estim <- cmpfun(f.estim)




# Load choice data from both experiments
dat1     <- read.table("choices_exp1.csv",  sep=",", header=T)
dat2     <- read.table("choices_exp2.csv",  sep=",", header=T)

# Remove erratic subjects (see 2nd paragraph of section 4; to replicate estimate individual parameters)
dropped  <- read.table("dropped_subjects_section4paragraph2.csv", sep=",", header=F)[,1]
dat1     <- dat1[(dat1$sid %in% dropped) == F, ]
dat2     <- dat2[(dat2$sid %in% dropped) == F, ]


# Estimating Finite Model with nc=3 Types using start values from aggregate estimation
cat("\n\n\n================== ESTIMATING SESSION 1 ================================\n")
pe1 <- f.estim(dat1, "svExp1_1cl.csv", "sid", nc=3, "choice_x", c("s_x","r_x","q","v"), c("s_y","r_y","q","v"), "self_x", "other_x", "self_y", "other_y", saem=T, loglik1=-5472.314, maxemiter=30, diffcrit=1e-3,sortrow=1)

cat("\n\n\n================== ESTIMATING SESSION 2 ================================\n")
pe2 <- f.estim(dat2, "svExp2_1cl.csv", "sid", nc=3, "choice_x", c("s_x","r_x","q","v"), c("s_y","r_y","q","v"), "self_x", "other_x", "self_y", "other_y", saem=T, loglik1=-4540.739, maxemiter=30, diffcrit=1e-3,sortrow=1)

cat("\n\nTo extract the individual ex-post probabilities of type membership type pe1$tau or pe2$tau.\n")


# Computing Robust Standard Errors
cat("\n\n\n================== ROBUST STANDARD ERRORS: SESSION 1 ===================\n")
rse1 <- f.clusterse(pe1, dat1, "sid", "choice_x", c("s_x","r_x","q","v"), c("s_y","r_y","q","v"), "self_x", "other_x", "self_y", "other_y", 3, sortrow=1)
cat("\n\n\n================== ROBUST STANDARD ERRORS: SESSION 2 ===================\n")
rse2 <- f.clusterse(pe2, dat2, "sid", "choice_x", c("s_x","r_x","q","v"), c("s_y","r_y","q","v"), "self_x", "other_x", "self_y", "other_y", 3, sortrow=1)

# Save the R workspace
save.image("finmix_est_3cl.RData")
