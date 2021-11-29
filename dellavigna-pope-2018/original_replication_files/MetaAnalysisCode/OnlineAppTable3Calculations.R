### This script contains the calculations underlying the implied effort under the
### average probability weights for p=0.01 and p=0.5 at the bottom of online appendix table 3.

rm(list=ls())

######################## Functions for the various calculations ##############################
# Simple function to check that value of gamma, k, s and p gives the actual level of effort
effort_fn = function(gamma,k,s,p) {((s+p^theta)/k)^(1/gamma)}
# Function that calculates the implied effort under for a given probability weight for p=0.01 or 0.5.
effort_pweight = function(pweight, p, theta){
  ((s + pweight*p^theta)/k)^(1/(gamma))
}

# Actual levels of effort
e0 = 1521 # Effort under zero piece rate treatment
e001 = 2029 # Effort under 1 cent piece rate
e01 = 2175 # Effort under 10 cent piece rate

# Average probability weights from the meta-analysis
pweight01 = 0.060 # Average p-weight for p=0.01
pweight50 = 0.452 # Average p-weight for p=0.5
p01 = 1 # Reward (in dollars) that occurs with probability p=0.01
p50 = 0.02 # Reward (in dollars) that occurs with probability p=0.5

### Functions that solve gamma, k and s (for a given value of theta)
# This is a system of 3 equations with 3 unknowns, but due to its non-linear nature,
# we need to solve numerically for one of the variables (using the "uniroot" function).
find_s = function(gamma, k){
  k*(e0^gamma)
}
find_k = function(gamma){
  (0.01^theta)/( (e001^gamma)-(e0^gamma) )
}

gamma_eqn = function(gamma){
  (e01^gamma) - (10^theta)*(e001^gamma) + (10^theta-1)*(e0^gamma)
}

################## Assuming theta = 1 (linear value function) ################
theta=1

gamma_root = uniroot(f=gamma_eqn, interval = c(1,50))
gamma = gamma_root$root
k = find_k(gamma=gamma)
s = find_s(gamma=gamma, k=k)

### Check that it fits the moments perfectly
effort_fn(gamma=gamma, k=k, s=s, p=0)
effort_fn(gamma=gamma, k=k, s=s, p=0.01)
effort_fn(gamma=gamma, k=k, s=s, p=0.1)
gamma; k; s

### Check 4-cent fit
effort_fn(gamma=gamma, k=k, s=s, p=0.04) # Actual is 2132, predicted is 2115.68

### Predict probabilistic treatments ###
effort_pweight(pweight=pweight01, p=p01, theta=theta) ### 2141.728
effort_pweight(pweight=pweight50, p=p50, theta=theta) ### 2022.83

################## Assuming theta = 0.88 (as in TK) ################
theta=0.88

gamma_root = uniroot(f=gamma_eqn, interval = c(1,50))
gamma = gamma_root$root
k = find_k(gamma=gamma)
s = find_s(gamma=gamma, k=k)

### Check that it fits the moments perfectly
effort_fn(gamma=gamma, k=k, s=s, p=0)
effort_fn(gamma=gamma, k=k, s=s, p=0.01)
effort_fn(gamma=gamma, k=k, s=s, p=0.1)
gamma; k; s

### Check 4-cent fit
effort_fn(gamma=gamma, k=k, s=s, p=0.04) # Actual is 2132, predicted is 2115.68

### Predict probabilistic treatments ###
effort_pweight(pweight=pweight01, p=p01, theta=theta) ### 2117.073
effort_pweight(pweight=pweight50, p=p50, theta=theta) ### 2016.233



################## Assuming theta = 0.7 ################
theta=0.7

gamma_root = uniroot(f=gamma_eqn, interval = c(1,50))
gamma = gamma_root$root
k = find_k(gamma=gamma)
s = find_s(gamma=gamma, k=k)

### Check that it fits the moments perfectly
effort_fn(gamma=gamma, k=k, s=s, p=0)
effort_fn(gamma=gamma, k=k, s=s, p=0.01)
effort_fn(gamma=gamma, k=k, s=s, p=0.1)
gamma; k; s

### Check 4-cent fit
effort_fn(gamma=gamma, k=k, s=s, p=0.04) # Actual is 2132, predicted is 2115.667

### Predict probabilistic treatments ###
effort_pweight(pweight=pweight01, p=p01, theta=theta) ### 2065.186
effort_pweight(pweight=pweight50, p=p50, theta=theta) ### 2002.185



# The remaining calculations assume values of theta not shown in the paper may possibly be
# of interest to the exceptionally curious reader.
################## Assuming theta = 0.6 ################
theta=0.6

gamma_root = uniroot(f=gamma_eqn, interval = c(1,50))
gamma = gamma_root$root
k = find_k(gamma=gamma)
s = find_s(gamma=gamma, k=k)

### Check that it fits the moments perfectly
effort_fn(gamma=gamma, k=k, s=s, p=0)
effort_fn(gamma=gamma, k=k, s=s, p=0.01)
effort_fn(gamma=gamma, k=k, s=s, p=0.1)
gamma; k; s

### Check 4-cent fit
effort_fn(gamma=gamma, k=k, s=s, p=0.04) # Actual is 2132, predicted is 2115.643

### Predict probabilistic treatments ###
effort_pweight(pweight=pweight01, p=p01, theta=theta) ### 2023.881
effort_pweight(pweight=pweight50, p=p50, theta=theta) ### 1991.855


################## Assuming theta = 0.5 ################
theta=0.5

gamma_root = uniroot(f=gamma_eqn, interval = c(1,50))
gamma = gamma_root$root
k = find_k(gamma=gamma)
s = find_s(gamma=gamma, k=k)

### Check that it fits the moments perfectly
effort_fn(gamma=gamma, k=k, s=s, p=0)
effort_fn(gamma=gamma, k=k, s=s, p=0.01)
effort_fn(gamma=gamma, k=k, s=s, p=0.1)
gamma; k; s

### Check 4-cent fit
effort_fn(gamma=gamma, k=k, s=s, p=0.04) # Actual is 2132, predicted is 2115.585

### Predict probabilistic treatments ###
effort_pweight(pweight=pweight01, p=p01, theta=theta) ### 1967.771
effort_pweight(pweight=pweight50, p=p50, theta=theta) ### 1975.24

################## Assuming theta = 0.4 ################
theta=0.4

gamma_root = uniroot(f=gamma_eqn, interval = c(1,50))
gamma = gamma_root$root
k = find_k(gamma=gamma)
s = find_s(gamma=gamma, k=k)

### Check that it fits the moments perfectly
effort_fn(gamma=gamma, k=k, s=s, p=0)
effort_fn(gamma=gamma, k=k, s=s, p=0.01)
effort_fn(gamma=gamma, k=k, s=s, p=0.1)
gamma; k; s

### Check 4-cent fit
effort_fn(gamma=gamma, k=k, s=s, p=0.04) # Actual is 2132, predicted is 2115.446

### Predict probabilistic treatments ###
effort_pweight(pweight=pweight01, p=p01, theta=theta) ### 1888.897
effort_pweight(pweight=pweight50, p=p50, theta=theta) ### 1952.564

################## Assuming theta = 0.3 ################
theta=0.3

gamma_root = uniroot(f=gamma_eqn, interval = c(1,50))
gamma = gamma_root$root
k = find_k(gamma=gamma)
s = find_s(gamma=gamma, k=k)

### Check that it fits the moments perfectly
effort_fn(gamma=gamma, k=k, s=s, p=0)
effort_fn(gamma=gamma, k=k, s=s, p=0.01)
effort_fn(gamma=gamma, k=k, s=s, p=0.1)
gamma; k; s

### Check 4-cent fit
effort_fn(gamma=gamma, k=k, s=s, p=0.04) # Actual is 2132, predicted is 2115.116

### Predict probabilistic treatments ###
effort_pweight(pweight=pweight01, p=p01, theta=theta) ### 1779.35
effort_pweight(pweight=pweight50, p=p50, theta=theta) ### 1917.417