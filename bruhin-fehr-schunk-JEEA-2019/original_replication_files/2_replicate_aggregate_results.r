# If you use this code or parts of it, please cite the following paper
# Bruhin A., E. Fehr, D. Schunk (2018): "The Many Faces of Human Sociality: Uncovering the Distribution and Stability of Social Preferences",
# Journal of the European Economic Association, forthcoming.
# ===========================================================================================================================================

rm(list=ls()) # Clear the memory
source("99_estim_functions_aggregate_individual.r") # Load the required functions for performing the estimations
set.seed(51)


# SESSION 1
cat("\n---- SESSION 1 -------------\n")

dat     <- read.table("choices_exp1.csv",  sep=",", header=T) # Load choice data
dropped <- read.table("dropped_subjects_section4paragraph2.csv", sep=",", header=F)[,1] # dropped subjects (see 2nd paragraph of section 4; to replicate estimate individual parameters)
dat     <- dat[(dat$sid %in% dropped) == F, ] # Eliminate dropped subjects

cat("\nEstimation Results without Clustered Standard Errors:\n\n")
pe.exp1 <- f.estim(dat, "choice_x", c("s_x","r_x","q","v"), c("s_y","r_y","q","v"), "self_x", "other_x", "self_y", "other_y", silent=F, comp.se=T) # Aggregate estimation

cat("\nEstimation Results with Clustered Standard Errors:\n")
clust.exp1 <- f.clusterse(pe.exp1, dat, "choice_x", "self_x", "other_x", "self_y", "other_y", c("s_x","r_x","q","v"), c("s_y","r_y","q","v"), "sid")

write.table(pe.exp1$out[,1],"svExp1_1cl.csv",col.names=F, row.names=F) # Save point estimates as start values for finite mixture model



# SESSION 2
cat("\n---- SESSION 2 -------------\n")

dat     <- read.table("choices_exp2.csv",  sep=",", header=T) # Load choice data
dropped <- read.table("dropped_subjects_section4paragraph2.csv", sep=",", header=F)[,1] # dropped subjects (see 2nd paragraph of section 4; to replicate estimate individual parameters)
dat     <- dat[(dat$sid %in% dropped) == F, ] # Eliminate dropped subjects

cat("\nEstimation Results without Clustered Standard Errors:\n\n")
pe.exp2 <- f.estim(dat, "choice_x", c("s_x","r_x","q","v"), c("s_y","r_y","q","v"), "self_x", "other_x", "self_y", "other_y", silent=F, comp.se=T) # Aggregate estimation

cat("\nEstimation Results with Clustered Standard Errors:\n")
clust.exp2 <- f.clusterse(pe.exp2, dat, "choice_x", "self_x", "other_x", "self_y", "other_y", c("s_x","r_x","q","v"), c("s_y","r_y","q","v"), "sid")

write.table(pe.exp2$out[,1],"svExp2_1cl.csv",col.names=F, row.names=F) # Save point estimates as start values for finite mixture model
