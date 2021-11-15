# If you use this code or parts of it, please cite the following paper
# Bruhin A., E. Fehr, D. Schunk (2018): "The Many Faces of Human Sociality: Uncovering the Distribution and Stability of Social Preferences",
# Journal of the European Economic Association, forthcoming.
# ===========================================================================================================================================

rm(list=ls()) # Clear the memory
source("99_estim_functions_aggregate_individual.r") # Load the required functions for performing the estimations
set.seed(51)

d.exp1 <- read.table("choices_exp1.csv", sep=",", header=T)
d.exp2 <- read.table("choices_exp2.csv", sep=",", header=T)

cat("\n\nSESSION 1\n")
iest1 <- f.indivest(d.exp1) # Estimate on Session 1 Data
cat("\n\nSESSION 2\n")
iest2 <- f.indivest(d.exp2) # Estimate on Session 2 Data

# Identify erratic subjects with parameter estiamates outside the identifiable range and drop them (see 2nd paragraph of section 4)
drop        <- f.droperratic(iest1, iest2)
sid.drop    <- as.numeric(names(drop)[drop>0]) # IDs of subjects that need to be dropped
iest1.valid <- iest1[(row.names(iest1) %in% names(drop)[drop>0]) == F,]
iest2.valid <- iest2[(row.names(iest2) %in% names(drop)[drop>0]) == F,]
write.table(sid.drop, "dropped_subjects_section4paragraph2.csv", sep=",", col.names=F, row.names=F) # Create file with IDs of dropped subjects

cat("\n\n\n---- SUMMARY STATS OF INDIVIDUAL ESTIMATES: SESSION 1 ----\n")
print(summary(iest1.valid))
cat("\nStandard Deviations\n")
print(apply(iest1.valid,2,sd))

cat("\n\n\n---- SUMMARY STATS OF INDIVIDUAL ESTIMATES: SESSION 2 ----\n")
print(summary(iest2.valid))
cat("\nStandard Deviations\n")
print(apply(iest2.valid,2,sd))

# Write the results into two csv files
out1 <- cbind(as.numeric(rownames(iest1.valid)), iest1.valid)
out2 <- cbind(as.numeric(rownames(iest2.valid)), iest2.valid)
colnames(out1)[1] <- "sid"
colnames(out2)[1] <- "sid"

write.table(out1, "individual_est_exp1.csv", col.names=T, row.names=F, sep=",")
write.table(out2, "individual_est_exp2.csv", col.names=T, row.names=F, sep=",")
