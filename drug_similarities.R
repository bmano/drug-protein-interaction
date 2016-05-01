library("ChemmineR") 
library("ChemmineOB")
setwd('/Users/bia/Desktop/bioinformatics2/drug-protein-interaction/') 

# load data
drugs_id <- read.table("data_project/Drug_ids.txt", header = FALSE)
# SMILES format into SDF format
convertFormatFile("SMI","SDF","data_project/Drug_smiles.txt","data_project/drugs.sdf")
sdfset <- read.SDFset("data_project/drugs.sdf")

# verifying sdfset content
header(sdfset[1:4])
atomcount(sdfset[1:4])

# computing fingerprints using 4 different methods
# FP2, FP3, FP4, MACCS
fingp_2 <- fingerprintOB(sdfset, "FP2")
fingp_3 <- fingerprintOB(sdfset, "FP3")
fingp_4 <- fingerprintOB(sdfset, "FP4")
fingp_m <- fingerprintOB(sdfset, "MACCS")

# computing fingerprint-based Tanimoto Kernels 
fp2_kernel <- sapply(cid(fingp_2), function(x) fpSim(x=fingp_2[x], fingp_2, sorted=FALSE))
fp3_kernel <- sapply(cid(fingp_3), function(x) fpSim(x=fingp_3[x], fingp_3, sorted=FALSE))
fp4_kernel <- sapply(cid(fingp_4), function(x) fpSim(x=fingp_4[x], fingp_4, sorted=FALSE))
fpm_kernel <- sapply(cid(fingp_m), function(x) fpSim(x=fingp_m[x], fingp_m, sorted=FALSE))


