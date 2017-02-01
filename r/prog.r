##############################################################
#                                                            #
# Sample program for running KC-RBM. Note that results may   #
# slightly differ from those in the paper, e.g. due to using #
# different Ensembl releases.                                #
#                                                            #
##############################################################

source("../../kcrbm/r/kcrbm.r")
source("../../kcrbm/r/kcrbm_tools.r")

# retrieve ensembl data
edata <- preprocessEnsemblData(dataset="mmusculus")
#save(edata,file="ensembl_data.RData")

#############################
#                           #
#      MuLV default         #
#                           #
#############################

# load MuLV insertion data
idata <- get(load("../data/MuLV_mm10.RData"))

# map insertions to genes/transcripts
# Setting map_to = "transcripts" will use transcript, intron, and 
# exon related information. Note that this is much slower than 
# setting map_to = "genes".
t1 <- Sys.time()
iset <- kcrbm(
	edata = edata,
	idata = idata,
	rules = "MuLV", # parameter values used for MuLV in NAR paper
	map_to = "transcripts"
)
t2 <- Sys.time()
print(difftime(t2,t1))

# Use kcrbm to map to genes: returns only the 
# ensembl ID of the target gene, and no transcript, intron, exon related
# information.
# NOTE: mapping to genes is much (!) faster than mapping to transcripts
t1 <- Sys.time()
iset.gene <- kcrbm(
	edata = edata,
	idata = idata,
	rules = "MuLV", # parameter values used for MuLV in NAR paper
	map_to = "genes"
)
t2 <- Sys.time()
print(difftime(t2,t1))

# The above mapping results in multiple targets per insertion. Now 
# select a single target per insertion: for each insertion, select 
# that gene that was most commonly targeted across insertions.
iset.single <- selectSingleTargets(iset.gene,type="nearest") # "ctg"
# retrieve the commonly targeted genes (CTGs) for this insertion to single-gene mapping
ctgs <- get_ctgs(iset.single)
plot_ctgs(ctgs[1:20])


#############################
#                           #
#       SB defaults         #
#                           #
#############################

# load the insertion data
# Note that for demonstration purposes MuLV insertion data is used with SB defaults
idata <- get(load("../data/MuLV_mm10.RData"))

# map insertions to genes/transcripts
t1 <- Sys.time()
iset <- kcrbm(
	edata = edata,
	idata = idata,
	rules = "SB" # parameter values used for SB in NAR paper; uses "TA"-motif for correction
)
t2 <- Sys.time()
print(difftime(t2,t1))

# The above mapping results in multiple targets per insertion. Now 
# select a single target per insertion: for each insertion, select 
# that gene that was most commonly targeted across insertions.
iset.single <- selectSingleTargets(iset,type="ctg")
# retrieve the commonly targeted genes (CTGs) for this insertion to single-gene mapping
ctgs <- get_ctgs(iset.single)
# plot the top 20 in a barplot
plot_ctgs(ctgs[1:20])


#############################
#                           #
#   use other than default  #
#   parameter values        #
#                           #
#############################

# load MuLV insertion data
idata <- get(load("../data/MuLV_mm10.RData"))

# map insertions to genes/transcripts
t1 <- Sys.time()
iset <- kcrbm(
	edata = edata,
	idata = idata,
	windows = rep(1000,4), # window sizes
	scale = 20000, # GKC scale
	orientation_homogeneity = 0.6 # required insertional homogeneity of a cluster
)
t2 <- Sys.time()
print(difftime(t2,t1))

# The above mapping results in multiple targets per insertion. Now 
# select a single target per insertion: for each insertion, select 
# that gene that was most commonly targeted across insertions.
iset.single <- selectSingleTargets(iset,type="nearest")
# retrieve the commonly targeted genes (CTGs) for this insertion to single-gene mapping
ctgs <- get_ctgs(iset.single)
# plot the top 20 in a barplot
plot_ctgs(ctgs[1:20])






