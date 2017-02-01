kcrbm <- function(
	edata, # ensembl data
	idata, # insertion data
	rules = NULL, # specify the use of default trans_rules for various model systems
	windows = NULL, # window sizes
	scale = NULL, # GKC scale
	orientation_homogeneity = NULL, # required insertional homogeneity of a cluster
	map_to = "genes" # either "genes" or "transcripts"
) {
	if (map_to == "genes") {
		kcrbm.gene(
			edata, # ensembl data
			idata, # insertion data
			rules = rules, # specify the use of default trans_rules for various model systems
			windows = windows, # window sizes
			scale = scale, # GKC scale
			orientation_homogeneity = orientation_homogeneity # required insertional homogeneity of a cluster
		)
	} else if (map_to == "transcripts") {
		kcrbm.trans(
			edata, # ensembl data
			idata, # insertion data
			rules = rules, # specify the use of default trans_rules for various model systems
			windows = windows, # window sizes
			scale = scales, # GKC scale
			orientation_homogeneity = orientation_homogeneity # required insertional homogeneity of a cluster
		)
	} else {
		stop("Argument 'map_to' must be set to one of \"genes\" or \"transcripts\"")
	}
}

# faster version of kcrbm() below, which only returns target gene IDs for each insertion,
# and the corresponding distances to the gene start site (fractional distance
# in case the insertion lies within a gene)
kcrbm.gene <- function(
	edata, # ensembl data
	idata, # insertion data
	rules = NULL, # specify the use of default trans_rules for various model systems
	windows = NULL, # window sizes
	scale = NULL, # GKC scale
	orientation_homogeneity = NULL # required insertional homogeneity of a cluster
) {


	if (!("gstart" %in% colnames(edata))) {
		edata$gstart <- edata$start
		edata$gend <- edata$end
	}
	if (!("gss" %in% colnames(edata))) {
		edata$gss <- edata$gstart
		edata$gss[edata$strand == -1] <- edata$gend[edata$strand == -1]
		edata$gts <- edata$gend
		edata$gts[edata$strand == -1] <- edata$gstart[edata$strand == -1]
	}
	if (!("geneid" %in% colnames(edata))) {
		edata$geneid <- edata$ensid
	}
		
	# Sort the insertions, just to be sure
	idata <- idata[order(idata$chr,idata$base,idata$ori),]
	
	if (!("id" %in% colnames(idata))) {
		idata$id <- 1:nrow(idata)
	}
	idata$internal_id <- 1:nrow(idata)
	
	if (!is.null(rules)) {
		if (rules %in% c("MuLV","MMTV")) {
			if (is.null(windows)) windows <- list(us=20000,ua=120000,ds=40000,da=5000)
			if (is.null(scale)) scale <- 10000
			if (is.null(orientation_homogeneity)) orientation_homogeneity <- 0.75
		} else if (rules %in% c("SB")) {
			if (is.null(windows)) windows <- list(us=20000,ua=10000,ds=25000,da=5000)
			if (is.null(scale)) scale <- 2000
			if (is.null(orientation_homogeneity)) orientation_homogeneity <- 0.75		
		}
	} else {
		if (
			is.null(windows) || 
			is.null(scale) || 
			is.null(orientation_homogeneity)
		) {
			stop(paste(
				"If argument rules is not provided,",
				"arguments windows, scale, and",
				"orientation_homogeneity need to be provided.\n"
			))
		} else {			
			if (!(
				is.numeric(orientation_homogeneity) &&
				0.5 < orientation_homogeneity && 
				orientation_homogeneity <= 1.0
			)) {
				stop("For parameter orientation_homogeneity is required: 0.5 < orientation_homogeneity <= 1.0\n")
			}
			
			if (!(
				is.numeric(scale) &&
				scale >= 0
			)) {
				stop("For parameter scale is required: scale >= 0.\n")
			}
			
			if (!(
				is.numeric(windows) &&
				length(windows) == 4
			)) {
				stop("RBM requires a windows parameter vector of length 4 (numeric).\n")
			}
			
			if (length(windows) > 1) {
				windows <- as.list(windows)
				names(windows) <- c("us","ua","ds","da")
			}			
		}
	}
	cat(sprintf("\nRunning KC-RBM with the following settings:\n"))
	cat(sprintf("\twindows:   (%s)\n",paste(paste(names(windows),windows,sep=" = "),collapse=", ")))
	cat(sprintf("\tscale:   %i\n", scale))
	cat(sprintf("\torientation_homogeneity:   %g\n\n", orientation_homogeneity))
	
	idata <- preprocess(idata,scale,orientation_homogeneity)

	edata <- edata[,which(colnames(edata) %in% c("chr","gstart","gend","strand","geneid","gss","gts"))]
	edata <- unique(edata,MARGIN=1)

	chrs <- sort(intersect(as.numeric(unique(edata$chr)),as.numeric(unique(idata$chr))))
	targets <- character(0)
	inss <- numeric(0)
	ds <- numeric(0)
	for (chr in chrs) {
		cat(sprintf("chr %s: ", chr))

		# all genes (transcripts/exons) on chromosome chr
		edata.chr <- edata[edata$chr == chr,]
		idata.chr <- idata[idata$chr == chr,]
		
		# for all insertions on chromosome chr ...
		inss.chr <- lapply(1:nrow(edata.chr), function(i) {
#			print(i)
			if (i %% 1000 == 0) {cat(sprintf("%i ",i))}
			if (edata.chr$strand[i] == 1) {
				I <- which(
					edata.chr$gstart[i] - windows$us <= idata.chr$clusterpeak & 
					idata.chr$clusterpeak <= edata.chr$gend[i] + windows$ds &
					idata.chr$clusterorientation == 1
				)
				if (length(I) > 0) {
					J <- findInterval(
						x = idata.chr$base[I],
						vec = c(edata.chr$gstart[i],edata.chr$gend[i]),
						rightmost.closed = TRUE
					)
					ds <- c(
						(idata.chr$base[I] - edata.chr$gstart[i])[J==0],
						((idata.chr$base[I] - edata.chr$gstart[i]) / 
							(edata.chr$gend[i] - edata.chr$gstart[i]))[J==1],
						(idata.chr$base[I] - edata.chr$gstart[i])[J==2]
					)
					IDs <- idata.chr$internal_id[I]
				} else {
					IDs <- NULL
					ds <- NULL
				}
				I <- which(
					edata.chr$gstart[i] - windows$ua <= idata.chr$clusterpeak & 
					idata.chr$clusterpeak <= edata.chr$gend[i] + windows$da &
					idata.chr$clusterorientation == -1
				)
				if (length(I) > 0) {
					J <- findInterval(
						x = idata.chr$base[I],
						vec = c(edata.chr$gstart[i],edata.chr$gend[i]),
						rightmost.closed = TRUE
					)
					da <- c(
						(idata.chr$base[I] - edata.chr$gstart[i])[J==0],
						((idata.chr$base[I] - edata.chr$gstart[i]) / 
							(edata.chr$gend[i] - edata.chr$gstart[i]))[J==1],
						(idata.chr$base[I] - edata.chr$gstart[i])[J==2]
					)
					IDa <- idata.chr$internal_id[I]
				} else {
					IDa <- NULL
					da <- NULL
				}
			} else {
				I <- which(
					edata.chr$gstart[i] - windows$ds <= idata.chr$clusterpeak & 
					idata.chr$clusterpeak <= edata.chr$gend[i] + windows$us &
					idata.chr$clusterorientation == -1
				)
				if (length(I) > 0) {
					J <- findInterval(
						x = idata.chr$base[I],
						vec = c(edata.chr$gstart[i],edata.chr$gend[i]),
						rightmost.closed = TRUE
					)
					ds <- c(
						(idata.chr$base[I] - edata.chr$gend[i])[J==0],
						((idata.chr$base[I] - edata.chr$gend[i]) / 
							(edata.chr$gend[i] - edata.chr$gstart[i]))[J==1],
						(idata.chr$base[I] - edata.chr$gend[i])[J==2]
					)
					IDs <- idata.chr$internal_id[I]
				} else {
					IDs <- NULL
					ds <- NULL
				}
				I <- which(
					edata.chr$gstart[i] - windows$da <= idata.chr$clusterpeak & 
					idata.chr$clusterpeak <= edata.chr$gend[i] + windows$ua &
					idata.chr$clusterorientation == 1
				)
				if (length(I) > 0) {	
					J <- findInterval(
						x = idata.chr$base[I],
						vec = c(edata.chr$gstart[i],edata.chr$gend[i]),
						rightmost.closed = TRUE
					)
					da <- c(
						(idata.chr$base[I] - edata.chr$gend[i])[J==0],
						((idata.chr$base[I] - edata.chr$gend[i]) / 
							(edata.chr$gend[i] - edata.chr$gstart[i]))[J==1],
						(idata.chr$base[I] - edata.chr$gend[i])[J==2]
					)
					IDa <- idata.chr$internal_id[I]
				} else {
					IDa <- NULL
					da <- NULL
				}
			}
			list(ID = c(IDs,IDa), d = c(ds,da))
		})
		ds.chr <- lapply(inss.chr,function(elem){elem$d})
		inss.chr <- lapply(inss.chr,function(elem){elem$ID})
		targets.chr <- unlist(lapply(1:nrow(edata.chr),function(i){rep(edata.chr$geneid[i],length(inss.chr[[i]]))}))
		inss.chr <- unlist(inss.chr)
		ds.chr <- unlist(ds.chr)
		targets <- c(targets,targets.chr)
		inss <- c(inss,inss.chr)
		ds <- c(ds,ds.chr)
		cat(sprintf("%i\n",nrow(edata.chr)))
	}
	iset <- idata[inss,]
	iset$ensid <- targets
	iset$d2gss <- ds
	iset$internal_id <- NULL
	idata$internal_id <- NULL
	notmapped <- idata[which(!(idata$id %in% iset$id)),]
	notmapped$ensid <- rep(NA,nrow(notmapped))
	notmapped$d2gss <- rep(NA,nrow(notmapped))
	iset <- rbind(iset,notmapped)
	iset <- iset[order(iset$chr,iset$base,iset$ori,iset$ensid),]
	mask <- is.na(iset$ensid)
	narows <- iset[mask,]
	iset <- iset[!mask,]

	I <- match(iset$ensid,edata$geneid)
	STRAND <- edata$strand[I]
	SENSE <- iset$ori == edata$strand[I]
	WITHIN <- abs(iset$d2gss) < 1
	UPSTREAM <- !WITHIN & sign(STRAND) == -sign(iset$d2gss)
	DOWNSTREAM <- !WITHIN & sign(STRAND) == sign(iset$d2gss)
	SA <- c("a","s")[as.numeric(SENSE)+1]
	UDI <- rep(NA,nrow(iset))
	UDI[WITHIN] <- "i"
	UDI[UPSTREAM] <- "u"
	UDI[DOWNSTREAM] <- "d"
	narows$mechanism <- rep(NA,nrow(narows))
	iset$mechanism <- sprintf("%s%s",UDI,SA)
	narows$d2gts <- rep(NA,nrow(narows))
	iset$d2gts <- rep(NA,nrow(iset))
	iset$d2gts[WITHIN] <- 1 - iset$d2gss[WITHIN]
	iset$d2gts <- iset$base - edata$gts[I]

	iset <- rbind(iset,narows)
	iset <- iset[order(iset$chr,iset$base,iset$ori,iset$ensid),]
	iset	
}

kcrbm.trans <- function(
	edata, # ensembl data
	idata, # insertion data
	rules = NULL, # specify the use of default trans_rules for various model systems
	windows = NULL, # window sizes
	scale = NULL, # GKC scale
	orientation_homogeneity = NULL # required insertional homogeneity of a cluster
) {

	# This program maps each insertion to all putative target genes.
	# The output is a data.frame with unique (insertion, gene ID) pairs for each row.
	# ===============================================================================
	#
	# INPUT
	# (Note: positions (in basepairs) are always specified relative to the start
	# position of the chromosome.
	#
	# edata (data.frame) 
	#		contains genome information from Ensembl (using 
	# 		Biomart), sorted first on chromosome, then on bp, with the 
	# 		following fields:
	# 			$chr: (numeric) chromosome number, with X and Y indicated by numbers 
	#         		(chrom X --> chrom 20; chrom Y --> chrom 21)
	# 			$gstart: (numeric) start position of gene (bp)
	# 			$gend: (numeric) end position of gene (bp)
	# 			$tstart: (numeric) start position of transcript (bp)
	# 			$tend: (numeric) end position of transcript (bp)
	# 			$gss: (numeric) orientation-relative start position of gene (bp)
	# 			$gts: (numeric) orientation-relative end position of gene (bp)
	# 			$tss: (numeric) orientation-relative start position of transcript (bp)
	# 			$tts: (numeric) orientation-relative end position of transcript (bp)
	# 			$css: (numeric) orientation-relative start position of coding region (bp)
	# 			$cts: (numeric) orientation-relative end position of coding region (bp)
	#			$strand: (numeric) on which strand lies the transcript (+1/-1)
	# 			$geneid: (character) Ensembl geneID for the transcript
	#			$transid: (character) Ensembl transcriptID for the transcript
	# 			$estart: (numeric) start position of exon (bp)
	# 			$eend: (numeric) end position of exon (bp)
	# 			$cstart: (numeric) start position of coding sequence, per exon (bp)
	# 			$cend: (numeric) end position of coding sequence, per exon (bp)
	#			$cmin: (numeric) start position of coding sequence, per transcript (bp)
	#			$cmax: (numeric) end position of coding sequence, per transcript (bp)
	#			$utr5start: (numeric) start position of 5' UTR
	#			$utr5end: (numeric) end position of 5' UTR
	#			$utr3start: (numeric) start position of 3' UTR
	#			$utr3end: (numeric) end position of 3' UTR
	#
	# idata (data.frame): 
	#		contains insertion information, first sorted on chromosome, 
	# 		then on bp, with at least the following fields for each 
	#		insertion:
	#	 		$chr (numeric): chromosome number, with X and Y indicated by numbers
	#		  	$base (numeric): insertion locus(bp)
	#		  	$ori (numeric): insertion orientation, sense or antisense (+1/-1)
	#		NB if idata has a column "id", then it has to be a unique identifier for
	#		an insertion. If no such column is present, a column "id" will be added.
	# 
	# {optional} trans_rules (character): 
	#		trans_rules to be used, one of:
	#			"MuLV" (Murine Leukaemia Virus)
	#			"SB" (Sleeping Beauty transposon)
	#			"MMTV" (Mouse Mammary Tumor Virus)
	# 		Note: If argument trans_rules is not provided, arguments windows, type, scale, and
	#		orientation_homogeneity need to be provided.
	#
	# {optional} windows (numeric): 
	# 		IF type == "windows": 
	#			4 window sizes to be used for mapping
	#		ELSE: 
	#			1 maximum window size to be used for nearest-gene or nearest-
	#			transcript mapping
	#
	# {optional} scale (numeric):
	# 		If scale > 0, then insertions in idata will be clustered by 
	#		assigning all insertions to the nearest peak, using CIMPL with 
	# 		scale = scale. The clustered loci will then be mapped according to 
	# 		the trans_rules, and the resulting targets for a cluster will then be assigned 
	# 		to the insertions in idata composing that cluster. Note: All distances 
	# 		returned will be relative to the insertion cluster means.
	#
	# {optional} orientation_homogeneity (0.5 < numeric <= 1.0):
	#		If at least orientation_homogeneity insertions in a cluster have orientation ori,
	#		then all insertions in this cluster will be treated as having orientation
	# 		ori.
	#
	# OUTPUT
	# insertions_trans (data.frame): 
	#		A data.frame with all fields in the original idata, and with additional 
	#		fields indicating for each insertion (or insertion cluster):
	#			$transid: (character) to which transcript ID(s) this insertion maps,
	#   		     	including all possible transID with all possible mechanisms
	#			$d2tss: (numeric) distance to TSS
	#			$d2tts: (numeric) distance to TTS
	#			$dfrac_trans: (numeric) fractional distance w.r.t. transcript length, 
	#				relative to TSS (only if insertion is within a transcript).
	#			$d2css: (numeric) distance to coding start
	#			$d2cts: (numeric) distance to coding end
	#			$dfrac_cds: (numeric) fractional distance w.r.t. coding region length, 
	#				relative to coding start (only if insertion is within a coding region).
	#  			$ensid: (character) to which gene ID(s) the transcript ID(s) belong
	#			$d2gss: (numeric) distance to gene start
	#			$d2gts: (numeric) distance to gene end
	#			$dfrac_gene: (numeric) fractional distance w.r.t. gene length, 
	#				relative to gene start (only if insertion is within a gene).
	#	  		$distance: (numeric) the distance between the insertion and the TSS
	#        	of the mapped transcripts(bp, not absolute value)
	#		        	(distance = insertion_locus - TSS)
	#	  		$mechanism: (cell) by which mechanism the insertion is mapped to a transcript
	#
	# NOTE: The output data.frame has row-wise unique (insertion, gene ID) pairs.
	#
	# DATE      14 April 2010
	# AUTHOR    Johann de Jong 
	# -------------------------------------------------------------------------
		
	# Sort the insertions, just to be sure
	idata <- idata[order(idata$chr,idata$base,idata$ori),]
	
	if (!("id" %in% colnames(idata))) {
		idata$id <- 1:nrow(idata)
	}
	
	if (!is.null(rules)) {
		if (rules %in% c("MuLV","MMTV")) {
			windows <- list(us=20000,ua=120000,ds=40000,da=5000)
			scale <- 10000
			orientation_homogeneity <- 0.75
		} else if (rules %in% c("SB")) {
			windows <- list(us=20000,ua=10000,ds=25000,da=5000)
			scale <- 2000
			orientation_homogeneity <- 0.75		
		}
	} else {
		if (
			is.null(windows) || 
			is.null(scale) || 
			is.null(orientation_homogeneity)
		) {
			stop(paste(
				"If argument rules is not provided,",
				"arguments windows, scale, and",
				"orientation_homogeneity need to be provided.\n"
			))
		} else {			
			if (!(
				is.numeric(orientation_homogeneity) &&
				0.5 < orientation_homogeneity && 
				orientation_homogeneity <= 1.0
			)) {
				stop("For parameter orientation_homogeneity is required: 0.5 < orientation_homogeneity <= 1.0\n")
			}
			
			if (!(
				is.numeric(scale) &&
				scale >= 0
			)) {
				stop("For parameter scale is required: scale >= 0.\n")
			}
			
			if (!(
				is.numeric(windows) &&
				length(windows) == 4
			)) {
				stop("RBM requires a windows parameter vector of length 4 (numeric).\n")
			}
			
			if (length(windows) > 1) {
				windows <- as.list(windows)
				names(windows) <- c("us","ua","ds","da")
			}			
		}
	}
	
	idata <- preprocess(idata,scale,orientation_homogeneity)
	classes <- sapply(1:ncol(idata),function(i) {class(idata[,i])})
	
	nms <- c(
		colnames(idata),
		"ensid",
		"transid",
		"d2tss",
		"d2tts",
		"dfrac_trans",
		"d2css",
		"d2cts",
		"dfrac_cds",
		"d2gss",
		"d2gts",
		"dfrac_gene",
		"mechanism",
		"gene_mechanism"
	)
	classes <- c(
		classes,
		"character",
		"character",
		"character",
		"character",
		"character",
		"character",
		"character",
		"character",
		"numeric",
		"numeric",
		"numeric",
		"character",
		"character"
	)
	names(classes) <- nms

	# lists of indices of the insertions/genes on different chromosomes	
	chrs <- sort(intersect(as.numeric(unique(edata$chr)),as.numeric(unique(idata$chr))))
	
	targets <- character(0)
	for (chr in chrs) {
		cat(sprintf("chr %s: ", chr))

		# all insertions on chromosome chr
		edata.chr <- edata[edata$chr == chr,]
		# all genes (transcripts/exons) on chromosome chr
		idata.chr <- idata[idata$chr == chr,]
		
		# for all insertions on chromosome chr ...
		targets.chr <- unlist(sapply(1:dim(idata.chr)[1], function(i) {
#			print(i)
			if (i %% 1000 == 0) {cat(sprintf("%i ",i))}
			map_insertion( # ... find the targets
				edata.chr,
				ins = idata.chr[i,,drop=FALSE],
				windows
			)
		}))
		targets <- c(targets,targets.chr)
		cat(sprintf("%i\n",nrow(idata.chr)))
	}
	
	# targets is now a vector, reshape it to a dataframe
	targets <- as.data.frame(
		t(matrix(
			targets,
			nrow = length(nms),
			ncol = length(targets) / length(nms)
		)),
		stringsAsFactors=FALSE
	)
	names(targets) <- nms
	
	# all columns in the dataframe are still of class character,
	# cast them to their original class
	for (i in 1:ncol(targets)) {
		col <- which(colnames(targets) == names(classes)[i])
		w <- withWarnings(class(targets[,col]) <- classes[i])
		cat(sprintf("Warning: casting column %s: %s\n",nms[i],w$warnings))
	}
	targets <- targets[order(targets$chr,targets$base,targets$ori,targets$ensid),]
	targets
}	

###### MAIN MAPPING FUNCTION ######

# Maps insertion 'ins' to genes 'edata.chr' within maximum windows 'windows'.
map_insertion <- function(
	edata.chr, # all TSSs within window distance from ins
	ins, # insertion
	windows
) {
	# limit the search space for target genes
	max_u_window <- max(windows$us,windows$ua)
	max_d_window <- max(windows$ds,windows$da)
	mask <- rep(NA,dim(edata.chr)[1])
	mask[edata.chr$strand == 1] <- # gene-centered
		(edata.chr$tss - max_u_window <= ins$clusterpeak &
		ins$clusterpeak <= edata.chr$tts + max_d_window)	[edata.chr$strand == 1]
	mask[edata.chr$strand == -1] <- # gene-centered
		(edata.chr$tts - max_d_window <= ins$clusterpeak &
		ins$clusterpeak <= edata.chr$tss + max_u_window)	[edata.chr$strand == -1]
	edata.window <- edata.chr[mask,,drop=FALSE]
	gene_rules <- as.data.frame(matrix(NA, nrow = nrow(edata.window), ncol = 0))
	trans_rules <- as.data.frame(matrix(NA, nrow = nrow(edata.window), ncol = 0))
	
#	edata.window <- edata.window[edata.window$transid %in% c("ENSMUST00000064976","ENSMUST00000027059"),]
	
	if (nrow(edata.window) > 0) { 
	
		# A FEW CHARACTERIZATIONS
	
		# REGARDING THE INSERTION
		# which transcripts have this insertion between their start and end position
		within_transcript <- 
			edata.window$tstart <= ins$base & ins$base <= edata.window$tend
		# which genes have this insertion between their start and end position
		within_gene <- 
			edata.window$gstart <= ins$base & ins$base <= edata.window$gend
		# which transcripts have this insertion within one of their exons
		within_exon <- 
			edata.window$estart <= ins$base & ins$base <= edata.window$eend
		II <- split(1:nrow(edata.window),edata.window$transid)
		within_intron <- rep(NA,nrow(edata.window))
		within_intron[unlist(II)] <- 
			unlist(sapply(II, function(I) { # for all exons in a transcript
				if (any(within_exon[I])) { # if this insertion is in any exon in this transcript ...
					within_transcript[I] & rep(FALSE,length(I)) # ... then it is not in an intron
				} else {
					within_transcript[I] & rep(TRUE,length(I)) # ... else it is
				}
			}))
		# which transcripts have this insertion within one of their coding regions
		within_cds <- 
			edata.window$cmin <= ins$base & ins$base <= edata.window$cmax
		within_cds[is.na(within_cds)] <- FALSE
		# distances w.r.t. coding start within this window
		edata.window$d2css <- 
			ins$base - edata.window$css
		# distances w.r.t. coding end within this window
		edata.window$d2cts <- 
			ins$base - edata.window$cts
		# distances w.r.t. TSSs within this window
		edata.window$d2tss <- 
			ins$base - edata.window$tss
		# distances w.r.t. TTSs within this window
		edata.window$d2tts <- 
			ins$base - edata.window$tts
		# distances of cluster peaks w.r.t. TSSs within this window
		edata.window$peak2tss <- 
			ins$clusterpeak - edata.window$tss
		# distances of cluster peaks w.r.t. TTSs within this window
		edata.window$peak2tts <- 
			ins$clusterpeak - edata.window$tts
		# set within transcript distances to fractional distance from TSS
		edata.window$dfrac_trans <- rep(NA,length(within_transcript))
		edata.window$dfrac_trans[within_transcript] <-
			(ins$base - edata.window$tss[within_transcript]) /
			(edata.window$tts[within_transcript] - edata.window$tss[within_transcript])
		# set within coding sequence distances to fractional distance from CSS
		edata.window$dfrac_cds <- rep(NA,length(within_cds))
		edata.window$dfrac_cds[which(within_cds)] <-
			(ins$base - edata.window$css[which(within_cds)]) /
			(edata.window$cts[which(within_cds)] - edata.window$css[which(within_cds)])
		# distances w.r.t. gene start sites within this window
		edata.window$d2gss <- 
			ins$base - edata.window$gss
		# distances w.r.t. gene termination sites within this window
		edata.window$d2gts <- 
			ins$base - edata.window$gts
		# set within transcript distances to fractional distance from genes
		edata.window$dfrac_gene <- rep(NA,length(within_gene))
		edata.window$dfrac_gene[within_gene] <-
			(ins$base - edata.window$gss[within_gene]) /
			(edata.window$gts[within_gene] - edata.window$gss[within_gene])
		# which are sense genes and which antisense w.r.t. this insertion?
		sense <- sign(edata.window$strand) == sign(ins$ori)
		antisense <- !sense
		# w.r.t. which transcripts is this insertion upstream or downstream?
		upstream <- sign(edata.window$d2tss) == -sign(edata.window$strand)
		downstream <- !upstream & ! within_transcript
		# w.r.t. which genes is this insertion upstream or downstream?
		upstream_gene <- sign(edata.window$d2gss) == -sign(edata.window$strand)
		downstream_gene <- !upstream_gene & !within_gene
		# which transcripts have this insertion before their first coding exon
		upstream_cds <- sign(edata.window$d2css) == -sign(edata.window$strand)
		upstream_cds[is.na(upstream_cds)] <- FALSE
		downstream_cds <- sign(edata.window$d2cts) == sign(edata.window$strand)
		downstream_cds[is.na(downstream_cds)] <- FALSE

		# REGARDING THE INSERTION CLUSTER PEAK
		# which transcripts have this insertion cluster peak between their start and end position
		cluster_within_transcript <- 
			edata.window$tstart <= ins$clusterpeak & ins$clusterpeak <= edata.window$tend
		# which transcripts are sense and which antisense w.r.t. this insertion cluster?
		clusterupstream <- sign(edata.window$peak2tss) == -sign(edata.window$strand)
		clusterdownstream <- !clusterupstream
		# w.r.t. which transcripts is this insertion cluster upstream or downstream?
		clustersense <- sign(edata.window$strand) == sign(ins$clusterorientation)
		clusterantisense <- !clustersense
		# which transcripts are within either one of the windows w.r.t. this insertion cluster peak?
		cluster_within_mapping_range <- 
			clusterupstream &
			clustersense &
			abs(edata.window$peak2tss) < windows$us |
			clusterupstream &
			clusterantisense &
			abs(edata.window$peak2tss) < windows$ua |
			clusterdownstream &
			clustersense &
			abs(edata.window$peak2tts) < windows$ds |
			clusterdownstream &
			clusterantisense &
			abs(edata.window$peak2tts) < windows$da | 
			cluster_within_transcript

		# SPECIFYING THE RULES
		# us -> upstream and sense
		trans_rules$us <- 
			!within_transcript &
			upstream & 
			sense & 
			cluster_within_mapping_range 
		# ua -> upstream and antisense
		trans_rules$ua <- 
			!within_transcript &
			upstream & 
			antisense & 
			cluster_within_mapping_range 
		# ds -> downstream and sense
		trans_rules$ds <- 
			!within_transcript &
			downstream & 
			sense & 
			cluster_within_mapping_range 
		# da -> downstream and antisense
		trans_rules$da <- 
			!within_transcript &
			downstream & 
			antisense & 
			cluster_within_mapping_range 
		# cds-exon-s
		trans_rules$cds_exon_s <- 
			within_transcript &
			within_cds &
			within_exon &		
			sense & 
			cluster_within_mapping_range
		# cds-exon-a
		trans_rules$cds_exon_a <- 
			within_transcript &
			within_cds &
			within_exon &		
			antisense & 
			cluster_within_mapping_range
		# cds-intron-s
		trans_rules$cds_intron_s <- 
			within_transcript &
			within_cds &		
			within_intron &		
			sense & 
			cluster_within_mapping_range
		# cds-intron-a
		trans_rules$cds_intron_a <- 
			within_transcript &
			within_cds &
			within_intron &		
			antisense & 
			cluster_within_mapping_range
		# utr5-exon-s
		trans_rules$utr5_exon_s <- 
			within_transcript &
			upstream_cds &
			within_exon &		
			sense & 
			cluster_within_mapping_range
		# utr5-exon-a
		trans_rules$utr5_exon_a <- 
			within_transcript &
			upstream_cds &
			within_exon &		
			antisense & 
			cluster_within_mapping_range
		# utr5-intron-s
		trans_rules$utr5_intron_s <- 
			within_transcript &
			upstream_cds &		
			within_intron &		
			sense & 
			cluster_within_mapping_range
		# utr5-intron-a
		trans_rules$utr5_intron_a <- 
			within_transcript &
			upstream_cds &
			within_intron &		
			antisense & 
			cluster_within_mapping_range
		# utr3-exon-s
		trans_rules$utr3_exon_s <- 
			within_transcript &
			downstream_cds &
			within_exon &		
			sense & 
			cluster_within_mapping_range
		# utr3-exon-a
		trans_rules$utr3_exon_a <- 
			within_transcript &
			downstream_cds &
			within_exon &		
			antisense & 
			cluster_within_mapping_range
		# utr3-intron-s
		trans_rules$utr3_intron_s <- 
			within_transcript &
			downstream_cds &		
			within_intron &		
			sense & 
			cluster_within_mapping_range
		# utr3-intron-a
		trans_rules$utr3_intron_a <- 
			within_transcript &
			downstream_cds &
			within_intron &		
			antisense & 
			cluster_within_mapping_range
		# within transcript - other - sense
		trans_rules$other_within_s <- 
			within_transcript & 
			sense & 
			cluster_within_mapping_range & 
			apply(trans_rules,1,function(row){!any(row)})
		# within transcript - other - antisense
		trans_rules$other_within_a <- 
			within_transcript & 
			antisense & 
			cluster_within_mapping_range & 
			apply(trans_rules,1,function(row){!any(row)})
		
		# us -> upstream and sense
		gene_rules$us <- 
			!within_gene &
			upstream_gene & 
			sense & 
			cluster_within_mapping_range 
		# ua -> upstream and antisense
		gene_rules$ua <- 
			!within_gene &
			upstream_gene & 
			antisense & 
			cluster_within_mapping_range 
		# ds -> downstream and sense
		gene_rules$ds <- 
			!within_gene &
			downstream_gene & 
			sense & 
			cluster_within_mapping_range 
		# da -> downstream and antisense
		gene_rules$da <- 
			!within_gene &
			downstream_gene & 
			antisense & 
			cluster_within_mapping_range 
		# is
		gene_rules$is <- 
			within_gene &
			sense & 
			cluster_within_mapping_range
		# ia
		gene_rules$ia <- 
			within_gene &
			antisense & 
			cluster_within_mapping_range
	}
	
	# APPLYING THE trans_rules	
	# (exons are not relevant anymore)
	apply_rules(edata.window,ins,trans_rules,gene_rules)
}

###### OTHER FUNCTIONS ######

apply_rules <- function(edata.window,ins,trans_rules,gene_rules) {
	ncols <- 13
	ins <- unlist(ins)
	if (nrow(edata.window) == 0) {
		return(c(ins,rep(NA,ncols)))
	} else {
		# remove empty rows (they are no targets)
		mask <- apply(trans_rules,1,any)
		edata.target <- edata.window[mask,]
		trans_rules <- trans_rules[mask,]
		gene_rules <- gene_rules[mask,]
		if (nrow(edata.target) == 0) {
			return(c(ins,rep(NA,ncols)))			
		} else {
	
			# remove non unique rows
			I <- sapply(split(1:nrow(edata.target),edata.target$transid), function(I){I[1]})
			edata.target <- edata.target[I,]
			trans_rules <- trans_rules[I,]
			gene_rules <- gene_rules[I,]
	
			II <- split(1:nrow(edata.target),edata.target$geneid)
			m <- sapply(II, function(I) {
				c(
					edata.target$geneid[I[1]],
					paste(edata.target$transid[I],collapse="|"),
					paste(edata.target$d2tss[I],collapse="|"),
					paste(edata.target$d2tts[I],collapse="|"),
					paste(edata.target$dfrac_trans[I],collapse="|"),
					paste(edata.target$d2css[I],collapse="|"),
					paste(edata.target$d2cts[I],collapse="|"),
					paste(edata.target$dfrac_cds[I],collapse="|"),
					edata.target$d2gss[I[1]],
					edata.target$d2gts[I[1]],
					paste(edata.target$dfrac_gene[I[1]],collapse="|"),
					paste(sapply(I,function(i){
						colnames(trans_rules)[as.logical(trans_rules[i,])]}),collapse="|"),
					colnames(gene_rules)[as.logical(gene_rules[I[1],])]
				)
			})
			nrows <- length(m)/ncols
			m <- rbind(
				matrix(rep(ins,nrows),nrow=length(ins),ncol=nrows),
				m
			)
			return(m)
		} 
	}
}

# evaluates an expression and collects possible warnings
withWarnings <- function(expression) {
	warnings <- NULL
	warningHandler <- function(warning) {
		warnings <<- c(warnings, list(warning))
		invokeRestart("muffleWarning")
	}
	list(
		value = withCallingHandlers(expression, warning = warningHandler), 
		warnings = warnings
	)
} 

preprocess <- function(idata,scale,orientation_homogeneity) {
	if (scale == 0) {
		idata$clusterpeak <- idata$base
		idata$clusterorientation <- idata$ori
		idata$ids <- idata$id
		idata$n_insertions <- rep(1,nrow(idata))
	} else {
	
		library("BSgenome.Mmusculus.UCSC.mm9")
		library(cimpl)

		names(idata)[which(names(idata) == "base")] <- "location"
		idata$chr <- sprintf("%i",idata$chr)
		idata$chr[idata$chr == "20"] <- "X"
		idata$chr[idata$chr == "21"] <- "Y"
		idata$chr <- sprintf("chr%s",idata$chr)
		cimpl <- doCimplAnalysis(
			data = idata,
			scales = scale,
			n_iterations = 1, # 1000
			system = "MuLV",
			BSgenome = Mmusculus,
			verbose = FALSE
		)
		names(idata)[which(names(idata) == "location")] <- "base"
		idata$chr <- gsub("chr","",idata$chr)
		idata$chr[idata$chr == "X"] <- "20"
		idata$chr[idata$chr == "Y"] <- "21"
		idata$chr <- as.numeric(idata$chr)
	
		cimpl@chromosomes <- gsub("chr","",cimpl@chromosomes)
		cimpl@chromosomes[cimpl@chromosomes == "X"] <- "20"
		cimpl@chromosomes[cimpl@chromosomes == "Y"] <- "21"
	
		ncol <- 4
		idata_cimpl <- matrix(NA,nrow=0,ncol=ncol)
		colnames(idata_cimpl) <- c("clusterpeak","clusterorientation","ids","n_insertions")
		for (i in 1:length(unique(idata$chr))) {
			chr <- unique(idata$chr)[i]
			idata.chr <- idata[idata$chr == chr,]
			chr.idx <- which(cimpl@chromosomes == as.character(chr))
			kw.idx <- which(cimpl@scales == scale)
			co <- cimpl@cimplObjects[[chr.idx]][[kw.idx]]
			peaks <- sapply(1:dim(idata.chr)[1],function(j){ # map each insertion to the nearest peak
				which.min(abs(co@peaks$x - idata.chr$base[j]))
			})
			L <- split(1:length(peaks),peaks)
			idata_cimpl.chr <- unlist(sapply(L,function(l) {
				sense <- which(idata.chr$ori[l] == 1)
				antisense <- which(idata.chr$ori[l] == -1)
				sense_fraction <- length(sense) / length(l)
				antisense_fraction <- length(antisense) / length(l)
				ori <- if (sense_fraction >= orientation_homogeneity) { # should we cluster the insertions as sense?
					peaklocus <- mean(idata.chr$base[l])
					rep(c(
						peaklocus,
						1,
						paste(idata.chr$id[l],collapse="|"),
						length(l)
					),length(l))
				} else if (antisense_fraction >= orientation_homogeneity) { # should we cluster the insertions as antisense?
					peaklocus <- mean(idata.chr$base[l])
					rep(c(
						peaklocus,
						-1,
						paste(idata.chr$id[l],collapse="|"),
						length(l)
					),length(l))
				} else { # cluster is orientation-wise too heterogeneous: cluster orientations separately
					peaklocus_sense <- mean(idata.chr$base[l[sense]])
					peaklocus_antisense <- mean(idata.chr$base[l[antisense]])
					lsense <- length(sense)
					lantisense <- length(antisense)
					res <- rep(NA,length(l)*ncol)
					res[(sense-1)*ncol+1] <- rep(peaklocus_sense,lsense)
					res[(sense-1)*ncol+2] <- rep(1,lsense)
					res[(sense-1)*ncol+3] <- rep(paste(idata.chr$id[l[sense]],collapse="|"),lsense)
					res[(sense-1)*ncol+4] <- rep(lsense,lsense)
					res[(antisense-1)*ncol+1] <- rep(peaklocus_antisense,lantisense)
					res[(antisense-1)*ncol+2] <- rep(-1,lantisense)
					res[(antisense-1)*ncol+3] <- rep(paste(idata.chr$id[l[antisense]],collapse="|"),lantisense)
					res[(antisense-1)*ncol+4] <- rep(lantisense,lantisense)
					res
				}
			}))
			idata_cimpl.chr <- t(matrix(idata_cimpl.chr,nrow=ncol,ncol=length(idata_cimpl.chr)/ncol))
			idata_cimpl <- rbind(idata_cimpl,idata_cimpl.chr)
		}
		idata_cimpl <- as.data.frame(idata_cimpl,stringsAsFactors=FALSE)
		idata_cimpl[,1] <- as.numeric(idata_cimpl[,1])
		idata_cimpl[,2] <- round(as.numeric(idata_cimpl[,2]))
		idata_cimpl[,4] <- as.numeric(idata_cimpl[,4])
		idata <- cbind(idata,idata_cimpl)
	
		# first sort on chromosome number
		idata <- idata[sort(idata$chr,index.return=TRUE)$ix,]

		# then sort on start bp no. for each chromosome separately
		idata$base <- as.numeric(idata$base)
		for (chr in 1:21) {
			mask <- idata$chr == chr
			idata[mask,] <- idata[mask,][sort(as.numeric(idata$base[mask]),index.return=TRUE)$ix,]
		}
		idata$id <- 1:(dim(idata)[1])
	}
	idata
}

