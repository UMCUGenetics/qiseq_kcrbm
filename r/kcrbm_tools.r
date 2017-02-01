# DESCRIPTION: Set of KC-RBM utility functions for pre- and postprocessing
# DEPENDS: biomaRt, BSgenome package for genome of choice
# ===============================================================================
# DATE      14 April 2010
# AUTHOR    Johann de Jong 
# LICENSE	GPL version 3 or newer
# -------------------------------------------------------------------------------




# INPUT - A data.frame ensembl_data with the following columns
# chr			(chromosome: strictly numeric, e.g. mmusculus Y = 20 and X = 21)
# gstart			(gene start)
# gend			(gene end)
# tstart			(transcript start)
# tend			(transcript end)
# strand			(strand)
# geneid			(ensembl gene id)
# transid		(ensembl transcript id)
# estart			(exon start)
# eend			(exon end)
# cstart			(coding start)
# cend 			(coding end)
# utr5start 		(utr5 start)
# utr5end 		(utr5 end)
# utr3start 		(utr3 start)
# utr3end 		(utr3 end)
#
# If ensembl_data is not provided, data will be downloaded through biomart, and the argument
# 'dataset' should be provided.
#
# OUTPUT - ensembl_data data.frame for input to kcrbm.r
preprocessEnsemblData <- function(ensembl_data,dataset="mmusculus") {
	
	if (missing(ensembl_data)) {
		# retrieve biomart data
		getloci <- function() {
			library(biomaRt)
			mart <- useMart(
				"ENSEMBL_MART_ENSEMBL", 
				dataset = sprintf("%s_gene_ensembl",dataset), 
				host = "www.ensembl.org"
			)
			cols <- c("chr", "tstart", "tend", "strand", "geneid", "transid", "estart","eend", 
				"cstart", "cend","utr5start","utr5end","utr3start","utr3end","gstart","gend")
			attributes <- c("chromosome_name","transcript_start","transcript_end","strand",
				"ensembl_gene_id","ensembl_transcript_id","exon_chrom_start","exon_chrom_end",
				"cds_start","cds_end","5_utr_start","5_utr_end","3_utr_start","3_utr_end",
				"start_position","end_position")
			base <- which(attributes == "exon_chrom_start")
			loci <- getBM(
				attributes = attributes, 
				mart = mart
			)
			w <- getOption("warn")
			options(warn = -1) # suppress warnings
			colnames(loci) <- cols
			chr <- max(as.numeric(loci$chr),na.rm=TRUE)
			loci$chr <- gsub("Y",chr+2,gsub("X",chr+1,loci$chr))
			loci$chr <- as.numeric(loci$chr)
			loci <- loci[!is.na(loci$chr),]
			options(warn = w)

			# sort
			loci <- loci[order(loci$chr,loci$estart),]
			loci
		}
		ensembl_data <- getloci()
	}

	ensembl_data <- ensembl_data[order(ensembl_data$chr,ensembl_data$estart),]
	
	# add fields tss and tts (for the transcription start and termination sites)
	mask <- ensembl_data$strand == 1
	ensembl_data$tss[mask] <- ensembl_data$tstart[mask]
	ensembl_data$tss[!mask] <- ensembl_data$tend[!mask]
	ensembl_data$tts[mask] <- ensembl_data$tend[mask]
	ensembl_data$tts[!mask] <- ensembl_data$tstart[!mask]
	ensembl_data$gss[mask] <- ensembl_data$gstart[mask]
	ensembl_data$gss[!mask] <- ensembl_data$gend[!mask]
	ensembl_data$gts[mask] <- ensembl_data$gend[mask]
	ensembl_data$gts[!mask] <- ensembl_data$gstart[!mask]
	ensembl_data$ess[mask] <- ensembl_data$estart[mask]
	ensembl_data$ess[!mask] <- ensembl_data$eend[!mask]
	ensembl_data$ets[mask] <- ensembl_data$eend[mask]
	ensembl_data$ets[!mask] <- ensembl_data$estart[!mask]
	
	# add fields for minimum cds start and maximum cds end
	w <- getOption("warn")
	options(warn = -1) # suppress warnings
	I <- split(1:length(ensembl_data$transid),ensembl_data$transid)
	ensembl_data$cmin[unlist(I)] <- unlist(sapply(I,function(i) {
		rep(min(ensembl_data$cstart[i],na.rm=TRUE),length(i))
	}))
	ensembl_data$cmin[is.infinite(ensembl_data$cmin)] <- NA
	I <- split(1:length(ensembl_data$transid),ensembl_data$transid)
	ensembl_data$cmax[unlist(I)] <- unlist(sapply(I,function(i) {
		rep(max(ensembl_data$cend[i],na.rm=TRUE),length(i))
	}))
	ensembl_data$cmax[is.infinite(ensembl_data$cmax)] <- NA
	options(warn = w)
		
	ensembl_data$css[mask] <- ensembl_data$cmin[mask]
	ensembl_data$css[!mask] <- ensembl_data$cmax[!mask]
	ensembl_data$cts[mask] <- ensembl_data$cmax[mask]
	ensembl_data$cts[!mask] <- ensembl_data$cmin[!mask]

	ensembl_data
}

# get the lengths of the chromosomes
get_chrl <- function(genome="mm9") {
	if (genome == "mm9") {
		library("BSgenome.Mmusculus.UCSC.mm9")
		GENOME <- Mmusculus
	}
	chrs <- 1:which(seqnames(GENOME) == "chrY")
	as.numeric(seqlengths(GENOME)[chrs])
}


# retrieve CTGs from kcrbm output
get_ctgs <- function(iset) {
	iset <- iset[!is.na(iset$ensid),]
	I <- sort(iset$ensid,index.return=TRUE)$ix
	iset <- iset[I,]
	if (!("weight" %in% colnames(iset))) {
		iset$weight <- rep(1,nrow(iset))
	}
	scores <- aggregate(iset$weight,list(iset$ensid),FUN=sum)
	ensids <- scores[,1]
	scores <- scores[,2]
	names(scores) <- ensids
	scores <- sort(scores,decreasing=TRUE)
}

# return a barplot of kcrbm CTG output
plot_ctgs <- function(
	ctgs, 
	file = NULL,
	main = "Most commonly targeted genes",
	ylab = "no. of times targeted",
	cex.lab = 1.2,
	cex.axis = 1.0
) {
	if (any(grepl("ENSMUSG",names(ctgs)))) {
		names(ctgs) <- ens2sym(names(ctgs))
	}

	if (!is.null(file)) {
		pdf(file, width = 7, height = 5)
	} 
	dmar <- 3
	margins <- par("mar")
	margins[1] <- margins[1] + dmar
	par(mar=margins)
	b <- barplot(
		ctgs,
		main = main,
		names.arg = names(ctgs),
		ylab = ylab,
		las = 2,
		legend.text = FALSE,
		cex.lab = cex.lab,
		cex.axis = cex.axis
	)
	margins[1] <- margins[1] - dmar
	par(mar=margins)
	if (!is.null(file)) {
		dev.off()
	} else {
		b
	}
}

getEnsemblGenes <- function() {
	library(biomaRt)
	getBM(
		attributes = c(
			"external_gene_name", 
			"ensembl_gene_id"
		),
#		mart = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
		mart = useMart("ENSEMBL_MART_ENSEMBL", dataset = "mmusculus_gene_ensembl", host = "www.ensembl.org")
	)
}

# generate a random set of insertions
random_idata <- function(ninsertions,chrl) {
	iset <- lapply(as.numeric(names(chrl)), function(chr) {
		n <- round(ninsertions*chrl[chr]/sum(chrl))
		data.frame(
			chr = rep(chr,n),
			base = sample(1:chrl[chr],n,replace=TRUE),
			ori = sample(c(-1,1),n,replace=TRUE)
		)
	})
	iset <- do.call("rbind",iset)
	iset
}

nulldist <- function(idata,edata,N,windows,scale,orientation_homogeneity,type,GENOME) {
	chrl <- get_chrl(GENOME)
	ensids <- sort(unique(edata$geneid))
	tmain <- Sys.time()
	CTGS <- sapply(1:N, function(n) {
		cat(sprintf("--- Permutation %i of %i ---\n", n, N))
		t <- Sys.time()
		iset <- fast.kcrbm(
			edata = edata,
			idata = random_idata(nrow(idata),chrl),
			windows = windows,
			scale = scale,
			orientation_homogeneity = orientation_homogeneity
		)
		ctgs <- get_ctgs(selectSingleTargets(iset,type=type))
		nottargeted <- ensids[!(ensids %in% names(ctgs))]
		nms <- c(names(ctgs),nottargeted)
		ctgs <- c(ctgs,rep(0,length(nottargeted)))
		names(ctgs) <- nms
		ctgs <- ctgs[order(names(ctgs))]
		cat(sprintf("Time: %g secs\n",difftime(Sys.time(),t,units="secs")))
		ctgs
	})
	cat(sprintf("Total time: %g secs\n",difftime(Sys.time(),tmain,units="secs")))
	cat("\n")
	rownames(CTGS) <- ensids
	colnames(CTGS) <- 1:N
	CTGS
}

# From kcrbm output, select single targets for each insertion.
# This implements ngKC-RBM and ctgKC-RBM as a post-processing step of KC-RBM.
selectSingleTargets <- function(
	iset,
	type = NULL
) {
	iset <- iset[!is.na(iset$ensid),]
	if (!is.null(type) && type == "nearest") {
		# select only the nearest
		iset$id <- as.numeric(iset$id)
		if ("d2tss" %in% colnames(iset)) {
			ds <- strsplit(iset$d2tss,"\\|")
			ds <- sapply(ds, function(d) {median(abs(as.numeric(d)))})
		
			I <- unlist(sapply(split(1:length(iset$id),iset$id), function(ID) {
				within_trans <- !grepl("NA",iset$dfrac_trans[ID])
				if (any(within_trans)) {
					ID[within_trans]
				} else {
					ID[which.min(ds[ID])]
				}
			}))
		} else { # assume only d2gss is present (as generated by function kcrbm.genes)
			I <- unlist(sapply(split(1:length(iset$id),iset$id), function(ID) {
				within_gene <- 0 <= iset$d2gss[ID] & iset$d2gss[ID] <= 1
				if (any(within_gene)) {
					ID[within_gene]
				} else {
					ID[which.min(abs(iset$d2gss[ID]))]
				}
			}))
		}
		
		cnames <- colnames(iset)
		rnames <- rownames(iset)
		iset <- iset[I,]
		colnames(iset) <- cnames
		rownames(iset) <- rnames[I]
	} else if (!is.null(type) && type == "ctg") {
		# select only the most commonly targeted
		
		# make a list of CTGs
		ensidCol <- which(colnames(iset) == "ensid")
		ctgs <- sort(table(iset$ensid,useNA="no"),decreasing=TRUE)
		
		# make the selection
		iset$id <- as.numeric(iset$id)
		I <- unlist(sapply(split(1:length(iset$id),iset$id), function(ID) {
			if (length(ID) > 1) {
				within_trans <- !grepl("NA",iset$dfrac_trans[ID])
				if (any(within_trans)) {
					ID[within_trans]
				} else {
					ID[which.max(ctgs[match(iset$ensid[ID],names(ctgs))])]
				}
			} else {
				ID
			}
		}))
		
		cnames <- colnames(iset)
		rnames <- rownames(iset)
		iset <- iset[I,]
		colnames(iset) <- cnames
		rownames(iset) <- rnames[I]
		
	}
	
	iset
}

# find nearest gene
ngm <- function(idata,edata,type="gss") {
	edata$gss <- rep(NA,nrow(edata))
	mask <- edata$start < edata$end
	edata$gss[mask] <- edata$start[mask]
	mask <- edata$start > edata$end
	edata$gss[mask] <- edata$end[mask]
	edata$gts <- rep(NA,nrow(edata))
	mask <- edata$start < edata$end
	edata$gts[mask] <- edata$end[mask]
	mask <- edata$start > edata$end
	edata$gts[mask] <- edata$start[mask]
	if (type == "gss") { # take the nearest gene start site
		ensids <- rep(NA,nrow(idata))
		I <- split(1:nrow(idata),idata$chr)
		chrs <- as.numeric(names(I))
		for (i in 1:length(I)) {
			edata.chr <- edata[edata$chr == chrs[i],]
			ensids[I[[i]]] <- sapply(I[[i]], function(j) {
				paste(edata.chr$ensid[which.min(abs(edata.chr$gss - idata$base[j]))],collapse="|")
			})
		}
	} else if (type == "ngm") { # take the nearest gene
		ensids <- rep(NA,nrow(idata))
		I <- split(1:nrow(idata),idata$chr)
		chrs <- as.numeric(names(I))
		for (i in 1:length(I)) {
			edata.chr <- edata[edata$chr == chrs[i],]
			ensids[I[[i]]] <- sapply(I[[i]], function(j) {
				mask <- edata.chr$start < idata$base[j] & idata$base[j] <= edata.chr$end
				if (any(mask)) {
					paste(edata.chr$ensid[mask],collapse="|")
				} else {
					v <- abs(c(edata.chr$gss - idata$base[j],edata.chr$gts - idata$base[j]))
					paste(edata.chr$ensid[((which.min(v) - 1) %% nrow(edata.chr)) + 1],collapse="|")
				}
			})
		}
	}
	idata$ensid <- ensids
	expand(idata,which(colnames(idata)=="ensid"),"\\|")
}

expand <- function(data, colno, sep) {
	DF <- is.data.frame(data)
	data[data == ""] <- NA

	L <- sapply(strsplit(data[,colno[1]],sep),function(s) {max(length(s),1)})
	I <- unlist(sapply(1:length(L),function(i){rep(i,L[i])}))
	
	newdata <- matrix(NA,ncol=ncol(data),nrow=length(I))
	for (i in 1:ncol(data)) {
		if (i %in% colno) {
			newdata[,i] <- unlist(strsplit(data[,i],sep))
		} else {
			newdata[,i] <- data[I,i]
		}
	}

	if (DF) {
		newdata <- as.data.frame(newdata,stringsAsFactors=FALSE)
		for (i in 1:dim(data)[2]) {
			newdata[,i] <- as(newdata[,i], class(data[,i]))
		}
	}
	colnames(newdata) <- colnames(data)
	newdata
}

ens2sym <- function(ensids) {
	mapper <- getEnsemblGenes()
	v <- unlist(sapply(ensids, function(ensid) {
		paste(unique(mapper$external_gene_name[ensid == mapper$ensembl_gene_id]), collapse = "|")
	}))
	v[v == "NA"] <- NA
	v
}

sym2ens <- function(symbols) {
	mapper <- getEnsemblGenes()
	v <- unlist(sapply(symbols, function(symbol) {
		paste(unique(mapper$ensembl_gene_id[symbol == mapper$external_gene_name]), collapse = "|")
	}))
	v[v == "NA"] <- NA
	v
}





