condenseType.single <- function(sequences) {
	
	getNewIDs <- function(.unique.list) {
		
		seqType <- strsplit(.unique.list[[1]][1], '-')[[1]][1]
		variantMatchIDs <- .unique.list %>%
			sapply(function(x) strsplit(x, '-')[1:length(x)]) %>%
			sapply(function(x) sapply(x, "[[", 2)) %>%
			sapply(function(x) paste0(x, sep='', collapse=',')) %>%
			sapply(function(x) paste(seqType, x, sep='-')) %>%
			unname()
		
		variantCounts <- unname(sapply(.unique.list, length))
		
		newIDs <- paste(variantCounts, variantMatchIDs, sep = '#_')
		
		newIDs
		
	}	
	
	tmp <- apply(sequences, 1, function(x) strsplit(x[2], '')[[1]]) %>%
		t()
	rownames(tmp) <- sequences$seqID
	
	seqs.distMat <- tmp %>%
		ape::as.DNAbin() %>%
		ape::dist.dna(model = 'raw', pairwise = T, as.matrix = T)
		
	avg.len <- mean(sapply(sequences$seqs, function(x) stringr::str_count(x) - stringr::str_count(x, pattern = '-')))
	
	seqs.distMat <- seqs.distMat * avg.len
	
	if(sum(seqs.distMat) == 0) {
		seqType <- strsplit(sequences$seqID[1], '-')[[1]][1]
		sequenceIDs <- sequences$seqID %>%
			sapply(strsplit, '-') %>%
			sapply('[[', 2) 
		newID <- paste0(nrow(sequences), '#_', seqType, '-', paste(sequenceIDs, sep='', collapse=','))
		
		data.frame(
			seqID = newID,
			seqs = sequences$seqs[1],
			stringsAsFactors = FALSE
		)
		
	} else {
		variantSets <- apply(seqs.distMat, 1, function(x) colnames(seqs.distMat)[which(x == 0)])
		unique.list <- unique(variantSets)
		if(length(unique.list) == nrow(sequences)) {
			newIDs <- paste0("1#_", unique.list)
		} else {	
			newIDs <- getNewIDs(unique.list)
		}
		newSeqs <- unique.list %>%
			sapply("[[", 1) %>%
			sapply(function(x) sequences$seqs[sequences$seqID %in% x])
	
		data.frame(
			seqID = newIDs,
			seqs = newSeqs,
			stringsAsFactors = FALSE
		)

	}
	
}

combineCondense.multi <- function(sequences, bind = FALSE) {
	if(!bind) seq_set <- dplyr::bind_rows(sequences)
	tmp <- apply(sequences, 1, function(x) strsplit(x[2], '')[[1]]) %>%
		t()
	rownames(tmp) <- sequences$seqID
	
	seqs.distMat <- tmp %>%
		ape::as.DNAbin() %>%
		ape::dist.dna(model = 'raw', pairwise = T, as.matrix = T)
		
	avg.len <- mean(sapply(sequences$seqs, function(x) stringr::str_count(x) - stringr::str_count(x, pattern = '-')))
	
	seqs.distMat <- seqs.distMat * avg.len
		
	getNewIDs <- function(.unique.list) {
	
		varCounts <- .unique.list %>%
			strsplit('#') %>%
			sapply('[[', 1) %>%
			as.numeric()
		if(length(varCounts) == 1) { return (.unique.list) }
		
		totalCount <- sum(varCounts)
		
		newNames <- .unique.list %>%
			paste(collapse = '+', sep = '') %>%
			paste(sum(varCounts), ., sep = '#_') %>%
			unname()
			
		newNames
	}
		
	unique.list <- seqs.distMat %>%
		apply(1, function(x) colnames(seqs.distMat)[which(x == 0)]) %>%
		unique()
		
	newIDs <- unique.list %>%
		sapply(getNewIDs)
	
	newSeqs <- unique.list %>%
		sapply("[[", 1) %>%
		sapply(function(x) sequences$seqs[sequences$seqID %in% x])
	
	data.frame(
		seqID = newIDs,
		seqs = newSeqs,
		stringsAsFactors = FALSE
	)
	
 }

splitMultifPECS <- function(seq_set) {

	car_groups <- seq_set$seqID %>%
		strsplit('-') %>%
		sapply(function(x) x[1]) %>%
    		lapply(function(y) unlist(strsplit(y, '\\.'))) %>%
    		sapply(function(z) paste(z[1:(length(z) - 1)], sep = '', collapse='.'))
	seq_set$factorID <- as.factor(car_groups)
	
	seq_set.split <- split(seq_set[, 1:2], seq_set$factorID)
	seq_set.split
	
}

defuzz.CAR <- function(seq_set, isCondensed = TRUE, sizeLim = 5)
{
	if(!isCondensed) seq_set <- condenseType.single(seqs)
	
	varCounts <- seq_set$seqID %>%
		strsplit('#') %>%
		sapply('[[', 1) %>%
		as.numeric()
		
	if(!(max(varCounts) >= sizeLim) & (min(varCounts) == 1)) {
		return (seq_set)
	}
	
	tmp <- apply(seq_set, 1, function(x) strsplit(x[2], '')[[1]]) %>%
		t()
	rownames(tmp) <- seq_set$seqID
	
	seqs.distMat <- tmp %>%
		ape::as.DNAbin() %>%
		ape::dist.dna(model = 'raw', pairwise = T, as.matrix = T)
		
	avg.len <- mean(sapply(seq_set$seqs, function(x) stringr::str_count(x) - stringr::str_count(x, pattern = '-')))
	
	seqs.distMat <- round(seqs.distMat * avg.len, 0)
	
	seqs.distMat <- seqs.distMat[!(varCounts < sizeLim), !(varCounts > 1)]
	
	if(sum(varCounts >= sizeLim) == 1) {
		seqs.distMat <- seqs.distMat %>% t()
		rownames(seqs.distMat) <- seq_set$seqID[which(varCounts >= sizeLim)]
	}
	
	fuzz.list <- seqs.distMat %>%
		apply(1, function(x) colnames(seqs.distMat)[which(x == 1)])
		
	if(rlang::is_empty(fuzz.list)) { return( seq_set ) }
	
	fuzz.list.edit <- fuzz.list %>%
		unlist() %>%
		unname()
	
	seq_set <- seq_set[!(seq_set$seqID %in% fuzz.list.edit), ]

	
	 getFuzzIDs <- function(x) {
		tmp <- sapply(sapply(x, function(i) strsplit(i, '-')), '[', 2) %>%
			unname()
		tmp
	}
	fuzz.ids <- sapply(fuzz.list, getFuzzIDs)
	
	
	if(is.null(names(fuzz.list))) {
		tmp.fuzz.names <- data.frame(
			base = colnames(fuzz.list),
			fuzz = paste(paste0(fuzz.ids, 'f'), collapse = ',', sep =''),
			stringsAsFactors = F,
			row.names = NULL
		)
	} else {
		tmp.fuzz.names <- data.frame(
			base = names(fuzz.ids),
			fuzz = sapply(fuzz.ids, function(x) paste(paste0(x, 'f'), collapse=',', sep='')),
			stringsAsFactors = F,
			row.names = NULL
		)
	}
	
	
	getFuzzNames <- function(fuzzIDs.df) {

		baseName <- fuzzIDs.df[1]
		baseSuf <- strsplit(baseName, '#')[[1]][2]
		baseCount <- as.numeric(strsplit(baseName, '#')[[1]][1])
		
		newCount <- baseCount + length(strsplit(fuzzIDs.df[2], ',')[[1]])
		newBase <- paste0(newCount, '#', baseSuf)
		tmp <- c(newBase, fuzzIDs.df[2])
		
		newID <- paste(tmp, collapse = ',', sep='')
		
		newID

	}
	
	tmp.fuzz.names$newNames <- apply(tmp.fuzz.names, 1, getFuzzNames)
	
	seq_set$seqID[seq_set$seqID %in% tmp.fuzz.names$base] <- tmp.fuzz.names$newNames
		
	seq_set
}

fPECS <- function(condSeqs, iCAD, cellsPerAli, totalAli = 8, iCAD.denom = 1000000)
{

	varCounts = lapply(condSeqs, function(x) x$seqID %>%
					strsplit('#') %>%
					sapply('[[', 1) %>%
					as.numeric()
					)
	cellCounts = sapply(varCounts, length)
	
	numInf = iCAD/iCAD.denom * cellsPerAli
	fPVE = cellCounts / numInf
	
	avg.fPVE = sum(fPVE) / totalAli
	sd.fPVE = sqrt((1/(totalAli - 1)) * (sum((fPVE - avg.fPVE)^2) + (totalAli - length(fPVE)) * avg.fPVE^2))
	iCAD = iCAD
	cellsPerAli = cellsPerAli
	totalAli = totalAli
	.iCAD.denom = iCAD.denom

	value <- list(
		varCounts = varCounts,
		cellCounts = cellCounts,
		numInf = numInf,
		fPVE = fPVE,
		avg.fPVE = avg.fPVE,
		sd.fPVE = sd.fPVE,
		totalAli = totalAli,
		iCAD = iCAD,
		cellsPerAli = cellsPerAli,
		.iCAD.denom = iCAD.denom
	)
	
	attr(value, "class") <- "fPECS"
	
	value
	
}


multifPECS <- function(seqs.list, iCADs, cellsPerAli, totalAlis, iCAD.denom) {
	seq_along(seqs.list) %>%
		lapply(
			function(x) fPECS(
				seqs.list[[x]],
				iCAD = iCADs[x], 
				cellsPerAli = cellsPerAli[x],
				totalAli = totalAlis[x],
				iCAD.denom = iCAD.denom[x]
			)
		) -> multifPECS.list
		
	names(multifPECS.list) <- names(seqs.list)
	multifPECS.df <- data.frame(
		avg.fPVE.list = sapply(multifPECS.list, function(x) x$avg.fPVE),
		sd.fPVE.list = sapply(multifPECS.list, function(x) x$sd.fPVE),
		alis.fPVE.list = sapply(multifPECS.list, function(x) x$totalAli)
	)	
		
	fPECS.glb.avg <- sum(multifPECS.df$alis.fPVE.list * multifPECS.df$avg.fPVE.list / multifPECS.df$sd.fPVE.list^2) /
		sum(multifPECS.df$alis.fPVE.list / multifPECS.df$sd.fPVE.list^2)

	fPECS.glb.sd <- 1 / sqrt(sum(multifPECS.df$alis.fPVE.list / multifPECS.df$sd.fPVE.list^2))

	value <- list(
		fPECStypes = multifPECS.list,
		glb.avg = fPECS.glb.avg,
		glb.sd = fPECS.glb.sd		
	)
	
	attr(value, "class") <- 'multifPECS'
	
	value
	
}



multifPECS.condenseType.single <- function(seqs.list) {
	lapply(seqs.list, function(x) lapply(x, condenseType.single))
}

multifPECS.defuzz.CAR <- function(seqs.list) {
	lapply(seqs.list, function(x) lapply(x, defuzz.CAR))
}



print.fPECS <- function(x, ...) {

	cat("Avg fPECS", signif(x$avg.fPVE, 3), "+/-", signif(x$sd.fPVE, 3),
		"in CAR-SGS experiment with", round(x$numInf, 0), 
		"infected cells per aliquot in", x$totalAli, "aliquots.")
	
}

print.multifPECS <- function(x, ...) {
	print(x$fPECStypes)
	
	cat("Global Average for multi-type fPECS", signif(x$glb.avg, 3), "+/-", signif(x$glb.sd))
}

# writefPECSsheet <- function(fPECS, he.cutoff = 20) {
	
# }

readSeqs <- function(inPath, .format="fasta") {
	s <- Biostrings::readDNAStringSet(inPath, format=.format) %>%
			as.data.frame() %>%
			tibble::rownames_to_column(var = "seqID") %>%
			dplyr::rename(seqs = x)
			
	s
			
	
}

sepFas <- function(seq_set) {

	seq_types <- seq_set$seqID %>%
		strsplit('-') %>%
		sapply(function(x) x[1])%>%
		unique()
		
	factor_IDs <- seq_set$seqID %>%
		sapply(function(x) strsplit(x, '-')) %>%
		sapply(function(x) which(seq_types %in% x[1])) %>%
		unname()
		
	seq_set$factorID <- seq_types[factor_IDs]
	
	set.Separated <- split(seq_set[,1:2], as.factor(seq_set$factorID))
	
	set.Separated
	
}

writeSeqs <- function(seq_set, outpath, writeNames = FALSE) {

	seq_set$seqs %>%
		as("XStringSet") %>%
		setNames(seq_set$seqID) %>%
		Biostrings::writeXStringSet(filepath=outpath)
	
	if(writeNames) {
		secondaryOut <- paste0(outpath, "_names.txt")
		write(seq_set$seqID, secondaryOut)
	}
	nrow(seq_set)
	
}
