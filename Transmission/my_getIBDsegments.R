my_getIBDsegments = function (ped.genotypes, parameters, number.cores = 1, minimum.snps = 20, 
          minimum.length.bp = 50000, error = 0.001) 
{
  library(Rcpp)
  library(foreach)
  library(doParallel)
  #library(isoRelate)
  sourceCpp("/nfsscratch/Karim/DELGEME/isoRelate/Codes/isolatePairs.cpp")
  #sourceCpp("/media/Data/Data/Documents_Karim/DELGEME/Deconvolution/Codes/calculateViterbi.cpp")
  #source('/media/Data/Data/Documents_Karim/DELGEME/Deconvolution/Codes/IBDTable.R')
  
  if (!is.list(ped.genotypes) | length(ped.genotypes) != 2) 
    stop("'ped.genotypes' must be a named list containing 2 objects: 'pedigree' and 'genotypes'")
  if (any(names(ped.genotypes) != c("pedigree", "genotypes"))) 
    stop("'ped.genotypes' must be a named list containing 'pedigree' and 'genotypes'")
  pedigree <- ped.genotypes[["pedigree"]]
  genotypes <- ped.genotypes[["genotypes"]]
  if (!is.data.frame(pedigree)) 
    stop("'ped.genotypes' has incorrect format - 'pedigree' is not a data.frame")
  if (!is.data.frame(genotypes)) 
    stop("'ped.genotypes' has incorrect format - 'genotypes' is not a data.frame")
  if (ncol(pedigree) != 6) 
    stop("'ped.genotypes' has incorrect format - 'pedigree' must have 6 columns: fid, iid, pid, mid, moi and aff")
  if (any(colnames(pedigree) != c("fid", "iid", "pid", "mid", 
                                  "moi", "aff"))) 
    stop("'ped.genotypes' has incorrect format - 'pedigree' must have columns labelled: fid, iid, pid, mid, moi and aff")
  if (ncol(genotypes) <= 6) 
    stop("'ped.genotypes' has incorrect format - minimum of 2 isolates required for analysis")
  if (nrow(genotypes) <= 1) 
    stop("'ped.genotypes' has incorrect format - too few SNPs for analysis")
  if (any(colnames(genotypes)[1:5] != c("chr", "snp_id", "pos_M", 
                                        "pos_bp", "freq"))) 
    stop("'ped.genotypes' has incorrect format - 'genotypes' must have columns labelled: chr, snp_id, pos_M, pos_bp and freq")
  if (!is.data.frame(parameters)) 
    stop("'parameters' has incorrect format - must be a data.frame")
  if (ncol(parameters) != 8) 
    stop("'parameters' has incorrect format - must have 8 columns labelled: fid1, iid1, fid2, iid2, m, ibd0, ibd1 and ibd2")
  if (any(colnames(parameters) != c("fid1", "iid1", "fid2", 
                                    "iid2", "m", "ibd0", "ibd1", "ibd2"))) 
    stop("'parameters' has incorrect format - must have 8 columns labelled: fid1, iid1, fid2, iid2, m, ibd0, ibd1 and ibd2")
  if (!is.vector(number.cores)) 
    stop("'number.cores' has incorrect format - must be a vector")
  if (!is.numeric(number.cores)) 
    stop("'number.cores' has incorrect format - must be numeric")
  if (length(number.cores) != 1) 
    stop("'number.cores' has incorrect format - must be a single numeric value")
  if (number.cores < 1) 
    stop("'number.cores' has incorrect format - must be >= 1")
  if (!is.vector(minimum.snps)) 
    stop("'minimum.snps' has incorrect format - must be a vector")
  if (!is.numeric(minimum.snps)) 
    stop("'minimum.snps' has incorrect format - must be numeric")
  if (length(minimum.snps) != 1) 
    stop("'minimum.snps' has incorrect format - must be a single numeric value")
  if (minimum.snps < 0) 
    stop("'minimum.snps' has incorrect format - must be >= 0")
  if (!is.vector(minimum.length.bp)) 
    stop("'minimum.length.bp' has incorrect format - must be a vector")
  if (!is.numeric(minimum.length.bp)) 
    stop("'minimum.length.bp' has incorrect format - must be numeric")
  if (length(minimum.length.bp) != 1) 
    stop("'minimum.length.bp' has incorrect format - must be a single numeric value")
  if (minimum.length.bp < 0) 
    stop("'minimum.length.bp' has incorrect format - must be >= 0")
  if (!is.vector(error)) 
    stop("'error' has incorrect format - must be a vector")
  if (!is.numeric(error)) 
    stop("'error' has incorrect format - must be numeric")
  if (length(error) != 1) 
    stop("'error' has incorrect format - must be a single numeric value")
  if (error < 0 | error > 1) 
    stop("'error' has incorrect format - must be between 0 and 1 (inclusive)")
  isolate.pairs <- isolatePair(pedigree[, 1], pedigree[, 2])
  number.pairs <- 1:nrow(isolate.pairs)
  pair.quantiles <- unique(round(quantile(number.pairs, probs = seq(0, 
                                                                    0.9999, 0.01))))
  number.quantiles <- length(pair.quantiles)
  pb <- txtProgressBar(min = 0, max = number.quantiles, style = 3)
  cl = makeCluster(number.cores)
  doParallel::registerDoParallel(cl)  #doParallel::registerDoParallel(cores = number.cores)
  chromosome <- unique(as.character(genotypes[, "chr"]))
  start <- 1
  ibd.segments <- NULL
  for (quantile.group in 1:number.quantiles) {
    if (number.quantiles == nrow(isolate.pairs)) 
      pair.group <- quantile.group
    if (number.quantiles < nrow(isolate.pairs) & quantile.group != 
        max(number.quantiles)) 
      pair.group <- start:(start + pair.quantiles[quantile.group + 
                                                    1] - pair.quantiles[quantile.group] - 1)
    if (number.quantiles < nrow(isolate.pairs) & quantile.group == 
        max(number.quantiles)) 
      pair.group <- start:nrow(isolate.pairs)
    ibd.segments.0 <- foreach::foreach(pair.i = pair.group, .combine = "rbind", .packages=c('Rcpp'), .noexport = c("IBDLabel","roundDecimal","emissionProbHH","emissionProbHD","emissionProbDD","transitionProbHH","transitionProbHD","transitionProbDD","genotypeErrorH","genotypeErrorD",
                                                                                                                   "trueGenotypes","emissionProbMissingGeno","calculateAlpha","calculateScale","calculateBeta","calculateGamma","calculateLogLikelihood")) %dopar% 
    {
        sourceCpp("/nfsscratch/Karim/DELGEME/isoRelate/Codes/calculateViterbi.cpp")
        source('/nfsscratch/Karim/DELGEME/isoRelate/Codes/IBDTable.R')
                                         fid.1 <- as.character(isolate.pairs[pair.i, 1])
                                         iid.1 <- as.character(isolate.pairs[pair.i, 2])
                                         fid.2 <- as.character(isolate.pairs[pair.i, 3])
                                         iid.2 <- as.character(isolate.pairs[pair.i, 4])
                                         gender.1 <- pedigree[pedigree[, "fid"] == fid.1 & 
                                                                pedigree[, "iid"] == iid.1, "moi"]
                                         gender.2 <- pedigree[pedigree[, "fid"] == fid.2 & 
                                                                pedigree[, "iid"] == iid.2, "moi"]
                                         if (gender.1 == 2 & gender.2 == 2) {
                                           number.states = 3
                                         }
                                         else number.states = 2
                                         meiosis <- as.numeric(parameters[parameters[, "fid1"] == 
                                                                            fid.1 & parameters[, "fid2"] == fid.2 & parameters[, 
                                                                                                                               "iid1"] == iid.1 & parameters[, "iid2"] == iid.2, 
                                                                          "m"])
                                         initial.prob <- as.numeric(parameters[parameters[, 
                                                                                          "fid1"] == fid.1 & parameters[, "fid2"] == fid.2 & 
                                                                                 parameters[, "iid1"] == iid.1 & parameters[, 
                                                                                                                            "iid2"] == iid.2, c("ibd0", "ibd1", "ibd2")])
                                         if (initial.prob[1] == 1) {
                                           initial.prob[1] <- 0.999
                                           initial.prob[2] <- 0.001
                                         }
                                         if (initial.prob[2] == 1) {
                                           initial.prob[1] <- 0.001
                                           initial.prob[2] <- 0.999
                                         }
                                         ibd.table.2 <- NULL
                                         for (chrom in chromosome) {
                                           pop.allele.freqs <- genotypes[genotypes[, "chr"] == 
                                                                           chrom, "freq"]
                                           pair.genotypes <- cbind(genotypes[genotypes[, 
                                                                                       "chr"] == chrom, paste(fid.1, iid.1, sep = "/")], 
                                                                   genotypes[genotypes[, "chr"] == chrom, paste(fid.2, 
                                                                                                                iid.2, sep = "/")])
                                           positions.m <- genotypes[genotypes[, "chr"] == 
                                                                      chrom, "pos_M"]
                                           positions.bp <- genotypes[genotypes[, "chr"] == 
                                                                       chrom, "pos_bp"]
                                           chromosomes <- as.character(genotypes[genotypes[, 
                                                                                           "chr"] == chrom, "chr"])
                                           markers <- as.character(genotypes[genotypes[, 
                                                                                       "chr"] == chrom, "snp_id"])
                                           number.snps <- length(positions.m)
                                           viterbi <- calculateViterbi(number.states, initial.prob, 
                                                                       meiosis, number.snps, pair.genotypes, pop.allele.freqs, 
                                                                       positions.m, error, gender.1, gender.2)
                                           ibd.results <- cbind(fid.1, iid.1, fid.2, iid.2, 
                                                                1:number.snps, chromosomes, markers, positions.m, 
                                                                positions.bp, viterbi)
                                           colnames(ibd.results) <- c("fid1", "iid1", "fid2", 
                                                                      "iid2", "markerNo", "chr", "marker", "pos.m", 
                                                                      "pos.bp", "viterbi")
                                           ibd.table.1 <- IBDTable(ibd.results)
                                           if (length(ibd.table.1) != 0) {
                                             ibd.table.2 <- rbind(ibd.table.2, ibd.table.1[as.numeric(ibd.table.1[, 
                                                                                                                  "number.snps"]) >= minimum.snps & as.numeric(ibd.table.1[, 
                                                                                                                                                                           "length.bp"]) >= minimum.length.bp, ])
                                           }
                                         }
                                         ibd.table.2
                                       }
    ibd.segments <- rbind(ibd.segments, ibd.segments.0)
    if (number.quantiles < nrow(isolate.pairs) & quantile.group != 
        max(number.quantiles)) 
      start <- start + pair.quantiles[quantile.group + 
                                        1] - pair.quantiles[quantile.group]
    setTxtProgressBar(pb, quantile.group)
  }
  close(pb)
  rownames(ibd.segments) <- NULL
  ibd.segments <- data.frame(ibd.segments)
  if (nrow(ibd.segments) > 0) {
    for (i in 1:7) ibd.segments[, i] <- as.character(ibd.segments[, 
                                                                  i])
    for (i in 8:15) ibd.segments[, i] <- as.numeric(as.character(ibd.segments[, 
                                                                              i]))
    colnames(ibd.segments) <- c("fid1", "iid1", "fid2", 
                                "iid2", "chr", "start_snp", "end_snp", "start_position_bp", 
                                "end_position_bp", "start_position_M", "end_position_M", 
                                "number_snps", "length_bp", "length_M", "ibd_status")
    number.pairs.ibd <- length(unique(paste(ibd.segments[, 
                                                         1], ibd.segments[, 2], ibd.segments[, 3], ibd.segments[, 
                                                                                                                4])))
    cat(paste(number.pairs.ibd, "pairs inferred IBD\n"))
    cat(paste(nrow(ibd.segments), "IBD segments detected\n"))
  }
  else {
    cat("0 pairs inferred IBD\n")
    cat("0 IBD segments detected\n")
  }
  return(ibd.segments)
}
