my_getIBDparameters = function (ped.genotypes, number.cores = 1) 
{
  library(Rcpp)
  library(foreach)
  library(doParallel)
  #library(isoRelate)
  sourceCpp("/nfsscratch/Karim/DELGEME/isoRelate/Codes/isolatePairs.cpp")
  #sourceCpp("/Users/karimmane/Documents/Alfred/my_IBDparameters.cpp")
  #source('/media/Data/Data/Documents_Karim/DELGEME/Deconvolution/Codes/my_IBDparameters.R')
  
  
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
  if (!is.vector(number.cores)) 
    stop("'number.cores' has incorrect format - must be a vector")
  if (!is.numeric(number.cores)) 
    stop("'number.cores' has incorrect format - must be numeric")
  if (length(number.cores) != 1) 
    stop("'number.cores' has incorrect format - must be a single numeric value")
  if (number.cores < 1) 
    stop("'number.cores' has incorrect format - must be >= 1")
  isolate.pairs <- isolatePair(pedigree[, 1], pedigree[, 2])
  number.pairs <- 1:nrow(isolate.pairs)
  pair.quantiles <- unique(round(quantile(number.pairs, probs = seq(0, 0.9999, 0.01))))
  number.quantiles <- length(pair.quantiles)
  pb <- txtProgressBar(min = 0, max = number.quantiles, style = 3)
  cl = makeCluster(number.cores)
  doParallel::registerDoParallel(cl) #doParallel::registerDoParallel(cores = number.cores)
  start <- 1
  ibd.estimates <- NULL
  for (quantile.group in 1:number.quantiles) 
  {
    if (number.quantiles == nrow(isolate.pairs)) 
      pair.group <- quantile.group
    if (number.quantiles < nrow(isolate.pairs) & quantile.group != max(number.quantiles)) 
      pair.group <- start:(start + pair.quantiles[quantile.group + 1] - pair.quantiles[quantile.group] - 1)
    if (number.quantiles < nrow(isolate.pairs) & quantile.group == max(number.quantiles)) 
      pair.group <- start:nrow(isolate.pairs)
    ibd.estimates.0 <- foreach::foreach(pair_i = pair.group, .combine = "rbind", .packages=c('Rcpp'), .noexport = c("bVectorHH","bVectorHD","bVectorDD","AmatrixHH","AmatrixHD","AmatrixDD")) %dopar% 
    {
        source('/nfsscratch/Karim/DELGEME/isoRelate/Codes/my_IBDparameters.R')
      fid.1 <- as.character(isolate.pairs[pair_i, 1])
      iid.1 <- as.character(isolate.pairs[pair_i, 2])
      fid.2 <- as.character(isolate.pairs[pair_i, 3])
      iid.2 <- as.character(isolate.pairs[pair_i, 4])
      gender.1 <- pedigree[pedigree[, "fid"] == fid.1 & pedigree[, "iid"] == iid.1, "moi"]
      gender.2 <- pedigree[pedigree[, "fid"] == fid.2 & pedigree[, "iid"] == iid.2, "moi"]
      pair.genotypes <- cbind(genotypes[, paste(fid.1, iid.1, sep = "/")], genotypes[, paste(fid.2, iid.2, sep = "/")])
      pop.allele.freqs <- genotypes[, "freq"]
      ibd.estimates <- cbind(fid.1, iid.1, fid.2, iid.2, IBDparameters(pair.genotypes, pop.allele.freqs, gender.1, gender.2))
      colnames(ibd.estimates) <- c("fid1", "iid1", "fid2", "iid2", "m", "ibd0", "ibd1", "ibd2")
      ibd.estimates
    }
    ibd.estimates <- rbind(ibd.estimates, ibd.estimates.0)
    if (number.quantiles < nrow(isolate.pairs) & quantile.group != max(number.quantiles)) 
      start <- start + pair.quantiles[quantile.group + 1] - pair.quantiles[quantile.group]
    setTxtProgressBar(pb, quantile.group)
  }
  stopCluster(cl)
  close(pb)
  ibd.estimates <- data.frame(ibd.estimates)
  for (i in 1:4) ibd.estimates[, i] <- as.character(ibd.estimates[, 
                                                                  i])
  for (i in 5:8) ibd.estimates[, i] <- as.numeric(as.character(ibd.estimates[, 
                                                                             i]))
  return(ibd.estimates)
}
