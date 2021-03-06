
# this function is a replacement for the 'isoRelate::getGenotypes'
# set the parameter 'maf', 'isolate.max.missing', 'snp.max.missing' to fit your input data
# the imput file 'ped.map' is a list that contains the ped and maf file for you data. The 5th column of the ped file is the MOI for 
    #each sample. Use the phased data to create the ped file. 1->allele A, 2->allele B, 0->missing data

getGenotypes = function(ped.map, reference.ped.map=NULL, maf=0.01, isolate.max.missing=0.2, snp.max.missing=0.2, chromosomes=NULL, input.map.distance="cM", reference.map.distance="cM")
{
    source('/nfsscratch/Karim/DELGEME/isoRelate/Codes/haplotypeToGenotype.R')
    source('/nfsscratch/Karim/DELGEME/isoRelate/Codes/calculatePopAlleleFreq.R') 
    source('/nfsscratch/Karim/DELGEME/isoRelate/Codes/calculateMissingness.R')
    
    # check input PED and MAP files
    if (!is.list(ped.map) | length(ped.map) != 2) stop ("'ped.map' must be a list containing 2 objects: 'PED' and 'MAP'")
    input.ped <- as.data.frame(ped.map[[1]])
    input.map <- as.data.frame(ped.map[[2]])
    if (!is.data.frame(input.ped)) stop ("'ped.map' has incorrect format - PED is not a data.frame")
    if (!is.data.frame(input.map)) stop ("'ped.map' has incorrect format - MAP is not a data.frame")
    
    # check the PED and MAP files have the same number of SNPs
    if (ncol(input.ped) != (2*nrow(input.map)+6)) stop ("'ped.map' has incorrect format - PED and MAP must have the same number of SNPs")
    
    # check the MAP file has 4 coloumns
    if (ncol(input.map) != 4) stop ("'ped.map' has incorrect format - MAP must have four columns")
    colnames(input.map) <- c("chr", "snp_id", "pos_M", "pos_bp")
    if (!is.numeric(input.map[,"pos_M"])) stop ("'ped.map' has incorrect format - genetic-map positions in MAP file are non-numeric")
    if (!is.numeric(input.map[,"pos_bp"])) stop ("'ped.map' has incorrect format - base-pair positions in MAP file are non-numeric")
    if (is.factor(input.map[,"chr"]))
        input.map[,"chr"] <- as.character(input.map[,"chr"])
    if (is.factor(input.map[,"snp_id"]))
        input.map[,"snp_id"] <- as.character(input.map[,"snp_id"])
    
    
    # check PED for factors
    if (is.factor(input.ped[,1]))
        input.ped[,1] <- as.character(input.ped[,1])
    if (is.factor(input.ped[,2]))
        input.ped[,2] <- as.character(input.ped[,2])
    
    
    # check reference data
    if (!is.null(reference.ped.map)) {
        if (!is.list(reference.ped.map) | length(reference.ped.map) != 2) stop ("'reference.ped.map' must be a list containing 2 objects: 'PED' and 'MAP'")
        reference.ped <- reference.ped.map[[1]]
        reference.map <- reference.ped.map[[2]]
        if (!is.data.frame(reference.ped)) stop ("'reference.ped.map' has incorrect format - PED is not a data.frame")
        if (!is.data.frame(reference.map)) stop ("'reference.ped.map' has incorrect format - MAP is not a data.frame")
        
        # check the PED and MAP files have the same number of SNPs
        if (ncol(reference.ped) != (2*nrow(reference.map)+6)) stop ("'reference.ped.map' has incorrect format - PED and MAP must have the same number of SNPs")
        
        # check the MAP file has 4 coloumns
        if (ncol(reference.map) != 4) stop ("'reference.ped.map' has incorrect format - MAP must have four columns")
        colnames(reference.map) <- c("chr", "snp_id", "pos_M", "pos_bp")
        if (!is.numeric(reference.map[,"pos_M"])) stop ("'reference.ped.map' has incorrect format - genetic-map positions in reference MAP file are non-numeric")
        if (!is.numeric(reference.map[,"pos_bp"])) stop ("'reference.ped.map' has incorrect format - base-pair positions in reference MAP file are non-numeric")
        if (is.factor(reference.map[,"chr"]))
            reference.map[,"chr"] <- as.character(reference.map[,"chr"])
        if (is.factor(reference.map[,"snp_id"]))
            reference.map[,"snp_id"] <- as.character(reference.map[,"snp_id"])
        
        # check PED for factors
        if (is.factor(reference.ped[,1]))
            reference.ped[,1] <- as.character(reference.ped[,1])
        if (is.factor(reference.ped[,2]))
            reference.ped[,2] <- as.character(reference.ped[,2])
    }
    
    # check maf
    if (!is.vector(maf)) stop ("'maf' has incorrect format - must be a vector")
    if (!is.numeric(maf)) stop ("'maf' has incorrect format - must be numeric")
    if (length(maf) != 1) stop ("'maf' has incorrect format - must be a single numeric value")
    if (maf > 1 | maf < 0) stop ("'maf' has incorrect format - must be between 0 and 1")
    
    # check isolate.max.missing
    if (!is.vector(isolate.max.missing)) stop ("'isolate.max.missing' has incorrect format - must be a vector")
    if (!is.numeric(isolate.max.missing)) stop ("'isolate.max.missing' has incorrect format - must be numeric")
    if (length(isolate.max.missing) != 1) stop ("'isolate.max.missing' has incorrect format - must be a single numeric value")
    if (isolate.max.missing > 1 | isolate.max.missing < 0) stop ("'isolate.max.missing' has incorrect format - must be between 0 and 1 (inclusive)")
    
    # check snp.max.missing
    if (!is.vector(snp.max.missing)) stop ("'snp.max.missing' has incorrect format - must be a vector")
    if (!is.numeric(snp.max.missing)) stop ("'snp.max.missing' has incorrect format - must be numeric")
    if (length(snp.max.missing) != 1) stop ("'snp.max.missing' has incorrect format - must be a single numeric value")
    if (snp.max.missing > 1 | snp.max.missing < 0) stop ("'snp.max.missing' has incorrect format - must be between 0 and 1 (inclusive)")
    
    # check chromosomes
    if (!is.null(chromosomes)) {
        if (!is.vector(chromosomes)) stop ("'chromosomes' has incorrect format - must be a vector")
        if(!all(chromosomes %in% input.map[,"chr"]))
            stop(paste0("chromosome ",paste0(chromosomes[!(chromosomes %in% input.map[,"chr"])]," not in 'ped.map'\n")))
        if (!is.null(reference.ped.map)) {
            if(!all(chromosomes %in% reference.map[,"chr"]))
                stop(paste0("chromosome ",paste0(chromosomes[!(chromosomes %in% reference.map[,"chr"])]," not in 'reference.ped.map'\n")))
        }
    } else
        chromosomes <- unique(as.character(input.map[,"chr"]))
    
    # check input map distance
    if (!is.vector(input.map.distance)) stop ("'input.map.distance' has incorrect format - must be a vector")
    if (!is.character(input.map.distance)) stop ("'input.map.distance' has incorrect format - must be either character M or cM")
    if (length(input.map.distance) != 1) stop ("'input.map.distance' has incorrect format - must be either character M or cM")
    if (input.map.distance != "M" & input.map.distance != "cM")
        stop ("'input.map.distance' has incorrect format - must be either M or cM")
    if (input.map.distance == "cM") {
        input.map[,"pos_M"] <- input.map[,"pos_M"]/100
    }
    
    # check reference map distance
    if (!is.vector(reference.map.distance)) stop ("'reference.map.distance' has incorrect format - must be a vector")
    if (!is.character(reference.map.distance)) stop ("'reference.map.distance' has incorrect format - must be either character M or cM")
    if (length(reference.map.distance) != 1) stop ("'reference.map.distance' has incorrect format - must be either character M or cM")
    if (reference.map.distance != "M" & reference.map.distance != "cM")
        stop ("'reference.map.distance' has incorrect format - must be either M or cM")
    if (reference.map.distance == "cM" & !is.null(reference.ped.map)) {
        reference.map[,"pos_M"] <- reference.map[,"pos_M"]/100
    }
    
    # begin data filtering
    
    # create new isolate IDs from PED FIDs and IIDs
    isolate.names <- paste(input.ped[,1], input.ped[,2], sep="/")
    if (any(duplicated(isolate.names)))
        stop ("duplicate sample IDs found")
    
    # merge input data with reference data
    if (!is.null(reference.ped.map)) {
        input.map.v1      <- cbind(1:nrow(input.map), input.map)
        reference.map.v1  <- cbind(1:nrow(reference.map), reference.map)
        input.map.v1      <- merge(input.map.v1, reference.map.v1, by.x="snp_id", by.y="snp_id")
        if (nrow(input.map.v1) == 0)
            stop ("no SNPs remaining after merging 'ped.map' and 'reference.ped.map'")
        input.map.v1  <- input.map.v1[order(input.map.v1[,"1:nrow(input.map)"]),]
        if (!is.null(chromosomes))
            input.map.v1 <- input.map.v1[input.map.v1[,"chr.x"] %in% chromosomes,]
        if (nrow(input.map.v1) == 0)
            stop ("no SNPs remaining after merging 'ped.map' and 'reference.ped.map' for selected chromosomes")
        input.ped.columns <- c(1:6, 2*input.map.v1[,"1:nrow(input.map)"] + 5, 2*input.map.v1[,"1:nrow(input.map)"] + 6)
        input.ped.columns <- input.ped.columns[order(input.ped.columns)]
        input.ped.v1      <- input.ped[,input.ped.columns]
        reference.ped.columns  <- c(1:6, 2*input.map.v1[,"1:nrow(reference.map)"] + 5, 2*input.map.v1[,"1:nrow(reference.map)"] + 6)
        reference.ped.columns  <- reference.ped.columns[order(reference.ped.columns)]
        reference.ped.v1       <- reference.ped[,reference.ped.columns]
        input.map.v2           <- input.map.v1[,c("chr.x", "snp_id", "pos_M.x", "pos_bp.x")]
        colnames(input.map.v2) <- c("chr", "snp_id", "pos_M","pos_bp")
    } else {
        if (!is.null(chromosomes)) {
            input.map.v1 <- cbind(1:nrow(input.map), input.map)
            input.map.v1 <- input.map.v1[input.map.v1[,"chr"] %in% chromosomes,]
            if (nrow(input.map.v1) == 0)
                stop ("no SNPs remaining after subsetting 'ped.map'by selected chromosomes")
            input.ped.columns <- c(1:6, 2*input.map.v1[,"1:nrow(input.map)"] + 5, 2*input.map.v1[,"1:nrow(input.map)"] + 6)
            input.ped.columns <- input.ped.columns[order(input.ped.columns)]
            input.ped.v1      <- input.ped[,input.ped.columns]
            input.map.v2      <- input.map.v1[,c("chr", "snp_id", "pos_M", "pos_bp")]
        } else {
            input.map.v2 <- input.map
            input.ped.v1 <- input.ped
        }
    }
    
    # call genotypes
    input.matrix        <- as.matrix(input.ped.v1[,7:ncol(input.ped.v1)])
    input.genders       <- input.ped.v1[,5]
    input.genotypes.v0  <- cbind(input.map.v2, haplotypeToGenotype(input.matrix, input.genders))
    if (!is.null(reference.ped.map)) {
        reference.matrix       <- as.matrix(reference.ped.v1[,7:ncol(reference.ped.v1)])
        reference.genders      <- reference.ped.v1[,5]
        reference.genotypes.v0 <- cbind(input.map.v2, haplotypeToGenotype(reference.matrix, reference.genders))
    }
    
    
    # calculate allele frequencies form reference data
    if (is.null(reference.ped.map)) {
        pop.allele.freq    <- calculatePopAlleleFreq(as.matrix(input.genotypes.v0[,5:ncol(input.genotypes.v0)]), input.ped.v1[,5])
        input.genotypes.v1 <- cbind(input.genotypes.v0[,c(1:4)],pop.allele.freq,input.genotypes.v0[,5:ncol(input.genotypes.v0)])
    } else {
        pop.allele.freq    <- calculatePopAlleleFreq(as.matrix(reference.genotypes.v0[,5:ncol(reference.genotypes.v0)]), reference.ped.v1[,5])
        input.genotypes.v1 <- cbind(input.genotypes.v0[,c(1:4)],pop.allele.freq,input.genotypes.v0[,c(5:ncol(input.genotypes.v0))])
    }
    colnames(input.genotypes.v1) <- c("chr", "snp_id", "pos_M","pos_bp", "freq", isolate.names)
    cat(paste("Begin filtering of ",length(isolate.names)," isolates and ",nrow(input.genotypes.v1)," SNPs...\n",sep=""))
    
    
    # remove SNPs with low population MAF
    input.genotypes.v2 <- subset(input.genotypes.v1, pop.allele.freq <= (1-maf) & pop.allele.freq >= maf) #removing SNPs with AF>0.99 and AF<0.01
    if (nrow(input.genotypes.v2) == 0)
        stop("0 SNPs remain after MAF removal")
    cat(paste(nrow(input.genotypes.v2)," SNPs remain after MAF removal...\n",sep=""))
    
    
    # remove snps with high missingness
    snp.missingness    <- calculateMissingness(as.matrix(t(input.genotypes.v2[,6:ncol(input.genotypes.v2)])))
    input.genotypes.v3 <- input.genotypes.v2[snp.missingness <= snp.max.missing,]
    if (nrow(input.genotypes.v3) == 0)
        stop("0 SNPs remain after missingness removal")
    cat(paste(nrow(input.genotypes.v3)," SNPs remain after missingness removal...\n",sep=""))
    
    
    # remove samples with high missingness
    isolate.missingness <- round(calculateMissingness(as.matrix(input.genotypes.v3[,6:ncol(input.genotypes.v3)])),digits=3)
    if (length(isolate.names[isolate.missingness > isolate.max.missing]) > 0) {
        my.remove <- isolate.names[isolate.missingness > isolate.max.missing]
        warning("isolates removed due to genotype missingness: ",paste(my.remove, collapse=", "))
        sample.keep        <- input.ped.v1[isolate.missingness <= isolate.max.missing,1:6]
        input.genotypes.v4 <- input.genotypes.v3[,c(1:5, which(isolate.missingness <= isolate.max.missing) + 5)]
        if(nrow(sample.keep) < 1) stop(paste("All isolates removed with missingness > ",isolate.max.missing*100,"%. No isolates remaining.",sep=""))
    } else {
        sample.keep        <- input.ped.v1[,1:6]
        input.genotypes.v4 <- input.genotypes.v3
    }
    colnames(sample.keep) <- c("fid", "iid", "pid", "mid", "moi", "aff")
    if ((ncol(input.genotypes.v4)-5) == 0) {
        stop("0 samples remain after missingness removal")
    }
    cat(paste(ncol(input.genotypes.v4)-5," isolates remain after missingness removal...\n",sep=""))
    
    
    return.genotypes <- list(sample.keep, input.genotypes.v4)
    names(return.genotypes) <- c("pedigree", "genotypes")
    return(return.genotypes)
    
}
