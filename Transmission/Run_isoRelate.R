#!/usr/bin/env Rscript

runIsorelate = function()
{
  library(data.table)
  library(isoRelate)
  library(tictoc)
  library(Rcpp)
  
  arguments = commandArgs(trailingOnly = TRUE)
  pedFile = as.character(arguments[1])
  mapFile = as.character(arguments[2])
  outDir = as.character(arguments[3])
  pathToRcodes = as.character(arguments[4])
  groupID = as.character(arguments[5])
  print(arguments)
  
  if(!file.exists(pedFile))
    stop(pedFile, "\nNo such file or directory")
  if(!file.exists(mapFile))
    stop(mapFile, "\nNo such file or directory")
  if(!file.exists(groupID))
    stop(groupID, "\nNo such file or directory")
  if(!dir.exists(outDir))
    system(sprintf("mkdir -p %s", outDir))
  if(!dir.exists(pathToRcodes))
    stop(pathToRcodes, "\nNo such file or directory")
  
  source(paste0(pathToRcodes,'/','getGenotypes.R'))
  sourceCpp(paste0(pathToRcodes,'/',"my_IBDparameters.cpp"))
  source(paste0(pathToRcodes,'/',"my_getIBDparameters.R"))
  source(paste0(pathToRcodes,'/',"my_getIBDsegments.R"))
  
  ped = fread(pedFile, header = FALSE)
  map = fread(mapFile, header = FALSE)
  group = fread(groupID, header = TRUE)
  ped.map = list(ped, map)
  
  tic("Time to reformat the data for isoRelate")
  return.genotypes = getGenotypes(ped.map = ped.map, reference.ped.map = NULL, maf = 0.01, chromosomes = NULL, input.map.distance = "cM", reference.map.distance = "cM")
  #write.table(return.genotypes[[1]], file=paste0(outDir,'/',"isoRelateGenotypesFormat.map"), sep = "\t", col.names = TRUE, row.names = FALSE,quote = FALSE)
  #write.table(return.genotypes[[2]], file=paste0(outDir,'/',"isoRelateGenotypesFormat.ped"), sep = "\t", col.names = TRUE, row.names = FALSE,quote = FALSE)
  toc()
  
  tic("Time to get IBD parameters")
  my_parameters = my_getIBDparameters(ped.genotypes = return.genotypes, number.cores = 50)
  #write.table(my_parameters, file=paste0(outDir,'/',"IBDparameters.txt"), sep = "\t", col.names = TRUE, row.names = FALSE,quote = FALSE)
  toc()
  
  # tic("Time to get IBD parameters")
  # my_parameters = fread(ibdParam)
  # toc()
  
  tic("Time to get IBD segments")
  geneLength = 20000
  snpNumber = 15
  
  minLength = paste0(outDir, '/', 'IBDsegment_',geneLength,'_',snpNumber,'.txt')  #snpNumber[j]
  my_ibd = my_getIBDsegments(ped.genotypes = return.genotypes,parameters = my_parameters, number.cores = 50,minimum.snps = snpNumber,minimum.length.bp = geneLength,error = 0.001)
  write.table(my_ibd, file=minLength, sep = "\t", col.names = TRUE, row.names = FALSE,quote = FALSE)
  toc()
  
  tic("Time to get IBD matrix")
  my_matrix = getIBDmatrix(ped.genotypes = return.genotypes, ibd.segments = my_ibd)
  toc()
  
  tic("Time to get iR")
  my_iR = getIBDiR(ped.genotypes = return.genotypes, ibd.matrix = my_matrix, groups = group)
  minLength = paste0(outDir, '/', 'IBD_iR.txt')  
  write.table(my_iR, file=minLength, sep = "\t", col.names = TRUE, row.names = FALSE,quote = FALSE)
  toc()
}


#----------------------------------- calling the main function -----
runIsorelate()






