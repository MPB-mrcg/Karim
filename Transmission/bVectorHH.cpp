#include <Rcpp.h>
using namespace Rcpp;

// Below are the functions used to perform the method-of-moments parameter
// estimation as in PLINK, Purcell et al. (2007). This method solves
// Ax=b for x in order to estimate omega values. The matrices A and b are
// created in this script.


//----- 2 haploid chromosomes, b and A -----

//' The vector b in the equation Ax=b for 2 haploid chromosomes
//' @param genotypes An integer matrix of genotype calls for a pair of isolates. Each coloumn represents and isolate and each row represents a SNP.
// [[Rcpp::export]]
IntegerVector bVectorHH(IntegerMatrix genotypes){
  IntegerVector Z(2);
  int IBS_0, IBS_1 = 0;
  int number_snps = genotypes.nrow();
  int number_snps_1 = 0;

  for(int i = 0; i < number_snps; ++i){ 
    if(genotypes(i,0) != -1 && genotypes(i,1) != -1){
      if(genotypes(i,0) == genotypes(i,1)){
        IBS_1 += 1;
      }
      number_snps_1 += 1;
    }
  }
  IBS_0 = number_snps_1 - IBS_1;
  Z[0] = IBS_0;
  Z[1] = IBS_1;
  return(Z);
}