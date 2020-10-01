#include <Rcpp.h>
using namespace Rcpp;

//----- 1 haploid chromosome and 1 diploid chromosome, b and A -----

//' The vector b in the equation Ax=b for 1 haploid chromosome and 1 diploid chromosome
//' @param genotypes An integer matrix of genotype calls for a pair of isolates. Each coloumn represents and isolate and each row represents a SNP.
// [[Rcpp::export]]
IntegerVector bVectorHD(IntegerMatrix genotypes){
  IntegerVector Z(2);
  int IBS_0, IBS_1 = 0;
  int number_snps = genotypes.nrow();
  int number_snps_1 = 0;

  for(int i = 0; i < number_snps; ++i){ 
    if(genotypes(i,0) != -1 && genotypes(i,1) != -1){
      if(genotypes(i,0) == genotypes(i,1) || genotypes(i,0) == 1 || genotypes(i,1) == 1){
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