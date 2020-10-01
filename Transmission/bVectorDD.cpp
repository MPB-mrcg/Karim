#include <Rcpp.h>
using namespace Rcpp;


//----- 2 diploid chromosomes, b and A -----

//' The vector b in the equation Ax=b for 2 diploid chromosomes
//' @param genotypes An integer matrix of genotype calls for a pair of isolates. Each coloumn represents and isolate and each row represents a SNP.
// [[Rcpp::export]]
IntegerVector bVectorDD(IntegerMatrix genotypes){
  IntegerVector Z(3);
  int IBS_0 = 0, IBS_1 = 0, IBS_2 = 0;
  int number_snps = genotypes.nrow();

  for(int i = 0; i < number_snps; ++i){ 
    if(genotypes(i,0) != -1 && genotypes(i,1) != -1){
      if(genotypes(i,0) == genotypes(i,1)){
        IBS_2 += 1;
      }
      if((genotypes(i,0) == 1 && genotypes(i,1) == 0) || (genotypes(i,0) == 1 && genotypes(i,1) == 2) ||
         (genotypes(i,0) == 0 && genotypes(i,1) == 1) || (genotypes(i,0) == 2 && genotypes(i,1) == 1)){
        IBS_1 += 1;
      }
      if((genotypes(i,0) == 0 && genotypes(i,1) == 2) || (genotypes(i,0) == 2 && genotypes(i,1) == 0)){
        IBS_0 += 1;
      }
    }
  }
  Z[0] = IBS_0;
  Z[1] = IBS_1;
  Z[2] = IBS_2;
  return(Z);
}