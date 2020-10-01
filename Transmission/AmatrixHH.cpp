#include <Rcpp.h>
using namespace Rcpp;

//' The matrix A in the equation Ax=b for 2 haploid chromosomes
//' @param pop_allele_freqs A numeric vector of population allele frequencies for each SNP
//' @param genotypes An integer matrix of genotype calls for a pair of isolates. Each coloumn represents and isolate and each row represents a SNP.
// [[Rcpp::export]]
NumericMatrix AmatrixHH(NumericVector pop_allele_freqs,IntegerMatrix genotypes){
  NumericMatrix A(2,2);
  int number_snps = pop_allele_freqs.size();
  double i0_z0 = 0, i0_z1 = 0, i1_z0 = 0, i1_z1 = 0;
  double p, q;

  for(int i = 0; i < number_snps; ++i){
    if(genotypes(i,0) != -1 && genotypes(i,1) != -1){
      p = pop_allele_freqs[i]; 
      q = 1 - p;
      i0_z0 += 2*p*q;
      i1_z0 += pow(p,2)+pow(q,2);
      i1_z1 += 1;
    }
  }
  A(0,0) = i0_z0;
  A(0,1) = i1_z0;
  A(1,0) = i0_z1;
  A(1,1) = i1_z1;
  
  return(A);
}