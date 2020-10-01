calculateMissingness = function(genotypes) 
{
    proportion_missing=vector('numeric',dim(genotypes)[2])
    number_snps = dim(genotypes)[2]
    number_isolates = dim(genotypes)[1]
    number_snps_1 = dim(genotypes)[1]
    
    for (i in 1:number_snps) 
    {
        #print(i)
        number_missing = 0.0
        for (j in 1:number_isolates) {
            if(genotypes[j,i] == -1 ) number_missing = number_missing+1;
        }
        proportion_missing[i] = number_missing/number_snps_1;
    }
    return(proportion_missing)
}
#tic();ggg = calculateMissingness(g);toc()