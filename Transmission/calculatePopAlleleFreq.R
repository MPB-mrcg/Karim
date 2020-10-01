calculatePopAlleleFreq = function(genotypes, moi)
{
    pop_allele_freqs = vector('numeric',dim(genotypes)[1])
    number_isolates = dim(genotypes)[2]
    number_snps = dim(genotypes)[1]
    
    for (t in 1:number_snps) 
    {
        #print(t)
        A = 0
        B = 0 
        for (i in 1:number_isolates) 
        {
            if (genotypes[t,i] == 0 && moi[i] == 1)  A=A+1
            if (genotypes[t,i] == 0 && moi[i] == 2)  A=A+2
            if (genotypes[t,i] == 1 && moi[i] == 2){ A=A+1; B=B+1 }
            if (genotypes[t,i] == 2 && moi[i] == 1)  B=B+1
            if (genotypes[t,i] == 2 && moi[i] == 2)  B=B+1
        }
        if (A + B == 0) pop_allele_freqs[t] = -1
        else pop_allele_freqs[t] = A/(A+B)
    }
    return(pop_allele_freqs)
}
#tic();gg = calculatePopAlleleFreq(g, input.genders);toc()