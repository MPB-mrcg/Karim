haplotypeToGenotype = function(haplotypes, moi)
{
    #haplotypes[haplotypes==0]=-1
    #haplotypes[haplotypes==1]=0
    #haplotypes[haplotypes=='.']=0
    #genotypes = matrix(-9, dim(haplotypes)[2]/2, dim(haplotypes)[1])
    genotypes = matrix(-9, dim(haplotypes)[1], dim(haplotypes)[2]/2)
    for(i in 1:dim(haplotypes)[1])
    {
        #print(i)
        j=1
        k=1
        while(j <= dim(haplotypes)[2])
        {
            if(haplotypes[i,j] == 1 && haplotypes[i,(j+1)] == 1) genotypes[i,k] = 0
            if(haplotypes[i,j] == 2 && haplotypes[i,(j+1)] == 2) genotypes[i,k] = 2
            if(haplotypes[i,j] == 0 && haplotypes[i,(j+1)] == 0) genotypes[i,k] = -1
            if((moi[i] == 1) && ((haplotypes[i,j] == 1 && haplotypes[i,(j+1)] == 2) || (haplotypes[i,j] == 2 && haplotypes[i,(j+1)] == 1))) genotypes[i,k] = -1
            if((moi[i] == 2) && ((haplotypes[i,j] == 1 && haplotypes[i,(j+1)] == 2) || (haplotypes[i,j] == 2 && haplotypes[i,(j+1)] == 1))) genotypes[i,k] = 1
            if(haplotypes[i,j] != 0 && haplotypes[i,j] != 1 && haplotypes[i,j] != 2) genotypes[i,k] = -1
            if(haplotypes[i,(j+1)] != 0 && haplotypes[i,(j+1)] != 1 && haplotypes[i,(j+1)] != 2) genotypes[i,k] = -1
            j=j+2
            k=k+1
        }
    }
    return(t(genotypes))
}