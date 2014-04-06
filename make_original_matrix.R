### Function that makes a population of organisms with one chromosome, with a given number of loci, with always two alleles per locus!


make_initial_genotype_matrix<-function  (
                                        popsize,
                                        n_loci,
                                        n_HTs,                              # Number of haplotypes
                                        allelefrequencies="norm",           # choose from "norm","unif"
                                        allele_exp=0.5,                     # expected value for allele frequencies
                                        allele_stdev=0.05,                  # standard deviation for allele frequencies in the case of normal distr.
                                        HT_lambda=1                         # lambda for frequencies of haplotypes. High lambda's give equal distr. Lower can lead to loss of haplotypes!
                                        )
{
  
  ### first define a function that calculates normal values, but truncates between 0 and 1
  rtnorm <- function(n, mean, sd, a = 0, b = 1)
  {
  qnorm(runif(n, pnorm(a, mean, sd), pnorm(b, mean, sd)), mean, sd)
  }
  
  if(n_loci>26){stop("maximal number of loci is 26")}
  #create vectors of alleles for each locus

  alleles = matrix(0, nrow = n_loci, ncol = 2)
  ### generate loci
  for (i in 1:n_loci)
  {
    alleles[i,] = c(toupper(letters)[i], letters[i])
  }
  
  
### First determine the allele frequencies for all n_loci loci
### Next, make a matrix with the haplotypes
  if(tolower(allelefrequencies)=="norm")  
  {  
    allelefreqmatrix<-rtnorm(n_loci,allele_exp,allele_stdev)   
    randomHTs<-matrix(0,ncol=n_loci,nrow=n_HTs)
    for(i in 1:n_HTs)
    {
      for(j in 1:n_loci)
      {
        randomHTs[i,j]<-sample(alleles[j,],1,replace=T,prob=c(allelefreqmatrix[j],1-allelefreqmatrix[j]))
      }
    }
  }      
  
  ### make a deterministic matrix for the uniformly distributed ones
  if(tolower(allelefrequencies)=="unif")  
  {
    numberoffirstallele   <-ceiling(   allele_exp *n_HTs)
    numberofsecondallele  <-ceiling((1-allele_exp)*n_HTs)
    firstallele<-matrix(alleles[,1],ncol=numberoffirstallele,nrow=n_loci)
    secondallele<-matrix(alleles[,2],ncol=numberofsecondallele,nrow=n_loci)
    alleles<-cbind(firstallele,secondallele)
    randomHTs<-matrix(0,ncol=n_loci,nrow=n_HTs)
    #print(randomHTs)
    for(i in 1:n_loci)
    {
      allelevector<-sample(alleles[i,],n_HTs,replace=F)
      randomHTs[,i]<-allelevector
    }
  }
  
#print(randomHTs)

  
### Now determine haplotype frequencies for all n_HT HTs
  HTfreqmatrix<-rpois(n_HTs,HT_lambda)
  
### Now sample from these random HTs with a poisson distribution
  populationmatrix<-matrix(0,ncol=n_loci,nrow=popsize)
  for (i in 1:popsize)
  {
    populationmatrix[i,]<-randomHTs[sample(1:n_HTs,1,replace=T,HTfreqmatrix),]
  }
  
  
### Return some stats
allelefrequencies<-NULL
for(i in 1:n_loci)
{
  allelefrequency<-sum(populationmatrix[,i]==unique(populationmatrix[,i])[1])/popsize
  allelefrequencies<-c(allelefrequencies,allelefrequency)
}

haplotypefrequencies<-NULL
uniquehaplotypes<-unique(populationmatrix)
for (i in 1:nrow(uniquehaplotypes))
{
  counter=0
  haplostring<-paste(uniquehaplotypes[i,],collapse="")
  for (j in 1:nrow(populationmatrix))
  {
    vec<-paste(populationmatrix[j,],collapse="")
    #print(vec)
    if (vec == haplostring){counter=counter+1}
  }
  haplotypefrequencies<-c(haplotypefrequencies,counter)
}

haplotypefrequencies<-cbind(uniquehaplotypes,haplotypefrequencies)


cat("allelefrequencies:\n")
print(allelefrequencies)
cat("haplotypefrequencies:\n")
print(haplotypefrequencies)
  
return(populationmatrix)
}
mat<-make_initial_genotype_matrix(popsize=1000,n_loci=12,n_HTs=15,allelefrequencies="unif",allele_exp=0.5,allele_stdev=0.1,HT_lambda=1)


#mat<-make_initial_genotype_matrix(1000,5,10,"norm",0.5,0.05)
