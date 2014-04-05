### This is a function that takes in two chromosomes, in the form of two vectors of the same length.
### With the probability of recombination per site, p_recomb, it will form a new chromosome, which is the output of this function.
##################################################################################################################################
 

recombine<-function(maternalchromosome,paternalchromosome,p_recomb)
{
    ### Define a function that can split a vector into subvectors, based on the coordinates of the breaks. 
    ### It takes as an input x, the vector it should cut, and pos, a vector with coordinates to cut at. 
    ### It outputs a list with all the new vectors, in the original order. 
    splitAt <- function(x, pos) unname(split(x, cumsum(seq_along(x) %in% pos)))

    ### The number of loci is the length of the chromosomes. The length of maternal and paternal chromosomes should be the same, by definition
    if(length(maternalchromosome) != length(paternalchromosome)){STOP("paternal and maternal chromosome are not of the same length")}
    Nloci=length(maternalchromosome)
    
    ### The number of breaksites should be the number of loci - 1    
    Nbreaksites=Nloci-1
    
    ### We will introduce breaks according to a poisson distribution. For this, we need to define a lambda.
    ### Smaller lambda means lower frequency of recombination. Very high lambda means complete scrambling. 
    ### Here, lambda is defined as the chance of recombination per breaksite, times the number of breaksites. 
    ### A lambda = 1 means that on average half of the chromosomes will have at least one break. 
    lambda<-p_recomb*Nbreaksites

    ### Now, put our chromosomes in a matrix (consider it analogous to a chromosome pair before crossing over starts)
    diploid<-as.matrix(rbind(maternalchromosome,paternalchromosome))
    
    ### Here pull from the poisson distribution the number of breaks we will introduce in this chromosome    
    nbreaks<-rpois(1,lambda)
    

    ### Randomly pick a number from {1,2} to choose the parental chromosome to start with
    randomsampler=sample(1:2,1)
    
    ### If the number of breaks is 0, then one of the chromosomes is inherited as it is. 
    ### This means we can just choose one, and go to the end of the function!
    ### In the other case, you go on after the else statement. 
    if(nbreaks==0){newchrom<-rbind(maternalchromosome,paternalchromosome)[randomsampler,]}else{
    
        ### Now make sure that we don't have more breaks than possible break sites. 
        if(nbreaks>Nbreaksites){nbreaks<-Nbreaksites}    
    
        ### Now let's make a vector with the sites at which our chromosome is going to break.
        ### For this, we sample from a uniform distribution between 2 and the end of your chromosome, 
        ### which is the number of breaksites +1
        ### The reason we use this interval, is because our splitAt function splits BEFORE the coordinates we give it!
        ### The reason we sort it, is because the splitAt function wants ordered coordinates. 
        breakcoordinates<-unique(sample(2:(Nbreaksites+1),nbreaks,replace=F))
        breakcoordinates<-sort(breakcoordinates)
        ### So the output of this, for chromosome abcdefghij should be a vector like 2,5,8 
        ### This chromosome will then break in the fragments a bcd efg hij 

        ### Now break mama and papa chromosomes according to the breakcoordinates, using the above defined function splitAt
        ### So the broken chromosome output is a list of fragments!
        brokenpaternal<-splitAt(paternalchromosome,breakcoordinates)
        brokenmaternal<-splitAt(maternalchromosome,breakcoordinates)
    
        ### Make an empty vector that is going to be our chromosome
        newchrom<-NULL
    
        ### And now stitch them together, and let them "grow" 
        ### You need nbreaks + 1 iterations, because n breaks make n+1 fragments. 
        for(chrombreak in 1:(nbreaks+1))
        {
            ### Remember the random sampler is just a random number between 1 and 2, to choose the chromosome to start with
            ### So of these two statements, only one is going to be executed at every permutation...
            if(randomsampler==1){newchrom<-c(newchrom,brokenpaternal[[chrombreak]])}
            if(randomsampler==2){newchrom<-c(newchrom,brokenmaternal[[chrombreak]])}
            
            ### And the next permutation, the randomsampler will always have changed from 1 to 2 or from 2 to 1
            ### through this nice mathematical function that always turns a 1 in a 2 and other way around.
            randomsampler<-1+(randomsampler)%%2
        }

    }
    ### In the end, don't forget to return that nice little new chromosome you just made to your environment!
    return(newchrom)

}

