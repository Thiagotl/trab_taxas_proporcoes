# PROBABILITY DENSITY FUNCTION


dUGo<-function(x, mu=.5, sigma=1.2){
  
  fx<-mu*sigma * (x^-(sigma+1)) * exp(-mu*(x^(-sigma)-1))
  
  return(fx)
  
}

# CUMULATIVE DISTRIBUTION FUNCTION 
pUGo<-function(q, mu=.5, sigma=1.2){
  
  cdf<-exp(-mu*(q^(-sigma)-1))
  
  return(cdf)
  
}


pUGo(.25)
integrate(dUGo, 0, .25)

# QUANTILE FUNCTION
qUGo<-function(u, mu=.5, sigma=1.2)
{
  q<-((-log(u)/mu)+1)^(-1/sigma)
  
  return(q)
}

u=pUGo(.82)
u

qUGo(u)



# INVERSION METHOD FOR RANDOM GENERATION

rUGO <- function(n, mu=.5, sigma=1.2, tau=.5) {
  u <- runif(n)
 
  return(y)
}



