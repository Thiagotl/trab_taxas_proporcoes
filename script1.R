# Função Densidade

dUGo<-function(x, mu=.5, sigma=1.2){
  
  fx<-mu*sigma * (x^-(sigma+1)) * exp(-mu*(x^(-sigma)-1))
  
  return(fx)
  
}

# função acumulada

pUGo<-function(q, mu=.5, sigma=1.2){
  
  cdf<-exp(-mu*(q^(-sigma)-1))
  
  return(cdf)
  
}


pUGo(.25)
integrate(dUGo, 0, .25)

# Funcção quantílica


qUGo<-function(u, mu=.5, sigma=1.2)
{
  q<-((-log(u)/mu)+1)^(-1/sigma)
  
  return(q)
}

u=pUGo(.82)
u

qUGo(u)


# Método da geração de número aleatórios 


rUGO <- function(n, mu=.5, sigma=1.2) {
  u <- runif(n)
  y<-((-log(u)/mu)+1)^(-1/sigma)
  return(y)
}


# Função de Log-Verossimilhança

log_lik <- function(theta, x) {
  
  n <- length(x)
  mu<-theta[1]
  sigma<-theta[2]
  
  term1 <- n * log(mu)
  term2 <- n * log(sigma)
  term3 <- -(sigma + 1) * sum(log(x))
  term4 <- -mu * sum(x^(-sigma) - 1)
  
  ll <- term1 + term2 + term3 + term4
  return(ll)
}

# Exemplo da estimação dos parâmetros

set.seed(123)

x<-rUGO(100)

rUGO <- function(n, mu=.5, sigma=1.2) {
  u <- runif(n)
  y <- ((-log(u)/mu)+1)^(-1/sigma)
  return(y)
}

parametros_inicial <- c(0.5, 1.6)  # Valores iniciais para alpha e beta

optim_result <- optim(parametros_inicial, log_lik, x=x, hessian = T,
                      method="BFGS", 
                      control = list(fnscale=-1))


print(optim_result$par)




