
UG<-function (mu.link = "logit", sigma.link = "log") 
{
  mstats <- checklink("mu.link", "UG", substitute(mu.link), 
                      c("logit", "probit", "cloglog", "cauchit", "log", "own"))
  dstats <- checklink("sigma.link", "UW", substitute(sigma.link), 
                      c("inverse", "log", "identity", "own"))
  structure(list(family = c("UG", "Unit-Gamma"), 
                 parameters = list(mu = TRUE, sigma = TRUE), 
                 nopar = 2, 
                 type = "Continuous", 
                 mu.link = as.character(substitute(mu.link)), 
                 sigma.link = as.character(substitute(sigma.link)), 
                 
                 mu.linkfun = mstats$linkfun, 
                 sigma.linkfun = dstats$linkfun, 
                 
                 mu.linkinv = mstats$linkinv, 
                 sigma.linkinv = dstats$linkinv, 
                 
                 mu.dr = mstats$mu.eta, 
                 sigma.dr = dstats$mu.eta, 
                 
                 dldm = function(y, mu, sigma) {
                   dldm <- mu^(1/sigma)/((1-mu^(1/sigma))*mu^(1/sigma+1))*
                     (1+mu^(1/sigma)/((1-mu^(1/sigma))*sigma)*log(y))
                   dldm
                 }, 
                 d2ldm2 = function(y,mu, sigma) {
                   dldm <- mu^(1/sigma)/((1-mu^(1/sigma))*mu^(1/sigma+1))*
                     (1+mu^(1/sigma)/((1-mu^(1/sigma))*sigma)*log(y))
                   d2ldm2 <- -dldm * dldm
                   d2ldm2 <- ifelse(d2ldm2 < -1e-15, d2ldm2,-1e-15) 
                   d2ldm2
                 }, 
                 dldd = function(y, mu, sigma) {
                   dldd <- log(-log(y))-1/sigma*mu^(1/sigma)/(1-mu^(1/sigma))*log(mu)*
                     (1+mu^(1/sigma)/((1-mu^(1/sigma))*sigma*mu^(1/sigma))*log(y))-
                     log((1-mu^(1/sigma)))-digamma(sigma)
                   dldd
                 },
                 d2ldd2 = function(y,mu, sigma) {
                   dldd <- log(-log(y))-1/sigma*mu^(1/sigma)/(1-mu^(1/sigma))*log(mu)*
                     (1+mu^(1/sigma)/((1-mu^(1/sigma))*sigma*mu^(1/sigma))*log(y))-
                     log((1-mu^(1/sigma)))-digamma(sigma)
                   d2ldd2 = -dldd * dldd
                   d2ldd2 <- ifelse(d2ldd2 < -1e-15, d2ldd2,-1e-15)  
                   d2ldd2
                 },
                 d2ldmdd = function(y,mu, sigma) {
                   dldm <- mu^(1/sigma)/((1-mu^(1/sigma))*mu^(1/sigma+1))*
                     (1+mu^(1/sigma)/((1-mu^(1/sigma))*sigma)*log(y))
                   dldd <- log(-log(y))-1/sigma*mu^(1/sigma)/(1-mu^(1/sigma))*log(mu)*
                     (1+mu^(1/sigma)/((1-mu^(1/sigma))*sigma*mu^(1/sigma))*log(y))-
                     log((1-mu^(1/sigma)))-digamma(sigma)
                   d2ldmdd = -(dldm * dldd)
                   d2ldmdd<-ifelse(is.na(d2ldmdd)==TRUE,0,d2ldmdd)
                   d2ldmdd  
                 }, 
                 G.dev.incr = function(y, mu, sigma, w, ...) -2 * log(dUG(y=y, mu=mu, sigma=sigma)), 
                 rqres = expression(
                   rqres(pfun = "pUG",  type = "Continuous", y = y, mu = mu, sigma = sigma)
                 ),
                 
                 mu.initial = expression(mu <- rep(mean(y),length(y))),   
                 sigma.initial = expression(sigma<- rep(2, length(y))),
                 mu.valid = function(mu) all(mu > 0 & mu < 1), 
                 sigma.valid = function(sigma)  all(sigma > 0),
                 y.valid = function(y) all(y > 0 &  y < 1)
  ), 
  class = c("gamlss.family", "family"))
}

#------------------------------------------------------------------------------------------
# density function
dUG<-function(y, mu = 0.7, sigma = 2.1, log = FALSE)
{
  if (any(mu <= 0) | any(mu >= 1)) stop(paste("mu must be between 0 and 1", "\n", ""))
  if (any(sigma <= 0)) stop(paste("sigma must be positive", "\n", "")) 
  if (any(y <= 0) | any(y >= 1)) stop(paste("x must be between 0 and 1", "\n", ""))
  
  fy1 <- 1/y*dgamma(-log(y),sigma,mu^(1/sigma)/(1-mu^(1/sigma)))
  if(log==FALSE) fy<-fy1 else fy<-log(fy1)
  fy
}
# integrate(dUG,0,1) # checking the pdf
#------------------------------------------------------------------------------------------
# cumulative distribution function
pUG<-function(q, mu = 0.7, sigma = 2.1, lower.tail = TRUE, log.p = FALSE){
  if (any(mu <= 0) | any(mu >= 1)) stop(paste("mu must be between 0 and 1", "\n", ""))
  if (any(sigma < 0))  stop(paste("sigma must be positive", "\n", "")) 
  if (any(q <= 0) | any(q >= 1)) stop(paste("x must be between 0 and 1", "\n", ""))
  cdf1<-  1-pgamma(-log(q),sigma,mu^(1/sigma)/(1-mu^(1/sigma)))
  if(lower.tail==TRUE) cdf<-cdf1 else cdf<- 1-cdf1
  if(log.p==FALSE) cdf<- cdf else cdf<- log(cdf)
  cdf
}
# pUG(.5)
# integrate(dUG,0,.5) # checking the cdf with the pdf
#------------------------------------------------------------------------------------------
# quantile function
qUG<-function(u,mu,sigma)
{
  q<- exp(-qgamma(1-u,sigma,mu^(1/sigma)/(1-mu^(1/sigma))))
  q
}	
# u=pUG(.5)
# qUG(u,mu=.7,sigma=2.1) # checking the qf with the cdf
#------------------------------------------------------------------------------------------
# inversion method for randon generation
rUG<-function(n,mu,sigma)
{
  u<- runif(n)
  y<- qUG(u,mu =mu, sigma =sigma)
  y
}

# # Checking the results
library(gamlss)
set.seed(10)
n<-1000
# Case 1: without regressors
mu_true<-.7
sigma_true<-7
mu_result<-sigma_result<-c()
logit_link<-make.link("logit")
log_link<-make.link("log")
for (i in 1:100) {
  y<-rUG(n,mu_true,sigma_true)
  fit1<-gamlss(y~1, family="UG", trace = F)
  mu_result[i]<-logit_link$linkinv(fit1$mu.coefficients)
  sigma_result[i]<-log_link$linkinv(fit1$sigma.coefficients)
}

result1<- matrix(c(mu_true, mean(mu_result),
            sigma_true, mean(sigma_result)),2,2)
colnames(result1)<-c("mu","sigma")
rownames(result1)<-c("true value","mean")
print(round(result1,2))

# # Checking the results
set.seed(10)
n<-1000
# Case 2: with regressors
X<-runif(n)
logit_link<-make.link("logit")
log_link<-make.link("log")
b1<-.7
b2<-3
mu_true<-logit_link$linkinv(b1+b2*X)
g1<-.5
g2<-1.5
sigma_true<-log_link$linkinv(g1+g2*X)
R<-100
mu_result<-sigma_result<-matrix(NA,R,2)
for (i in 1:R) {
  y<-rUG(n,mu_true,sigma_true)
  fit1<-gamlss(y~X,sigma.formula =~ X, family=UG(), trace = F)
  mu_result[i,]<-fit1$mu.coefficients
  sigma_result[i,]<-fit1$sigma.coefficients
}

true_values<-c(b1,b2, g1,g2)
mean_values<-c(apply(mu_result,2,mean),
               apply(sigma_result,2,mean))
b_values<-(true_values-mean_values)/true_values*100
eqm_values<-c(apply(mu_result,2,var),
              apply(sigma_result,2,var))+(true_values-mean_values)^2
result1<- cbind(true_values,
                mean_values,
                b_values,
                eqm_values
)
colnames(result1)<-c("true value","mean","bias","eqm")
rownames(result1)<-c("b1","b2","g1","g2")
print(round(result1,2))
