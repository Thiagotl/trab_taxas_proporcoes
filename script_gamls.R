
UGo <- expression(log(mu)+log(sigma)-(sigma + 1) * log(y)-mu * (y^(-sigma) - 1))

m1UGo<-D(UGo,"mu")
s1UGo<-D(UGo,"sigma")
ms2UGo<-D(m1UGo,"sigma")

UGo<-function (mu.link = "logit", sigma.link = "log")
{
  mstats <- checklink("mu.link", "UGo", substitute(mu.link),
                      c("logit", "probit", "cloglog", "cauchit", "log", "own"))
  dstats <- checklink("sigma.link", "UGo", substitute(sigma.link),
                      c("inverse", "log", "identity", "own"))
  structure(list(family = c("UGo", "Unit-Gompertz"),
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
                   tau=.5
                   dldm <- eval(m1UGo)
                   dldm
                 },
                 d2ldm2 = function(y,mu, sigma) {
                   tau=.5
                   dldm <- eval(m1UGo)
                   d2ldm2 <- -dldm * dldm
                   d2ldm2 <- ifelse(d2ldm2 < -1e-15, d2ldm2,-1e-15)
                   d2ldm2
                 },
                 dldd = function(y, mu, sigma) {
                   tau=.5
                   dldd <- eval(s1UGo)
                   dldd
                 },
                 d2ldd2 = function(y,mu, sigma) {
                   tau=.5
                   dldd <- eval(s1UGo)
                   d2ldd2 = -dldd * dldd
                   d2ldd2 <- ifelse(d2ldd2 < -1e-15, d2ldd2,-1e-15)
                   d2ldd2
                 },
                 d2ldmdd = function(y,mu, sigma) {
                   tau=.5
                   dldm <- eval(m1UGo)
                   dldd <- eval(s1UGo)
                   d2ldmdd = -(dldm * dldd)
                   d2ldmdd<-ifelse(is.na(d2ldmdd)==TRUE,0,d2ldmdd)
                   d2ldmdd
                 },
                 G.dev.incr = function(y, mu, sigma, w, ...) -2 * log(dUGo(y=y, mu=mu, sigma=sigma)),
                 rqres = expression(
                   rqres(pfun = "pUGo", type = "Continuous", y = y, mu = mu, sigma = sigma)
                 ),
                 mu.initial = expression(mu <- rep(mean(y),length(y))),
                 sigma.initial = expression(sigma<- rep(5, length(y))),
                 mu.valid = function(mu) all(mu > 0 & mu < 1),
                 sigma.valid = function(sigma) all(sigma > 0),
                 y.valid = function(y) all(y > 0 & y < 1)
  ),
  class = c("gamlss.family", "family"))
}

#-------------------------------------------------------------------------------

# density function
dUGo<-function(y, mu, sigma, log = FALSE)
{
  fx<-mu*sigma * (y^-(sigma+1)) * exp(-mu*(y^(-sigma)-1))
  
  return(fx)
}

#-------------------------------------------------------------------------------
# cumulative distribution function
pUGo<-function(q, mu, sigma, lower.tail = TRUE, log.p = FALSE)
{
  cdf<-exp(-mu*(q^(-sigma)-1))
  
  return(cdf)
}


#-------------------------------------------------------------------------------
# quantile function
qUGo<-function(u,mu, sigma)
{
  q<-((-log(u)/mu)+1)^(-1/sigma)
  
  return(q)
}

#-------------------------------------------------------------------------------
# inversion method for randon generation

rUGo<-function(n,mu, sigma){
  u <- runif(n)
  y<-((-log(u)/mu)+1)^(-1/sigma)
  return(y)
}

library(gamlss)

# Definir parâmetros globais
vn <- c(30, 70, 150, 300) # Tamanhos de amostra
logit_link <- make.link("logit")
log_link <- make.link("log")
b1 <- 0.7  # mu
b2 <- 0.3  # mu
g1 <- .7 # sigma
g2 <- .25  # sigma
R <- 1000  # Número de repetições

set.seed(10)

# Função auxiliar para calcular métricas
calculate_metrics <- function(mu_result, sigma_result, true_values) {
  mean_values <- c(
    apply(mu_result, 2, mean, na.rm = TRUE),
    apply(sigma_result, 2, mean, na.rm = TRUE)
  )
  bias_values <- (true_values - mean_values) / true_values * 100
  eqm_values <- c(
    apply(mu_result, 2, var, na.rm = TRUE),
    apply(sigma_result, 2, var, na.rm = TRUE)
  ) + (true_values - mean_values)^2
  
  result <- cbind(
    "True Value" = true_values,
    "Mean" = mean_values,
    "Bias (%)" = bias_values,
    "EQM" = eqm_values
  )
  rownames(result) <- c("b1", "b2", "g1", "g2")
  return(result)
}

# Inicialização de resultados finais
#final_results <- list()
bug_counter <- 0  
final_results <- data.frame()

for (n in vn) {
  
  mu_result <- matrix(NA, R, 2)  # Para b1, b2
  sigma_result <- matrix(NA, R, 2)  # Para g1, g2
  X <- runif(n)
  
  mu_true <- logit_link$linkinv(b1+b2*X)
  # mean(mu_true)
  sigma_true <- log_link$linkinv(g1+g2*X)
  # mean(sigma_true)
  pb <- txtProgressBar(min = 0, max = R, style = 3)
  
  i<-0
  while (i < R) {
    y <- rUGo(n, mu_true, sigma_true)  # Geração dos dados
    
    # Ajustar o modelo 
    fit1 <- try(
      gamlss(
        y ~ X, sigma.formula = ~ X,
        family = UGo(sigma.link = "log"),
        # sigma.start = start,
        c.crit = 0.001,
        n.cyc = 700,
        mu.step = 0.1,
        sigma.step = 0.1,
        trace = FALSE
      ),
      silent = TRUE
    )
    
    if (inherits(fit1, "try-error")) {
      bug_counter <- bug_counter + 1
      next
    }
    
    i<-i+1
    # print(i)
    #  coeficientes ajustados
    mu_result[i, ] <- fit1$mu.coefficients
    sigma_result[i, ] <- fit1$sigma.coefficients
    
    setTxtProgressBar(pb, i)
  }
  # print(sigma_result)
  close(pb)
  
  # Calcular métricas
  true_values <- c(b1, b2, g1, g2)
  result <- calculate_metrics(mu_result, sigma_result, true_values)
  result <- as.data.frame(result)
  result$Sample_Size <- n
  
  # result_file <- paste0("simulation_results_n_", n, ".txt")
  # write.table(round(result, 2), file = result_file, sep = "\t", col.names = NA, quote = FALSE)
  # final_results[[as.character(n)]] <- result
  
  
  final_results <- rbind(final_results, result)
  
  
  # Armazenar e exibir resultados
  #final_results[[as.character(n)]] <- result
  cat("\nTamanho da amostra:", n, "\n")
  print(round(result, 2))
  cat("N de erros no ajuste do modelo:", bug_counter, "\n")
}

# Exibir contadores de erros
# cat("\nNúmero total de erros no ajuste do modelo:", bug_counter, "\n")


