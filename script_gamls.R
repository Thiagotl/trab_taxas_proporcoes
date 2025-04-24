library(gamlss)


loglikUGo <- expression(
  log(mu) + log(sigma) - (sigma + 1) * log(y) - mu * (y^(-sigma) - 1)
)

dldm_expr <- D(loglikUGo, "mu")
dldd_expr <- D(loglikUGo, "sigma")


#### DEFINIÇÃO DA FAMÍLIA UGo PARA GAMLSS ----


UGo <- function(mu.link = "logit", sigma.link = "log") {
  mstats <- checklink("mu.link", "UGo", substitute(mu.link),
                      c("logit", "probit", "cloglog", "cauchit", "log", "own"))
  dstats <- checklink("sigma.link", "UGo", substitute(sigma.link),
                      c("inverse", "log", "identity", "own"))
  
  structure(list(
    family = c("UGo", "Unit-Gompertz"),
    parameters = list(mu = TRUE, sigma = TRUE),
    nopar = 2,
    type = "Continuous",
    control = list(trace = FALSE),
    
    mu.link = as.character(substitute(mu.link)),
    sigma.link = as.character(substitute(sigma.link)),
    
    mu.linkfun = mstats$linkfun,
    sigma.linkfun = dstats$linkfun,
    
    mu.linkinv = mstats$linkinv,
    sigma.linkinv = dstats$linkinv,
    
    mu.dr = mstats$mu.eta,
    sigma.dr = dstats$mu.eta,
    
    dldm = function(y, mu, sigma) eval(dldm_expr),
    
    d2ldm2 = function(y, mu, sigma) {
      dldm <- eval(dldm_expr)
      d2 <- -dldm^2
      ifelse(d2 < -1e-15, d2, -1e-15)
    },
    
    dldd = function(y, mu, sigma) eval(dldd_expr),
    
    d2ldd2 = function(y, mu, sigma) {
      dldd <- eval(dldd_expr)
      d2 <- -dldd^2
      ifelse(d2 < -1e-15, d2, -1e-15)
    },
    
    d2ldmdd = function(y, mu, sigma) {
      dldm <- eval(dldm_expr)
      dldd <- eval(dldd_expr)
      d2 <- -dldm * dldd
      ifelse(is.na(d2), 0, d2)
    },
    
    G.dev.incr = function(y, mu, sigma, ...) {
      -2 * dUGo(y = y, mu = mu, sigma = sigma, log = TRUE)
    },
    
    rqres = expression(
      rqres(pfun = "pUGo", type = "Continuous", y = y, mu = mu, sigma = sigma)
    ),
    
    mu.initial = expression(mu <- rep(0.5, length(y))),
    sigma.initial = expression(sigma <- rep(1, length(y))),
    
    mu.valid = function(mu) all(mu > 0 & mu < 1),
    sigma.valid = function(sigma) all(sigma > 0),
    y.valid = function(y) all(y > 0 & y < 1)
  ), class = c("gamlss.family", "family"))
}

# FUNÇÃO DENSIDADE
dUGo <- function(y, mu, sigma, log = FALSE) {
  if (any(y <= 0 | y >= 1)) {
    return(if (log) rep(-Inf, length(y)) else rep(0, length(y)))
  }
  
  log_fx <- log(mu) + log(sigma) - (sigma + 1) * log(y) - mu * (y^(-sigma) - 1)
  
  if (log) {
    return(log_fx)
  } else {
    return(exp(log_fx))
  }
}

# FUNÇÃO ACUMULADA (CDF)

# pUGo <- function(q, mu, sigma, lower.tail = TRUE, log.p = FALSE) {
#   # Verificação do suporte
#   if (any(q <= 0 | q >= 1)) {
#     if (log.p) {
#       return(ifelse(q <= 0, -Inf, ifelse(q >= 1, 0, NA)))
#     } else {
#       return(ifelse(q <= 0, 0, ifelse(q >= 1, 1, NA)))
#     }
#   }
#   
#   # Cálculo estável em log-space
#   log_cdf <- -mu * (q^(-sigma) - 1)
#   
#   # Aplica lower.tail e log.p
#   if (!lower.tail) {
#     log_cdf <- log(1 - exp(log_cdf))
#   }
#   
#   if (log.p) {
#     return(log_cdf)
#   } else {
#     return(exp(log_cdf))
#   }
# }

pUGo<-function(q, mu, sigma, lower.tail = TRUE, log.p = FALSE)
{
  cdf<-exp(-mu*(q^(-sigma)-1))

  return(cdf)
}

# QUANTIL

qUGo <- function(u, mu, sigma) {
  ((-log(u) / mu) + 1)^(-1 / sigma)
}


# GERAÇÃO DE AMOSTRAS

rUGo <- function(n, mu, sigma) {
  u <- runif(n, min = 1e-6, max = 1 - 1e-6)
  ((-log(u) / mu) + 1)^(-1 / sigma)
}


# SIMULAÇÃO E AJUSTE DO MODELO SEM REGRESSORES

set.seed(10)
n <- 150
mu_true <- 0.1
sigma_true <- 0.5
mu_result <- sigma_result <- numeric(100)

logit_link <- make.link("logit")
log_link <- make.link("log")

for (i in 1:100) {
  y <- rUGo(n, mu_true, sigma_true)

  fit <- tryCatch(
    gamlss(y ~ 1, family = UGo(), method = RS(), control = gamlss.control(n.cyc = 100), trace = FALSE),
    error = function(e) NULL
  )

  if (!is.null(fit) &&
      !is.null(fit$mu.coefficients) && !is.null(fit$sigma.coefficients) &&
      length(fit$mu.coefficients) == 1 && length(fit$sigma.coefficients) == 1) {

    mu_result[i] <- logit_link$linkinv(fit$mu.coefficients)
    sigma_result[i] <- log_link$linkinv(fit$sigma.coefficients)
  } else {
    mu_result[i] <- NA
    sigma_result[i] <- NA
  }
}


# RESULTADOS

result <- matrix(c(mu_true, mean(mu_result, na.rm = TRUE),
                   sigma_true, mean(sigma_result, na.rm = TRUE)), 2, 2)

colnames(result) <- c("mu", "sigma")
rownames(result) <- c("true value", "mean")

print(round(result, 3))


#### COM REGRESSORES---- 

#Definir parâmetros globais
vn <- c(30, 70, 150, 300) # Tamanhos de amostra
logit_link <- make.link("logit")
log_link <- make.link("log")
b1 <- 0.2  # mu
b2 <- 0.5  # mu
g1 <- 0.5 # sigma
g2 <- 0.8  # sigma
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
    y <- rUGo(n, mu_true, sigma_true)
    
    fit <- tryCatch(
      gamlss(
        y ~ X, sigma.formula = ~ X,
        family = UGo(sigma.link = "log"),
        
        method = CG(),
        control = gamlss.control(n.cyc = 700),
        mu.step = 0.1,
        sigma.step = 0.1,
        trace = FALSE
      ),
      
      error = function(e) NULL
    )
    
    if (is.null(fit) ||
        any(is.na(fit$mu.coefficients)) ||
        any(is.na(fit$sigma.coefficients))) {
      bug_counter <- bug_counter + 1
      next
    }
    
    i <- i + 1
    mu_result[i, ] <- fit$mu.coefficients
    sigma_result[i, ] <- fit$sigma.coefficients
    setTxtProgressBar(pb, i)
  }
  # print(sigma_result)
  close(pb)

  # Calcular métricas
  true_values <- c(b1, b2, g1, g2)
  result <- calculate_metrics(mu_result, sigma_result, true_values)
  result <- as.data.frame(result)
  result$Sample_Size <- n

  final_results <- rbind(final_results, result)


  # Armazenar e exibir resultados
  cat("\nTamanho da amostra:", n, "\n")
  print(round(result, 2))
  cat("N de erros no ajuste do modelo:", bug_counter, "\n")}

