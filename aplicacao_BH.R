####### Abordagem Bayesiana ########
install.packages("expertsurv")
library(expertsurv)
######## ABORDAGEM BAYSIANA PARA OS DADOS DE CONTAGEM ###########

# diretorio
caminho_dados <- "C:/Users/melga/OneDrive/Desktop/IC ANA JULIA/testes/Data_CPDO.txt"

# base bruta
dados <- read.table(caminho_dados, header = TRUE)

XI_2 <- dados                       
bronq_aguda <- dados$bronq_aguda   


# Supondo que NO seja a primeira coluna de XI_2
X_real <- as.matrix(XI_2[, "NO", drop = FALSE]) 
y_real <- as.numeric(bronq_aguda)

# DICA: Padronizar o NO ajuda muito na convergência Bayesiana
X_real_scale <- scale(X_real)


ajuste.bronquite.bayes <- function(y, X, niter = 20000) {
  library(nimble)
  
  n <- length(y)
  num_betas <- ncol(X)
  
  code <- nimbleCode({
    # Prioris
    alpha ~ dnorm(0, sd = 1)      # Intercepto
    beta_NO ~ dnorm(0, sd = 1)    # Efeito do NO
    phi ~ dunif(-0.99, 0.99)      # Memória temporal
    r_size ~ dunif(0, 50)         # Dispersão
    
    # Inicialização da recursão
    z[1] <- 0
    e[1] <- 0
    
    for (t in 1:n) {
      # Estrutura AR(1)
      z[t + 1] <- phi * (z[t] + e[t])
      
      # Preditor: Intercepto + NO * beta + Memória
      # Usamos X[t, 1] porque passamos apenas o NO na matriz
      log(mu[t]) <- alpha + beta_NO * X[t, 1] + z[t + 1]
      
      # Likelihood NegBin
      y[t] ~ dnegbin(prob = r_size / (r_size + mu[t]), size = r_size)
      
      # Resíduo de Pearson para alimentar o próximo passo
      var_t[t] <- mu[t] + (mu[t]^2 / r_size)
      e[t + 1] <- (y[t] - mu[t]) / sqrt(var_t[t] + 0.1)
    }
  })
  
  dados <- list(y = y, X = X)
  cts <- list(n = n)
  
  # Valores iniciais para o MCMC
  ini <- list(alpha = mean(log(y + 0.5)), beta_NO = 0, phi = 0.1, r_size = 1)
  
  samples <- nimbleMCMC(code = code, constants = cts, data = dados, inits = ini,
                        niter = niter, nburnin = floor(niter*0.3),
                        monitors = c("alpha", "beta_NO", "phi", "r_size"),
                        thin = 10) # Thinning ajuda a reduzir autocorrelação
  return(samples)
}


# Rodando o modelo
resultado_bronquite <- ajuste.bronquite.bayes(y_real, X_real)

# Resumo e Intervalos

summary(resultado_bronquite)
cred_int(resultado_bronquite, type = "HPD")
colMeans(resultado_bronquite)
