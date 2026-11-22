###### NOVA SIMULAÇAO NUCLEO FINAL BINOMIAL NEGATIVA #######

#install.packages(c("data.table", "nimble", "pkgbuild"))
library(nimble)
library(data.table)
library(parallel)

#pkgbuild::has_build_tools(debug = TRUE)

# 1. Configurações de Caminho e Cenários --------------------------------------
diretorio_base <- "C:/Users/melga/OneDrive/Desktop/IC ANA JULIA/testes"
phi_cenarios <- c(0.8)   # Adicionar os outros valores (-0.8, -0.3, 0.3) depois
n_replicas <- 1000
n_obs <- 200

# Número de processos paralelos (deixa 1 núcleo livre para o Windows)
n_cores <- max(1, detectCores(logical = FALSE) - 1)

set.seed(222)
X_fixo <- matrix(rnorm(n_obs * 2), ncol = 2)

# Parâmetros verdadeiros (os mesmos de antes)
alpha_real <- 1.0 #beta 0
betas_real <- c(0.5, -0.1)
r_real     <- 5.0 #alpha

# 2. Função de Geração de Dados (DGP - Lógica Davis/Co-autor) -----------------
gera_dados_glarma_davis <- function(n, phi, alpha, betas, X, r_disp) {
  y <- numeric(n)
  z <- numeric(n + 1)
  e <- numeric(n + 1)
  eta_fixo <- as.numeric(X %*% betas)
  
  for (t in 1:n) {
    # Recursão Davis: Z_t depende da persistência do estado + erro anterior
    z[t+1] <- phi * (z[t] + e[t])
    mu_t <- exp(alpha + eta_fixo[t] + z[t+1])
    y[t] <- rnbinom(1, size = r_disp, mu = mu_t)
    
    # Resíduo de Pearson conforme Davis (2003)
    var_t <- mu_t + (mu_t^2 / r_disp)
    e[t+1] <- (y[t] - mu_t) / sqrt(var_t + 0.1)
  }
  return(y)
}

# 3. O Motor do Bloco (NIMBLE compilado UMA vez por worker) -------------------
# Cada worker recebe um vetor de ids, compila o modelo uma única vez e depois
# só troca os dados y a cada réplica.
rodar_bloco <- function(ids, n, p_val, X_fixo, pasta_dest,
                        alpha_real, betas_real, r_real, worker_id) {
  library(nimble)
  library(data.table)
  
  log_arq <- file.path(pasta_dest, paste0("log_worker_", worker_id, ".txt"))
  log_msg <- function(...) {
    cat(format(Sys.time(), "[%d/%m %H:%M] "), ..., "\n",
        file = log_arq, append = TRUE)
  }
  
  # Mesma estrutura do modelo original
  code <- nimbleCode({
    alpha ~ dnorm(0, sd = 1)         # Priori menos restritiva
    beta[1] ~ dnorm(0, sd = 1)
    beta[2] ~ dnorm(0, sd = 1)
    r_size ~ dunif(0, 50)
    
    # Priori Uniforme: Essencial para não "puxar" o phi para o zero
    phi ~ dunif(-0.99, 0.99)
    
    z[1] <- 0
    e[1] <- 0
    for (t in 1:n) {
      # Recursão Exata Davis et al. (2003)
      z[t+1] <- phi * (z[t] + e[t])
      log(mu[t]) <- alpha + beta[1]*X_fixo[t,1] + beta[2]*X_fixo[t,2] + z[t+1]
      
      y[t] ~ dnegbin(prob = r_size / (r_size + mu[t]), size = r_size)
      
      # Resíduo de Pearson dinâmico
      var_t[t] <- mu[t] + (mu[t]^2 / r_size)
      e[t+1] <- (y[t] - mu[t]) / sqrt(var_t[t] + 0.1)
    }
  })
  
  inits <- list(alpha = 0, beta = c(0, 0), r_size = 2, phi = 0.0)
  
  # --- Compilação (acontece só aqui, uma vez por worker) ---
  set.seed(ids[1])
  y_ini <- gera_dados_glarma_davis(n, p_val, alpha_real, betas_real, X_fixo, r_real)
  
  modelo <- nimbleModel(code = code,
                        constants = list(n = n, X_fixo = X_fixo),
                        data = list(y = y_ini),
                        inits = inits)
  Cmodelo <- compileNimble(modelo)
  conf <- configureMCMC(modelo,
                        monitors = c("alpha", "beta", "phi", "r_size"))
  mcmc <- buildMCMC(conf)
  Cmcmc <- compileNimble(mcmc, project = modelo)
  log_msg("Modelo compilado. Réplicas deste worker:", length(ids))
  
  n_ok <- 0
  for (id in ids) {
    tryCatch({
      # Semente por réplica: o mesmo id gera o mesmo dataset (útil p/ comparar com a GLARMA)
      set.seed(id)
      y_sim <- gera_dados_glarma_davis(n, p_val, alpha_real, betas_real, X_fixo, r_real)
      
      # Troca os dados no modelo compilado e reinicia os valores iniciais
      Cmodelo$y <- y_sim
      Cmodelo$alpha <- inits$alpha
      Cmodelo$beta  <- inits$beta
      Cmodelo$r_size <- inits$r_size
      Cmodelo$phi   <- inits$phi
      Cmodelo$calculate()   # atualiza nós determinísticos (z, e, mu, var_t)
      
      samples <- runMCMC(Cmcmc,
                         niter = 6000, nburnin = 500,
                         progressBar = FALSE,
                         setSeed = id)
      
      # Escrita segura: grava em .tmp e renomeia (evita CSV corrompido se cair)
      nome_arq <- file.path(pasta_dest, paste0("rep_", sprintf("%03d", id), ".csv"))
      tmp_arq  <- paste0(nome_arq, ".tmp")
      fwrite(as.data.frame(samples), tmp_arq)
      file.rename(tmp_arq, nome_arq)
      
      n_ok <- n_ok + 1
      log_msg("Réplica", id, "salva com sucesso.")
    }, error = function(e) {
      log_msg("Erro na", id, ":", conditionMessage(e))
    })
    
    if (n_ok %% 10 == 0) gc()   # limpa a memória RAM
  }
  return(n_ok)
}

# 4. Pastas e sistema de retomada ---------------------------------------------
caminho_phi_08 <- file.path(diretorio_base, "phi_08")
if (!dir.exists(caminho_phi_08)) dir.create(caminho_phi_08)

# Remove restos de escrita interrompida
file.remove(list.files(caminho_phi_08, pattern = "\\.tmp$", full.names = TRUE))

arquivos_na_pasta <- list.files(caminho_phi_08, pattern = "^rep_.*\\.csv$")
ids_que_ja_existem <- as.numeric(gsub("rep_([0-9]+)\\.csv", "\\1", arquivos_na_pasta))
faltam_de_verdade <- setdiff(1:n_replicas, ids_que_ja_existem)

cat("Verificação concluída. Faltam exatamente", length(faltam_de_verdade), "réplicas.\n")

# 5. Rodar em paralelo --------------------------------------------------------
if (length(faltam_de_verdade) > 0) {
  n_cores <- min(n_cores, length(faltam_de_verdade))
  
  # Divide os ids de forma intercalada entre os workers
  blocos <- split(faltam_de_verdade, rep_len(1:n_cores, length(faltam_de_verdade)))
  
  cl <- makeCluster(n_cores)
  clusterExport(cl, c("gera_dados_glarma_davis", "rodar_bloco"))
  
  cat("Iniciando", n_cores, "workers. Acompanhe o progresso nos arquivos log_worker_*.txt em:\n",
      caminho_phi_08, "\n") #MUDAR AQUI
  
  # MUDAR OS PHIS AQUI: p_val é o phi (0.3)
  resultado <- clusterMap(
    cl, rodar_bloco,
    ids = blocos,
    worker_id = seq_along(blocos),
    MoreArgs = list(n = n_obs, p_val = 0.8, X_fixo = X_fixo, 
                    pasta_dest = caminho_phi_08,
                    alpha_real = alpha_real, betas_real = betas_real,
                    r_real = r_real),
    SIMPLIFY = FALSE
  )
  stopCluster(cl)
  
  cat("\nRéplicas salvas nesta execução:", sum(unlist(resultado)), "\n")
}

arquivos_na_pasta <- list.files(caminho_phi_08, pattern = "^rep_.*\\.csv$")
cat("\n\n================================================")
cat("\nTotal de réplicas na pasta:", length(arquivos_na_pasta), "de", n_replicas)
cat("\n================================================\n")

# 6. Consolidação (igual ao script original) ----------------------------------
arquivos_salvos <- list.files(caminho_phi_08, pattern = "^rep_.*\\.csv$", full.names = TRUE)
cat("Lendo", length(arquivos_salvos), "réplicas...\n")

ler_resumo <- function(arq) {
  dados <- fread(arq)
  return(colMeans(dados))
}

lista_medias <- lapply(arquivos_salvos, ler_resumo)
resultados_finais <- do.call(rbind, lista_medias)

valores_reais <- c(alpha = 1.0, "beta[1]" = 0.5, "beta[2]" = -0.1, phi = 0.8, r_size = 5.0)# Garante a mesma ordem das colunas dos resultados
valores_reais <- valores_reais[colnames(resultados_finais)]
estimativas_medias <- colMeans(resultados_finais)

tabela_paper <- data.frame(
  Real = valores_reais,
  Estimado = estimativas_medias,
  Vies = estimativas_medias - valores_reais,
  EQM = apply(resultados_finais, 2, var) + (estimativas_medias - valores_reais)^2
)

print(round(tabela_paper, 4))

tabela_phi08 <- tabela_paper
#tabela_phi03 <- tabela_paper

# 7. Diagnósticos de uma réplica ----------------------------------------------
amostras <- fread(file.path(caminho_phi_08, "rep_001.csv"))

plot(as.ts(amostras$phi), main = "Trace Plot de Phi", ylab = "Valor", col = "blue")
plot(as.ts(amostras))
plot(density(amostras$phi), main = "Densidade Posterior de Phi", col = "red", lwd = 2)
abline(v = 0.8, col = "black", lty = 2)   # Valor real (0.8 nesta pasta)
acf(amostras$phi, main = "Autocorrelação de Phi")

