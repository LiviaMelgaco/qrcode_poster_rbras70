# Modelos GLARMA, segundo teste IC

library(pacman)
p_load(glarma, rpart, rpart.plot,tidyr, ggplot2)

# diretorio
caminho_dados <- "C:/Users/melga/OneDrive/Desktop/IC ANA JULIA/Data_CPDO.txt"

# base bruta
dados <- read.table(caminho_dados, header = TRUE)

#matriz variavel resposta
y_bronq <- dados$bronq_aguda


# criando um for para fazer todas as matrizes de variaveis explicativas
#(tem que ter intercepto - coluna de 1's) 
poluentes <- c("CO", "PM10", "NO", "NO2", "NOX", "O3")
lista_matrizes_x <- list()

for (poluente in poluentes) {
  x <- model.matrix(~ dados[[poluente]] + sen12 + cos12 + sen6 + cos6 , data = dados)
  
  lista_matrizes_x[[poluente]] <- x
} 


###### CO (2,0,0)x(1,0,0) #####
modelo_co_ar1 <- glarma(y = y_bronq, 
                    X = lista_matrizes_x$CO, 
                    type = "Poi", 
                    method = "FS",          
                    residuals = "Pearson",  
                    phiLags = c(1,2,12))         

summary(modelo_co_ar1)
plot(modelo_co_ar1$residuals)
acf(modelo_co_ar1$residuals)
pacf(modelo_co_ar1$residuals)

plot(modelo_co_ar1)

##### CO (2,0,0) ####
modelo_co_ar2 <- glarma(y = y_bronq, 
                        X = lista_matrizes_x$CO, 
                        type = "Poi", 
                        method = "FS",          
                        residuals = "Pearson",  
                        phiLags = c(1,2))         

summary(modelo_co_ar2)
plot(modelo_co_ar2$residuals)
acf(modelo_co_ar2$residuals)
pacf(modelo_co_ar2$residuals)

plot(modelo_co_ar2)
##### PM10 ####
modelo_pm10 <- glarma(y = y_bronq, 
                        X = lista_matrizes_x$PM10, 
                        type = "Poi", 
                        method = "FS",          
                        residuals = "Pearson",  
                        phiLags = c(1)) 
summary(modelo_pm10)
plot(modelo_pm10$residuals)
acf(modelo_pm10$residuals)
pacf(modelo_pm10$residuals)

plot(modelo_pm10)

##### NO #####
modelo_no <- glarma(y = y_bronq,
                    X = lista_matrizes_x$NO,
                    type = "Poi", 
                    method = "FS",          
                    residuals = "Pearson",  
                    phiLags = c(1)) 
summary(modelo_no)
ts.plot(modelo_no$residuals)
acf(modelo_no$residuals)
pacf(modelo_no$residuals)

plot(modelo_no)

##### NO2 (2,0,0)x(1,0,0) #####
modelo_no2 <- glarma(y = y_bronq,
                    X = lista_matrizes_x$NO2,
                    type = "Poi", 
                    method = "FS",          
                    residuals = "Pearson",  
                    phiLags = c(1,2,12)) 
summary(modelo_no2)
plot(modelo_no2$residuals)
acf(modelo_no2$residuals)
pacf(modelo_no2$residuals)

##### NOX #####
modelo_nox <- glarma(y = y_bronq,
                     X = lista_matrizes_x$NOX,
                     type = "Poi", 
                     method = "FS",          
                     residuals = "Pearson",  
                     phiLags = c(1)) 
summary(modelo_nox)
plot(modelo_nox$residuals)
acf(modelo_nox$residuals)
pacf(modelo_nox$residuals)

##### O3 #####
modelo_o3 <- glarma(y = y_bronq,
                     X = lista_matrizes_x$O3,
                     type = "Poi", 
                     method = "FS",          
                     residuals = "Pearson",
                     phiLags = c(12),
                     thetaLags = c(1)) 
summary(modelo_o3)
plot(modelo_o3$residuals)
acf(modelo_o3$residuals)
pacf(modelo_o3$residuals)


###### GRÁFICOS E PAINEIS UNI #####

glarma_uni_poluentes <- list(modelo_co_ar1,modelo_pm10, modelo_no, modelo_no2,
                             modelo_nox, modelo_o3)

## TABELA DE AIC
tabela_uni_aic <- data.frame()

for (i in seq_along(poluentes)) {
  nome <- poluentes[[i]]
  modelo_glarma <- glarma_uni_poluentes[[i]]
  aic <- modelo_glarma$aic
  
  tabela_uni_aic <- rbind(tabela_uni_aic, 
                          data.frame(Poluente = nome, AIC = aic))
}

tabela_uni_aic
saveRDS(tabela_uni_aic, "tabela_uni_aic.rds")

## PAINEL DE RESIDUOS histograma e nuvem de pontos

paineis_resi_uni_glarma <- list()

for (i in seq_along(poluentes)) {
  nome <- poluentes[[i]]
  modelo_glarma <- glarma_uni_poluentes[[i]]
  
  df_residuos <- data.frame(
    indice = seq_along(modelo_glarma$residuals),
    residuo = modelo_glarma$residuals
  )
  
  p1 <- ggplot(df_residuos, aes(x = indice, y = residuo)) +
    geom_point() +
    geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
    labs(title = paste("Dispersão e Histograma dos Resíduos de", nome),
         x = "Observação",
         y = "Resíduo") +
    theme_minimal()
  
  p2 <- ggplot(df_residuos, aes(x = residuo)) +
    geom_histogram(bins = 20, fill = "steelblue", color = "white") +
    labs(title = "",
         x = "Resíduo",
         y = "Frequência") +
    theme_minimal()
  
  paineis_resi_uni_glarma[[nome]] <- p1 + p2
}


paineis_resi_uni_glarma[["O3"]]
hist(modelo_o3$residuals)


saveRDS(paineis_resi_uni_glarma, "paineis_resi_uni_glarma.rds")

## ACF E PACF RESIDUOS

acf_pacf_unires <- list()

for (i in seq_along(poluentes)) {
  nome <- poluentes[i]
  # vai guardar o modelo da vez
  modelo_poluente <- glarma_uni_poluentes[[i]]
  
  dados_grafico <- data.frame(residuos = modelo_poluente$residuals)
  
  #graficos dos resíduos
  
  acf <- ggAcf(modelo_poluente$residuals) +
    ggtitle(paste("ACF e PACF dos Resíduos GLARMA", nome)) +
    theme_minimal()
  
  pacf <- ggPacf(modelo_poluente$residuals) +
    ggtitle("") +
    theme_minimal()
  
  #arima <- forecast:::arima.string(modelo_poluente, padding = FALSE)
  #paineis
  
  acf_pacf_residuo <- acf + pacf 
  
  #listas que vao guardar as series
  acf_pacf_unires[[i]] <- acf_pacf_residuo

}

acf_pacf_unires[[1]]

saveRDS(acf_pacf_unires, "acf_pacf_unires.rds")

## GRÁFICOS DO MODELO

graficos_glarma_uni <- list()


for (i in seq_along(poluentes)) {
  nome <- poluentes[[i]]
  modelo_glarma <- glarma_uni_poluentes[[i]]
  
  df_plot <- data.frame(
    Time = seq_along(y_bronq),
    Observado = y_bronq,
    GLARMA = (modelo_glarma$mu) 
  )
  
  p <- ggplot(df_plot, aes(x = Time)) +
    geom_line(aes(y = Observado, color = "Observado", linetype = "Observado"), linewidth = 0.6) +
    geom_line(aes(y = GLARMA, color = "GLARMA", linetype = "GLARMA"), linewidth = 0.6) +
    scale_color_manual(name = "", values = c(Observado = "black", GLARMA = "red")) +
    scale_linetype_manual(name = "", values = c(Observado = "dashed", GLARMA = "solid")) +
    labs(title = paste("Observado vs GLARMA:", nome),
         x = "Tempo", y = "Contagem") +
    theme_minimal() +
    theme(legend.position = "top")
  
  graficos_glarma_uni[[nome]] <- p
}

graficos_glarma_uni[["NO"]]

saveRDS(graficos_glarma_uni, "glarmas_uni.rds")

##### TESTES MULTIVARIADO #####

summary(dados$bronq_aguda)
var(dados$bronq_aguda)

x <- model.matrix(~ NO + CO + sen12 + cos12 + sen6 + cos6,
                  data = dados)
#(2,0,0)x(0,1,1) -> usando a estimaçao de paremntros do arima
modelo_mult1 <- glarma(y = y_bronq,
                       X = x,
                       type = "NegBin",
                       method = "FS",
                       residuals =  "Pearson",
                       phiLags = c(1,6))
                       #thetaLags = c(12)) 
summary(modelo_mult1)

plot(modelo_mult1$residuals)
acf(modelo_mult1$residuals)
pacf(modelo_mult1$residuals)


modelo_mult3 <- glarma(y = y_bronq,
                       X = x,
                       type = "Poi", # comparando com poisson
                       method = "FS",
                       residuals =  "Pearson",
                       phiLags = c(1,6))
#thetaLags = c(12)) 
summary(modelo_mult3)

plot(modelo_mult1$residuals)
acf(modelo_mult3$residuals)
pacf(modelo_mult3$residuals)



################# ARVORE DE DECISAO VARIAVEIS ######
arvore_bronq <- rpart(
  formula = bronq_aguda ~ CO + PM10 + NO + NO2 + NOX + O3,
  data = dados,
  method = "anova"
)

rpart.plot(arvore_bronq, 
           type = 3, 
           extra = 101, 
           main = "Variáveis Mais Importantes para Internações")
summary(arvore_bronq)


# Converte o vetor nomeado em um data frame estruturado
var_importance_gl <- data.frame(
  Poluente = names(arvore_bronq$variable.importance),
  Importância = as.numeric(arvore_bronq$variable.importance)
)

saveRDS(var_importance_gl, "var_imp_gl.RDS")
######### MODELO MULTI COM ARVORE DE DECISÃO ######


x_tree <- model.matrix(~ NO + sen12 + cos12 + sen6 + cos6,
                  data = dados)
#(2,0,0)x(0,1,1)
modelo_mult2 <- glarma(y = y_bronq,
                       X = x_tree,
                       type = "Poi", #se eu mudo pra poisson, 
                       method = "FS",    # o p valor fica ok, mas meu aic sobre muito
                       residuals =  "Pearson",
                       phiLags = c(1,6))
                       
#thetaLags = c(12)) 
summary(modelo_mult2)

plot(modelo_mult2$residuals)
acf(modelo_mult2$residuals)
pacf(modelo_mult2$residuals)
plot.glarma(modelo_mult2)


######## graficos de residuos

df_multres <- data.frame(
  indice = seq_along(modelo_mult2$residuals),
  residuo = modelo_mult2$residuals
)

resmult1 <- ggplot(df_multres, aes(x = indice, y = residuo)) +
  geom_point() +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  labs(title = paste("Dispersão e Histograma dos Resíduos do Modelo Multivariado."),
       x = "Observação",
       y = "Resíduo") +
  theme_minimal()

resmult2 <- ggplot(df_multres, aes(x = residuo)) +
  geom_histogram(bins = 20, fill = "steelblue", color = "white") +
  labs(title = "",
       x = "Resíduo",
       y = "Frequência") +
  theme_minimal()

res_multi_glarma <- resmult1 + resmult2

saveRDS(res_multi_glarma, "res_mult_glarma.RDS")

######## com base na matriz de correlaçao
x4 <- model.matrix(~ O3 + NO + sen12 + cos12 + sen6 + cos6,
                      data = dados)

modelo_mult4 <- glarma(y = y_bronq,
                       X = x_tree,
                       type = "NegBin", #se eu mudo pra poisson, 
                       method = "FS",    # o p valor fica ok, mas meu aic sobre muito
                       residuals =  "Pearson",
                       phiLags = c(1,6))

#thetaLags = c(12)) 
summary(modelo_mult4)

plot(modelo_mult2$residuals)
acf(modelo_mult2$residuals)
pacf(modelo_mult2$residuals)
plot.glarma(modelo_mult4)

#os dados foram tirados da toposfera?

##### TABELAS E GRÁFICOS MULTI #####

## acrescentar aic do multi na tabela de aic

tbl_aic_geral <- rbind(tabela_uni_aic, data.frame(Poluente = "MULTI", AIC = modelo_mult2$aic))
tbl_aic_geral <- tbl_aic_geral[order(tbl_aic_geral$AIC), ]

saveRDS(tbl_aic_geral, "tbl_aic_geral.RDS")

## acrescentar modelo multi na lista de modelos uni para comparação de seus coeficientes

modelos_glarma_geral <- list(modelo_co_ar1,modelo_pm10, modelo_no, modelo_no2,
                             modelo_nox, modelo_o3, modelo_mult2)
# ACF E PACF
acf_mult <- ggAcf(modelo_mult2$residuals) +
  ggtitle(paste("ACF e PACF dos Resíduos GLARMA Multivariado")) +
  theme_minimal()

pacf_mult <- ggPacf(modelo_mult2$residuals) +
  ggtitle("") +
  theme_minimal()

acf_pacf_mult <- acf_mult + pacf_mult

saveRDS(acf_pacf_mult, "grf_acfpacf_mult.RDS")


# GRAFICOS DO MODELO: OBS VS MODEL

df_plot_mult <- data.frame(
    Time = seq_along(y_bronq),
    Observado = y_bronq,
    GLARMA = (modelo_mult2$mu) 
)

plot_mult <- ggplot(df_plot_mult, aes(x = Time)) +
          geom_line(aes(y = Observado, color = "Observado", linetype = "Observado"), linewidth = 0.6) +
          geom_line(aes(y = GLARMA, color = "GLARMA", linetype = "GLARMA"), linewidth = 0.6) +
          scale_color_manual(name = "", values = c(Observado = "black", GLARMA = "red")) +
          scale_linetype_manual(name = "", values = c(Observado = "dashed", GLARMA = "solid")) +
          labs(title = paste("Observado vs GLARMA Multivariado"),
               x = "Tempo", y = "Contagem") +
          theme_minimal() +
          theme(legend.position = "top")

saveRDS(plot_mult, "plot_mult.RDS")

###### TABELAS DOS MODELOS #####

resumo_mult <- summary(modelo_mult2)
coef_mult <- as.data.frame(resumo_mult$coefficients1)

saveRDS(coef_mult, "tbl_coef_mult.RDS")

resumo_uni <- summary(modelo_no)
coef_uni <- as.data.frame(resumo_uni$coefficients1)
rownames(coef_uni)[2] <- "NO"

saveRDS(coef_uni, "tbl_coef_uni.RDS")


###### GRAFICOS DE COMPARAÇAO DOS MODELOS #####
df_plot_comp <- data.frame(
  Time = seq_along(y_bronq),
  Observado = y_bronq,
  Uni = (modelo_no$mu),
  Multi = (modelo_mult2$mu) 
)

plot_comparacao <- ggplot(df_plot_comp, aes(x = Time)) +
        geom_line(aes(y = Observado, color = "Observado", linetype = "Observado"), linewidth = 0.6) +
        geom_line(aes(y = Uni, color = "Uni", linetype = "Uni"), linewidth = 0.6) +
        geom_line(aes(y = Multi, color = "Multi", linetype = "Multi"), linewidth = 0.9) +
        scale_color_manual(name = "", values = c(Observado = "black", Uni = "red", Multi = "blue")) +
        scale_linetype_manual(name = "", values = c(Observado = "dashed", Uni = "solid", Multi = "dotted")) +
        labs(title = paste("Observado vs Uni vs Multi"),
             x = "Tempo", y = "Contagem") +
        theme_minimal() +
        theme(legend.position = "top")

saveRDS(plot_comparacao, "plot_comparacao.RDS")
