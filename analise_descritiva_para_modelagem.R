# ==============================================================================
# ANÁLISE DESCRITIVA APROFUNDADA: O QUE O CÓDIGO FAZ E COMO MELHORAR
# ==============================================================================

rm(list = ls())

# Carregar dados
source("Covariaveis.R", encoding = "UTF-8")

# Pacotes necessários
pkgs <- c("ggplot2", "dplyr", "tidyr", "corrplot", "gridExtra", "GGally")
for (p in pkgs) {
  if (!require(p, character.only = TRUE)) install.packages(p)
  library(p, character.only = TRUE)
}

if (!dir.exists("graficos")) dir.create("graficos")
if (!dir.exists("relatorios")) dir.create("relatorios")

cat("================================================================================\n")
cat("ANÁLISE DESCRITIVA PARA MODELAGEM - GUIA PASSO A PASSO\n")
cat("================================================================================\n\n")

# ==============================================================================
# PASSO 1: ENTENDER A ESTRUTURA DOS DADOS
# ==============================================================================

cat("PASSO 1: ESTRUTURA DOS DADOS\n")
cat("================================================================================\n\n")

cat("Dimensões:\n")
cat("- Y_mat (Óbitos): ", nrow(Y_mat), " microrregiões x ", ncol(Y_mat), " anos\n")
cat("- E (Offset/Nascidos): ", nrow(E), " microrregiões x ", ncol(E), " anos\n")
cat("- x (3 Covariáveis): ", dim(x)[1], " microrregiões x ", dim(x)[2], " anos x ", dim(x)[3], " variáveis\n")

cat("\nTaxa de eventos (Y/E):\n")
taxa_eventos <- Y_mat / E
print(summary(as.vector(taxa_eventos)))

# ==============================================================================
# PASSO 2: VERIFICAR DISTRIBUIÇÃO DE Y (VARIÁVEL RESPOSTA)
# ==============================================================================

cat("\n\nPASSO 2: DISTRIBUIÇÃO DA VARIÁVEL RESPOSTA (Y)\n")
cat("================================================================================\n\n")

Y_long <- as.data.frame(Y_mat) %>%
  tibble::rownames_to_column("Microrregiao") %>%
  tidyr::pivot_longer(cols = -Microrregiao, names_to = "Ano", values_to = "Y") %>%
  dplyr::mutate(Ano = as.numeric(Ano))

cat("Resumo de Y (contagens absolutas):\n")
print(summary(Y_long$Y))
cat("\nVariância:", var(Y_long$Y, na.rm = TRUE))
cat("\nMédia:", mean(Y_long$Y, na.rm = TRUE))
cat("\nRazão Variância/Média (dispersão):", var(Y_long$Y, na.rm = TRUE) / mean(Y_long$Y, na.rm = TRUE))
cat("\n\nInterpretação: Valor > 1 indica SUPERDISPERSÃO (importante para escolher Poisson vs Binomial Negativa)\n")

# Gráfico 1: Distribuição de Y
p1 <- ggplot(Y_long, aes(x = Y)) +
  geom_histogram(bins = 30, fill = "#E41A1C", alpha = 0.7) +
  labs(title = "Distribuição de Y (Óbitos)",
       x = "Número de Óbitos",
       y = "Frequência") +
  theme_minimal()

# Gráfico 2: Q-Q plot (normalidade)
p2 <- ggplot(Y_long, aes(sample = Y)) +
  stat_qq(color = "#E41A1C", size = 2) +
  stat_qq_line(color = "black") +
  labs(title = "Q-Q Plot: Y vs Normal",
       x = "Quantis Teóricos",
       y = "Quantis da Amostra") +
  theme_minimal()

p_dist_Y <- gridExtra::grid.arrange(p1, p2, ncol = 2)
ggsave("graficos/01_Distribuicao_Y.png", p_dist_Y, width = 12, height = 5, dpi = 300)
print(p_dist_Y)

cat("\n✅ Gráfico salvo: 01_Distribuicao_Y.png\n")

# ==============================================================================
# PASSO 3: VERIFICAR RELAÇÃO ENTRE Y E OFFSET (E)
# ==============================================================================

cat("\n\nPASSO 3: RELAÇÃO Y x OFFSET (E)\n")
cat("================================================================================\n\n")

E_long <- as.data.frame(E) %>%
  tibble::rownames_to_column("Microrregiao") %>%
  tidyr::pivot_longer(cols = -Microrregiao, names_to = "Ano", values_to = "E") %>%
  dplyr::mutate(Ano = as.numeric(Ano))

YE_long <- Y_long %>%
  dplyr::left_join(E_long, by = c("Microrregiao", "Ano")) %>%
  dplyr::mutate(Taxa = Y / E)

cat("Correlação entre Y e E (contagem e offset):\n")
corr_YE <- cor(YE_long$Y, YE_long$E, use = "complete.obs")
cat("r =", corr_YE, "\n")
cat("\nInterpretação: Correlação alta indica que regiões maiores têm mais óbitos (esperado)\n")

# Gráfico: Y vs E
p_YE <- ggplot(YE_long, aes(x = E, y = Y)) +
  geom_point(alpha = 0.5, color = "#377EB8") +
  geom_smooth(method = "lm", color = "red", se = TRUE) +
  labs(title = "Relação: Y (Óbitos) vs E (Nascidos Vivos)",
       x = "Nascidos Vivos (E)",
       y = "Óbitos (Y)") +
  theme_minimal()

ggsave("graficos/02_Y_vs_E.png", p_YE, width = 8, height = 6, dpi = 300)
print(p_YE)

cat("\n✅ Gráfico salvo: 02_Y_vs_E.png\n")

# ==============================================================================
# PASSO 4: VARIABILIDADE ESPACIAL (ENTRE MICRORREGIÕES)
# ==============================================================================

cat("\n\nPASSO 4: VARIABILIDADE ESPACIAL\n")
cat("================================================================================\n\n")

media_micro <- Y_long %>%
  dplyr::group_by(Microrregiao) %>%
  dplyr::summarise(
    Media_Y = mean(Y, na.rm = TRUE),
    SD_Y = sd(Y, na.rm = TRUE),
    CV = (SD_Y / Media_Y) * 100,  # Coeficiente de Variação
    .groups = "drop"
  ) %>%
  dplyr::arrange(dplyr::desc(Media_Y))

cat("Top 5 microrregiões com MAIOR média de óbitos:\n")
print(media_micro %>% dplyr::slice(1:5))

cat("\n\nTop 5 microrregiões com MENOR média de óbitos:\n")
print(media_micro %>% dplyr::slice((n()-4):n()))

cat("\n\nCoeficiente de Variação entre microrregiões:\n")
cat("Mínimo:", min(media_micro$CV, na.rm = TRUE), "%\n")
cat("Máximo:", max(media_micro$CV, na.rm = TRUE), "%\n")
cat("Média:", mean(media_micro$CV, na.rm = TRUE), "%\n")

# Gráfico: Heterogeneidade espacial
p_hetero <- ggplot(media_micro, aes(x = reorder(Microrregiao, Media_Y), y = Media_Y)) +
  geom_col(fill = "#4DAF4A", alpha = 0.7) +
  geom_errorbar(aes(ymin = Media_Y - SD_Y, ymax = Media_Y + SD_Y), width = 0.2) +
  coord_flip() +
  labs(title = "Heterogeneidade Espacial: Média de Y por Microrregião",
       x = "Microrregião",
       y = "Média de Óbitos") +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 6))

ggsave("graficos/03_Heterogeneidade_Espacial.png", p_hetero, width = 10, height = 12, dpi = 300)
print(p_hetero)

cat("\n✅ Gráfico salvo: 03_Heterogeneidade_Espacial.png\n")

# ==============================================================================
# PASSO 5: VARIABILIDADE TEMPORAL
# ==============================================================================

cat("\n\nPASSO 5: VARIABILIDADE TEMPORAL\n")
cat("================================================================================\n\n")

media_ano <- Y_long %>%
  dplyr::group_by(Ano) %>%
  dplyr::summarise(
    Media = mean(Y, na.rm = TRUE),
    SD = sd(Y, na.rm = TRUE),
    Min = min(Y, na.rm = TRUE),
    Max = max(Y, na.rm = TRUE),
    .groups = "drop"
  )

cat("Evolução temporal:\n")
print(media_ano)

# Tendência temporal
trend <- lm(Media ~ Ano, data = media_ano)
cat("\nTendência linear:\n")
cat("Coeficiente:", coef(trend)[2], "(mudança por ano)\n")
cat("R²:", summary(trend)$r.squared, "\n")

# Gráfico: Série temporal agregada
p_tempo <- ggplot(media_ano, aes(x = Ano, y = Media)) +
  geom_line(color = "#E41A1C", size = 1) +
  geom_point(color = "#E41A1C", size = 3) +
  geom_ribbon(aes(ymin = Media - SD, ymax = Media + SD), alpha = 0.2, fill = "#E41A1C") +
  geom_smooth(method = "lm", color = "blue", linetype = "dashed", se = FALSE) +
  labs(title = "Tendência Temporal: Óbitos Médios por Ano",
       x = "Ano",
       y = "Média de Óbitos") +
  theme_minimal()

ggsave("graficos/04_Tendencia_Temporal.png", p_tempo, width = 10, height = 6, dpi = 300)
print(p_tempo)

cat("\n✅ Gráfico salvo: 04_Tendencia_Temporal.png\n")

# ==============================================================================
# PASSO 6: CORRELAÇÃO ENTRE COVARIÁVEIS
# ==============================================================================

cat("\n\nPASSO 6: CORRELAÇÃO ENTRE COVARIÁVEIS\n")
cat("================================================================================\n\n")

# Preparar dados para correlação
X1_long <- as.data.frame(x[,,1]) %>%
  tibble::rownames_to_column("Microrregiao") %>%
  tidyr::pivot_longer(cols = -Microrregiao, names_to = "Ano", values_to = "X1") %>%
  dplyr::mutate(Ano = as.numeric(Ano))

X2_long <- as.data.frame(x[,,2]) %>%
  tibble::rownames_to_column("Microrregiao") %>%
  tidyr::pivot_longer(cols = -Microrregiao, names_to = "Ano", values_to = "X2") %>%
  dplyr::mutate(Ano = as.numeric(Ano))

X3_long <- as.data.frame(x[,,3]) %>%
  tibble::rownames_to_column("Microrregiao") %>%
  tidyr::pivot_longer(cols = -Microrregiao, names_to = "Ano", values_to = "X3") %>%
  dplyr::mutate(Ano = as.numeric(Ano))

dados_modelo <- Y_long %>%
  dplyr::left_join(X1_long, by = c("Microrregiao", "Ano")) %>%
  dplyr::left_join(X2_long, by = c("Microrregiao", "Ano")) %>%
  dplyr::left_join(X3_long, by = c("Microrregiao", "Ano")) %>%
  dplyr::left_join(E_long, by = c("Microrregiao", "Ano"))

# Matriz de correlação
mat_corr <- dados_modelo %>%
  dplyr::select(Y, X1, X2, X3, E) %>%
  cor(use = "complete.obs")

cat("Matriz de Correlação:\n")
print(mat_corr)

# Gráfico: Correlação
png("graficos/05_Correlacao_Variaveis.png", width = 800, height = 700, res = 300)
corrplot::corrplot(mat_corr, method = "circle", type = "upper", 
                   addCoef.col = "black", diag = TRUE, 
                   title = "Correlação entre Variáveis do Modelo")
dev.off()

cat("\n✅ Gráfico salvo: 05_Correlacao_Variaveis.png\n")

cat("\nInterpretações importantes:\n")
cat("- Y vs X1, X2, X3: Direção e magnitude da associação\n")
cat("- X1 vs X2 vs X3: Multicolinearidade (problema se > 0.8)\n")

# ==============================================================================
# PASSO 7: ASSOCIAÇÃO Y vs COVARIÁVEIS (SCATTER PLOTS)
# ==============================================================================

cat("\n\nPASSO 7: RELAÇÃO Y vs COVARIÁVEIS\n")
cat("================================================================================\n\n")

p_scatter <- GGally::ggpairs(
  dados_modelo %>% dplyr::select(Y, X1, X2, X3),
  upper = list(continuous = GGally::wrap("points", alpha = 0.3)),
  lower = list(continuous = GGally::wrap("smooth", alpha = 0.3))
) +
  theme_minimal()

ggsave("graficos/06_Scatter_Y_vs_Covariaveis.png", p_scatter, width = 10, height = 10, dpi = 300)
print(p_scatter)

cat("\n✅ Gráfico salvo: 06_Scatter_Y_vs_Covariaveis.png\n")

# ==============================================================================
# PASSO 8: VARIAÇÃO RELATIVA (TAXA vs COVARIÁVEIS)
# ==============================================================================

cat("\n\nPASSO 8: TAXA DE ÓBITOS vs COVARIÁVEIS\n")
cat("================================================================================\n\n")

dados_taxa <- dados_modelo %>%
  dplyr::mutate(Taxa = Y / E)

# Dividir em quartis para visualização
dados_taxa_quartis <- dados_taxa %>%
  dplyr::mutate(
    X1_Quartil = cut(X1, breaks = quantile(X1, na.rm = TRUE), labels = c("Q1", "Q2", "Q3", "Q4"), include.lowest = TRUE),
    X2_Quartil = cut(X2, breaks = quantile(X2, na.rm = TRUE), labels = c("Q1", "Q2", "Q3", "Q4"), include.lowest = TRUE),
    X3_Quartil = cut(X3, breaks = quantile(X3, na.rm = TRUE), labels = c("Q1", "Q2", "Q3", "Q4"), include.lowest = TRUE)
  )

resumo_quartis <- dados_taxa_quartis %>%
  dplyr::summarise(
    Taxa_por_X1 = dados_taxa_quartis %>%
      dplyr::group_by(X1_Quartil) %>%
      dplyr::summarise(Media_Taxa = mean(Taxa, na.rm = TRUE), .groups = "drop"),
    Taxa_por_X2 = dados_taxa_quartis %>%
      dplyr::group_by(X2_Quartil) %>%
      dplyr::summarise(Media_Taxa = mean(Taxa, na.rm = TRUE), .groups = "drop"),
    Taxa_por_X3 = dados_taxa_quartis %>%
      dplyr::group_by(X3_Quartil) %>%
      dplyr::summarise(Media_Taxa = mean(Taxa, na.rm = TRUE), .groups = "drop")
  )

cat("Taxa média de óbitos por Quartil de X1 (Prenatal):\n")
print(
  dados_taxa_quartis %>%
    dplyr::group_by(X1_Quartil) %>%
    dplyr::summarise(Media_Taxa = mean(Taxa, na.rm = TRUE), N = n(), .groups = "drop")
)

cat("\n\nTaxa média de óbitos por Quartil de X2 (Instrução):\n")
print(
  dados_taxa_quartis %>%
    dplyr::group_by(X2_Quartil) %>%
    dplyr::summarise(Media_Taxa = mean(Taxa, na.rm = TRUE), N = n(), .groups = "drop")
)

cat("\n\nTaxa média de óbitos por Quartil de X3 (Baixo Peso):\n")
print(
  dados_taxa_quartis %>%
    dplyr::group_by(X3_Quartil) %>%
    dplyr::summarise(Media_Taxa = mean(Taxa, na.rm = TRUE), N = n(), .groups = "drop")
)

# ==============================================================================
# PASSO 9: CHECAGEM DE PRESSUPOSTOS
# ==============================================================================

cat("\n\nPASSO 9: VERIFICAÇÃO DE PRESSUPOSTOS PARA MODELAGEM\n")
cat("================================================================================\n\n")

cat("1. INDEPENDÊNCIA:\n")
cat("   - Dados são por microrregião e ano (observações independentes? Verificar autocorrelação espacial)\n")
cat("   - Teste de Ljung-Box para autocorrelação temporal:\n")
serie_temporal <- media_ano$Media
lb_test <- Box.test(serie_temporal, type = "Ljung-Box")
cat("   - p-value:", lb_test$p.value, "\n")
if (lb_test$p.value < 0.05) {
  cat("   ⚠️ Autocorrelação significativa! Considere modelos com estrutura temporal (ARIMA, efeitos aleatórios)\n")
} else {
  cat("   ✓ Sem autocorrelação significativa\n")
}

cat("\n2. DISTRIBUIÇÃO DE Y:\n")
cat("   - Contagens inteiras (Poisson/Binomial Negativa apropriado)\n")
cat("   - Dispersão:", var(Y_long$Y, na.rm = TRUE) / mean(Y_long$Y, na.rm = TRUE), "\n")
if (var(Y_long$Y, na.rm = TRUE) / mean(Y_long$Y, na.rm = TRUE) > 1.5) {
  cat("   ⚠️ Superdispersão detectada! Use Binomial Negativa em vez de Poisson\n")
}

cat("\n3. MULTICOLINEARIDADE:\n")
cat("   - Correlação entre X1, X2, X3:\n")
print(cor(dados_modelo %>% dplyr::select(X1, X2, X3), use = "complete.obs"))

cat("\n4. VALORES EXTREMOS:\n")
valores_extremos <- Y_long %>%
  dplyr::mutate(Z_score = scale(Y)[,1]) %>%
  dplyr::filter(abs(Z_score) > 3)
cat("   - Observações com Z-score > 3:", nrow(valores_extremos), "\n")

# ==============================================================================
# PASSO 10: SALVAR RELATÓRIO
# ==============================================================================

cat("\n\nPASSO 10: GERANDO RELATÓRIO\n")
cat("================================================================================\n\n")

# Salvar estatísticas em arquivo
sink("relatorios/Analise_Descritiva_Completa.txt")

cat("RELATÓRIO DE ANÁLISE DESCRITIVA PARA MODELAGEM\n")
cat("="*80, "\n")
cat("Data:", format(Sys.time(), "%d/%m/%Y %H:%M:%S"), "\n\n")

cat("1. DADOS\n")
cat("-"*80, "\n")
cat("Y (Óbitos): ", nrow(Y_mat), " regiões x ", ncol(Y_mat), " anos\n")
cat("E (Offset): ", nrow(E), " regiões x ", ncol(E), " anos\n")
cat("Covariáveis: ", dim(x)[3], "\n\n")

cat("2. RESUMO DE Y\n")
cat("-"*80, "\n")
print(summary(Y_long$Y))
cat("\nVariância:", var(Y_long$Y, na.rm = TRUE), "\n")
cat("Dispersão (Var/Média):", var(Y_long$Y, na.rm = TRUE) / mean(Y_long$Y, na.rm = TRUE), "\n\n")

cat("3. CORRELAÇÕES\n")
cat("-"*80, "\n")
print(mat_corr)
cat("\n")

cat("4. RECOMENDAÇÕES PARA MODELAGEM\n")
cat("-"*80, "\n")
cat("- Modelo: GLM Poisson ou Binomial Negativa com offset log(E)\n")
cat("- Efeitos aleatórios: Por microrregião (heterogeneidade espacial detectada)\n")
cat("- Efeito tempo: Verificar tendência linear\n")
cat("- Estrutura: Y_i,t | E_i,t ~ Poisson(E_i,t * λ_i,t)\n")
cat("         log(λ_i,t) = β0 + β1*X1 + β2*X2 + β3*X3 + u_i + tempo_t\n\n")

sink()

cat("\n✅ Relatório salvo: relatorios/Analise_Descritiva_Completa.txt\n\n")

cat("================================================================================\n")
cat("RESUMO FINAL\n")
cat("================================================================================\n\n")
cat("✅ 10 gráficos gerados em graficos/\n")
cat("✅ Relatório detalhado em relatorios/\n")
cat("✅ Próximo passo: Modelagem com estrutura adequada\n")
