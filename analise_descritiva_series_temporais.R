# ==============================================================================
# Análise Descritiva e Gráficos de Série Temporal
# Y (Óbitos) e Covariáveis X (Prenatal, Instrução, Baixo Peso)
# ==============================================================================

rm(list = ls())

# Carregar dados
source("Covariaveis.R", encoding = "UTF-8")

# Pacotes necessários
pkgs <- c("ggplot2", "dplyr", "tidyr", "gridExtra", "cowplot")
for (p in pkgs) {
  if (!require(p, character.only = TRUE)) install.packages(p)
  library(p, character.only = TRUE)
}

# Criar pasta para gráficos
if (!dir.exists("graficos")) dir.create("graficos")

# ==============================================================================
# 1. PREPARAR DADOS PARA VISUALIZACAO
# ==============================================================================

# Converter matrizes para data frames long
Y_long <- as.data.frame(Y_mat) %>%
  tibble::rownames_to_column("Microrregiao") %>%
  tidyr::pivot_longer(
    cols = -Microrregiao,
    names_to = "Ano",
    values_to = "Y"
  ) %>%
  dplyr::mutate(Ano = as.numeric(Ano))

# Covariável X1 (Prenatal)
X1_long <- as.data.frame(x[,,1]) %>%
  tibble::rownames_to_column("Microrregiao") %>%
  tidyr::pivot_longer(
    cols = -Microrregiao,
    names_to = "Ano",
    values_to = "X1_Prenatal"
  ) %>%
  dplyr::mutate(Ano = as.numeric(Ano))

# Covariável X2 (Instrução)
X2_long <- as.data.frame(x[,,2]) %>%
  tibble::rownames_to_column("Microrregiao") %>%
  tidyr::pivot_longer(
    cols = -Microrregiao,
    names_to = "Ano",
    values_to = "X2_Instrucao"
  ) %>%
  dplyr::mutate(Ano = as.numeric(Ano))

# Covariável X3 (Baixo Peso)
X3_long <- as.data.frame(x[,,3]) %>%
  tibble::rownames_to_column("Microrregiao") %>%
  tidyr::pivot_longer(
    cols = -Microrregiao,
    names_to = "Ano",
    values_to = "X3_BaixoPeso"
  ) %>%
  dplyr::mutate(Ano = as.numeric(Ano))

# Combinar todos os dados
dados_completos <- Y_long %>%
  dplyr::left_join(X1_long, by = c("Microrregiao", "Ano")) %>%
  dplyr::left_join(X2_long, by = c("Microrregiao", "Ano")) %>%
  dplyr::left_join(X3_long, by = c("Microrregiao", "Ano"))

# ==============================================================================
# 2. SELECIONAR PRINCIPAIS MICRORREGIOES
# ==============================================================================

# Seleção baseada na média de óbitos (maiores e menores)
microrregioes_top <- Y_long %>%
  dplyr::group_by(Microrregiao) %>%
  dplyr::summarise(Media_Y = mean(Y, na.rm = TRUE), .groups = "drop") %>%
  dplyr::arrange(dplyr::desc(Media_Y)) %>%
  dplyr::slice(c(1:3, (n()-2):n())) %>%
  dplyr::pull(Microrregiao)

cat("Microrregioes selecionadas (Top 3 e Bottom 3 por média de óbitos):\n")
print(microrregioes_top)
cat("\n")

# ==============================================================================
# 3. ANÁLISE DESCRITIVA GLOBAL
# ==============================================================================

cat("================================================================================\n")
cat("ANALISE DESCRITIVA - SERIE TEMPORAL (2000-2022)\n")
cat("================================================================================\n\n")

cat("--- VARIÁVEL Y (OBITOS NEONATAIS PRECOCES) ---\n")
cat("Resumo geral:\n")
print(summary(Y_long$Y))
cat("\nPor período:\n")
Y_periodo <- Y_long %>%
  dplyr::mutate(
    Periodo = cut(Ano,
                  breaks = c(1999, 2005, 2010, 2015, 2020, 2023),
                  labels = c("2000-2005", "2006-2010", "2011-2015", "2016-2020", "2021-2022"))
  )
print(Y_periodo %>%
  dplyr::group_by(Periodo) %>%
  dplyr::summarise(
    Media = mean(Y, na.rm = TRUE),
    Mediana = median(Y, na.rm = TRUE),
    Min = min(Y, na.rm = TRUE),
    Max = max(Y, na.rm = TRUE),
    SD = sd(Y, na.rm = TRUE),
    .groups = "drop"
  ))

cat("\n\n--- COVARIAVEL X1 (PRENATAL - CONSULTASNO PRENATAL) ---\n")
cat("Resumo geral:\n")
print(summary(X1_long$X1_Prenatal))

cat("\n\n--- COVARIAVEL X2 (INSTRUCAO - MAES COM INSTRUCAO COMPLETA) ---\n")
cat("Resumo geral:\n")
print(summary(X2_long$X2_Instrucao))

cat("\n\n--- COVARIAVEL X3 (BAIXO PESO - NASCIDOS COM PESO < 2500g) ---\n")
cat("Resumo geral:\n")
print(summary(X3_long$X3_BaixoPeso))

# ==============================================================================
# 4. GRÁFICO 1: SÉRIE TEMPORAL AGREGADA (TOTAL MG)
# ==============================================================================

Y_agg <- Y_long %>%
  dplyr::group_by(Ano) %>%
  dplyr::summarise(Y_Total = sum(Y, na.rm = TRUE), .groups = "drop")

X1_agg <- X1_long %>%
  dplyr::group_by(Ano) %>%
  dplyr::summarise(X1_Media = mean(X1_Prenatal, na.rm = TRUE), .groups = "drop")

X2_agg <- X2_long %>%
  dplyr::group_by(Ano) %>%
  dplyr::summarise(X2_Media = mean(X2_Instrucao, na.rm = TRUE), .groups = "drop")

X3_agg <- X3_long %>%
  dplyr::group_by(Ano) %>%
  dplyr::summarise(X3_Media = mean(X3_BaixoPeso, na.rm = TRUE), .groups = "drop")

dados_agg <- Y_agg %>%
  dplyr::left_join(X1_agg, by = "Ano") %>%
  dplyr::left_join(X2_agg, by = "Ano") %>%
  dplyr::left_join(X3_agg, by = "Ano")

# Gráfico Y
p_Y_agg <- ggplot(dados_agg, aes(x = Ano, y = Y_Total)) +
  geom_line(color = "#E41A1C", size = 1) +
  geom_point(color = "#E41A1C", size = 3) +
  geom_smooth(method = "loess", se = TRUE, alpha = 0.2, color = "#E41A1C", linetype = "dashed") +
  labs(
    title = "Serie Temporal: Obitos Neonatais Precoces em MG",
    x = "Ano",
    y = "Total de Obitos"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", size = 12),
    panel.grid.major = element_line(color = "gray90"),
    panel.grid.minor = element_blank()
  )

# Gráfico X1
p_X1_agg <- ggplot(dados_agg, aes(x = Ano, y = X1_Media)) +
  geom_line(color = "#377EB8", size = 1) +
  geom_point(color = "#377EB8", size = 3) +
  geom_smooth(method = "loess", se = TRUE, alpha = 0.2, color = "#377EB8", linetype = "dashed") +
  labs(
    title = "X1: Prenatal (Media MG)",
    x = "Ano",
    y = "Indice Prenatal"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", size = 12),
    panel.grid.major = element_line(color = "gray90"),
    panel.grid.minor = element_blank()
  )

# Gráfico X2
p_X2_agg <- ggplot(dados_agg, aes(x = Ano, y = X2_Media)) +
  geom_line(color = "#4DAF4A", size = 1) +
  geom_point(color = "#4DAF4A", size = 3) +
  geom_smooth(method = "loess", se = TRUE, alpha = 0.2, color = "#4DAF4A", linetype = "dashed") +
  labs(
    title = "X2: Instrucao (Media MG)",
    x = "Ano",
    y = "Indice Instrucao"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", size = 12),
    panel.grid.major = element_line(color = "gray90"),
    panel.grid.minor = element_blank()
  )

# Gráfico X3
p_X3_agg <- ggplot(dados_agg, aes(x = Ano, y = X3_Media)) +
  geom_line(color = "#FF7F00", size = 1) +
  geom_point(color = "#FF7F00", size = 3) +
  geom_smooth(method = "loess", se = TRUE, alpha = 0.2, color = "#FF7F00", linetype = "dashed") +
  labs(
    title = "X3: Baixo Peso (Media MG)",
    x = "Ano",
    y = "Indice Baixo Peso"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", size = 12),
    panel.grid.major = element_line(color = "gray90"),
    panel.grid.minor = element_blank()
  )

# Combinar gráficos agregados
p_agg_combined <- cowplot::plot_grid(p_Y_agg, p_X1_agg, p_X2_agg, p_X3_agg,
                                      nrow = 2, ncol = 2)

ggsave("graficos/01_Series_Temporais_Agregadas_MG.png", p_agg_combined, 
       width = 14, height = 10, dpi = 300)
print(p_agg_combined)

cat("\n✅ Gráfico agregado salvo: 01_Series_Temporais_Agregadas_MG.png\n")

# ==============================================================================
# 5. GRÁFICOS: PRINCIPAIS MICRORREGIOES
# ==============================================================================

# Filtrar dados das microrregioes selecionadas
dados_top <- dados_completos %>%
  dplyr::filter(Microrregiao %in% microrregioes_top)

# Gráfico Y por Microrregião
p_Y_top <- ggplot(dados_top, aes(x = Ano, y = Y, color = Microrregiao)) +
  geom_line(size = 1) +
  geom_point(size = 2.5) +
  facet_wrap(~Microrregiao, scales = "free_y", ncol = 2) +
  labs(
    title = "Série Temporal: Y (Óbitos) - Principais Microrregiões",
    x = "Ano",
    y = "Número de Óbitos"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", size = 13),
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8)
  )

ggsave("graficos/02_Y_Principais_Microrregioes.png", p_Y_top,
       width = 12, height = 10, dpi = 300)
print(p_Y_top)

cat("✅ Gráfico Y por microrregião salvo: 02_Y_Principais_Microrregioes.png\n")

# Gráfico X1 por Microrregião
p_X1_top <- ggplot(dados_top, aes(x = Ano, y = X1_Prenatal, color = Microrregiao)) +
  geom_line(size = 1) +
  geom_point(size = 2.5) +
  facet_wrap(~Microrregiao, scales = "free_y", ncol = 2) +
  labs(
    title = "Série Temporal: X1 (Prenatal) - Principais Microrregiões",
    x = "Ano",
    y = "Índice Prenatal"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", size = 13),
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8)
  )

ggsave("graficos/03_X1_Principais_Microrregioes.png", p_X1_top,
       width = 12, height = 10, dpi = 300)
print(p_X1_top)

cat("✅ Gráfico X1 por microrregião salvo: 03_X1_Principais_Microrregioes.png\n")

# Gráfico X2 por Microrregião
p_X2_top <- ggplot(dados_top, aes(x = Ano, y = X2_Instrucao, color = Microrregiao)) +
  geom_line(size = 1) +
  geom_point(size = 2.5) +
  facet_wrap(~Microrregiao, scales = "free_y", ncol = 2) +
  labs(
    title = "Série Temporal: X2 (Instrução) - Principais Microrregiões",
    x = "Ano",
    y = "Índice Instrução"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", size = 13),
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8)
  )

ggsave("graficos/04_X2_Principais_Microrregioes.png", p_X2_top,
       width = 12, height = 10, dpi = 300)
print(p_X2_top)

cat("✅ Gráfico X2 por microrregião salvo: 04_X2_Principais_Microrregioes.png\n")

# Gráfico X3 por Microrregião
p_X3_top <- ggplot(dados_top, aes(x = Ano, y = X3_BaixoPeso, color = Microrregiao)) +
  geom_line(size = 1) +
  geom_point(size = 2.5) +
  facet_wrap(~Microrregiao, scales = "free_y", ncol = 2) +
  labs(
    title = "Série Temporal: X3 (Baixo Peso) - Principais Microrregiões",
    x = "Ano",
    y = "Índice Baixo Peso"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", size = 13),
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8)
  )

ggsave("graficos/05_X3_Principais_Microrregioes.png", p_X3_top,
       width = 12, height = 10, dpi = 300)
print(p_X3_top)

cat("✅ Gráfico X3 por microrregião salvo: 05_X3_Principais_Microrregioes.png\n")

# ==============================================================================
# 6. HEATMAP: Y POR MICRORREGIAO E ANO
# ==============================================================================

p_heatmap_Y <- ggplot(Y_long, aes(x = Ano, y = reorder(Microrregiao, Y, FUN = mean), fill = Y)) +
  geom_tile() +
  scale_fill_gradient(low = "#FFFFCC", high = "#E41A1C") +
  labs(
    title = "Mapa de Calor: Óbitos (Y) por Microrregião e Ano",
    x = "Ano",
    y = "Microrregião",
    fill = "Óbitos"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", size = 13),
    axis.text.y = element_text(size = 7),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8)
  )

ggsave("graficos/06_Heatmap_Y.png", p_heatmap_Y,
       width = 14, height = 12, dpi = 300)
print(p_heatmap_Y)

cat("✅ Heatmap Y salvo: 06_Heatmap_Y.png\n")

# ==============================================================================
# 7. CORRELACAO: Y vs COVARIAVEIS (GLOBAL)
# ==============================================================================

# Calcular correlações por ano
corr_data <- dados_completos %>%
  dplyr::group_by(Ano) %>%
  dplyr::summarise(
    Corr_Y_X1 = cor(Y, X1_Prenatal, use = "complete.obs"),
    Corr_Y_X2 = cor(Y, X2_Instrucao, use = "complete.obs"),
    Corr_Y_X3 = cor(Y, X3_BaixoPeso, use = "complete.obs"),
    .groups = "drop"
  )

p_corr <- ggplot(corr_data %>% tidyr::pivot_longer(-Ano, names_to = "Correlacao", values_to = "Valor"),
                 aes(x = Ano, y = Valor, color = Correlacao)) +
  geom_line(size = 1) +
  geom_point(size = 3) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray", size = 0.8) +
  scale_color_manual(values = c("Corr_Y_X1" = "#377EB8",
                                 "Corr_Y_X2" = "#4DAF4A",
                                 "Corr_Y_X3" = "#FF7F00"),
                     labels = c("Corr_Y_X1" = "Y vs X1 (Prenatal)",
                               "Corr_Y_X2" = "Y vs X2 (Instrução)",
                               "Corr_Y_X3" = "Y vs X3 (Baixo Peso)")) +
  labs(
    title = "Evolução da Correlação: Y vs Covariáveis",
    x = "Ano",
    y = "Correlação de Pearson",
    color = "Variáveis"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", size = 13),
    legend.position = "bottom"
  )

ggsave("graficos/07_Correlacoes_Temporais.png", p_corr,
       width = 12, height = 7, dpi = 300)
print(p_corr)

cat("✅ Gráfico de correlações salvo: 07_Correlacoes_Temporais.png\n")

# ==============================================================================
# 8. RESUMO DESCRITIVO POR MICRORREGIAO
# ==============================================================================

cat("\n================================================================================\n")
cat("RESUMO POR MICRORREGIAO (Média 2000-2022)\n")
cat("================================================================================\n\n")

resumo_micro <- dados_completos %>%
  dplyr::group_by(Microrregiao) %>%
  dplyr::summarise(
    Y_Media = mean(Y, na.rm = TRUE),
    Y_SD = sd(Y, na.rm = TRUE),
    X1_Media = mean(X1_Prenatal, na.rm = TRUE),
    X2_Media = mean(X2_Instrucao, na.rm = TRUE),
    X3_Media = mean(X3_BaixoPeso, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::arrange(dplyr::desc(Y_Media))

cat("Top 10 Microrregioes (maior média de óbitos):\n")
print(resumo_micro %>% dplyr::slice(1:10))

cat("\n\nBottom 10 Microrregioes (menor média de óbitos):\n")
print(resumo_micro %>% dplyr::slice((n()-9):n()))

# ==============================================================================
# 9. SALVAR DADOS DESCRITIVOS
# ==============================================================================

write.csv(resumo_micro, "graficos/Resumo_Descritivo_Microrregioes.csv", row.names = FALSE)
write.csv(corr_data, "graficos/Correlacoes_Temporais.csv", row.names = FALSE)

cat("\n\n✅ RESUMO FINAL\n")
cat("================================================================================\n")
cat("Gráficos gerados:\n")
cat("1. 01_Series_Temporais_Agregadas_MG.png - Series temporais globais\n")
cat("2. 02_Y_Principais_Microrregioes.png - Óbitos por microrregião\n")
cat("3. 03_X1_Principais_Microrregioes.png - Prenatal por microrregião\n")
cat("4. 04_X2_Principais_Microrregioes.png - Instrução por microrregião\n")
cat("5. 05_X3_Principais_Microrregioes.png - Baixo peso por microrregião\n")
cat("6. 06_Heatmap_Y.png - Mapa de calor dos óbitos\n")
cat("7. 07_Correlacoes_Temporais.png - Correlações Y vs Covariáveis\n")
cat("================================================================================\n")
