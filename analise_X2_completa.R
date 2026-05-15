# ==============================================================================
# Painel Completo: Série Temporal de X2 (Instrução) - Todas as Microrregiões
# ==============================================================================

rm(list = ls())

# Carregar dados
source("Covariaveis.R", encoding = "UTF-8")

# Pacotes necessários
pkgs <- c("ggplot2", "dplyr", "tidyr", "gridExtra")
for (p in pkgs) {
  if (!require(p, character.only = TRUE)) install.packages(p)
  library(p, character.only = TRUE)
}

# Criar pasta para gráficos
if (!dir.exists("graficos")) dir.create("graficos")

# ==============================================================================
# 1. PREPARAR DADOS DE X2
# ==============================================================================

# Covariável X2 (Instrução)
X2_long <- as.data.frame(x[,,2]) %>%
  tibble::rownames_to_column("Microrregiao") %>%
  tidyr::pivot_longer(
    cols = -Microrregiao,
    names_to = "Ano",
    values_to = "X2_Instrucao"
  ) %>%
  dplyr::mutate(Ano = as.numeric(Ano))

# ==============================================================================
# 2. GRÁFICO PAINEL: TODAS AS 75 MICRORREGIÕES
# ==============================================================================

p_X2_all <- ggplot(X2_long, aes(x = Ano, y = X2_Instrucao, color = Microrregiao)) +
  geom_line(linewidth = 0.7, alpha = 0.8) +
  geom_point(size = 1.5) +
  facet_wrap(~Microrregiao, scales = "free_y", ncol = 8) +
  labs(
    title = "Série Temporal: X2 (Instrução) - Mães sem Ensino Superior",
    subtitle = "Proporção de nascidos vivos cuja mãe não possui ensino superior (2000-2022)",
    x = "Ano",
    y = "Proporção [0-1]",
    caption = "Fonte: DATASUS"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(size = 11, color = "gray50"),
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
    axis.text.y = element_text(size = 7),
    strip.text = element_text(size = 8, face = "bold"),
    panel.grid.major = element_line(color = "gray90", linewidth = 0.3),
    panel.grid.minor = element_blank()
  )

ggsave("graficos/X2_Series_Temporais_Todas_Microrregioes.png", p_X2_all,
       width = 20, height = 24, dpi = 300)
print(p_X2_all)

cat("✅ Painel completo X2 salvo: X2_Series_Temporais_Todas_Microrregioes.png\n")

# ==============================================================================
# 3. ANÁLISE DESCRITIVA DE X2
# ==============================================================================

cat("\n================================================================================\n")
cat("ANALISE DESCRITIVA DE X2 (INSTRUCAO - MAES SEM ENSINO SUPERIOR)\n")
cat("================================================================================\n\n")

cat("Resumo Geral (2000-2022):\n")
print(summary(X2_long$X2_Instrucao))

cat("\n\nPor Período:\n")
X2_periodo <- X2_long %>%
  dplyr::mutate(
    Periodo = cut(Ano,
                  breaks = c(1999, 2005, 2010, 2015, 2020, 2023),
                  labels = c("2000-2005", "2006-2010", "2011-2015", "2016-2020", "2021-2022"))
  )
print(X2_periodo %>%
  dplyr::group_by(Periodo) %>%
  dplyr::summarise(
    Media = mean(X2_Instrucao, na.rm = TRUE),
    Mediana = median(X2_Instrucao, na.rm = TRUE),
    Min = min(X2_Instrucao, na.rm = TRUE),
    Max = max(X2_Instrucao, na.rm = TRUE),
    SD = sd(X2_Instrucao, na.rm = TRUE),
    .groups = "drop"
  ))

# ==============================================================================
# 4. RANKING: TOP 10 E BOTTOM 10 MICRORREGIÕES
# ==============================================================================

cat("\n\nRanking de Microrregiões por Média de X2 (2000-2022):\n")

ranking_X2 <- X2_long %>%
  dplyr::group_by(Microrregiao) %>%
  dplyr::summarise(
    Media_X2 = mean(X2_Instrucao, na.rm = TRUE),
    Min_X2 = min(X2_Instrucao, na.rm = TRUE),
    Max_X2 = max(X2_Instrucao, na.rm = TRUE),
    SD_X2 = sd(X2_Instrucao, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::arrange(dplyr::desc(Media_X2))

cat("\nTop 10 (MAIOR proporção de mães sem ensino superior):\n")
print(ranking_X2 %>% dplyr::slice(1:10))

cat("\n\nBottom 10 (MENOR proporção de mães sem ensino superior):\n")
print(ranking_X2 %>% dplyr::slice((n()-9):n()))

# ==============================================================================
# 5. EVOLUÇÃO TEMPORAL: MUDANÇAS DE 2000 PARA 2022
# ==============================================================================

cat("\n\n================================================================================\n")
cat("EVOLUCAO TEMPORAL: MUDANCAS DE 2000 PARA 2022\n")
cat("================================================================================\n\n")

X2_evolucao <- X2_long %>%
  tidyr::pivot_wider(names_from = Ano, values_from = X2_Instrucao) %>%
  dplyr::mutate(
    Mudanca_Absoluta = `2022` - `2000`,
    Mudanca_Relativa = (((`2022` - `2000`) / `2000`) * 100)
  ) %>%
  dplyr::select(Microrregiao, `2000`, `2022`, Mudanca_Absoluta, Mudanca_Relativa) %>%
  dplyr::arrange(Mudanca_Absoluta)

cat("Maiores REDUÇÕES (Melhorias):\n")
print(X2_evolucao %>% dplyr::slice(1:10))

cat("\n\nMaiores AUMENTOS (Pioras):\n")
print(X2_evolucao %>% dplyr::slice((n()-9):n()))

# ==============================================================================
# 6. SALVAR RESULTADOS
# ==============================================================================

write.csv(ranking_X2, "graficos/X2_Ranking_Microrregioes.csv", row.names = FALSE)
write.csv(X2_evolucao, "graficos/X2_Evolucao_2000_2022.csv", row.names = FALSE)

cat("\n\n✅ RESUMO FINAL\n")
cat("================================================================================\n")
cat("Arquivos gerados:\n")
cat("1. X2_Series_Temporais_Todas_Microrregioes.png - Painel com 75 microrregiões\n")
cat("2. X2_Ranking_Microrregioes.csv - Ranking por média de X2\n")
cat("3. X2_Evolucao_2000_2022.csv - Mudanças de 2000 para 2022\n")
cat("================================================================================\n")
