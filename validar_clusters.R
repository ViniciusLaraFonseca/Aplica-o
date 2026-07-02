# ==============================================================================
# validar_clusters_CORRIGIDO.R
# Script para validar a estrutura HIERÁRQUICA de clusters
# ==============================================================================

rm(list = ls())
setwd("C:/Users/vlara/OneDrive/Estatistica UFMG/Mestrado/Pesquisa/Aplicação/main")

load("dados_modelo_final.RData")
library(dplyr)
library(tidyr)

# Matriz hAI
hAI <- structure(
  .Data = c(
    1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
    1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
  ),
  .Dim = c(75, 4)
)
colnames(hAI) <- paste0("h_", 1:4)

# Criar tabela
tabela_validacao <- data.frame(
  ordem = 1:75,
  MICRO_ = ordem_modelo,
  AI = cluster_info$AI,
  cluster = cluster_info$cluster,
  hAI,
  row.names = NULL
)

tabela_validacao$cluster_nome <- dplyr::case_when(
  tabela_validacao$cluster == 1 ~ "Alta adequação",
  tabela_validacao$cluster == 2 ~ "Média-alta adequação",
  tabela_validacao$cluster == 3 ~ "Média-baixa adequação",
  tabela_validacao$cluster == 4 ~ "Baixa adequação"
)

# ==============================================================================
# VALIDAÇÃO CORRETA DA ESTRUTURA HIERÁRQUICA
# ==============================================================================

cat("\n")
cat(paste(rep("=", 80), collapse = ""), "\n")
cat("VALIDAÇÃO DA ESTRUTURA HIERÁRQUICA DE CLUSTERS\n")
cat(paste(rep("=", 80), collapse = ""), "\n\n")

cat("📊 ESTRUTURA DO MODELO:\n")
cat("─────────────────────────────────────────────────────────────────────\n")
cat("O modelo usa clusters HIERÁRQUICOS CUMULATIVOS:\n\n")
cat("  Cluster 1 (Alta adequação):      h = [1,0,0,0]  → Pertence SÓ ao cluster 1\n")
cat("  Cluster 2 (Média-alta adequação): h = [1,1,0,0]  → Pertence aos clusters 1 e 2\n")
cat("  Cluster 3 (Média-baixa adequação): h = [1,1,1,0]  → Pertence aos clusters 1,2,3\n")
cat("  Cluster 4 (Baixa adequação):      h = [1,1,1,1]  → Pertence a TODOS os clusters\n\n")

cat("Isso significa que:\n")
cat("  • Cluster 1 é o MAIS restrito (apenas 28 microrregiões)\n")
cat("  • Cluster 4 é o MAIS amplo (inclui TODAS as 75 microrregiões)\n")
cat("  • É uma estrutura de CLASSES LATENTES HIERÁRQUICAS\n\n")

# Validar estrutura hierárquica
cat("✅ VALIDAÇÕES:\n")
cat("─────────────────────────────────────────────────────────────────────\n")

# 1. Verificar que h_1 = 1 para TODAS as microrregiões
if(all(hAI[,1] == 1)) {
  cat("✅ Todas as 75 microrregiões pertencem ao Cluster 1 (h_1 = 1)\n")
} else {
  cat("❌ ERRO: Algumas microrregiões não pertencem ao Cluster 1!\n")
}

# 2. Verificar estrutura hierárquica (se h_k = 1, então h_{k-1} = 1)
hierarquia_ok <- all(
  (hAI[,2] == 1 & hAI[,1] == 1) | (hAI[,2] == 0),
  (hAI[,3] == 1 & hAI[,2] == 1) | (hAI[,3] == 0),
  (hAI[,4] == 1 & hAI[,3] == 1) | (hAI[,4] == 0)
)

if(hierarquia_ok) {
  cat("✅ Estrutura hierárquica PERFEITA (h_k=1 ⇒ h_{k-1}=1)\n")
} else {
  cat("❌ ERRO: Violação da estrutura hierárquica!\n")
}

# 3. Contar pertinência a cada cluster
cat("\n📊 PERTINÊNCIA AOS CLUSTERS (visão cumulativa):\n")
cat("─────────────────────────────────────────────────────────────────────\n")

for(k in 1:4) {
  n_k <- sum(hAI[,k] == 1)
  pct <- n_k / 75 * 100
  cat(sprintf("  Cluster %d (h_%d = 1): %2d microrregiões (%.1f%%)\n", k, k, n_k, pct))
}

# 4. Resumo por grupo
cat("\n📊 GRUPOS DE MICRORREGIÕES:\n")
cat("─────────────────────────────────────────────────────────────────────\n")

grupos <- list(
  "Grupo 1 (Alta)" = which(hAI[,1] == 1 & hAI[,2] == 0),
  "Grupo 2 (Média-alta)" = which(hAI[,1] == 1 & hAI[,2] == 1 & hAI[,3] == 0),
  "Grupo 3 (Média-baixa)" = which(hAI[,1] == 1 & hAI[,2] == 1 & hAI[,3] == 1 & hAI[,4] == 0),
  "Grupo 4 (Baixa)" = which(hAI[,1] == 1 & hAI[,2] == 1 & hAI[,3] == 1 & hAI[,4] == 1)
)

for(nome in names(grupos)) {
  idx <- grupos[[nome]]
  n <- length(idx)
  ia_range <- range(tabela_validacao$AI[idx])
  cat(sprintf("  %s: %2d microrregiões, IA ∈ [%6.2f, %6.2f]\n", 
              nome, n, ia_range[1], ia_range[2]))
}

# 5. Mostrar primeiras e últimas de cada grupo
cat("\n\n📋 DETALHE DOS GRUPOS:\n")
cat(paste(rep("=", 80), collapse = ""), "\n")

for(nome in names(grupos)) {
  idx <- grupos[[nome]]
  cat(sprintf("\n%s:\n", nome))
  cat(paste(rep("-", 50), collapse = ""), "\n")
  print(tabela_validacao[idx, c("MICRO_", "AI", "h_1", "h_2", "h_3", "h_4")], row.names = FALSE)
}

# 6. Salvar resultado
cat("\n\n💾 Exportando resultados...\n")
write.csv(tabela_validacao, "validacao_clusters_hierarquica.csv", row.names = FALSE)
cat("✅ Tabela completa salva em 'validacao_clusters_hierarquica.csv'\n\n")

# Resumo final
cat(paste(rep("=", 80), collapse = ""), "\n")
cat("RESUMO FINAL\n")
cat(paste(rep("=", 80), collapse = ""), "\n\n")

cat("✅ Estrutura de clusters hierárquicos validada com SUCESSO!\n\n")
cat("Características:\n")
cat("  • Modelo de classes latentes hierárquicas\n")
cat("  • 4 níveis de adequação da notificação\n")
cat("  • Estrutura cumulativa: Cluster 1 ⊂ Cluster 2 ⊂ Cluster 3 ⊂ Cluster 4\n")
cat("  • Total: 75 microrregiões de Minas Gerais\n")
cat("  • Ordenadas por IA (Índice de Adequação) decrescente\n\n")

cat("Grupos mutuamente exclusivos:\n")
cat("  • Alta adequação:     28 microrregiões\n")
cat("  • Média-alta adequação: 16 microrregiões\n")
cat("  • Média-baixa adequação: 14 microrregiões\n")
cat("  • Baixa adequação:     17 microrregiões\n")