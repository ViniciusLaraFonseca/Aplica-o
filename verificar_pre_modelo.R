# Verificar tudo antes de rodar o modelo
load("dados_modelo_final.RData")

cat("Verificações pré-modelo:\n")
cat("1. Dimensões:\n")
cat("  Y:", dim(Y_mat), "\n")
cat("  E:", dim(E), "\n")
cat("  x:", dim(x), "\n")
cat("  ordem_modelo:", length(ordem_modelo), "\n\n")

cat("2. Estrutura hierárquica:\n")
cat("  Clusters únicos:", unique(cluster_info$cluster), "\n")
cat("  Distribuição:", table(cluster_info$cluster), "\n\n")

cat("3. Ordens consistentes:\n")
cat("  rownames(Y) == ordem_modelo:", all(rownames(Y_mat) == ordem_modelo), "\n")
cat("  rownames(E) == ordem_modelo:", all(rownames(E) == ordem_modelo), "\n")
cat("  dimnames(x)[[1]] == ordem_modelo:", all(dimnames(x)[[1]] == ordem_modelo), "\n\n")

cat("✅ Pronto para o modelo!\n")