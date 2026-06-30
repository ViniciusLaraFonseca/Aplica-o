# ==============================================================================
# ajuste_TODOS_clusters_gamma_COMPLETO_vs_SIMPLIFICADO_FINAL.R
# Modelo com TODOS os clusters simultâneos (75 microrregiões)
# Estrutura hierárquica cumulativa via gamma[k,t]
# COM TODOS OS GRÁFICOS CORRIGIDOS
# ==============================================================================

inicio_global <- Sys.time()

setwd("C:/Users/vlara/OneDrive/Estatistica UFMG/Mestrado/Pesquisa/Aplicação/main")

pkgs <- c("nimble", "coda", "dplyr", "ggplot2", "tidyr", 
          "readr", "stringr", "tibble", "scales")
for (pkg in pkgs) {
  if (!require(pkg, character.only = TRUE)) install.packages(pkg)
  library(pkg, character.only = TRUE)
}

# ==============================================================================
# 1. CARREGAR DADOS
# ==============================================================================

load("dados_modelo_final.RData")

E_norm <- E / mean(E)
n_times <- ncol(Y_mat)
p_full <- dim(x)[3]
K <- 4
N_regions <- 75

h_mat <- structure(
  .Data = c(
    1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
    1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
  ),
  .Dim = c(75, 4)
)

region_names <- rownames(Y_mat)
anos_label <- colnames(Y_mat)

# Grupos
grupo <- rep(NA, N_regions)
for(i in 1:N_regions) {
  if(h_mat[i,1] == 1 && h_mat[i,2] == 0) grupo[i] <- 1
  else if(h_mat[i,1] == 1 && h_mat[i,2] == 1 && h_mat[i,3] == 0) grupo[i] <- 2
  else if(h_mat[i,1] == 1 && h_mat[i,2] == 1 && h_mat[i,3] == 1 && h_mat[i,4] == 0) grupo[i] <- 3
  else if(all(h_mat[i,] == 1)) grupo[i] <- 4
}

x_full <- x
x_simpl <- x[, , c(1, 3), drop = FALSE]

cat(sprintf("\n✅ Modelo com TODAS as %d microrregiões e K=%d clusters\n", N_regions, K))
cat(sprintf("   Cluster 1: %d | Cluster 2: %d | Cluster 3: %d | Cluster 4: %d\n",
            sum(grupo==1), sum(grupo==2), sum(grupo==3), sum(grupo==4)))

# ==============================================================================
# 2. DADOS DE ADJACÊNCIA
# ==============================================================================

data_adj <- list(
  num = c(6,5,5,4,5,4,3,7,3,7,5,6,3,4,4,4,8,6,4,8,6,4,7,3,6,9,5,7,2,3,2,4,5,2,7,8,3,5,9,5,8,5,4,2,6,2,9,3,9,6,6,7,5,6,6,5,2,7,5,3,9,5,5,4,4,5,6,4,4,6,5,6,3,9,2),
  adj = c(
    2,3,5,12,17,18, 1,4,5,6,17, 1,5,7,12,19, 2,6,17,31,
    1,2,3,6,7, 2,4,5,7, 3,5,6, 9,12,15,19,21,30,58,
    8,15,19, 11,18,27,28,39,58,74, 10,17,18,20,28,
    1,3,8,18,19,58, 14,17,20, 13,20,23,42, 8,9,21,36,
    22,26,45,68, 1,2,4,11,13,18,20,31, 1,10,11,12,17,58,
    3,8,9,12, 11,13,14,17,23,28,67,73, 8,15,27,30,36,58,
    16,26,42,45, 14,20,32,35,42,67,73, 57,64,71,
    36,41,49,51,63,66, 16,22,33,45,47,62,65,68,70,
    10,21,36,39,58, 10,11,20,39,55,67,74, 41,48, 8,21,58,
    4,17, 23,35,55,67, 26,44,65,68,70, 40,54,
    23,32,42,45,47,55,56, 15,21,25,27,37,39,51,66,
    36,51,71, 40,43,52,54,72, 10,27,28,36,53,60,61,66,74,
    34,38,46,52,54, 25,29,48,49,51,54,64,72,
    14,22,23,35,45, 38,52,59,72, 33,68,
    16,22,26,35,42,47, 40,52,
    26,35,45,49,50,56,61,62,74, 29,41,54,
    25,41,47,50,59,61,63,69,72, 47,49,52,59,62,70,
    25,36,37,41,64,71, 38,40,43,46,50,59,70,
    39,61,63,66,69, 34,38,40,41,48,72,
    28,32,35,56,67,74, 35,47,55,61,74, 24,71,
    8,10,12,18,21,27,30, 43,49,50,52,72, 39,61,74,
    39,47,49,53,56,60,69,74,75, 26,47,50,65,70,
    25,49,53,66,69, 24,41,51,71, 26,33,62,70,
    25,36,39,53,63, 20,23,28,32,55,73, 16,26,33,44,
    49,53,61,63, 26,33,50,52,62,65, 24,37,51,57,64,
    38,41,43,49,54,59, 20,23,67,
    10,28,39,47,55,56,61,60,75, 61,74
  )
)

n_adj <- length(data_adj$adj)

# ==============================================================================
# 3. SELEÇÃO DE 3 REGIÕES POR CLUSTER PARA GRÁFICOS
# ==============================================================================

set.seed(123)

regioes_por_cluster <- list()
for(g in 1:4) {
  idx_grupo <- which(grupo == g)
  
  if(g == 1) {
    bh <- grep("BELO_HORIZONTE", region_names)[1]
    outras <- sample(setdiff(idx_grupo, bh), min(2, length(idx_grupo)-1))
    regioes_por_cluster[[g]] <- c(bh, outras)
  } else if(g == 4) {
    aracuai <- grep("ARACUAI", region_names)[1]
    outras <- sample(setdiff(idx_grupo, aracuai), min(2, length(idx_grupo)-1))
    regioes_por_cluster[[g]] <- c(aracuai, outras)
  } else {
    regioes_por_cluster[[g]] <- sample(idx_grupo, min(3, length(idx_grupo)))
  }
}

# Lista plana: 12 regiões (3 por cluster)
all_regions <- unlist(regioes_por_cluster)

cat("\nRegiões selecionadas para gráficos (3 por cluster):\n")
for(g in 1:4) {
  cat(sprintf("  Cluster %d: %s\n", g, 
              paste(region_names[regioes_por_cluster[[g]]], collapse = ", ")))
}

# ==============================================================================
# 4. FUNÇÕES AUXILIARES
# ==============================================================================

make_label <- function(names) {
  sapply(names, function(n) {
    n_clean <- gsub("_", " ", n)
    if (nchar(n_clean) > 30) paste0(substr(n_clean, 1, 27), "...") else n_clean
  })
}

safe_hpd <- function(sv) {
  if (var(sv) < 1e-12) return(c(NA_real_, NA_real_))
  as.numeric(HPDinterval(as.mcmc(sv), prob = 0.95))
}

safe_gelman <- function(obj) {
  tryCatch(gelman.diag(obj, autoburnin = FALSE)$psrf[, 1],
           error = function(e) rep(NA_real_, nvar(obj)))
}

cores_cluster <- c("Cluster 1" = "#1f77b4", "Cluster 2" = "#ff7f0e", 
                   "Cluster 3" = "#2ca02c", "Cluster 4" = "#d62728")

# ==============================================================================
# 5. FUNÇÃO DE GRÁFICOS (CORRIGIDA - 3 colunas, facetas separadas)
# ==============================================================================

# ==============================================================================
# FUNÇÃO DE GRÁFICOS CORRIGIDA - Nomes reais nas colunas
# ==============================================================================

# ==============================================================================
# FUNÇÃO DE GRÁFICOS CORRIGIDA - USANDO facet_wrap (como no código que funciona)
# ==============================================================================

criar_graficos_corrigido <- function(model_type, scenario_dir, x_data, p, cov_names) {
  
  cat(sprintf("\n📊 Criando gráficos para modelo %s...\n", model_type))
  
  samples <- readRDS(file.path(scenario_dir, "samples.rds"))
  samples_mat <- as.matrix(samples)
  
  beta_cols <- grep("^beta\\[", colnames(samples_mat), value = TRUE)
  
  # =======================================================================
  # PAINEL 1: Lambda[i,t] - USANDO facet_wrap (como no código que funciona)
  # =======================================================================
  
  cat("  [1/8] Painel lambda (facet_wrap)...\n")
  
  # Criar data.frame com todas as regiões selecionadas
  lambda_panel <- do.call(rbind, lapply(all_regions, function(i) {
    do.call(rbind, lapply(1:n_times, function(t) {
      nm <- paste0("lambda[", i, ", ", t, "]")
      sv <- samples_mat[, nm]
      hpd <- safe_hpd(sv)
      data.frame(
        Region = i,
        Nome = make_label(region_names[i]),
        Cluster = factor(paste("Cluster", grupo[i]), levels = paste("Cluster", 1:4)),
        Time = t,
        Mean = mean(sv),
        Lower = hpd[1],
        Upper = hpd[2],
        stringsAsFactors = FALSE
      )
    }))
  }))
  
  # Ordenar por cluster para organização visual
  lambda_panel$Region <- factor(lambda_panel$Region, 
                                levels = all_regions[order(grupo[all_regions])])
  
  ggsave(
    file.path(scenario_dir, "painel_lambda_por_cluster.png"),
    ggplot(lambda_panel, aes(x = Time)) +
      geom_ribbon(aes(ymin = Lower, ymax = Upper), fill = "darkorange", alpha = 0.3) +
      geom_line(aes(y = Mean), color = "darkorange", linewidth = 0.8) +
      facet_wrap(~ Region, ncol = 3, scales = "free_y",
                 labeller = labeller(Region = function(x) {
                   paste(make_label(region_names[as.numeric(x)]), 
                         "- Cluster", grupo[as.numeric(x)])
                 })) +
      scale_x_continuous(breaks = seq(1, n_times, by = 5), 
                         labels = anos_label[seq(1, n_times, by = 5)]) +
      theme_bw(base_size = 11) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
            strip.text = element_text(size = 9)) +
      labs(title = sprintf("λ[i,t] - %s", model_type),
           subtitle = "3 regiões selecionadas por cluster",
           x = "Ano", y = expression(lambda[it])),
    width = 14, height = 10
  )
  
  # =======================================================================
  # PAINEL 2: Mu[i,t] - USANDO facet_wrap
  # =======================================================================
  
  cat("  [2/8] Painel mu (facet_wrap)...\n")
  
  mu_panel <- do.call(rbind, lapply(all_regions, function(i) {
    do.call(rbind, lapply(1:n_times, function(t) {
      ldraws <- samples_mat[, paste0("lambda[", i, ", ", t, "]")]
      eps_t  <- samples_mat[, paste0("epsilon[", i, ", ", t, "]")]
      bdraws <- samples_mat[, beta_cols, drop = FALSE]
      
      # Cálculo vetorizado
      lin_pred <- as.vector(bdraws %*% as.numeric(x_data[i, t, ]))
      mu <- ldraws * E_norm[i, t] * eps_t * exp(lin_pred)
      
      hpd <- safe_hpd(mu)
      data.frame(
        Region = i,
        Nome = make_label(region_names[i]),
        Cluster = factor(paste("Cluster", grupo[i]), levels = paste("Cluster", 1:4)),
        Time = t,
        Mean = mean(mu),
        Lower = hpd[1],
        Upper = hpd[2],
        Y_obs = Y_mat[i, t],
        stringsAsFactors = FALSE
      )
    }))
  }))
  
  mu_panel$Region <- factor(mu_panel$Region, 
                            levels = all_regions[order(grupo[all_regions])])
  
  ggsave(
    file.path(scenario_dir, "painel_mu_por_cluster.png"),
    ggplot(mu_panel, aes(x = Time)) +
      geom_ribbon(aes(ymin = Lower, ymax = Upper), fill = "steelblue", alpha = 0.3) +
      geom_line(aes(y = Mean), color = "steelblue", linewidth = 0.8) +
      geom_point(aes(y = Y_obs), shape = 21, fill = "white", color = "black", size = 1.2) +
      facet_wrap(~ Region, ncol = 3, scales = "free_y",
                 labeller = labeller(Region = function(x) {
                   paste(make_label(region_names[as.numeric(x)]), 
                         "- Cluster", grupo[as.numeric(x)])
                 })) +
      scale_x_continuous(breaks = seq(1, n_times, by = 5), 
                         labels = anos_label[seq(1, n_times, by = 5)]) +
      theme_bw(base_size = 11) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
            strip.text = element_text(size = 9)) +
      labs(title = sprintf("μ[i,t] - %s", model_type),
           subtitle = "Pontos = Y observado | 3 regiões selecionadas por cluster",
           x = "Ano", y = expression(mu[it])),
    width = 14, height = 10
  )
  
  # =======================================================================
  # PAINEL 3: Theta[i,t] - USANDO facet_wrap
  # =======================================================================
  
  cat("  [3/8] Painel theta (facet_wrap)...\n")
  
  theta_panel <- do.call(rbind, lapply(all_regions, function(i) {
    do.call(rbind, lapply(1:n_times, function(t) {
      ldraws <- samples_mat[, paste0("lambda[", i, ", ", t, "]")]
      bdraws <- samples_mat[, beta_cols, drop = FALSE]
      
      # Cálculo vetorizado
      lin_pred <- as.vector(bdraws %*% as.numeric(x_data[i, t, ]))
      theta <- ldraws * exp(lin_pred)
      
      hpd <- safe_hpd(theta)
      data.frame(
        Region = i,
        Nome = make_label(region_names[i]),
        Cluster = factor(paste("Cluster", grupo[i]), levels = paste("Cluster", 1:4)),
        Time = t,
        Mean = mean(theta),
        Lower = hpd[1],
        Upper = hpd[2],
        stringsAsFactors = FALSE
      )
    }))
  }))
  
  theta_panel$Region <- factor(theta_panel$Region, 
                               levels = all_regions[order(grupo[all_regions])])
  
  ggsave(
    file.path(scenario_dir, "painel_theta_por_cluster.png"),
    ggplot(theta_panel, aes(x = Time)) +
      geom_ribbon(aes(ymin = Lower, ymax = Upper), fill = "darkgreen", alpha = 0.3) +
      geom_line(aes(y = Mean), color = "darkgreen", linewidth = 0.8) +
      facet_wrap(~ Region, ncol = 3, scales = "free_y",
                 labeller = labeller(Region = function(x) {
                   paste(make_label(region_names[as.numeric(x)]), 
                         "- Cluster", grupo[as.numeric(x)])
                 })) +
      scale_x_continuous(breaks = seq(1, n_times, by = 5), 
                         labels = anos_label[seq(1, n_times, by = 5)]) +
      theme_bw(base_size = 11) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
            strip.text = element_text(size = 9)) +
      labs(title = sprintf("θ[i,t] (risco relativo) - %s", model_type),
           subtitle = "3 regiões selecionadas por cluster",
           x = "Ano", y = expression(theta[it])),
    width = 14, height = 10
  )
  
  # =======================================================================
  # PAINEL 4: Gamma[k,t] - facet_wrap
  # =======================================================================
  
  cat("  [4/8] Painel gamma...\n")
  
  gamma_panel <- do.call(rbind, lapply(1:K, function(k) {
    do.call(rbind, lapply(1:n_times, function(t) {
      nm <- paste0("gamma[", k, ", ", t, "]")
      sv <- samples_mat[, nm]
      hpd <- safe_hpd(sv)
      data.frame(
        Cluster = factor(paste("Cluster", k), levels = paste("Cluster", 1:4)),
        Time = t,
        Mean = mean(sv),
        Lower = hpd[1],
        Upper = hpd[2],
        stringsAsFactors = FALSE
      )
    }))
  }))
  
  ggsave(
    file.path(scenario_dir, "painel_gamma_temporal.png"),
    ggplot(gamma_panel, aes(x = Time)) +
      geom_ribbon(aes(ymin = Lower, ymax = Upper), fill = "steelblue", alpha = 0.3) +
      geom_line(aes(y = Mean), color = "steelblue", linewidth = 0.8) +
      facet_wrap(~ Cluster, ncol = 2, scales = "free_y") +
      scale_x_continuous(breaks = seq(1, n_times, by = 5), 
                         labels = anos_label[seq(1, n_times, by = 5)]) +
      theme_bw(base_size = 11) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8)) +
      labs(title = sprintf("γ[k,t] - %s", model_type),
           subtitle = expression(gamma[1][t] ~ "~ U(0, 0.05) | Demais: U(0, 1 - soma acumulada)"),
           x = "Ano", y = expression(gamma[k][t])),
    width = 12, height = 8
  )
  
  # =======================================================================
  # PAINEL 5: Epsilon médio por cluster
  # =======================================================================
  
  cat("  [5/8] Painel epsilon...\n")
  
  epsilon_cluster <- do.call(rbind, lapply(1:K, function(g) {
    idx_grupo <- which(grupo == g)
    do.call(rbind, lapply(1:n_times, function(t) {
      eps_draws <- sapply(idx_grupo, function(i) {
        samples_mat[, paste0("epsilon[", i, ", ", t, "]")]
      })
      eps_mean_draws <- rowMeans(eps_draws)
      hpd <- safe_hpd(eps_mean_draws)
      data.frame(
        Cluster = factor(paste("Cluster", g), levels = paste("Cluster", 1:4)),
        Time = t,
        Mean = mean(eps_mean_draws),
        Lower = hpd[1],
        Upper = hpd[2],
        stringsAsFactors = FALSE
      )
    }))
  }))
  
  ggsave(
    file.path(scenario_dir, "painel_epsilon_por_cluster.png"),
    ggplot(epsilon_cluster, aes(x = Time)) +
      geom_ribbon(aes(ymin = Lower, ymax = Upper), fill = "darkorange", alpha = 0.3) +
      geom_line(aes(y = Mean), color = "darkorange", linewidth = 0.8) +
      facet_wrap(~ Cluster, ncol = 2, scales = "free_y") +
      scale_x_continuous(breaks = seq(1, n_times, by = 5), 
                         labels = anos_label[seq(1, n_times, by = 5)]) +
      theme_bw(base_size = 11) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8)) +
      labs(title = sprintf("Epsilon médio por Cluster - %s", model_type),
           subtitle = expression(epsilon ~ "= 1 -" ~ Sigma ~ "h[i,k]" ~ gamma[k][t]),
           x = "Ano", y = expression(epsilon[t])),
    width = 12, height = 8
  )
  
  # =======================================================================
  # PAINEL 6: Beta
  # =======================================================================
  
  cat("  [6/8] Painel beta...\n")
  
  beta_names <- paste0("beta[", seq_len(p), "]")
  beta_summary <- do.call(rbind, lapply(beta_names, function(nm) {
    sv <- samples_mat[, nm]
    hpd <- safe_hpd(sv)
    data.frame(
      Parameter = nm,
      Covariate = cov_names[which(beta_names == nm)],
      Mean = mean(sv),
      Lower = hpd[1],
      Upper = hpd[2],
      stringsAsFactors = FALSE
    )
  }))
  
  ggsave(
    file.path(scenario_dir, "beta_posterior.png"),
    ggplot(beta_summary, aes(x = Covariate, y = Mean)) +
      geom_point(size = 4, color = "steelblue") +
      geom_errorbar(aes(ymin = Lower, ymax = Upper), 
                    width = 0.2, color = "steelblue", linewidth = 1) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
      theme_bw(base_size = 13) +
      labs(title = sprintf("Coeficientes beta - %s", model_type), y = expression(beta)),
    width = 8, height = 5
  )
  
  # =======================================================================
  # PAINEL 7: Efeito espacial s[i]
  # =======================================================================
  
  cat("  [7/8] Painel efeito espacial...\n")
  
  s_names <- paste0("s[", 1:N_regions, "]")
  s_summary <- do.call(rbind, lapply(1:N_regions, function(i) {
    sv <- samples_mat[, s_names[i]]
    hpd <- safe_hpd(sv)
    data.frame(
      Region = i,
      Nome = region_names[i],
      Cluster = factor(paste("Cluster", grupo[i]), levels = paste("Cluster", 1:4)),
      Mean = mean(sv),
      Lower = hpd[1],
      Upper = hpd[2],
      stringsAsFactors = FALSE
    )
  }))
  s_summary <- s_summary[order(s_summary$Mean), ]
  s_summary$Nome <- factor(s_summary$Nome, levels = s_summary$Nome)
  
  ggsave(
    file.path(scenario_dir, "s_posterior_por_cluster.png"),
    ggplot(s_summary, aes(x = Nome, y = Mean, color = Cluster)) +
      geom_point(size = 1.5) +
      geom_errorbar(aes(ymin = Lower, ymax = Upper), width = 0.4, linewidth = 0.3) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
      scale_color_manual(values = cores_cluster) +
      theme_bw() +
      theme(
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 5),
        legend.position = "bottom"
      ) +
      labs(title = sprintf("Efeito espacial s[i] - %s", model_type),
           y = "s[i]", x = ""),
    width = 18, height = 8
  )
  
  # =======================================================================
  # PAINEL 8: TABELA DE DIAGNÓSTICO (ESS e Rhat)
  # =======================================================================
  
  cat("  [8/8] Tabela de diagnóstico ESS e Rhat...\n")
  
  # Criar tabela de diagnóstico
  mcmc_list_full <- mcmc.list(lapply(1:2, function(ch) as.mcmc(samples[[ch]])))
  
  # Selecionar parâmetros principais para diagnóstico
  params_diag <- c(
    beta_names,
    paste0("gamma[", rep(1:K, each = 3), ", ", rep(c(1, round(n_times/2), n_times), times = K), "]"),
    "sigma_s"
  )
  params_diag <- params_diag[params_diag %in% colnames(samples_mat)]
  
  diag_table <- do.call(rbind, lapply(params_diag, function(nm) {
    sv <- samples_mat[, nm]
    data.frame(
      Parameter = nm,
      ESS = as.numeric(effectiveSize(mcmc_list_full[, nm])),
      Rhat = safe_gelman(mcmc_list_full[, nm]),
      stringsAsFactors = FALSE
    )
  }))
  
  # Adicionar informações de lambda para regiões selecionadas
  for (i in all_regions) {
    for (t in c(1, round(n_times/2), n_times)) {
      nm <- paste0("lambda[", i, ", ", t, "]")
      if (nm %in% colnames(samples_mat)) {
        diag_table <- rbind(diag_table, data.frame(
          Parameter = paste0(nm, " (", region_names[i], ")"),
          ESS = as.numeric(effectiveSize(mcmc_list_full[, nm])),
          Rhat = safe_gelman(mcmc_list_full[, nm]),
          stringsAsFactors = FALSE
        ))
      }
    }
  }
  
  # Salvar tabela
  write_csv(diag_table, file.path(scenario_dir, "diagnostico_ess_rhat.csv"))
  
  # Criar gráfico da tabela
  diag_table_plot <- diag_table %>%
    mutate(
      Parameter = gsub("\\[", " [", Parameter),
      Parameter = gsub(", ", ", ", Parameter),
      ESS_cat = case_when(
        ESS < 100 ~ "Ruim (<100)",
        ESS < 500 ~ "Razoável (100-500)",
        ESS < 1000 ~ "Bom (500-1000)",
        TRUE ~ "Excelente (>1000)"
      ),
      Rhat_cat = case_when(
        Rhat > 1.1 ~ "Preocupante (>1.1)",
        Rhat > 1.05 ~ "Aceitável (1.05-1.1)",
        TRUE ~ "Bom (<1.05)"
      )
    )
  
  # Tabela como gráfico
  ggsave(
    file.path(scenario_dir, "diagnostico_ess_rhat.png"),
    ggplot(diag_table_plot, aes(x = Parameter, y = ESS, fill = ESS_cat)) +
      geom_col() +
      geom_hline(yintercept = c(100, 500, 1000), linetype = "dashed", color = "red", alpha = 0.5) +
      scale_fill_manual(values = c("Ruim (<100)" = "red",
                                   "Razoável (100-500)" = "orange",
                                   "Bom (500-1000)" = "yellow",
                                   "Excelente (>1000)" = "green")) +
      theme_bw(base_size = 9) +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7),
            legend.position = "bottom") +
      labs(title = sprintf("ESS (Effective Sample Size) - %s", model_type),
           subtitle = "Linhas vermelhas: 100, 500, 1000",
           x = "", y = "ESS"),
    width = 14, height = 8
  )
  
  # Tabela Rhat
  ggsave(
    file.path(scenario_dir, "diagnostico_rhat.png"),
    ggplot(diag_table_plot, aes(x = Parameter, y = Rhat, color = Rhat_cat)) +
      geom_point(size = 3) +
      geom_hline(yintercept = c(1.05, 1.1), linetype = "dashed", color = "red", alpha = 0.5) +
      scale_color_manual(values = c("Preocupante (>1.1)" = "red",
                                    "Aceitável (1.05-1.1)" = "orange",
                                    "Bom (<1.05)" = "green")) +
      theme_bw(base_size = 9) +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7),
            legend.position = "bottom") +
      labs(title = sprintf("Rhat (Gelman-Rubin) - %s", model_type),
           subtitle = "Linhas vermelhas: 1.05 e 1.1",
           x = "", y = expression(hat(R))),
    width = 14, height = 8
  )
  
  # Também salvar como tabela CSV legível
  write_csv(diag_table, file.path(scenario_dir, "diagnostico_completo.csv"))
  
  cat(sprintf("✅ Gráficos corrigidos para %s concluídos!\n", model_type))
}

# ==============================================================================
# 6. FUNÇÃO PRINCIPAL DO MODELO
# ==============================================================================

run_model <- function(model_type, output_dir_base) {
  
  cat(sprintf("\n%s=== MODELO %s ===%s\n", 
              paste(rep("=", 40), collapse = ""), model_type, 
              paste(rep("=", 40), collapse = "")))
  
  if (model_type == "COMPLETO") {
    x_data <- x_full; p <- 3
    cov_names <- c("prenatal", "instrucao", "baixo_peso")
  } else {
    x_data <- x_simpl; p <- 2
    cov_names <- c("prenatal", "baixo_peso")
  }
  
  a0 <- 1.0; b0 <- 1.0; w <- 0.85
  a_unif <- 0.0; b_unif <- 0.05
  
  constants <- list(
    n_regions = N_regions, n_times = n_times, p = p, K = K,
    a0 = a0, b0 = b0, w = w,
    a_unif = a_unif, b_unif = b_unif,
    adj = data_adj$adj, num = data_adj$num, weights = rep(1, n_adj), n_adj = n_adj
  )
  
  data_nimble <- list(Y = Y_mat, E = E_norm, x = x_data, h = h_mat)
  
  code <- nimbleCode({
    for (j in 1:p) { beta[j] ~ dnorm(0, sd = 1) }
    
    for (t in 1:n_times) {
      gamma[1, t] ~ dunif(a_unif, b_unif)
      for (k in 2:K) {
        gamma[k, t] ~ dunif(0, 1 - sum(gamma[1:(k-1), t]))
      }
    }
    
    for (i in 1:n_regions) {
      for (t in 1:n_times) {
        epsilon[i, t] <- 1 - inprod(h[i, 1:K], gamma[1:K, t])
      }
    }
    
    sigma_s ~ T(dt(0, 1, 1), 0, )
    tau_s <- 1 / (sigma_s^2)
    s[1:n_regions] ~ dcar_normal(adj[1:n_adj], weights[1:n_adj], num[1:n_regions], tau_s, zero_mean = 1)
    
    for (i in 1:n_regions) {
      for (t in 1:n_times) { lambda[i, t] ~ dgamma(a0, b0) }
    }
    
    for (i in 1:n_regions) {
      for (t in 1:n_times) {
        log_mu[i, t] <- log(lambda[i, t]) + log(E[i, t]) + log(epsilon[i, t]) +
          inprod(beta[1:p], x[i, t, 1:p]) + s[i]
        mu[i, t] <- exp(log_mu[i, t])
        Y[i, t] ~ dpois(mu[i, t])
        logLik_Y[i, t] <- dpois(Y[i, t], mu[i, t], log = TRUE)
      }
    }
  })
  
  # FFBS
  ffbs_lambda <- nimbleFunction(
    contains = sampler_BASE,
    setup = function(model, mvSaved, target, control) {
      n_times <- control$n_times; p <- control$p
      a0 <- control$a0; b0 <- control$b0; w <- control$w
      region_i <- control$region_i
      at_buf <- nimNumeric(n_times + 1, 0); bt_buf <- nimNumeric(n_times + 1, 0)
      calcNodes <- model$getDependencies(target, self = FALSE)
      targetNodes <- model$expandNodeNames(target)
      setupOutputs(at_buf, bt_buf)
    },
    run = function() {
      declare(t, integer()); declare(t_idx, integer())
      declare(t_back, integer()); declare(k, integer())
      declare(prod_val, double()); declare(g_it, double()); declare(nu, double())
      
      at_buf[1] <<- a0; bt_buf[1] <<- b0
      
      for (t in 1:n_times) {
        prod_val <- 0.0
        for (k in 1:p) prod_val <- prod_val + model$x[region_i, t, k] * model$beta[k]
        g_it <- model$E[region_i, t] * model$epsilon[region_i, t] * 
          exp(prod_val + model$s[region_i])
        at_buf[t + 1] <<- w * at_buf[t] + model$Y[region_i, t]
        bt_buf[t + 1] <<- w * bt_buf[t] + g_it
      }
      
      model$lambda[region_i, n_times] <<- rgamma(1, shape = at_buf[n_times + 1], 
                                                 rate = max(bt_buf[n_times + 1], 1e-10))
      for (t_idx in 1:(n_times - 1)) {
        t_back <- n_times - t_idx
        nu <- rgamma(1, shape = (1 - w) * at_buf[t_back + 1],
                     rate = max(bt_buf[t_back + 1], 1e-10))
        model$lambda[region_i, t_back] <<- nu + w * model$lambda[region_i, t_back + 1]
      }
      
      model$calculate(calcNodes)
      copy(from = model, to = mvSaved, row = 1, nodes = targetNodes, logProb = TRUE)
    },
    methods = list(reset = function() {})
  )
  
  # Inicializações
  set.seed(123)
  lambda_init <- matrix(1.0, nrow = N_regions, ncol = n_times)
  gamma_init <- matrix(0, nrow = K, ncol = n_times)
  gamma_init[1, ] <- runif(n_times, 0.01, 0.03)
  for (k in 2:K) gamma_init[k, ] <- runif(n_times, 0, 0.01)
  
  inits_1 <- list(lambda = lambda_init, beta = rep(0, p), gamma = gamma_init,
                  sigma_s = 0.5, s = rep(0, N_regions))
  inits_2 <- list(lambda = matrix(rgamma(N_regions * n_times, a0, b0), N_regions, n_times),
                  beta = rnorm(p, 0, 0.3), gamma = gamma_init * 0.5,
                  sigma_s = 1.0, s = rep(0, N_regions))
  
  # Compilar
  cat(sprintf("[%s] Compilando...\n", model_type))
  model <- nimbleModel(code = code, constants = constants, data = data_nimble, 
                       inits = inits_1, check = FALSE)
  Cmodel <- compileNimble(model)
  
  conf <- configureMCMC(model)
  conf$removeSamplers("lambda")
  for (i in seq_len(N_regions)) {
    conf$addSampler(target = paste0("lambda[", i, ", 1:", n_times, "]"),
                    type = ffbs_lambda,
                    control = list(n_times = n_times, p = p, a0 = a0, b0 = b0, w = w, region_i = i))
  }
  conf$removeSamplers("gamma")
  for (k in 1:K) for (t in 1:n_times) {
    conf$addSampler(target = paste0("gamma[", k, ", ", t, "]"), type = "slice")
  }
  
  conf$addMonitors(c("beta", "gamma", "epsilon", "lambda", "logLik_Y", "s", "sigma_s"))
  
  Rmcmc <- buildMCMC(conf)
  Cmcmc <- compileNimble(Rmcmc, project = model)
  
  # Executar
  niter <- 50000; nburnin <- 10000; nchains <- 2; thin <- 10
  
  cat(sprintf("[%s] MCMC: %d iter, %d burnin, %d chains, thin=%d\n",
              model_type, niter, nburnin, nchains, thin))
  
  t_inicio <- Sys.time()
  samples <- runMCMC(Cmcmc, niter = niter, nburnin = nburnin, nchains = nchains, thin = thin,
                     inits = list(inits_1, inits_2), samplesAsCodaMCMC = TRUE,
                     summary = FALSE, WAIC = FALSE)
  t_fim <- Sys.time()
  
  cat(sprintf("[%s] Tempo: %.1f min\n", model_type, difftime(t_fim, t_inicio, units = "mins")))
  
  scenario_dir <- file.path(output_dir_base, model_type)
  dir.create(scenario_dir, recursive = TRUE, showWarnings = FALSE)
  saveRDS(samples, file.path(scenario_dir, "samples.rds"))
  
  samples_mat <- as.matrix(samples)
  mcmc_list_full <- mcmc.list(lapply(seq_len(nchains), function(ch) as.mcmc(samples[[ch]])))
  
  # Beta
  beta_names <- paste0("beta[", seq_len(p), "]")
  beta_summary <- do.call(rbind, lapply(beta_names, function(nm) {
    sv <- samples_mat[, nm]; hpd <- safe_hpd(sv)
    data.frame(Model = model_type, Parameter = nm, Covariate = cov_names[which(beta_names == nm)],
               Mean = mean(sv), SD = sd(sv), HPD_Lower = hpd[1], HPD_Upper = hpd[2],
               ESS = as.numeric(effectiveSize(mcmc_list_full[, nm])),
               Rhat = safe_gelman(mcmc_list_full[, nm]), stringsAsFactors = FALSE)
  }))
  write_csv(beta_summary, file.path(scenario_dir, "beta_summary.csv"))
  
  # Gamma
  gamma_summary <- do.call(rbind, lapply(1:K, function(k) {
    do.call(rbind, lapply(1:n_times, function(t) {
      nm <- paste0("gamma[", k, ", ", t, "]")
      sv <- samples_mat[, nm]
      data.frame(Cluster = k, Time = t, Mean = mean(sv), SD = sd(sv),
                 stringsAsFactors = FALSE)
    }))
  }))
  write_csv(gamma_summary, file.path(scenario_dir, "gamma_summary.csv"))
  
  # Epsilon
  epsilon_summary <- do.call(rbind, lapply(all_regions, function(i) {
    do.call(rbind, lapply(1:n_times, function(t) {
      nm <- paste0("epsilon[", i, ", ", t, "]")
      sv <- samples_mat[, nm]; hpd <- safe_hpd(sv)
      data.frame(Region = i, Nome = region_names[i], Cluster = grupo[i],
                 Time = t, Mean = mean(sv), Lower = hpd[1], Upper = hpd[2],
                 stringsAsFactors = FALSE)
    }))
  }))
  write_csv(epsilon_summary, file.path(scenario_dir, "epsilon_summary.csv"))
  
  # WAIC
  loglik_names <- grep("logLik_Y", colnames(samples_mat), value = TRUE)
  lm_mat <- samples_mat[, loglik_names, drop = FALSE]
  lppd <- sum(apply(lm_mat, 2, function(x) { mx <- max(x); mx + log(mean(exp(x - mx))) }))
  p_waic <- sum(apply(lm_mat, 2, var))
  waic <- -2 * (lppd - p_waic)
  LPML <- sum(log(1 / apply(lm_mat, 2, function(x) mean(exp(-x)))))
  
  criteria <- data.frame(Model = model_type, WAIC = waic, LPML = LPML, 
                         lppd = lppd, pWAIC = p_waic, stringsAsFactors = FALSE)
  write_csv(criteria, file.path(scenario_dir, "criteria.csv"))
  
  cat(sprintf("[%s] WAIC = %.2f | LPML = %.2f\n", model_type, waic, LPML))
  
  # Gráficos
  criar_graficos(model_type, scenario_dir, x_data, p, cov_names)
  
  data.frame(Model = model_type, p = p, WAIC = waic, LPML = LPML,
             ESS_beta_min = min(beta_summary$ESS, na.rm = TRUE),
             Rhat_beta_max = max(beta_summary$Rhat, na.rm = TRUE),
             tempo_min = as.numeric(difftime(t_fim, t_inicio, units = "mins")),
             stringsAsFactors = FALSE)
}

# ==============================================================================
# 7. EXECUÇÃO
# ==============================================================================

output_dir <- "C:/Users/vlara/OneDrive/Estatistica UFMG/Mestrado/Pesquisa/Aplicação/resultados_TODOS_clusters_gamma"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

resultados <- list()
resultados[["COMPLETO"]] <- run_model("COMPLETO", output_dir)
resultados[["SIMPLIFICADO"]] <- run_model("SIMPLIFICADO", output_dir)

# ==============================================================================
# 8. CONSOLIDAÇÃO
# ==============================================================================

resumo <- bind_rows(resultados)
write_csv(resumo, file.path(output_dir, "resumo_comparativo.csv"))

cat("\n")
cat(paste(rep("=", 80), collapse = ""), "\n")
cat("RESUMO COMPARATIVO - TODOS OS CLUSTERS\n")
cat(paste(rep("=", 80), collapse = ""), "\n\n")
print(resumo)

# ==============================================================================
# 9. GRÁFICOS COMPARATIVOS
# ==============================================================================

cat("\n📊 Criando gráficos comparativos...\n")

# Beta comparativo
beta_C <- read_csv(file.path(output_dir, "COMPLETO", "beta_summary.csv"), show_col_types = FALSE)
beta_S <- read_csv(file.path(output_dir, "SIMPLIFICADO", "beta_summary.csv"), show_col_types = FALSE)
beta_all <- bind_rows(beta_C, beta_S)

ggsave(
  file.path(output_dir, "beta_COMPARATIVO.png"),
  ggplot(beta_all, aes(x = Covariate, y = Mean, color = Model)) +
    geom_point(position = position_dodge(width = 0.5), size = 4) +
    geom_errorbar(aes(ymin = HPD_Lower, ymax = HPD_Upper), 
                  position = position_dodge(width = 0.5), width = 0.3, linewidth = 1) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
    scale_color_manual(values = c("COMPLETO" = "steelblue", "SIMPLIFICADO" = "darkorange")) +
    theme_bw(base_size = 13) + theme(legend.position = "bottom") +
    labs(title = "Comparação dos Coeficientes beta", y = expression(beta)),
  width = 10, height = 6
)

# Gamma comparativo (facetas separadas por cluster)
gamma_C <- read_csv(file.path(output_dir, "COMPLETO", "gamma_summary.csv"), show_col_types = FALSE)
gamma_S <- read_csv(file.path(output_dir, "SIMPLIFICADO", "gamma_summary.csv"), show_col_types = FALSE)
gamma_all <- bind_rows(
  gamma_C %>% mutate(Model = "COMPLETO"),
  gamma_S %>% mutate(Model = "SIMPLIFICADO")
)
gamma_all$Cluster <- factor(paste("Cluster", gamma_all$Cluster), 
                            levels = paste("Cluster", 1:4))

ggsave(
  file.path(output_dir, "gamma_COMPARATIVO.png"),
  ggplot(gamma_all, aes(x = Time, y = Mean, color = Model)) +
    geom_line(linewidth = 0.8) +
    scale_color_manual(values = c("COMPLETO" = "steelblue", "SIMPLIFICADO" = "darkorange")) +
    facet_wrap(~ Cluster, scales = "free_y", ncol = 2) +
    scale_x_continuous(breaks = seq(1, n_times, by = 5), 
                       labels = anos_label[seq(1, n_times, by = 5)]) +
    theme_bw(base_size = 11) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
          legend.position = "bottom") +
    labs(title = "Comparação gamma[k,t]: COMPLETO vs SIMPLIFICADO",
         x = "Ano", y = expression(gamma[k][t])),
  width = 12, height = 8
)

# WAIC comparativo
ggsave(
  file.path(output_dir, "WAIC_COMPARATIVO.png"),
  ggplot(resumo, aes(x = Model, y = WAIC, fill = Model)) +
    geom_col(width = 0.5) +
    geom_text(aes(label = round(WAIC, 1)), vjust = -0.5, size = 5) +
    scale_fill_manual(values = c("COMPLETO" = "steelblue", "SIMPLIFICADO" = "darkorange")) +
    theme_bw(base_size = 13) + theme(legend.position = "none") +
    labs(title = "WAIC: COMPLETO vs SIMPLIFICADO", 
         subtitle = paste("ΔWAIC =", round(resumo$WAIC[1] - resumo$WAIC[2], 1)),
         y = "WAIC"),
  width = 6, height = 6
)

cat(sprintf("\nΔWAIC = %.2f\n", resumo$WAIC[1] - resumo$WAIC[2]))
cat("\n⏱️ Tempo total:", format(Sys.time() - inicio_global), "\n")