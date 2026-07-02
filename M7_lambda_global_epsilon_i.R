# ==============================================================================
# M7: Lambda[t] GLOBAL + Epsilon[i] FIXO via Gamma[k] Clássico (Uniforme)
# COMPLETO vs SIMPLIFICADO - SCRIPT ROBUSTO E OTIMIZADO
# ==============================================================================

inicio_global <- Sys.time()

setwd("C:/Users/vlara/OneDrive/Estatistica UFMG/Mestrado/Pesquisa/Aplicação/main")

pkgs <- c("nimble", "coda", "dplyr", "ggplot2", "tidyr", "readr", "stringr", "tibble", "scales")
for (pkg in pkgs) {
  if (!require(pkg, character.only = TRUE)) install.packages(pkg)
  library(pkg, character.only = TRUE)
}

# ==============================================================================
# 1. CARREGAR DADOS E ESTRUTURAR COVARIÁVEIS
# ==============================================================================
load("dados_modelo_final.RData")

E_norm    <- E / mean(E)
n_times   <- ncol(Y_mat)
K         <- 4
N_regions <- 75

# Configuração dos Cenários de Covariáveis
x_full  <- x
x_simpl <- x[, , c(1, 3), drop = FALSE]

h_mat <- structure(
  .Data = c(
    1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
    1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
  ), .Dim = c(75, 4)
)

region_names <- rownames(Y_mat)
anos_label   <- colnames(Y_mat)

grupo <- rep(NA, N_regions)
for (i in 1:N_regions) {
  if      (h_mat[i,1]==1 & h_mat[i,2]==0)                                       grupo[i] <- 1
  else if (h_mat[i,1]==1 & h_mat[i,2]==1 & h_mat[i,3]==0)                       grupo[i] <- 2
  else if (h_mat[i,1]==1 & h_mat[i,2]==1 & h_mat[i,3]==1 & h_mat[i,4]==0)      grupo[i] <- 3
  else if (all(h_mat[i,] == 1))                                                 grupo[i] <- 4
}

# Amostragem de microrregiões para monitoramento visual (3 por cluster)
set.seed(123)
regioes_por_cluster <- list()
for (g in 1:4) {
  idx_grupo <- which(grupo == g)
  if (g == 1) {
    bh     <- grep("BELO_HORIZONTE", region_names)[1]
    outras <- sample(setdiff(idx_grupo, bh), min(2, length(idx_grupo) - 1))
    regioes_por_cluster[[g]] <- c(bh, outras)
  } else if (g == 4) {
    aracuai <- grep("ARACUAI", region_names)[1]
    outras  <- sample(setdiff(idx_grupo, aracuai), min(2, length(idx_grupo) - 1))
    regioes_por_cluster[[g]] <- c(aracuai, outras)
  } else {
    regioes_por_cluster[[g]] <- sample(idx_grupo, min(3, length(idx_grupo)))
  }
}
all_regions <- unlist(regioes_por_cluster)

data_adj <- list(
  num = c(6,5,5,4,5,4,3,7,3,7,5,6,3,4,4,4,8,6,4,8,6,4,7,3,6,9,5,7,2,3,2,4,5,2,7,8,3,5,9,5,8,5,4,2,6,2,9,3,9,6,6,7,5,6,6,5,2,7,5,3,9,5,5,4,4,5,6,4,4,6,5,6,3,9,2),
  adj = c(2,3,5,12,17,18, 1,4,5,6,17, 1,5,7,12,19, 2,6,17,31,1,2,3,6,7, 2,4,5,7, 3,5,6, 9,12,15,19,21,30,58,8,15,19, 11,18,27,28,39,58,74, 10,17,18,20,28,1,3,8,18,19,58, 14,17,20, 13,20,23,42, 8,9,21,36,22,26,45,68, 1,2,4,11,13,18,20,31, 1,10,11,12,17,58,3,8,9,12, 11,13,14,17,23,28,67,73, 8,15,27,30,36,58,16,26,42,45, 14,20,32,35,42,67,73, 57,64,71,36,41,49,51,63,66, 16,22,33,45,47,62,65,68,70,10,21,36,39,58, 10,11,20,39,55,67,74, 41,48, 8,21,58,4,17, 23,35,55,67, 26,44,65,68,70, 40,54,23,32,42,45,47,55,56, 15,21,25,27,37,39,51,66,36,51,71, 40,43,52,54,72, 10,27,28,36,53,60,61,66,74,34,38,46,52,54, 25,29,48,49,51,54,64,72,14,22,23,35,45, 38,52,59,72, 33,68,16,22,26,35,42,47, 40,52,26,35,45,49,50,56,61,62,74, 29,41,54,25,41,47,50,59,61,63,69,72, 47,49,52,59,62,70,25,36,37,41,64,71, 38,40,43,46,50,59,70,39,61,63,66,69, 34,38,40,41,48,72,28,32,35,56,67,74, 35,47,55,61,74, 24,71,8,10,12,18,21,27,30, 43,49,50,52,72, 39,61,74,39,47,49,53,56,60,69,74,75, 26,47,50,65,70,25,49,53,66,69, 24,41,51,71, 26,33,62,70,25,36,39,53,63, 20,23,28,32,55,73, 16,26,33,44,49,53,61,63, 26,33,50,52,62,65, 24,37,51,57,64,38,41,43,49,54,59, 20,23,67,10,28,39,47,55,56,61,60,75, 61,74)
)
n_adj <- length(data_adj$adj)

safe_hpd <- function(sv) {
  if (var(sv) < 1e-12) return(c(NA_real_, NA_real_))
  as.numeric(HPDinterval(as.mcmc(sv), prob = 0.95))
}
safe_gelman <- function(obj) {
  tryCatch(gelman.diag(obj, autoburnin = FALSE)$psrf[, 1], error = function(e) rep(NA_real_, nvar(obj)))
}
make_label <- function(nms) {
  sapply(nms, function(n) {
    nc <- gsub("_", " ", n)
    if (nchar(nc) > 30) paste0(substr(nc, 1, 27), "...") else nc
  })
}
cores_cluster <- c("Cluster 1" = "#1f77b4", "Cluster 2" = "#ff7f0e", "Cluster 3" = "#2ca02c", "Cluster 4" = "#d62728")

ordered_idx    <- all_regions[order(grupo[all_regions])]
ordered_labels <- unique(paste(make_label(region_names[ordered_idx]), "- Cluster", grupo[ordered_idx]))

# ==============================================================================
# FUNÇÃO DA ENGENHARIA DE EXECUÇÃO - M7 CLÁSSICO
# ==============================================================================
run_model <- function(model_type, output_dir_base) {
  
  cat(sprintf("\n%s=== INICIANDO EXECUÇÃO DO CENÁRIO M7 %s ===%s\n", 
              paste(rep("=", 20), collapse = ""), model_type, paste(rep("=", 20), collapse = "")))
  
  if (model_type == "COMPLETO") {
    x_data <- x_full; p <- 3
    cov_names <- c("prenatal", "instrucao", "baixo_peso")
  } else {
    x_data <- x_simpl; p <- 2
    cov_names <- c("prenatal", "baixo_peso")
  }
  
  a0 <- 1.0; b0 <- 1.0; w <- 0.85; a_unif <- 0.0; b_unif <- 0.05
  
  constants <- list(
    n_regions = N_regions, n_times = n_times, p = p, K = K,
    a0 = a0, b0 = b0, w = w, a_unif = a_unif, b_unif = b_unif,
    adj = data_adj$adj, num = data_adj$num, weights = rep(1, n_adj), n_adj = n_adj
  )
  data_nimble <- list(Y = Y_mat, E = E_norm, x = x_data, h = h_mat)
  
  # Código Nimble - ESTRUTURA CLÁSSICA DO GAMMA (Uniforme com restrição)
  code_M7 <- nimbleCode({
    for (j in 1:p) { beta[j] ~ dnorm(0, sd = 1) }
    
    # Gamma Estático via Uniforme
    gamma[1] ~ dunif(a_unif, b_unif)
    for (k in 2:K) {
      gamma[k] ~ dunif(0, 1 - sum(gamma[1:(k-1)]))
    }
    
    # Epsilon Fixo
    for (i in 1:n_regions) {
      epsilon[i] <- 1 - inprod(h[i, 1:K], gamma[1:K])
    }
    
    sigma_s ~ T(dt(0, 1, 1), 0, )
    tau_s   <- 1 / (sigma_s^2)
    s[1:n_regions] ~ dcar_normal(adj[1:n_adj], weights[1:n_adj], num[1:n_regions], tau_s, zero_mean = 1)
    
    for (t in 1:n_times) { lambda[t] ~ dgamma(a0, b0) }
    
    for (i in 1:n_regions) {
      for (t in 1:n_times) {
        log_mu[i, t] <- log(lambda[t]) + log(E[i, t]) + log(epsilon[i]) + inprod(beta[1:p], x[i, t, 1:p]) + s[i]
        mu[i, t]       <- exp(log_mu[i, t])
        Y[i, t]        ~ dpois(mu[i, t])
        logLik_Y[i, t] <- dpois(Y[i, t], mu[i, t], log = TRUE)
      }
    }
  })
  
  ffbs_lambda_global_M7 <- nimbleFunction(
    contains = sampler_BASE,
    setup = function(model, mvSaved, target, control) {
      n_times <- control$n_times; n_regions <- control$n_regions; p <- control$p
      a0 <- control$a0; b0 <- control$b0; w <- control$w
      at_buf <- nimNumeric(n_times + 1, 0); bt_buf <- nimNumeric(n_times + 1, 0)
      calcNodes <- model$getDependencies(target, self = FALSE)
      targetNodes <- model$expandNodeNames(target)
      setupOutputs(at_buf, bt_buf)
    },
    run = function() {
      declare(t, integer()); declare(t_idx, integer()); declare(t_back, integer()); declare(i, integer()); declare(k, integer())
      declare(sum_Y, double()); declare(sum_g, double()); declare(g_it, double()); declare(lin, double()); declare(nu, double())
      
      at_buf[1] <<- a0; bt_buf[1] <<- b0
      for (t in 1:n_times) {
        sum_Y <- 0.0; sum_g <- 0.0
        for (i in 1:n_regions) {
          sum_Y <- sum_Y + model$Y[i, t]
          lin <- 0.0
          for (k in 1:p) lin <- lin + model$x[i, t, k] * model$beta[k]
          g_it <- model$E[i, t] * model$epsilon[i] * exp(lin + model$s[i])
          sum_g <- sum_g + g_it
        }
        at_buf[t + 1] <<- w * at_buf[t] + sum_Y
        bt_buf[t + 1] <<- w * bt_buf[t] + sum_g
      }
      model$lambda[n_times] <<- rgamma(1, shape = at_buf[n_times + 1], rate = max(bt_buf[n_times + 1], 1e-10))
      for (t_idx in 1:(n_times - 1)) {
        t_back <- n_times - t_idx
        nu <- rgamma(1, shape = (1 - w) * at_buf[t_back + 1], rate = max(bt_buf[t_back + 1], 1e-10))
        model$lambda[t_back] <<- nu + w * model$lambda[t_back + 1]
      }
      model$calculate(calcNodes)
      copy(from = model, to = mvSaved, row = 1, nodes = targetNodes, logProb = TRUE)
    },
    methods = list(reset = function() {})
  )
  
  set.seed(42)
  gamma_init_vec <- c(0.02, 0.005, 0.005, 0.005)
  
  inits_1 <- list(lambda = rep(1.0, n_times), beta = rep(0, p), gamma = gamma_init_vec, sigma_s = 0.5, s = rep(0, N_regions))
  inits_2 <- list(lambda = rgamma(n_times, a0, b0), beta = rnorm(p, 0, 0.3), gamma = gamma_init_vec * 0.5, sigma_s = 1.0, s = rep(0, N_regions))
  
  model  <- nimbleModel(code = code_M7, constants = constants, data = data_nimble, inits = inits_1, check = FALSE)
  Cmodel <- compileNimble(model)
  
  conf <- configureMCMC(model)
  conf$removeSamplers("lambda")
  conf$addSampler(target = paste0("lambda[1:", n_times, "]"), type = ffbs_lambda_global_M7, 
                  control = list(n_times = n_times, n_regions = N_regions, p = p, a0 = a0, b0 = b0, w = w))
  
  # Voltando aos samplers clássicos de Slice para o Gamma Uniforme
  conf$removeSamplers("gamma")
  for (k in 1:K) conf$addSampler(target = paste0("gamma[", k, "]"), type = "slice")
  
  conf$addMonitors(c("beta", "gamma", "epsilon", "lambda", "logLik_Y", "s", "sigma_s"))
  
  Rmcmc <- buildMCMC(conf)
  Cmcmc <- compileNimble(Rmcmc, project = model)
  
  niter <- 50000; nburnin <- 10000; nchains <- 2; thin <- 10
  
  t_ini <- Sys.time()
  samples <- runMCMC(Cmcmc, niter = niter, nburnin = nburnin, nchains = nchains, thin = thin,
                     inits = list(inits_1, inits_2), samplesAsCodaMCMC = TRUE, summary = FALSE, WAIC = FALSE)
  t_fim <- Sys.time()
  tempo_min <- as.numeric(difftime(t_fim, t_ini, units = "mins"))
  
  scenario_dir <- file.path(output_dir_base, model_type)
  dir.create(scenario_dir, recursive = TRUE, showWarnings = FALSE)
  saveRDS(samples, file.path(scenario_dir, "samples.rds"))
  
  samples_mat <- as.matrix(samples)
  mcmc_list_full <- mcmc.list(lapply(1:nchains, function(ch) as.mcmc(samples[[ch]])))
  
  # ============================================================================
  # PROCESSAMENTO DE SUMÁRIOS
  # ============================================================================
  beta_names <- paste0("beta[", seq_len(p), "]")
  beta_summary <- do.call(rbind, lapply(beta_names, function(nm) {
    sv <- samples_mat[, nm]; hpd <- safe_hpd(sv)
    data.frame(Model = model_type, Parameter = nm, Covariate = cov_names[which(beta_names == nm)],
               Mean = mean(sv), SD = sd(sv), HPD_Lower = hpd[1], HPD_Upper = hpd[2],
               ESS = as.numeric(effectiveSize(mcmc_list_full[, nm])), Rhat = safe_gelman(mcmc_list_full[, nm]), stringsAsFactors = FALSE)
  }))
  write_csv(beta_summary, file.path(scenario_dir, "beta_summary.csv"))
  
  lambda_names <- paste0("lambda[", 1:n_times, "]")
  lambda_summary <- do.call(rbind, lapply(1:n_times, function(t) {
    nm <- lambda_names[t]; sv <- samples_mat[, nm]; hpd <- safe_hpd(sv)
    data.frame(Time = t, Ano = anos_label[t], Mean = mean(sv), SD = sd(sv), HPD_Lower = hpd[1], HPD_Upper = hpd[2],
               ESS = as.numeric(effectiveSize(mcmc_list_full[, nm])), Rhat = safe_gelman(mcmc_list_full[, nm]), stringsAsFactors = FALSE)
  }))
  write_csv(lambda_summary, file.path(scenario_dir, "lambda_global_summary.csv"))
  
  gamma_names_k <- paste0("gamma[", 1:K, "]")
  gamma_summary <- do.call(rbind, lapply(1:K, function(k) {
    nm <- gamma_names_k[k]; sv <- samples_mat[, nm]; hpd <- safe_hpd(sv)
    data.frame(Cluster = k, Mean = mean(sv), SD = sd(sv), HPD_Lower = hpd[1], HPD_Upper = hpd[2],
               ESS = as.numeric(effectiveSize(mcmc_list_full[, nm])), Rhat = safe_gelman(mcmc_list_full[, nm]), stringsAsFactors = FALSE)
  }))
  write_csv(gamma_summary, file.path(scenario_dir, "gamma_summary.csv"))
  
  s_names <- paste0("s[", 1:N_regions, "]")
  s_summary <- do.call(rbind, lapply(1:N_regions, function(i) {
    nm <- s_names[i]; sv <- samples_mat[, nm]; hpd <- safe_hpd(sv)
    data.frame(Region_Index = i, Region = region_names[i], Cluster = grupo[i], Mean = mean(sv), SD = sd(sv), HPD_Lower = hpd[1], HPD_Upper = hpd[2], stringsAsFactors = FALSE)
  }))
  write_csv(s_summary, file.path(scenario_dir, "spatial_effects_summary.csv"))
  
  # ============================================================================
  # SUBSTITUA A PARTIR DAQUI:
  # ============================================================================
  
  epsilon_summary_all <- do.call(rbind, lapply(1:N_regions, function(i) {
    nm <- paste0("epsilon[", i, "]"); sv <- samples_mat[, nm]; hpd <- safe_hpd(sv) # <-- Espaço removido aqui
    data.frame(Region_Index = i, Region = region_names[i], Cluster = grupo[i], Mean = mean(sv), SD = sd(sv), HPD_Lower = hpd[1], HPD_Upper = hpd[2], stringsAsFactors = FALSE)
  }))
  write_csv(epsilon_summary_all, file.path(scenario_dir, "epsilon_summary.csv"))
  
  # DIAGNÓSTICO ACF
  cat(sprintf("[%s] Calculando diagnósticos estendidos de autocorrelação (ACF)...\n", model_type))
  all_params_diag <- c(beta_names, gamma_names_k, "sigma_s", lambda_names)
  all_params_diag <- all_params_diag[all_params_diag %in% colnames(samples_mat)]
  
  acf_results <- do.call(rbind, lapply(all_params_diag, function(nm) {
    ac   <- acf(samples_mat[, nm], lag.max = 200, plot = FALSE)
    lags <- as.vector(ac$lag[-1]); acfs <- as.vector(ac$acf[-1])
    ess_v  <- tryCatch(as.numeric(effectiveSize(mcmc_list_full[, nm])), error = function(e) NA_real_)
    rhat_v <- tryCatch(safe_gelman(mcmc_list_full[, nm]),               error = function(e) NA_real_)
    tibble(Model = model_type, Parameter = nm, ESS = ess_v, Rhat = rhat_v,
           lag_0.10 = { v <- lags[which(abs(acfs) < 0.10)[1]]; ifelse(is.na(v), Inf, v) },
           lag_0.05 = { v <- lags[which(abs(acfs) < 0.05)[1]]; ifelse(is.na(v), Inf, v) },
           acf_lag1 = acfs[1])
  }))
  write_csv(acf_results, file.path(scenario_dir, "acf_diagnostics.csv"))
  
  # WAIC / LPML
  loglik_names <- grep("^logLik_Y", colnames(samples_mat), value = TRUE)
  lm_mat <- samples_mat[, loglik_names, drop = FALSE]
  lppd   <- sum(apply(lm_mat, 2, function(x) { mx <- max(x); mx + log(mean(exp(x - mx))) }))
  p_waic <- sum(apply(lm_mat, 2, var))
  waic   <- -2 * (lppd - p_waic)
  LPML   <- sum(log(1 / apply(lm_mat, 2, function(x) mean(exp(-x)))))
  
  # ============================================================================
  # GERAÇÃO DE GRÁFICOS
  # ============================================================================
  p_lam <- ggplot(lambda_summary, aes(x = Time, y = Mean)) +
    geom_line(color = "blue", linewidth = 1) + geom_ribbon(aes(ymin = HPD_Lower, ymax = HPD_Upper), fill = "blue", alpha = 0.15) +
    scale_x_continuous(breaks = seq(1, n_times, by = 5), labels = anos_label[seq(1, n_times, by = 5)]) +
    theme_bw() + labs(title = paste("Trajetória Lambda Global (M7) -", model_type), x = "Ano", y = expression(lambda[t]))
  ggsave(file.path(scenario_dir, "plot_lambda_global.png"), p_lam, width = 8, height = 5)
  
  mu_panel <- do.call(rbind, lapply(all_regions, function(i) {
    label_str <- paste(make_label(region_names[i]), "- Cluster", grupo[i])
    eps_i     <- samples_mat[, paste0("epsilon[", i, "]")] # <-- Espaço removido aqui
    s_i       <- samples_mat[, paste0("s[", i, "]")]       # <-- Espaço removido aqui
    do.call(rbind, lapply(1:n_times, function(t) {
      ldraws   <- samples_mat[, paste0("lambda[", t, "]")] # <-- Espaço removido aqui
      bdraws   <- samples_mat[, beta_names, drop = FALSE]
      lin_pred <- as.vector(bdraws %*% as.numeric(x_data[i, t, ]))
      mu       <- ldraws * E_norm[i, t] * eps_i * exp(lin_pred + s_i)
      hpd      <- safe_hpd(mu)
      data.frame(Region = i, FacetLabel = label_str, Time = t, Mean = mean(mu), Lower = hpd[1], Upper = hpd[2], Y_obs = Y_mat[i, t], stringsAsFactors = FALSE)
    }))
  }))
  mu_panel$FacetLabel <- factor(mu_panel$FacetLabel, levels = ordered_labels)
  
  p_mu <- ggplot(mu_panel, aes(x = Time)) +
    geom_ribbon(aes(ymin = Lower, ymax = Upper), fill = "steelblue", alpha = 0.3) +
    geom_line(aes(y = Mean), color = "steelblue", linewidth = 0.8) +
    geom_point(aes(y = Y_obs), shape = 21, fill = "white", color = "black", size = 1.2) +
    facet_wrap(~ FacetLabel, ncol = 3, scales = "free_y") +
    scale_x_continuous(breaks = seq(1, n_times, by = 5), labels = anos_label[seq(1, n_times, by = 5)]) +
    theme_bw(base_size = 10) + theme(axis.text.x = element_text(angle = 45, hjust = 1), strip.text = element_text(size = 9)) +
    labs(title = paste("M7 — mu[i,t]: Esperados vs Observados (", model_type, ")"), x = "Ano", y = expression(mu[it]))
  ggsave(file.path(scenario_dir, "painel_mu_por_cluster.png"), p_mu, width = 14, height = 10)
  
  eps_ord <- epsilon_summary_all[order(epsilon_summary_all$Mean), ]
  eps_ord$Nome        <- factor(eps_ord$Region, levels = eps_ord$Region)
  eps_ord$Cluster_lab <- factor(paste("Cluster", eps_ord$Cluster), levels = paste("Cluster", 1:4))
  
  p_eps_reg <- ggplot(eps_ord, aes(x = Nome, y = Mean, color = Cluster_lab)) +
    geom_point(size = 1.5) + geom_errorbar(aes(ymin = HPD_Lower, ymax = HPD_Upper), width = 0.4, linewidth = 0.3) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "red") + scale_color_manual(values = cores_cluster) +
    theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 5), legend.position = "bottom") +
    labs(title = paste("M7 — Epsilon[i] por Microrregião (", model_type, ")"), y = expression(epsilon[i]), x = "")
  ggsave(file.path(scenario_dir, "painel_epsilon_por_regiao.png"), p_eps_reg, width = 18, height = 8)
  
  p_eps_box <- ggplot(eps_ord, aes(x = Cluster_lab, y = Mean, fill = Cluster_lab)) +
    geom_boxplot(alpha = 0.6) + scale_fill_manual(values = cores_cluster) +
    theme_bw(base_size = 12) + theme(legend.position = "none") +
    labs(title = paste("M7 — Distribuição de Epsilon[i] por Cluster (", model_type, ")"), x = "Cluster", y = expression(epsilon[i]))
  ggsave(file.path(scenario_dir, "painel_epsilon_por_cluster.png"), p_eps_box, width = 8, height = 6)
  
  params_diag_tbl <- c(beta_names, "sigma_s", gamma_names_k, lambda_names)
  diag_table <- do.call(rbind, lapply(params_diag_tbl, function(nm) {
    data.frame(Parameter = nm, ESS  = as.numeric(effectiveSize(mcmc_list_full[, nm])), Rhat = safe_gelman(mcmc_list_full[, nm]), stringsAsFactors = FALSE)
  }))
  diag_plot <- diag_table %>% mutate(
    ESS_cat  = case_when(ESS < 100 ~ "Ruim (<100)", ESS < 500 ~ "Razoável (100-500)", ESS < 1000 ~ "Bom (500-1000)", TRUE ~ "Excelente (>1000)"),
    Rhat_cat = case_when(Rhat > 1.1 ~ "Preocupante (>1.1)", Rhat > 1.05 ~ "Aceitável (1.05-1.1)", TRUE ~ "Bom (<1.05)")
  )
  p_ess <- ggplot(diag_plot, aes(x = Parameter, y = ESS, fill = ESS_cat)) +
    geom_col() + geom_hline(yintercept = c(100, 500, 1000), linetype = "dashed", color = "red", alpha = 0.5) +
    scale_fill_manual(values = c("Ruim (<100)" = "red", "Razoável (100-500)" = "orange", "Bom (500-1000)" = "yellow", "Excelente (>1000)" = "green")) +
    theme_bw(base_size = 9) + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7), legend.position = "bottom") +
    labs(title = paste("M7 — ESS por Parâmetro (", model_type, ")"), x = "", y = "ESS")
  ggsave(file.path(scenario_dir, "diagnostico_ess_rhat.png"), p_ess, width = 14, height = 8)
  
  return(data.frame(Model = model_type, p = p, WAIC = waic, LPML = LPML,
                    ESS_beta_min = min(beta_summary$ESS, na.rm = TRUE),
                    Rhat_beta_max = max(beta_summary$Rhat, na.rm = TRUE),
                    tempo_min = tempo_min, stringsAsFactors = FALSE))
} # <-- FINAL DA FUNÇÃO run_model

# ==============================================================================
# EXECUÇÃO E COMPARAÇÃO ENTRE CENÁRIOS
# ==============================================================================
output_dir <- "C:/Users/vlara/OneDrive/Estatistica UFMG/Mestrado/Pesquisa/Aplicação/resultados_M7_lambda_global_epsilon_i"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

resultados <- list()
resultados[["COMPLETO"]]     <- run_model("COMPLETO", output_dir)
resultados[["SIMPLIFICADO"]] <- run_model("SIMPLIFICADO", output_dir)

resumo <- bind_rows(resultados)
write_csv(resumo, file.path(output_dir, "resumo_comparativo.csv"))
print(resumo)

# ==============================================================================
# GRÁFICOS COMPARATIVOS GLOBAIS (SALVOS NA RAIZ)
# ==============================================================================
beta_C <- read_csv(file.path(output_dir, "COMPLETO", "beta_summary.csv"), show_col_types = FALSE)
beta_S <- read_csv(file.path(output_dir, "SIMPLIFICADO", "beta_summary.csv"), show_col_types = FALSE)
beta_all <- bind_rows(beta_C, beta_S)

ggsave(file.path(output_dir, "beta_COMPARATIVO.png"),
       ggplot(beta_all, aes(x = Covariate, y = Mean, color = Model)) +
         geom_point(position = position_dodge(width = 0.5), size = 4) +
         geom_errorbar(aes(ymin = HPD_Lower, ymax = HPD_Upper), position = position_dodge(width = 0.5), width = 0.3, linewidth = 1) +
         geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
         scale_color_manual(values = c("COMPLETO" = "steelblue", "SIMPLIFICADO" = "darkorange")) +
         theme_bw(base_size = 13) + theme(legend.position = "bottom") +
         labs(title = "Comparação Coeficientes beta - M7", y = expression(beta)),
       width = 10, height = 6)

gamma_C <- read_csv(file.path(output_dir, "COMPLETO", "gamma_summary.csv"), show_col_types = FALSE)
gamma_S <- read_csv(file.path(output_dir, "SIMPLIFICADO", "gamma_summary.csv"), show_col_types = FALSE)
gamma_all <- bind_rows(gamma_C %>% mutate(Model = "COMPLETO"), gamma_S %>% mutate(Model = "SIMPLIFICADO"))
gamma_all$Cluster_lab <- factor(paste("Cluster", gamma_all$Cluster), levels = paste("Cluster", 1:4))

ggsave(file.path(output_dir, "gamma_COMPARATIVO.png"),
       ggplot(gamma_all, aes(x = Cluster_lab, y = Mean, fill = Model)) +
         geom_col(position = position_dodge(width = 0.7), width = 0.6) +
         geom_errorbar(aes(ymin = HPD_Lower, ymax = HPD_Upper), position = position_dodge(width = 0.7), width = 0.2, linewidth = 0.8) +
         scale_fill_manual(values = c("COMPLETO" = "steelblue", "SIMPLIFICADO" = "darkorange")) +
         theme_bw(base_size = 13) + theme(legend.position = "bottom") +
         labs(title = "Comparação gamma[k] Fixo: COMPLETO vs SIMPLIFICADO - M7", x = "Cluster", y = expression(gamma[k])),
       width = 9, height = 6)

ggsave(file.path(output_dir, "WAIC_COMPARATIVO.png"),
       ggplot(resumo, aes(x = Model, y = WAIC, fill = Model)) +
         geom_col(width = 0.5) + geom_text(aes(label = round(WAIC, 1)), vjust = -0.5, size = 5) +
         scale_fill_manual(values = c("COMPLETO" = "steelblue", "SIMPLIFICADO" = "darkorange")) +
         theme_bw(base_size = 13) + theme(legend.position = "none") +
         labs(title = "WAIC: COMPLETO vs SIMPLIFICADO - M7", subtitle = paste("ΔWAIC =", round(resumo$WAIC[1] - resumo$WAIC[2], 1)), y = "WAIC"),
       width = 6, height = 6)