# ==============================================================================
# ajuste_cluster_COMPLETO_vs_SIMPLIFICADO_FINAL.R
# Versão COMPLETA com todos os diagnósticos e gráficos
# Pronto para Cluster 1 e Cluster 4
# ==============================================================================

inicio_global <- Sys.time()

setwd("C:/Users/vlara/OneDrive/Estatistica UFMG/Mestrado/Pesquisa/Aplicação/main")

pkgs <- c("nimble", "coda", "dplyr", "ggplot2", "tidyr", 
          "readr", "stringr", "tibble")
for (pkg in pkgs) {
  if (!require(pkg, character.only = TRUE)) install.packages(pkg)
  library(pkg, character.only = TRUE)
}

# ==============================================================================
# CONFIGURAÇÃO: Escolha o cluster aqui!
# ==============================================================================
CLUSTER_ESCOLHIDO <- 4  # Mude para 1 ou 4
# Cluster 1: Alta adequação (28 microrregiões) - epsilon ~ U(0.95, 1.00)
# Cluster 4: Baixa adequação (17 microrregiões) - epsilon ~ U(0.75, 0.85)

if (CLUSTER_ESCOLHIDO == 1) {
  EPS_LOWER <- 0.95
  EPS_UPPER <- 1.00
  CLUSTER_NOME <- "Cluster1_Alta"
} else if (CLUSTER_ESCOLHIDO == 4) {
  EPS_LOWER <- 0.75
  EPS_UPPER <- 0.85
  CLUSTER_NOME <- "Cluster4_Baixa"
} else {
  stop("Escolha CLUSTER_ESCOLHIDO = 1 ou 4")
}

cat(sprintf("\n🎯 Analisando %s\n", CLUSTER_NOME))
cat(sprintf("   epsilon[t] ~ U(%.2f, %.2f)\n\n", EPS_LOWER, EPS_UPPER))

# ==============================================================================
# 1. CARREGAR DADOS
# ==============================================================================

load("dados_modelo_final.RData")

E_norm <- E / mean(E)
n_times <- ncol(Y_mat)
p_full <- dim(x)[3]

# Matriz h
h_mat <- structure(
  .Data = c(
    1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
    1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
  ),
  .Dim = c(75, 4)
)

# Selecionar cluster
if (CLUSTER_ESCOLHIDO == 1) {
  cluster_idx <- which(h_mat[,1] == 1 & h_mat[,2] == 0)  # Só cluster 1
} else {
  cluster_idx <- which(h_mat[,1] == 1 & h_mat[,2] == 1 & 
                         h_mat[,3] == 1 & h_mat[,4] == 1)  # Todos clusters
}

N_cluster <- length(cluster_idx)

cat(sprintf("✅ %s: %d microrregiões\n", CLUSTER_NOME, N_cluster))

# Filtrar dados
Y_cluster <- Y_mat[cluster_idx, , drop = FALSE]
E_cluster <- E_norm[cluster_idx, , drop = FALSE]
x_full <- x[cluster_idx, , , drop = FALSE]
x_simpl <- x_full[, , c(1, 3), drop = FALSE]

region_names <- rownames(Y_cluster)
anos_label <- colnames(Y_mat)

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

# ==============================================================================
# 3. SUBGRAFO DE ADJACÊNCIA
# ==============================================================================

old_to_new <- setNames(seq_along(cluster_idx), cluster_idx)

num_sub <- integer(N_cluster)
adj_sub <- c()
weights_sub <- c()

for (i_new in seq_len(N_cluster)) {
  i_old <- cluster_idx[i_new]
  
  if (data_adj$num[i_old] > 0) {
    start_idx <- if (i_old == 1) 1 else sum(data_adj$num[1:(i_old-1)]) + 1
    end_idx   <- start_idx + data_adj$num[i_old] - 1
    neighbors_old <- data_adj$adj[start_idx:end_idx]
    
    neighbors_new <- old_to_new[as.character(neighbors_old)]
    neighbors_new <- neighbors_new[!is.na(neighbors_new)]
    
    num_sub[i_new] <- length(neighbors_new)
    if (length(neighbors_new) > 0) {
      adj_sub <- c(adj_sub, neighbors_new)
      weights_sub <- c(weights_sub, rep(1, length(neighbors_new)))
    }
  }
}

n_adj_sub <- length(adj_sub)

cat(sprintf("Subgrafo: %d arestas\n", n_adj_sub))

# ==============================================================================
# 4. REGIÕES PARA MONITORAR
# ==============================================================================

# Encontrar região de referência
if (CLUSTER_ESCOLHIDO == 1) {
  ref_idx <- grep("BELO_HORIZONTE", region_names)[1]
} else {
  ref_idx <- grep("ARACUAI", region_names)[1]
}

set.seed(42)
n_monitor <- min(6, N_cluster)
REGIONS_INTEREST <- c(ref_idx, sample(setdiff(seq_len(N_cluster), ref_idx), n_monitor - 1))
REGIONS_INTEREST <- sort(REGIONS_INTEREST)

cat("Regiões monitoradas:", paste(region_names[REGIONS_INTEREST], collapse = ", "), "\n\n")

# ==============================================================================
# 5. FUNÇÃO WORKER (COMPLETA)
# ==============================================================================

run_model <- function(model_type, output_dir_base) {
  
  cat(sprintf("\n%s=== MODELO %s ===%s\n", 
              paste(rep("=", 40), collapse = ""), model_type, 
              paste(rep("=", 40), collapse = "")))
  
  if (model_type == "COMPLETO") {
    x_data <- x_full
    p <- 3
    cov_names <- c("prenatal", "instrucao", "baixo_peso")
  } else {
    x_data <- x_simpl
    p <- 2
    cov_names <- c("prenatal", "baixo_peso")
  }
  
  a0 <- 1.0; b0 <- 1.0; w <- 0.85
  
  constants <- list(
    n_regions = N_cluster, n_times = n_times, p = p,
    a0 = a0, b0 = b0, w = w,
    eps_lower = EPS_LOWER, eps_upper = EPS_UPPER,
    adj = adj_sub, num = num_sub, weights = weights_sub, n_adj = n_adj_sub
  )
  
  data_nimble <- list(Y = Y_cluster, E = E_cluster, x = x_data)
  
  code <- nimbleCode({
    for (j in 1:p) { beta[j] ~ dnorm(0, sd = 1) }
    for (t in 1:n_times) { epsilon[t] ~ dunif(eps_lower, eps_upper) }
    sigma_s ~ T(dt(0, 1, 1), 0, )
    tau_s   <- 1 / (sigma_s^2)
    s[1:n_regions] ~ dcar_normal(
      adj[1:n_adj], weights[1:n_adj], num[1:n_regions], tau_s, zero_mean = 1
    )
    for (i in 1:n_regions) {
      for (t in 1:n_times) { lambda[i, t] ~ dgamma(a0, b0) }
    }
    for (i in 1:n_regions) {
      for (t in 1:n_times) {
        mu[i, t] <- lambda[i, t] * E[i, t] * epsilon[t] *
          exp(inprod(beta[1:p], x[i, t, 1:p]) + s[i])
        Y[i, t] ~ dpois(mu[i, t])
        logLik_Y[i, t] <- dpois(Y[i, t], mu[i, t], log = TRUE)
      }
    }
  })
  
  # FFBS para lambda
  ffbs_lambda <- nimbleFunction(
    contains = sampler_BASE,
    setup = function(model, mvSaved, target, control) {
      n_times  <- control$n_times; p <- control$p
      a0 <- control$a0; b0 <- control$b0; w <- control$w
      region_i <- control$region_i
      at_buf <- nimNumeric(n_times + 1, 0)
      bt_buf <- nimNumeric(n_times + 1, 0)
      calcNodes   <- model$getDependencies(target, self = FALSE)
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
        g_it <- model$E[region_i, t] * model$epsilon[t] * 
          exp(prod_val + model$s[region_i])
        at_buf[t + 1] <<- w * at_buf[t] + model$Y[region_i, t]
        bt_buf[t + 1] <<- w * bt_buf[t] + g_it
      }
      
      model$lambda[region_i, n_times] <<- rgamma(1,
                                                 shape = at_buf[n_times + 1], rate = max(bt_buf[n_times + 1], 1e-10))
      
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
  lambda_init <- matrix(1.0, nrow = N_cluster, ncol = n_times)
  
  inits_1 <- list(
    lambda  = lambda_init,
    beta    = rep(0, p),
    epsilon = rep((EPS_LOWER + EPS_UPPER) / 2, n_times),
    sigma_s = 0.5,
    s       = rep(0, N_cluster)
  )
  inits_2 <- list(
    lambda  = matrix(rgamma(N_cluster * n_times, a0, b0), N_cluster, n_times),
    beta    = rnorm(p, 0, 0.3),
    epsilon = runif(n_times, EPS_LOWER + 0.01, EPS_UPPER - 0.01),
    sigma_s = 1.0,
    s       = rep(0, N_cluster)
  )
  
  # Compilar
  cat(sprintf("[%s] Compilando...\n", model_type))
  model <- nimbleModel(code = code, constants = constants,
                       data = data_nimble, inits = inits_1, check = FALSE)
  Cmodel <- compileNimble(model)
  
  conf <- configureMCMC(model)
  
  conf$removeSamplers("lambda")
  for (i in seq_len(N_cluster)) {
    conf$addSampler(
      target = paste0("lambda[", i, ", 1:", n_times, "]"),
      type   = ffbs_lambda,
      control = list(n_times = n_times, p = p, a0 = a0, b0 = b0, w = w, region_i = i)
    )
  }
  
  conf$removeSamplers("epsilon")
  for (t in seq_len(n_times)) {
    conf$addSampler(target = paste0("epsilon[", t, "]"), type = "slice")
  }
  
  mu_monitors <- unlist(lapply(REGIONS_INTEREST,
                               function(r) paste0("mu[", r, ", ", seq_len(n_times), "]")))
  
  conf$addMonitors(c("beta", "epsilon", "lambda", "logLik_Y", "s", "sigma_s", "tau_s"))
  conf$addMonitors(mu_monitors)
  
  Rmcmc <- buildMCMC(conf)
  Cmcmc <- compileNimble(Rmcmc, project = model)
  
  # Executar
  niter   <- 50000
  nburnin <- 10000
  nchains <- 2
  thin    <- 10
  
  cat(sprintf("[%s] MCMC: %d iter, %d burnin, %d chains, thin=%d\n",
              model_type, niter, nburnin, nchains, thin))
  
  t_inicio <- Sys.time()
  
  samples <- runMCMC(
    Cmcmc, niter = niter, nburnin = nburnin, nchains = nchains, thin = thin,
    inits = list(inits_1, inits_2), samplesAsCodaMCMC = TRUE,
    summary = FALSE, WAIC = FALSE
  )
  
  t_fim <- Sys.time()
  cat(sprintf("[%s] Tempo: %.1f min\n", model_type, 
              difftime(t_fim, t_inicio, units = "mins")))
  
  # Salvar
  scenario_dir <- file.path(output_dir_base, model_type)
  dir.create(scenario_dir, recursive = TRUE, showWarnings = FALSE)
  saveRDS(samples, file.path(scenario_dir, "samples.rds"))
  
  # =========================================================================
  # PROCESSAMENTO COMPLETO DOS RESULTADOS
  # =========================================================================
  
  samples_mat    <- as.matrix(samples)
  mcmc_list_full <- mcmc.list(lapply(seq_len(nchains),
                                     function(ch) as.mcmc(samples[[ch]])))
  
  safe_hpd <- function(sv) {
    if (var(sv) < 1e-12) return(c(NA_real_, NA_real_))
    as.numeric(HPDinterval(as.mcmc(sv), prob = 0.95))
  }
  safe_gelman <- function(obj) {
    tryCatch(gelman.diag(obj, autoburnin = FALSE)$psrf[, 1],
             error = function(e) rep(NA_real_, nvar(obj)))
  }
  
  beta_names    <- paste0("beta[", seq_len(p), "]")
  epsilon_names <- paste0("epsilon[", seq_len(n_times), "]")
  lambda_names_all <- paste0("lambda[", rep(seq_len(N_cluster), each = n_times),
                             ", ", rep(seq_len(n_times), N_cluster), "]")
  s_names       <- paste0("s[", seq_len(N_cluster), "]")
  
  # ── BETA ─────────────────────────────────────────────────
  beta_summary <- do.call(rbind, lapply(beta_names, function(nm) {
    sv <- samples_mat[, nm]; hpd <- safe_hpd(sv)
    tibble(Model = model_type, Parameter = nm, 
           Covariate = cov_names[which(beta_names == nm)],
           Mean = mean(sv), SD = sd(sv),
           HPD_Lower = hpd[1], HPD_Upper = hpd[2],
           ESS  = as.numeric(effectiveSize(mcmc_list_full[, nm])),
           Rhat = safe_gelman(mcmc_list_full[, nm]))
  }))
  write_csv(beta_summary, file.path(scenario_dir, "beta_summary.csv"))
  
  # ── EPSILON ──────────────────────────────────────────────
  epsilon_summary <- do.call(rbind, lapply(seq_len(n_times), function(t) {
    nm <- epsilon_names[t]; sv <- samples_mat[, nm]; hpd <- safe_hpd(sv)
    tibble(Model = model_type, Time = t, Parameter = nm,
           Mean = mean(sv), SD = sd(sv),
           Lower = hpd[1], Upper = hpd[2],
           ESS  = as.numeric(effectiveSize(mcmc_list_full[, nm])),
           Rhat = safe_gelman(mcmc_list_full[, nm]))
  }))
  write_csv(epsilon_summary, file.path(scenario_dir, "epsilon_summary.csv"))
  
  # ── LAMBDA (regiões monitoradas) ─────────────────────────
  lambda_sel <- do.call(rbind, lapply(REGIONS_INTEREST, function(i) {
    do.call(rbind, lapply(1:n_times, function(t) {
      nm <- paste0("lambda[", i, ", ", t, "]")
      sv <- samples_mat[, nm]; hpd <- safe_hpd(sv)
      tibble(Model = model_type, Region = i, Nome = region_names[i], 
             Time = t, Mean = mean(sv), Lower = hpd[1], Upper = hpd[2],
             ESS = as.numeric(effectiveSize(mcmc_list_full[, nm])),
             Rhat = safe_gelman(mcmc_list_full[, nm]))
    }))
  }))
  write_csv(lambda_sel, file.path(scenario_dir, "lambda_selected.csv"))
  
  # ── S[i] (efeito espacial) ───────────────────────────────
  s_summary <- do.call(rbind, lapply(1:N_cluster, function(i) {
    nm <- s_names[i]; sv <- samples_mat[, nm]; hpd <- safe_hpd(sv)
    tibble(Model = model_type, Region = i, Nome = region_names[i],
           Mean = mean(sv), Lower = hpd[1], Upper = hpd[2],
           ESS = as.numeric(effectiveSize(mcmc_list_full[, nm])),
           Rhat = safe_gelman(mcmc_list_full[, nm]))
  }))
  write_csv(s_summary, file.path(scenario_dir, "s_summary.csv"))
  
  # ── THETA[i,t] e MU[i,t] ─────────────────────────────────
  beta_cols <- grep("^beta\\[", colnames(samples_mat), value = TRUE)
  
  theta_mu_list <- lapply(REGIONS_INTEREST, function(i) {
    do.call(rbind, lapply(1:n_times, function(t) {
      ldraws <- samples_mat[, paste0("lambda[", i, ", ", t, "]")]
      eps_t  <- samples_mat[, paste0("epsilon[", t, "]")]
      bdraws <- samples_mat[, beta_cols, drop = FALSE]
      
      # theta = lambda * exp(beta'x)
      theta <- ldraws * exp(as.vector(bdraws %*% as.numeric(x_data[i, t, ])))
      
      # mu = theta * E * epsilon
      mu <- theta * E_cluster[i, t] * eps_t
      
      theta_hpd <- safe_hpd(theta)
      mu_hpd <- safe_hpd(mu)
      
      # Criar duas linhas separadas e depois combinar
      theta_row <- data.frame(
        Region = i, Nome = region_names[i], Time = t, 
        Tipo = "theta", Mean = mean(theta), 
        Lower = theta_hpd[1], Upper = theta_hpd[2],
        Y_obs = NA_real_, stringsAsFactors = FALSE
      )
      
      mu_row <- data.frame(
        Region = i, Nome = region_names[i], Time = t,
        Tipo = "mu", Mean = mean(mu),
        Lower = mu_hpd[1], Upper = mu_hpd[2],
        Y_obs = Y_cluster[i, t], stringsAsFactors = FALSE
      )
      
      rbind(theta_row, mu_row)
    }))
  })
  
  theta_mu_df <- bind_rows(theta_mu_list)
  write_csv(theta_mu_df, file.path(scenario_dir, "theta_mu_summary.csv"))
  # ── WAIC e LPML ─────────────────────────────────────────
  loglik_names <- grep("logLik_Y", colnames(samples_mat), value = TRUE)
  lm_mat <- samples_mat[, loglik_names, drop = FALSE]
  lppd   <- sum(apply(lm_mat, 2, function(x) { mx <- max(x); mx + log(mean(exp(x - mx))) }))
  p_waic <- sum(apply(lm_mat, 2, var))
  waic   <- -2 * (lppd - p_waic)
  LPML   <- sum(log(1 / apply(lm_mat, 2, function(x) mean(exp(-x)))))
  
  criteria <- tibble(Model = model_type, WAIC = waic, LPML = LPML, 
                     lppd = lppd, pWAIC = p_waic)
  write_csv(criteria, file.path(scenario_dir, "criteria.csv"))
  
  cat(sprintf("[%s] WAIC = %.2f | LPML = %.2f\n", model_type, waic, LPML))
  
  # ── ACF DIAGNÓSTICOS ────────────────────────────────────
  params_diag <- c(beta_names, epsilon_names, "sigma_s")
  
  # Lambda ACF (apenas regiões monitoradas)
  lambda_monitor_names <- paste0("lambda[", 
                                 rep(REGIONS_INTEREST, each = n_times), ", ", 
                                 rep(1:n_times, length(REGIONS_INTEREST)), "]")
  
  all_params <- c(params_diag, lambda_monitor_names)
  all_params <- all_params[all_params %in% colnames(samples_mat)]
  
  acf_results <- do.call(rbind, lapply(all_params, function(nm) {
    ac   <- acf(samples_mat[, nm], lag.max = 200, plot = FALSE)
    lags <- as.vector(ac$lag[-1]); acfs <- as.vector(ac$acf[-1])
    ess_v  <- tryCatch(as.numeric(effectiveSize(mcmc_list_full[, nm])),
                       error = function(e) NA_real_)
    rhat_v <- tryCatch(safe_gelman(mcmc_list_full[, nm]),
                       error = function(e) NA_real_)
    tibble(Model = model_type, Parameter = nm, 
           ESS = ess_v, Rhat = rhat_v,
           lag_0.10 = { v <- lags[which(abs(acfs) < 0.10)[1]]; ifelse(is.na(v), Inf, v) },
           lag_0.05 = { v <- lags[which(abs(acfs) < 0.05)[1]]; ifelse(is.na(v), Inf, v) },
           acf_lag1 = acfs[1])
  }))
  write_csv(acf_results, file.path(scenario_dir, "acf_diagnostics.csv"))
  
  # =========================================================================
  # GRÁFICOS
  # =========================================================================
  
  make_label <- function(names) {
    sapply(names, function(n) {
      n_clean <- gsub("_", " ", n)
      if (nchar(n_clean) > 25) paste0(substr(n_clean, 1, 22), "...") else n_clean
    })
  }
  
  # Epsilon[t]
  ggsave(
    file.path(scenario_dir, "painel_epsilon.png"),
    ggplot(epsilon_summary, aes(x = Time)) +
      geom_ribbon(aes(ymin = Lower, ymax = Upper), fill = "darkorange", alpha = 0.3) +
      geom_line(aes(y = Mean), color = "darkorange", linewidth = 0.8) +
      geom_hline(yintercept = c(EPS_LOWER, EPS_UPPER), linetype = "dashed", color = "grey50") +
      scale_x_continuous(breaks = 1:n_times, labels = anos_label) +
      theme_bw(base_size = 12) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      labs(title = sprintf("epsilon[t] ~ U(%.2f, %.2f) - %s", EPS_LOWER, EPS_UPPER, model_type),
           x = "Ano", y = expression(epsilon[t])),
    width = 10, height = 5
  )
  
  # Beta
  ggsave(
    file.path(scenario_dir, "beta_posterior.png"),
    ggplot(beta_summary, aes(x = Covariate, y = Mean)) +
      geom_point(size = 3, color = "steelblue") +
      geom_errorbar(aes(ymin = HPD_Lower, ymax = HPD_Upper), width = 0.2, color = "steelblue") +
      geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
      theme_bw(base_size = 12) +
      labs(title = paste("Coeficientes beta -", model_type), y = expression(beta)),
    width = 8, height = 5
  )
  
  # Lambda[i,t] para regiões monitoradas
  lambda_sel$Label <- make_label(lambda_sel$Nome)
  
  ggsave(
    file.path(scenario_dir, "painel_lambda.png"),
    ggplot(lambda_sel, aes(x = Time)) +
      geom_ribbon(aes(ymin = Lower, ymax = Upper), fill = "grey70", alpha = 0.5) +
      geom_line(aes(y = Mean), color = "black", linewidth = 0.8) +
      facet_wrap(~ Label, scales = "free_y", ncol = 3) +
      scale_x_continuous(breaks = seq(1, n_times, by = 5), labels = anos_label[seq(1, n_times, by = 5)]) +
      theme_bw(base_size = 10) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      labs(title = sprintf("lambda[i,t] - %s", model_type), x = "Ano", y = expression(lambda[it])),
    width = 14, height = 8
  )
  
  # Mu[i,t]
  theta_mu_df$Label <- make_label(theta_mu_df$Nome)
  mu_df <- theta_mu_df %>% filter(Tipo == "mu")
  
  ggsave(
    file.path(scenario_dir, "painel_mu.png"),
    ggplot(mu_df, aes(x = Time)) +
      geom_ribbon(aes(ymin = Lower, ymax = Upper), fill = "steelblue", alpha = 0.3) +
      geom_line(aes(y = Mean), color = "steelblue", linewidth = 0.8) +
      geom_point(aes(y = Y_obs), shape = 21, fill = "white", color = "black", size = 2) +
      facet_wrap(~ Label, scales = "free_y", ncol = 3) +
      scale_x_continuous(breaks = seq(1, n_times, by = 5), labels = anos_label[seq(1, n_times, by = 5)]) +
      theme_bw(base_size = 10) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      labs(title = sprintf("mu[i,t] - %s", model_type), 
           subtitle = "Pontos = Y observado", x = "Ano", y = expression(mu[it])),
    width = 14, height = 8
  )
  
  # Theta[i,t]
  theta_df <- theta_mu_df %>% filter(Tipo == "theta")
  
  ggsave(
    file.path(scenario_dir, "painel_theta.png"),
    ggplot(theta_df, aes(x = Time)) +
      geom_ribbon(aes(ymin = Lower, ymax = Upper), fill = "darkgreen", alpha = 0.3) +
      geom_line(aes(y = Mean), color = "darkgreen", linewidth = 0.8) +
      facet_wrap(~ Label, scales = "free_y", ncol = 3) +
      scale_x_continuous(breaks = seq(1, n_times, by = 5), labels = anos_label[seq(1, n_times, by = 5)]) +
      theme_bw(base_size = 10) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      labs(title = sprintf("theta[i,t] (risco relativo) - %s", model_type), 
           x = "Ano", y = expression(theta[it])),
    width = 14, height = 8
  )
  
  # Efeito espacial s[i]
  ref_nome <- region_names[ref_idx]
  s_summary$Highlight <- ifelse(s_summary$Nome == ref_nome, "Ref", "Outras")
  
  ggsave(
    file.path(scenario_dir, "s_posterior.png"),
    ggplot(s_summary, aes(x = reorder(Nome, Mean), y = Mean, color = Highlight)) +
      geom_point(size = 1.5) +
      geom_errorbar(aes(ymin = Lower, ymax = Upper), width = 0.4, linewidth = 0.3) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
      scale_color_manual(values = c("Ref" = "red", "Outras" = "black")) +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6),
            legend.position = "none") +
      labs(title = sprintf("Efeito espacial s[i] - %s", model_type),
           subtitle = paste("Vermelho:", ref_nome), y = "s[i]", x = ""),
    width = max(10, N_cluster * 0.3), height = 6
  )
  
  # ACF de beta
  acf_beta_df <- do.call(rbind, lapply(beta_names, function(nm) {
    ac <- acf(samples_mat[, nm], lag.max = 100, plot = FALSE)
    tibble(Parameter = nm, Lag = as.vector(ac$lag[-1]), ACF = as.vector(ac$acf[-1]))
  }))
  
  ggsave(
    file.path(scenario_dir, "acf_beta.png"),
    ggplot(acf_beta_df, aes(x = Lag, y = ACF)) +
      geom_col(width = 0.6, fill = "grey50") +
      geom_hline(yintercept = c(-0.10, 0.10), linetype = "dashed", color = "blue", linewidth = 0.5) +
      geom_hline(yintercept = c(-0.05, 0.05), linetype = "dotted", color = "red", linewidth = 0.5) +
      facet_wrap(~ Parameter, scales = "free_y") +
      theme_bw(base_size = 11) +
      labs(title = paste("ACF beta -", model_type)),
    width = 10, height = 5
  )
  
  # ACF de epsilon
  acf_eps_df <- do.call(rbind, lapply(epsilon_names, function(nm) {
    ac <- acf(samples_mat[, nm], lag.max = 100, plot = FALSE)
    tibble(Time = as.integer(str_extract(nm, "\\d+")),
           Lag  = as.vector(ac$lag[-1]), ACF = as.vector(ac$acf[-1]))
  }))
  
  ggsave(
    file.path(scenario_dir, "acf_epsilon.png"),
    ggplot(acf_eps_df, aes(x = Lag, y = ACF)) +
      geom_col(width = 0.6, fill = "darkorange") +
      geom_hline(yintercept = c(-0.10, 0.10), linetype = "dashed", color = "blue", linewidth = 0.5) +
      geom_hline(yintercept = c(-0.05, 0.05), linetype = "dotted", color = "red", linewidth = 0.5) +
      facet_wrap(~ Time, scales = "free_y", ncol = 6) +
      theme_bw(base_size = 9) +
      labs(title = paste("ACF epsilon[t] -", model_type)),
    width = 14, height = 10
  )
  
  # ACF de lambda (apenas regiões monitoradas)
  acf_lambda_df <- do.call(rbind, lapply(REGIONS_INTEREST, function(i) {
    do.call(rbind, lapply(1:n_times, function(t) {
      nm <- paste0("lambda[", i, ", ", t, "]")
      if (!nm %in% colnames(samples_mat)) return(NULL)
      ac <- acf(samples_mat[, nm], lag.max = 100, plot = FALSE)
      tibble(Region = i, Nome = region_names[i], Time = t,
             Lag  = as.vector(ac$lag[-1]), ACF = as.vector(ac$acf[-1]))
    }))
  }))
  
  if (nrow(acf_lambda_df) > 0) {
    acf_lambda_df$Label <- make_label(acf_lambda_df$Nome)
    
    ggsave(
      file.path(scenario_dir, "acf_lambda.png"),
      ggplot(acf_lambda_df, aes(x = Lag, y = ACF)) +
        geom_col(width = 0.6, fill = "steelblue") +
        geom_hline(yintercept = c(-0.10, 0.10), linetype = "dashed", color = "blue", linewidth = 0.5) +
        geom_hline(yintercept = c(-0.05, 0.05), linetype = "dotted", color = "red", linewidth = 0.5) +
        facet_grid(Label ~ Time, scales = "free_y") +
        theme_bw(base_size = 7) +
        theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 5),
              strip.text.y = element_text(size = 6)) +
        labs(title = paste("ACF lambda[i,t] -", model_type)),
      width = 20, height = 3 * length(REGIONS_INTEREST)
    )
  }
  
  # Traceplots de beta
  cores_cadeia <- c("Cadeia 1" = "#2166AC", "Cadeia 2" = "#D6604D")
  df_trace <- do.call(rbind, lapply(seq_len(nchains), function(ch) {
    cm <- as.matrix(mcmc_list_full[[ch]])
    do.call(rbind, lapply(beta_names, function(nm) {
      vals <- cm[, nm]
      tibble(Iter = seq_along(vals), Value = vals,
             ErgMedia = cumsum(vals) / seq_along(vals),
             Parameter = nm, Cadeia = paste0("Cadeia ", ch))
    }))
  }))
  
  ggsave(
    file.path(scenario_dir, "traceplots_beta.png"),
    ggplot(df_trace, aes(x = Iter, color = Cadeia)) +
      geom_line(aes(y = Value), alpha = 0.25, linewidth = 0.20) +
      geom_line(aes(y = ErgMedia), alpha = 0.90, linewidth = 0.75) +
      scale_color_manual(values = cores_cadeia) +
      facet_wrap(~ Parameter, scales = "free_y") +
      theme_bw(base_size = 11) + theme(legend.position = "bottom") +
      labs(title = paste("Traceplots beta -", model_type), x = "Iteração", y = "Valor"),
    width = 10, height = 5
  )
  
  # Retorno resumido
  tibble(
    Model = model_type, p = p, WAIC = waic, LPML = LPML,
    ESS_beta_min    = min(beta_summary$ESS, na.rm = TRUE),
    ESS_epsilon_min = min(epsilon_summary$ESS, na.rm = TRUE),
    ESS_lambda_min  = min(lambda_sel$ESS, na.rm = TRUE),
    Rhat_beta_max    = max(beta_summary$Rhat, na.rm = TRUE),
    Rhat_epsilon_max = max(epsilon_summary$Rhat, na.rm = TRUE),
    Rhat_lambda_max  = max(lambda_sel$Rhat, na.rm = TRUE),
    lag_max_0.10 = max(acf_results$lag_0.10, na.rm = TRUE),
    lag_max_0.05 = max(acf_results$lag_0.05, na.rm = TRUE),
    tempo_min = as.numeric(difftime(t_fim, t_inicio, units = "mins"))
  )
}

# ==============================================================================
# 6. EXECUÇÃO
# ==============================================================================

output_dir <- file.path("C:/Users/vlara/OneDrive/Estatistica UFMG/Mestrado/Pesquisa/Aplicação",
                        paste0("resultados_", CLUSTER_NOME, "_COMPLETO_vs_SIMPLIFICADO"))
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

resultados <- list()
resultados[["COMPLETO"]] <- run_model("COMPLETO", output_dir)
resultados[["SIMPLIFICADO"]] <- run_model("SIMPLIFICADO", output_dir)

# ==============================================================================
# 7. CONSOLIDAÇÃO FINAL
# ==============================================================================

resumo <- bind_rows(resultados)
write_csv(resumo, file.path(output_dir, "resumo_comparativo.csv"))

cat("\n")
cat(paste(rep("=", 80), collapse = ""), "\n")
cat(sprintf("RESUMO COMPARATIVO - %s\n", CLUSTER_NOME))
cat(paste(rep("=", 80), collapse = ""), "\n\n")
print(resumo)

cat(sprintf("\nΔWAIC = %.2f\n", resumo$WAIC[1] - resumo$WAIC[2]))
cat("\n⏱️ Tempo total:", format(Sys.time() - inicio_global), "\n")
