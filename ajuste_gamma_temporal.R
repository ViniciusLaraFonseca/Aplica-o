# ==============================================================================
# ajuste_gamma_temporal.R
# Modelo Bayesiano Espaço-Estados — dados reais MG
# A = 75 microrregiões | T = 23 anos | K = 4 clusters | p = 3 covariáveis
#
# EXTENSÃO em relação a ajuste_lambda_shared.R:
#   gamma: VETOR K  →  MATRIZ K × T   (gamma[k, t])
#   epsilon: VETOR N  →  MATRIZ N × T  (epsilon[i, t])
#
#     epsilon[i, t] = 1 - sum_k h[i, k] * gamma[k, t]
#
#   Preditor:
#     log(mu[i,t]) = log(lambda[t]) + log(E[i,t]) + log(epsilon[i,t])
#                    + beta' x[i,t] + s[i]   (s[i] apenas no modelo espacial)
#
#   Priori de gamma por coluna-t (restrição de soma independente em cada t):
#     gamma[1, t] ~ dunif(a_unif, b_unif)
#     gamma[j, t] ~ dunif(0, 1 - sum(gamma[1:(j-1), t]))   j = 2..K
#
#   Amostrador: AF_slice separado por coluna-t de gamma (23 blocos de tamanho K).
#
#   FFBS: g[i,t] usa epsilon[i, t] (indexado em t).
# ==============================================================================

inicio_global <- Sys.time()

setwd("C:/Users/vlara/OneDrive/Estatistica UFMG/Mestrado/Pesquisa/PesquisaMestrado")

pkgs <- c("nimble", "coda", "parallel", "dplyr", "ggplot2",
          "tidyr", "readr", "stringr", "tibble")
for (pkg in pkgs) {
  if (!require(pkg, character.only = TRUE)) install.packages(pkg)
  library(pkg, character.only = TRUE)
}

Sys.setenv(OMP_NUM_THREADS = "1")
Sys.setenv(MKL_NUM_THREADS = "1")
if (requireNamespace("RhpcBLASctl", quietly = TRUE))
  RhpcBLASctl::blas_set_num_threads(1)

# ==============================================================================
# 1. DADOS
# ==============================================================================

# carrega Y_mat, E, x
source("_dataCaseStudy.r")   # carrega data$adj, data$num, data$sumNumNeigh, data$hAI

stopifnot(
  "Y_mat deve ter 75 linhas"       = nrow(Y_mat) == 75,
  "E deve ter 75 linhas"           = nrow(E)     == 75,
  "x deve ter dimensão correta"    = all(dim(x) == c(75, 23, 3)),
  "Y_mat e E devem ter mesma dim"  = identical(dim(Y_mat), dim(E)),
  "linhas de Y_mat e E alinhadas"  = identical(row.names(Y_mat), row.names(E)),
  "linhas de Y_mat e x alinhadas"  = identical(row.names(Y_mat), dimnames(x)[[1]]),
  "estrutura espacial consistente" = sum(data$num) == data$sumNumNeigh
)
cat("Verificacoes de consistencia: OK\n")

E_norm    <- E / mean(E)
N_regions <- nrow(Y_mat)
n_times   <- ncol(Y_mat)
p         <- dim(x)[3]
K         <- ncol(data$hAI)
cat(sprintf("N = %d | T = %d | p = %d | K = %d\n", N_regions, n_times, p, K))

adj_vec     <- as.integer(data$adj)
num_vec     <- as.integer(data$num)
n_adj_val   <- as.integer(data$sumNumNeigh)
weights_vec <- rep(1.0, n_adj_val)
h_mat       <- data$hAI   # 75 × 4

# ==============================================================================
# 2. REGIÕES DE INTERESSE
# ==============================================================================
cluster_ids <- apply(h_mat, 1, sum)

set.seed(42)
REGIONS_INTEREST <- unlist(lapply(1:K, function(cl) {
  regs <- which(cluster_ids == cl)
  if (length(regs) >= 3) sample(regs, 3) else regs
}))
cat("Regioes selecionadas:\n"); print(sort(REGIONS_INTEREST))

# ==============================================================================
# 3. CONSTANTES E DADOS NIMBLE
# ==============================================================================
constants_spatial <- list(
  n_regions = N_regions, n_times = n_times, p = p, K = K,
  h         = h_mat,
  mu_beta   = rep(0, p),
  a_unif    = 0.0, b_unif = 0.1,
  a0 = 1.0, b0 = 1.0, w = 0.9,
  adj     = adj_vec,    num     = num_vec,
  weights = weights_vec, n_adj  = n_adj_val
)
constants_nonspatial <- constants_spatial[
  setdiff(names(constants_spatial), c("adj", "num", "weights", "n_adj"))
]
data_nimble <- list(Y = Y_mat, E = E_norm, x = x)

# ==============================================================================
# 4. INICIALIZAÇÕES
# ==============================================================================
set.seed(123)
gamma_init1 <- matrix(c(0.05, 0.10, 0.10, 0.15), nrow = K, ncol = n_times)
gamma_init2 <- matrix(c(0.04, 0.09, 0.09, 0.14), nrow = K, ncol = n_times)

inits_list_spatial <- list(
  list(beta = rep(0, p),        gamma = gamma_init1,
       lambda = rep(1.0, n_times), sigma_s = 0.5, s = rep(0, N_regions)),
  list(beta = rnorm(p, 0, 0.3), gamma = gamma_init2,
       lambda = rgamma(n_times, 1, 1), sigma_s = 1.0, s = rep(0, N_regions))
)
inits_list_nonspatial <- list(
  list(beta = rep(0, p),        gamma = gamma_init1, lambda = rep(1.0, n_times)),
  list(beta = rnorm(p, 0, 0.3), gamma = gamma_init2, lambda = rgamma(n_times, 1, 1))
)

# ==============================================================================
# 5. FUNÇÃO WORKER
# ==============================================================================
run_model <- function(model_type, output_dir) {
  
  library(nimble); library(coda); library(dplyr)
  library(ggplot2); library(readr); library(stringr); library(tibble)
  
  # ── 5a. Código NIMBLE ────────────────────────────────────────────────────────
  code_spatial <- nimbleCode({
    for (j in 1:p) beta[j] ~ dnorm(mu_beta[j], sd = 1)
    
    for (t in 1:n_times) {
      gamma[1, t] ~ dunif(min = a_unif, max = b_unif)
      for (j in 2:K)
        gamma[j, t] ~ dunif(min = 0, max = (1 - sum(gamma[1:(j-1), t])))
    }
    
    sigma_s ~ T(dt(0, 1, 1), 0, )
    tau_s   <- 1 / (sigma_s^2)
    s[1:n_regions] ~ dcar_normal(
      adj[1:n_adj], weights[1:n_adj], num[1:n_regions], tau_s, zero_mean = 1
    )
    
    for (t in 1:n_times) lambda[t] ~ dgamma(1, 1)
    
    for (i in 1:n_regions) {
      for (t in 1:n_times) {
        epsilon[i, t] <- 1 - inprod(h[i, 1:K], gamma[1:K, t])
        log(mu[i, t]) <- log(lambda[t]) + log(E[i, t]) + log(epsilon[i, t]) +
          inprod(beta[1:p], x[i, t, 1:p]) + s[i]
        Y[i, t]        ~ dpois(mu[i, t])
        logLik_Y[i, t] <- dpois(Y[i, t], mu[i, t], log = TRUE)
      }
    }
  })
  
  code_nonspatial <- nimbleCode({
    for (j in 1:p) beta[j] ~ dnorm(mu_beta[j], sd = 1)
    
    for (t in 1:n_times) {
      gamma[1, t] ~ dunif(min = a_unif, max = b_unif)
      for (j in 2:K)
        gamma[j, t] ~ dunif(min = 0, max = (1 - sum(gamma[1:(j-1), t])))
    }
    
    for (t in 1:n_times) lambda[t] ~ dgamma(1, 1)
    
    for (i in 1:n_regions) {
      for (t in 1:n_times) {
        epsilon[i, t] <- 1 - inprod(h[i, 1:K], gamma[1:K, t])
        log(mu[i, t]) <- log(lambda[t]) + log(E[i, t]) + log(epsilon[i, t]) +
          inprod(beta[1:p], x[i, t, 1:p])
        Y[i, t]        ~ dpois(mu[i, t])
        logLik_Y[i, t] <- dpois(Y[i, t], mu[i, t], log = TRUE)
      }
    }
  })
  
  # ── 5b. FFBS espacial ────────────────────────────────────────────────────────
  ffbs_spatial <- nimbleFunction(
    contains = sampler_BASE,
    setup = function(model, mvSaved, target, control) {
      n_regions <- control$n_regions; n_times <- control$n_times
      p <- control$p; a0 <- control$a0; b0 <- control$b0; w <- control$w
      at_buf <- nimNumeric(n_times + 1, 0)
      bt_buf <- nimNumeric(n_times + 1, 0)
      calcNodes   <- model$getDependencies(target, self = FALSE)
      targetNodes <- model$expandNodeNames(target)
      setupOutputs(at_buf, bt_buf)
    },
    run = function() {
      declare(i, integer()); declare(t, integer()); declare(t_idx, integer())
      declare(t_back, integer()); declare(k, integer())
      declare(prod_val, double()); declare(g_it, double())
      declare(sum_Y_t, double()); declare(sum_g_t, double()); declare(nu, double())
      
      at_buf[1] <<- a0; bt_buf[1] <<- b0
      
      for (t in 1:n_times) {
        sum_Y_t <- 0; sum_g_t <- 0
        for (i in 1:n_regions) {
          sum_Y_t  <- sum_Y_t + model$Y[i, t]
          prod_val <- 0
          for (k in 1:p) prod_val <- prod_val + model$x[i, t, k] * model$beta[k]
          # epsilon[i, t] — indexado em t
          g_it    <- model$E[i, t] * model$epsilon[i, t] * exp(prod_val + model$s[i])
          sum_g_t <- sum_g_t + g_it
        }
        at_buf[t + 1] <<- w * at_buf[t] + sum_Y_t
        bt_buf[t + 1] <<- w * bt_buf[t] + sum_g_t
      }
      
      model$lambda[n_times] <<- rgamma(1, shape = at_buf[n_times + 1],
                                       rate  = bt_buf[n_times + 1])
      for (t_idx in 1:(n_times - 1)) {
        t_back <- n_times - t_idx
        nu <- rgamma(1, shape = (1 - w) * at_buf[t_back + 1],
                     rate  = bt_buf[t_back + 1])
        model$lambda[t_back] <<- nu + w * model$lambda[t_back + 1]
      }
      model$calculate(calcNodes)
      copy(from = model, to = mvSaved, row = 1, nodes = targetNodes, logProb = TRUE)
    },
    methods = list(reset = function() {})
  )
  
  # ── 5c. FFBS não-espacial ────────────────────────────────────────────────────
  ffbs_nonspatial <- nimbleFunction(
    contains = sampler_BASE,
    setup = function(model, mvSaved, target, control) {
      n_regions <- control$n_regions; n_times <- control$n_times
      p <- control$p; a0 <- control$a0; b0 <- control$b0; w <- control$w
      at_buf <- nimNumeric(n_times + 1, 0)
      bt_buf <- nimNumeric(n_times + 1, 0)
      calcNodes   <- model$getDependencies(target, self = FALSE)
      targetNodes <- model$expandNodeNames(target)
      setupOutputs(at_buf, bt_buf)
    },
    run = function() {
      declare(i, integer()); declare(t, integer()); declare(t_idx, integer())
      declare(t_back, integer()); declare(k, integer())
      declare(prod_val, double()); declare(g_it, double())
      declare(sum_Y_t, double()); declare(sum_g_t, double()); declare(nu, double())
      
      at_buf[1] <<- a0; bt_buf[1] <<- b0
      
      for (t in 1:n_times) {
        sum_Y_t <- 0; sum_g_t <- 0
        for (i in 1:n_regions) {
          sum_Y_t  <- sum_Y_t + model$Y[i, t]
          prod_val <- 0
          for (k in 1:p) prod_val <- prod_val + model$x[i, t, k] * model$beta[k]
          g_it    <- model$E[i, t] * model$epsilon[i, t] * exp(prod_val)
          sum_g_t <- sum_g_t + g_it
        }
        at_buf[t + 1] <<- w * at_buf[t] + sum_Y_t
        bt_buf[t + 1] <<- w * bt_buf[t] + sum_g_t
      }
      
      model$lambda[n_times] <<- rgamma(1, shape = at_buf[n_times + 1],
                                       rate  = bt_buf[n_times + 1])
      for (t_idx in 1:(n_times - 1)) {
        t_back <- n_times - t_idx
        nu <- rgamma(1, shape = (1 - w) * at_buf[t_back + 1],
                     rate  = bt_buf[t_back + 1])
        model$lambda[t_back] <<- nu + w * model$lambda[t_back + 1]
      }
      model$calculate(calcNodes)
      copy(from = model, to = mvSaved, row = 1, nodes = targetNodes, logProb = TRUE)
    },
    methods = list(reset = function() {})
  )
  
  # ── 5d. Selecionar objetos ───────────────────────────────────────────────────
  is_spatial <- (model_type == "spatial")
  model_code <- if (is_spatial) code_spatial      else code_nonspatial
  constants  <- if (is_spatial) constants_spatial  else constants_nonspatial
  inits_list <- if (is_spatial) inits_list_spatial else inits_list_nonspatial
  ffbs_fn    <- if (is_spatial) ffbs_spatial       else ffbs_nonspatial
  
  cat("\n=== Iniciando modelo:", model_type, "===\n")
  scenario_dir <- file.path(output_dir, model_type)
  dir.create(scenario_dir, recursive = TRUE, showWarnings = FALSE)
  
  # ── 5e. Construir e compilar ─────────────────────────────────────────────────
  model  <- nimbleModel(code = model_code, constants = constants,
                        data = data_nimble, inits = inits_list[[1]], check = FALSE)
  Cmodel <- compileNimble(model)
  
  conf <- configureMCMC(model)
  
  conf$removeSamplers("lambda")
  conf$addSampler(
    target  = "lambda", type = ffbs_fn,
    control = list(n_regions = N_regions, n_times = n_times, p = p,
                   a0 = constants$a0, b0 = constants$b0, w = constants$w)
  )
  
  # Um AF_slice por coluna-t de gamma (23 blocos de K parâmetros)
  conf$removeSamplers("gamma")
  for (t_idx in seq_len(n_times)) {
    conf$addSampler(
      target = paste0("gamma[", seq_len(K), ", ", t_idx, "]"),
      type   = "AF_slice"
    )
  }
  
  monitors_base <- c("beta", "gamma", "logLik_Y", "lambda")
  if (is_spatial) monitors_base <- c(monitors_base, "s", "sigma_s", "tau_s")
  conf$addMonitors(monitors_base)
  conf$printSamplers()
  
  Rmcmc <- buildMCMC(conf)
  Cmcmc <- compileNimble(Rmcmc, project = model)
  
  # ── 5f. MCMC ─────────────────────────────────────────────────────────────────
  niter <- 50000; nburnin <- 10000; nchains <- 2; thin <- 10
  
  cat(sprintf("[%s] niter=%d | nburnin=%d | thin=%d | cadeias=%d\n",
              model_type, niter, nburnin, thin, nchains))
  
  samples <- runMCMC(
    Cmcmc, niter = niter, nburnin = nburnin, nchains = nchains, thin = thin,
    inits = inits_list, samplesAsCodaMCMC = TRUE, summary = FALSE, WAIC = FALSE
  )
  saveRDS(samples, file.path(scenario_dir, "samples.rds"))
  cat("[", model_type, "] Amostras salvas.\n")
  
  samples_mat    <- as.matrix(samples)
  mcmc_list_full <- mcmc.list(lapply(seq_len(nchains),
                                     function(ch) as.mcmc(samples[[ch]])))
  rm(samples); gc()
  
  # ── 5g. Funções auxiliares ───────────────────────────────────────────────────
  safe_hpd <- function(sv) {
    if (var(sv) < 1e-12) return(c(NA_real_, NA_real_))
    as.numeric(HPDinterval(as.mcmc(sv), prob = 0.95))
  }
  safe_gelman <- function(obj) {
    tryCatch(gelman.diag(obj, autoburnin = FALSE)$psrf[, 1],
             error = function(e) rep(NA_real_, nvar(obj)))
  }
  
  beta_names   <- paste0("beta[",   seq_len(p),       "]")
  lambda_names <- paste0("lambda[", seq_len(n_times), "]")
  
  # Nomes "gamma[k, t]" na matriz de amostras
  gamma_names_mat <- outer(
    seq_len(K), seq_len(n_times),
    function(k, t) paste0("gamma[", k, ", ", t, "]")
  )  # K × n_times
  
  # ── 5h. Epsilon posterior — array (n_draw × N × T) ──────────────────────────
  n_draw        <- nrow(samples_mat)
  epsilon_draws <- array(NA_real_, dim = c(n_draw, N_regions, n_times))
  for (t in seq_len(n_times)) {
    g_t <- samples_mat[, gamma_names_mat[, t], drop = FALSE]  # n_draw × K
    epsilon_draws[, , t] <- 1 - g_t %*% t(h_mat)              # n_draw × N
  }
  
  # ── 5i. Sumário de gamma[k, t] ───────────────────────────────────────────────
  gamma_summary <- do.call(rbind, lapply(seq_len(K), function(k) {
    do.call(rbind, lapply(seq_len(n_times), function(t) {
      nm  <- gamma_names_mat[k, t]
      sv  <- samples_mat[, nm]; hpd <- safe_hpd(sv)
      tibble(k = k, Time = t, Parameter = nm,
             Mean = mean(sv), SD = sd(sv), Lower = hpd[1], Upper = hpd[2],
             ESS  = as.numeric(effectiveSize(mcmc_list_full[, nm])),
             Rhat = safe_gelman(mcmc_list_full[, nm]))
    }))
  }))
  write_csv(gamma_summary, file.path(scenario_dir, "gamma_summary.csv"))
  
  # Painel temporal de gamma[k, t]: facet por cluster k, eixo x = Tempo
  # (mesmo estilo de painel_lambda: ribbon + linha)
  ggsave(
    file.path(scenario_dir, "painel_gamma.png"),
    ggplot(gamma_summary, aes(x = Time)) +
      geom_ribbon(aes(ymin = Lower, ymax = Upper), fill = "steelblue", alpha = 0.3) +
      geom_line(aes(y = Mean), color = "steelblue", linewidth = 0.8) +
      facet_wrap(
        ~ k, ncol = 2, scales = "free_y",
        labeller = labeller(k = function(x) paste0("gamma[", x, ", t]"))
      ) +
      theme_bw(base_size = 11) +
      labs(title = paste("Trajetoria de gamma[k, t] (", model_type, ")"),
           x = "Tempo", y = expression(gamma[kt])),
    width = 10, height = 7
  )
  
  # ── 5j. Sumário de epsilon[i, t] por região e por cluster ────────────────────
  
  # Sumário completo (N × T) — salvo em CSV
  eps_full <- do.call(rbind, lapply(seq_len(N_regions), function(i) {
    do.call(rbind, lapply(seq_len(n_times), function(t) {
      sv  <- epsilon_draws[, i, t]; hpd <- safe_hpd(sv)
      tibble(Region = i, Cluster = cluster_ids[i], Time = t,
             Mean = mean(sv), Lower = hpd[1], Upper = hpd[2])
    }))
  }))
  write_csv(eps_full, file.path(scenario_dir, "epsilon_summary.csv"))
  
  # Agrega por cluster × tempo (pool de draws das regiões do cluster)
  eps_by_cluster <- do.call(rbind, lapply(seq_len(K), function(k) {
    regs_k <- which(cluster_ids == k)
    do.call(rbind, lapply(seq_len(n_times), function(t) {
      sv  <- as.vector(epsilon_draws[, regs_k, t])
      hpd <- safe_hpd(sv)
      tibble(k = k, Time = t, Mean = mean(sv), Lower = hpd[1], Upper = hpd[2])
    }))
  }))
  write_csv(eps_by_cluster, file.path(scenario_dir, "epsilon_cluster_summary.csv"))
  
  # Painel temporal de epsilon agregado por cluster k
  ggsave(
    file.path(scenario_dir, "painel_epsilon.png"),
    ggplot(eps_by_cluster, aes(x = Time)) +
      geom_ribbon(aes(ymin = Lower, ymax = Upper), fill = "darkorange", alpha = 0.3) +
      geom_line(aes(y = Mean), color = "darkorange", linewidth = 0.8) +
      facet_wrap(
        ~ k, ncol = 2, scales = "free_y",
        labeller = labeller(k = function(x) paste0("epsilon[", x, ", t]"))
      ) +
      theme_bw(base_size = 11) +
      labs(title = paste("Trajetoria de epsilon por cluster (", model_type, ") — K=4"),
           x = "Tempo", y = expression(epsilon[kt])),
    width = 10, height = 7
  )
  # ── 5k. Sumário de beta ───────────────────────────────────────────────────────
  beta_summary <- do.call(rbind, lapply(beta_names, function(nm) {
    sv  <- samples_mat[, nm]; hpd <- safe_hpd(sv)
    tibble(Parameter = nm, Mean = mean(sv), SD = sd(sv),
           HPD_Lower = hpd[1], HPD_Upper = hpd[2],
           ESS  = as.numeric(effectiveSize(mcmc_list_full[, nm])),
           Rhat = safe_gelman(mcmc_list_full[, nm]))
  }))
  write_csv(beta_summary, file.path(scenario_dir, "beta_summary.csv"))
  
  # ── 5l. tau_s e s[i] (apenas espacial) ───────────────────────────────────────
  ESS_tau <- NA_real_
  if (is_spatial) {
    tau_sv <- samples_mat[, "tau_s"]; hpd_t <- safe_hpd(tau_sv)
    tau_sum <- tibble(
      Parameter = "tau_s", Mean = mean(tau_sv), SD = sd(tau_sv),
      HPD_Lower = hpd_t[1], HPD_Upper = hpd_t[2],
      ESS  = as.numeric(effectiveSize(mcmc_list_full[, "tau_s"])),
      Rhat = safe_gelman(mcmc_list_full[, "tau_s"])
    )
    write_csv(tau_sum, file.path(scenario_dir, "tau_summary.csv"))
    ESS_tau <- tau_sum$ESS
    
    s_names   <- paste0("s[", seq_len(N_regions), "]")
    s_summary <- do.call(rbind, lapply(seq_len(N_regions), function(i) {
      sv  <- samples_mat[, s_names[i]]; hpd <- safe_hpd(sv)
      tibble(Region = i, Mean = mean(sv), SD = sd(sv),
             HPD_Lower = hpd[1], HPD_Upper = hpd[2],
             ESS = as.numeric(effectiveSize(mcmc_list_full[, s_names[i]])))
    }))
    write_csv(s_summary, file.path(scenario_dir, "s_summary.csv"))
    ggsave(
      file.path(scenario_dir, "s_posterior.png"),
      ggplot(s_summary, aes(x = Region, y = Mean)) +
        geom_point(size = 0.9) +
        geom_errorbar(aes(ymin = HPD_Lower, ymax = HPD_Upper),
                      width = 0.4, linewidth = 0.3) +
        geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
        theme_bw() +
        labs(title = paste("Efeito espacial s[i] (", model_type, ")"),
             y = "s[i]", x = "Regiao"),
      width = 10, height = 5
    )
  }
  
  # ── 5m. lambda[t] — painel temporal ──────────────────────────────────────────
  lambda_summary <- do.call(rbind, lapply(seq_len(n_times), function(t) {
    nm <- lambda_names[t]; sv <- samples_mat[, nm]; hpd <- safe_hpd(sv)
    tibble(Time = t, Mean = mean(sv), SD = sd(sv), Lower = hpd[1], Upper = hpd[2],
           ESS  = as.numeric(effectiveSize(mcmc_list_full[, nm])),
           Rhat = safe_gelman(mcmc_list_full[, nm]), model = model_type)
  }))
  write_csv(lambda_summary, file.path(scenario_dir, "lambda_summary.csv"))
  
  ggsave(
    file.path(scenario_dir, "painel_lambda.png"),
    ggplot(lambda_summary, aes(x = Time)) +
      geom_ribbon(aes(ymin = Lower, ymax = Upper), fill = "grey70", alpha = 0.5) +
      geom_line(aes(y = Mean), color = "black", linewidth = 0.8) +
      theme_bw(base_size = 12) +
      labs(title = paste("Lambda compartilhado lambda[t] (", model_type, ")"),
           x = "Tempo", y = expression(lambda[t])),
    width = 9, height = 5
  )
  
  # ── 5n. theta[i, t] e mu[i, t] ───────────────────────────────────────────────
  beta_cols <- grep("^beta\\[", colnames(samples_mat), value = TRUE)
  
  # Média de epsilon por região (sobre tempos) — para rótulo dos facets
  eps_mean_region <- eps_full %>%
    group_by(Region) %>%
    summarise(MeanEps = mean(Mean), .groups = "drop")
  make_label <- function(regs) {
    setNames(
      sprintf("Reg %d (C%d)\ne=%.3f", regs, cluster_ids[regs],
              eps_mean_region$MeanEps[match(regs, eps_mean_region$Region)]),
      as.character(regs)
    )
  }
  
  theta_summary <- do.call(rbind, lapply(REGIONS_INTEREST, function(i) {
    do.call(rbind, lapply(seq_len(n_times), function(t) {
      ldraws <- samples_mat[, lambda_names[t]]
      bdraws <- samples_mat[, beta_cols, drop = FALSE]
      theta  <- ldraws * exp(as.vector(bdraws %*% data_nimble$x[i, t, ]))
      hpd    <- safe_hpd(theta)
      tibble(Region = i, Time = t, Mean = mean(theta),
             Lower = hpd[1], Upper = hpd[2], model = model_type)
    }))
  }))
  write_csv(theta_summary, file.path(scenario_dir, "theta_selected.csv"))
  
  theta_regs <- sort(unique(theta_summary$Region))
  ggsave(
    file.path(scenario_dir, "painel_theta.png"),
    ggplot(theta_summary, aes(x = Time)) +
      geom_ribbon(aes(ymin = Lower, ymax = Upper), fill = "steelblue", alpha = 0.3) +
      geom_line(aes(y = Mean), color = "steelblue", linewidth = 0.8) +
      facet_wrap(~ Region, scales = "free_y", ncol = 3,
                 labeller = labeller(Region = make_label(theta_regs))) +
      theme_bw(base_size = 10) +
      labs(title = paste("Theta estimado (", model_type, ")"),
           x = "Tempo", y = expression(theta[it])),
    width = 12, height = 10
  )
  
  mu_summary <- do.call(rbind, lapply(REGIONS_INTEREST, function(i) {
    do.call(rbind, lapply(seq_len(n_times), function(t) {
      ldraws   <- samples_mat[, lambda_names[t]]
      bdraws   <- samples_mat[, beta_cols, drop = FALSE]
      eps_it   <- epsilon_draws[, i, t]          # epsilon[i, t]
      mu_draws <- ldraws * exp(as.vector(bdraws %*% data_nimble$x[i, t, ])) *
        data_nimble$E[i, t] * eps_it
      hpd <- safe_hpd(mu_draws)
      tibble(Region = i, Time = t, Mean = mean(mu_draws),
             Lower = hpd[1], Upper = hpd[2], model = model_type)
    }))
  }))
  write_csv(mu_summary, file.path(scenario_dir, "mu_selected.csv"))
  
  mu_regs <- sort(unique(mu_summary$Region))
  ggsave(
    file.path(scenario_dir, "painel_mu.png"),
    ggplot(mu_summary, aes(x = Time)) +
      geom_ribbon(aes(ymin = Lower, ymax = Upper), fill = "darkorange", alpha = 0.3) +
      geom_line(aes(y = Mean), color = "darkorange", linewidth = 0.8) +
      facet_wrap(~ Region, scales = "free_y", ncol = 3,
                 labeller = labeller(Region = make_label(mu_regs))) +
      theme_bw(base_size = 10) +
      labs(title = paste("Mu estimado (", model_type, ")"),
           x = "Tempo", y = expression(mu[it])),
    width = 12, height = 10
  )
  
  # ── 5o. WAIC e LPML ───────────────────────────────────────────────────────────
  loglik_names <- grep("logLik_Y", colnames(samples_mat), value = TRUE)
  waic <- NA_real_; LPML <- NA_real_
  if (length(loglik_names) > 0) {
    lm     <- samples_mat[, loglik_names, drop = FALSE]
    lppd   <- sum(apply(lm, 2, function(x) { mx <- max(x); mx + log(mean(exp(x - mx))) }))
    p_waic <- sum(apply(lm, 2, var))
    waic   <- -2 * (lppd - p_waic)
    LPML   <- sum(log(1 / apply(lm, 2, function(x) mean(exp(-x)))))
    write_csv(tibble(WAIC = waic, LPML = LPML, lppd = lppd, pWAIC = p_waic),
              file.path(scenario_dir, "criteria.csv"))
    cat(sprintf("[%s] WAIC = %.2f | LPML = %.2f\n", model_type, waic, LPML))
  }
  
  # ── 5p. Diagnósticos ACF ──────────────────────────────────────────────────────
  params_struct <- c(beta_names, lambda_names,
                     if (is_spatial) "tau_s" else character(0))
  all_diag <- c(params_struct, as.vector(gamma_names_mat))
  
  ESS_struct  <- effectiveSize(mcmc_list_full[, params_struct])
  Rhat_struct <- safe_gelman(mcmc_list_full[, params_struct])
  
  acf_results <- do.call(rbind, lapply(all_diag, function(nm) {
    ac   <- acf(samples_mat[, nm], lag.max = 200, plot = FALSE)
    lags <- as.vector(ac$lag[-1]); acfs <- as.vector(ac$acf[-1])
    ess_v  <- tryCatch(as.numeric(effectiveSize(mcmc_list_full[, nm])),
                       error = function(e) NA_real_)
    rhat_v <- tryCatch(safe_gelman(mcmc_list_full[, nm]),
                       error = function(e) NA_real_)
    tibble(Parameter = nm, ESS = ess_v, Rhat = rhat_v,
           lag_0.10 = { v <- lags[which(abs(acfs) < 0.10)[1]]; ifelse(is.na(v), Inf, v) },
           lag_0.05 = { v <- lags[which(abs(acfs) < 0.05)[1]]; ifelse(is.na(v), Inf, v) },
           acf_lag1 = acfs[1])
  }))
  write_csv(acf_results, file.path(scenario_dir, "acf_diagnostics.csv"))
  
  # ACF de beta (e tau_s)
  params_acf_beta <- c(beta_names, if (is_spatial) "tau_s" else character(0))
  acf_beta_df <- do.call(rbind, lapply(params_acf_beta, function(nm) {
    ac <- acf(samples_mat[, nm], lag.max = 100, plot = FALSE)
    tibble(Parameter = nm, Lag = as.vector(ac$lag[-1]), ACF = as.vector(ac$acf[-1]))
  }))
  ggsave(
    file.path(scenario_dir, "acf_beta.png"),
    ggplot(acf_beta_df, aes(x = Lag, y = ACF)) +
      geom_col(width = 0.6, fill = "grey50") +
      geom_hline(yintercept = c(-0.10, 0.10), linetype = "dashed",
                 color = "blue", linewidth = 0.5) +
      geom_hline(yintercept = c(-0.05, 0.05), linetype = "dotted",
                 color = "red",  linewidth = 0.5) +
      facet_wrap(~ Parameter, scales = "free_y") +
      theme_bw(base_size = 11) +
      labs(title    = paste("ACF de beta (", model_type, ")"),
           subtitle = "Azul: |0.10| | Vermelho: |0.05|"),
    width = 10, height = 5
  )
  
  # ACF de lambda[t]
  acf_lambda_df <- do.call(rbind, lapply(lambda_names, function(nm) {
    ac <- acf(samples_mat[, nm], lag.max = 100, plot = FALSE)
    tibble(Time = as.integer(str_extract(nm, "\\d+")),
           Lag  = as.vector(ac$lag[-1]), ACF = as.vector(ac$acf[-1]))
  }))
  ggsave(
    file.path(scenario_dir, "acf_lambda.png"),
    ggplot(acf_lambda_df, aes(x = Lag, y = ACF)) +
      geom_col(width = 0.6, fill = "steelblue") +
      geom_hline(yintercept = c(-0.10, 0.10), linetype = "dashed",
                 color = "blue", linewidth = 0.5) +
      geom_hline(yintercept = c(-0.05, 0.05), linetype = "dotted",
                 color = "red",  linewidth = 0.5) +
      facet_wrap(~ Time, scales = "free_y", ncol = 6,
                 labeller = label_bquote(lambda[.(Time)])) +
      theme_bw(base_size = 9) +
      labs(title    = paste("ACF de lambda[t] (", model_type, ")"),
           subtitle = "Azul: |0.10| | Vermelho: |0.05|"),
    width = 14, height = 10
  )
  
  # ACF de gamma[k, t] — um arquivo por cluster k
  # USA lapply (não for) para evitar captura tardia de variável de loop
  # em ambiente paralelo, o que causava "objeto 'k_plot' não encontrado"
  invisible(lapply(seq_len(K), function(k_val) {
    acf_g_df <- do.call(rbind, lapply(seq_len(n_times), function(t) {
      nm <- gamma_names_mat[k_val, t]
      ac <- acf(samples_mat[, nm], lag.max = 100, plot = FALSE)
      tibble(Time = t, Lag = as.vector(ac$lag[-1]), ACF = as.vector(ac$acf[-1]))
    }))
    # Força avaliação de k_val dentro do escopo local deste lapply
    k_local <- k_val
    p_acf <- ggplot(acf_g_df, aes(x = Lag, y = ACF)) +
      geom_col(width = 0.6, fill = "darkorange") +
      geom_hline(yintercept = c(-0.10, 0.10), linetype = "dashed",
                 color = "blue", linewidth = 0.5) +
      geom_hline(yintercept = c(-0.05, 0.05), linetype = "dotted",
                 color = "red",  linewidth = 0.5) +
      facet_wrap(~ Time, scales = "free_y", ncol = 6) +
      theme_bw(base_size = 9) +
      labs(title    = paste0("ACF de gamma[", k_local, ", t] (", model_type, ")"),
           subtitle = "Azul: |0.10| | Vermelho: |0.05|",
           x = "Lag", y = "ACF")
    ggsave(file.path(scenario_dir, sprintf("acf_gamma_k%d.png", k_local)),
           p_acf, width = 14, height = 10)
  }))
  
  # Traceplots de beta (e tau_s)
  cores_cadeia <- c("Cadeia 1" = "#2166AC", "Cadeia 2" = "#D6604D")
  df_trace <- do.call(rbind, lapply(seq_len(nchains), function(ch) {
    cm <- as.matrix(mcmc_list_full[[ch]])
    do.call(rbind, lapply(params_acf_beta, function(nm) {
      vals <- cm[, nm]
      tibble(Iter = seq_along(vals), Value = vals,
             ErgMedia = cumsum(vals) / seq_along(vals),
             Parameter = nm, Cadeia = paste0("Cadeia ", ch))
    }))
  }))
  ggsave(
    file.path(scenario_dir, "traceplots_beta.png"),
    ggplot(df_trace, aes(x = Iter, color = Cadeia)) +
      geom_line(aes(y = Value),    alpha = 0.25, linewidth = 0.20) +
      geom_line(aes(y = ErgMedia), alpha = 0.90, linewidth = 0.75) +
      scale_color_manual(values = cores_cadeia) +
      facet_wrap(~ Parameter, scales = "free_y") +
      theme_bw(base_size = 11) + theme(legend.position = "bottom") +
      labs(title    = paste("Traceplots beta (", model_type, ")"),
           subtitle = "Linha grossa = media ergodica | Linha fina = cadeia",
           x = "Iteracao (pos-burnin)", y = "Valor"),
    width = 10, height = max(5, 3 * ceiling(length(params_acf_beta) / 3))
  )
  
  # ── 5q. Retorno resumido ──────────────────────────────────────────────────────
  tibble(
    model = model_type, niter = niter, nburnin = nburnin, thin = thin,
    WAIC = waic, LPML = LPML,
    ESS_beta_min    = min(beta_summary$ESS,   na.rm = TRUE),
    ESS_gamma_min   = min(gamma_summary$ESS,  na.rm = TRUE),
    ESS_lambda_min  = min(lambda_summary$ESS, na.rm = TRUE),
    ESS_tau         = ESS_tau,
    Rhat_beta_max   = max(Rhat_struct[beta_names],   na.rm = TRUE),
    Rhat_lambda_max = max(Rhat_struct[lambda_names], na.rm = TRUE),
    Rhat_gamma_max  = max(gamma_summary$Rhat, na.rm = TRUE),
    lag_max_0.10    = max(acf_results$lag_0.10, na.rm = TRUE),
    lag_max_0.05    = max(acf_results$lag_0.05, na.rm = TRUE)
  )
}

# ==============================================================================
# 6. EXECUÇÃO PARALELA
# ==============================================================================
model_types <- c("spatial", "non_spatial")
n_cores     <- min(length(model_types), parallel::detectCores() - 1)
if (n_cores < 1) n_cores <- 1

output_dir <- file.path(
  "C:/Users/vlara/OneDrive/Estatistica UFMG/Mestrado/Pesquisa/Aplicação",
  "resultados_gamma_temporal"
)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

cl <- makeCluster(n_cores)
clusterExport(cl, c(
  "constants_spatial", "constants_nonspatial", "data_nimble",
  "inits_list_spatial", "inits_list_nonspatial",
  "N_regions", "n_times", "p", "K",
  "h_mat", "cluster_ids",
  "run_model", "output_dir", "REGIONS_INTEREST"
))
clusterEvalQ(cl, {
  library(nimble); library(coda); library(dplyr)
  library(ggplot2); library(readr); library(stringr); library(tibble)
  Sys.setenv(OMP_NUM_THREADS = "1"); Sys.setenv(MKL_NUM_THREADS = "1")
  if (requireNamespace("RhpcBLASctl", quietly = TRUE))
    RhpcBLASctl::blas_set_num_threads(1)
})

resultados <- parLapply(cl, model_types, function(m) run_model(m, output_dir))
stopCluster(cl)

# ==============================================================================
# 7. CONSOLIDAÇÃO E GRÁFICOS COMPARATIVOS
# ==============================================================================
resumo <- bind_rows(resultados)
write_csv(resumo, file.path(output_dir, "resumo_comparativo.csv"))
cat("\n=== RESUMO COMPARATIVO ===\n"); print(resumo)

# Comparativo lambda[t]
lambda_all <- lapply(model_types, function(m) {
  path <- file.path(output_dir, m, "lambda_summary.csv")
  if (!file.exists(path)) return(NULL)
  read_csv(path, show_col_types = FALSE)
}) |> bind_rows()

if (nrow(lambda_all) > 0) {
  ggsave(
    file.path(output_dir, "lambda_comparativo.png"),
    ggplot(lambda_all, aes(x = Time, y = Mean, color = model, fill = model)) +
      geom_ribbon(aes(ymin = Lower, ymax = Upper), alpha = 0.15, color = NA) +
      geom_line(linewidth = 0.9) +
      theme_bw(base_size = 12) + theme(legend.position = "bottom") +
      labs(title = "Comparacao lambda[t]: espacial vs. nao-espacial",
           x = "Tempo", y = expression(lambda[t]),
           color = "Modelo", fill = "Modelo"),
    width = 10, height = 5, dpi = 300
  )
}

# Comparativo gamma[k, t]
gamma_all <- lapply(model_types, function(m) {
  path <- file.path(output_dir, m, "gamma_summary.csv")
  if (!file.exists(path)) return(NULL)
  read_csv(path, show_col_types = FALSE) %>% mutate(model = m)
}) |> bind_rows()

if (nrow(gamma_all) > 0) {
  # Cria uma banda estática [0.0, 0.1] para a priori de gamma[1,t] no M4
  prior_ribbon <- tibble(
    Time    = seq_len(n_times),
    k       = 1L,
    a_lower = constants_spatial$a_unif, # 0.0
    b_upper = constants_spatial$b_unif  # 0.1
  )
  
  ggsave(
    file.path(output_dir, "gamma_comparativo.png"),
    ggplot(gamma_all, aes(x = Time, y = Mean, color = model, fill = model)) +
      geom_ribbon(aes(ymin = Lower, ymax = Upper), alpha = 0.12, color = NA) +
      geom_line(linewidth = 0.8) +
      # Banda da priori apenas para k=1 (âncora)
      geom_ribbon(
        data = prior_ribbon,
        aes(x = Time, ymin = a_lower, ymax = b_upper),
        inherit.aes = FALSE,
        fill = "grey40", alpha = 0.10, linetype = "dotted", color = "grey40"
      ) +
      facet_wrap(
        ~ k, ncol = 2, scales = "free_y",
        labeller = labeller(k = function(x) paste0("gamma[", x, ", t]"))
      ) +
      theme_bw(base_size = 12) + theme(legend.position = "bottom") +
      labs(title = "Comparacao gamma[k,t]: espacial vs. nao-espacial (K=4)",
           subtitle = "Banda cinza em k=1: suporte da priori estatica [0, 0.1]",
           x = "Tempo", y = expression(gamma[kt]),
           color = "Modelo", fill = "Modelo"),
    width = 10, height = 8, dpi = 300
  )
}

# Comparativos theta e mu
for (tipo in c("theta", "mu")) {
  all_df <- lapply(model_types, function(m) {
    path <- file.path(output_dir, m, paste0(tipo, "_selected.csv"))
    if (!file.exists(path)) return(NULL)
    read_csv(path, show_col_types = FALSE)
  }) |> bind_rows()
  if (nrow(all_df) > 0) {
    ggsave(
      file.path(output_dir, paste0(tipo, "_comparativo.png")),
      ggplot(all_df, aes(x = Time, y = Mean, color = model, fill = model)) +
        geom_ribbon(aes(ymin = Lower, ymax = Upper), alpha = 0.15, color = NA) +
        geom_line(linewidth = 0.8) +
        facet_wrap(~ Region, scales = "free_y", ncol = 3) +
        theme_bw(base_size = 12) + theme(legend.position = "bottom") +
        labs(title  = sprintf("Comparacao %s: espacial vs. nao-espacial", tipo),
             x = "Tempo", color = "Modelo", fill = "Modelo"),
      width = 14, height = 10, dpi = 300
    )
  }
}

# Comparativo epsilon por cluster
eps_cluster_all <- lapply(model_types, function(m) {
  path <- file.path(output_dir, m, "epsilon_cluster_summary.csv")
  if (!file.exists(path)) return(NULL)
  read_csv(path, show_col_types = FALSE) %>% mutate(model = m)
}) |> bind_rows()

if (nrow(eps_cluster_all) > 0) {
  ggsave(
    file.path(output_dir, "epsilon_comparativo.png"),
    ggplot(eps_cluster_all,
           aes(x = Time, y = Mean, color = model, fill = model)) +
      geom_ribbon(aes(ymin = Lower, ymax = Upper), alpha = 0.12, color = NA) +
      geom_line(linewidth = 0.8) +
      facet_wrap(
        ~ k, ncol = 2, scales = "free_y",
        labeller = labeller(k = function(x) paste0("epsilon[", x, ", t]"))
      ) +
      theme_bw(base_size = 12) + theme(legend.position = "bottom") +
      labs(title = "Comparacao epsilon por cluster: espacial vs. nao-espacial (K=4)",
           subtitle = "Banda = IC 95% HPD",
           x = "Tempo", y = expression(epsilon[kt]),
           color = "Modelo", fill = "Modelo"),
    width = 10, height = 8, dpi = 300
  )
}
cat("\nTempo total de execucao:\n")
print(Sys.time() - inicio_global)
