# ==============================================================================
# ajuste_intercepto_temporal_gamma_t.R
# Modelo Poisson log-linear com intercepto temporal i.i.d. e gamma variável no tempo
#
# log(mu[i,t]) = beta0t[t] + log(E[i,t]) + log(epsilon[i,t])
#                + beta[1:p] %*% x[i,t,1:p]  (+ s[i]  se espacial)
#
# epsilon[i,t] = 1 - sum_{k=1}^K h[i,k] * gamma[k,t]
#
# Priors:
#   beta0t[t]      ~ N(0, 1)  i.i.d. para cada t = 1,...,T
#   beta[j]        ~ N(0, 1)
#   gamma[1,t]     ~ Unif(0, 0.1)   para cada t
#   gamma[k,t]     ~ Unif(0, 1 - sum_{j<k} gamma[j,t]) para k=2,...,K e cada t
#   sigma_s        ~ Half-t(0, 1, 1)   [apenas espacial]
#   tau_s          = 1/sigma_s^2
#   s              ~ ICAR(tau_s)        [apenas espacial]
#
# Gráficos seguem o estilo de ajuste_gamma_temporal.R (painéis com ribbon, ACF separados, etc.)
# ==============================================================================

inicio_global <- Sys.time()

setwd("C:/Users/vlara/OneDrive/Estatistica UFMG/Mestrado/Pesquisa/PesquisaMestrado")

pkgs <- c("nimble", "coda", "parallel", "dplyr", "ggplot2",
          "tidyr", "readr", "stringr", "tibble")
for (pkg in pkgs) {
  if (!require(pkg, character.only = TRUE)) install.packages(pkg)
  library(pkg, character.only = TRUE)
}

Sys.setenv(OMP_NUM_THREADS = "1"); Sys.setenv(MKL_NUM_THREADS = "1")
if (requireNamespace("RhpcBLASctl", quietly = TRUE))
  RhpcBLASctl::blas_set_num_threads(1)

# ==============================================================================
# 1. DADOS
# ==============================================================================
source("_dataCaseStudy.r")

E_norm <- E / mean(E)

N_regions <- nrow(Y_mat)
n_times   <- ncol(Y_mat)
p         <- dim(x)[3]
K         <- ncol(data$hAI)

adj_vec     <- as.integer(data$adj)
num_vec     <- as.integer(data$num)
n_adj_val   <- as.integer(data$sumNumNeigh)
weights_vec <- rep(1.0, n_adj_val)
h_mat       <- data$hAI

stopifnot(
  identical(rownames(Y_mat), rownames(E)),
  identical(dim(Y_mat), dim(E)),
  identical(dim(Y_mat)[1:2], dim(x)[1:2]),
  sum(num_vec) == n_adj_val
)
cat(sprintf("N = %d | T = %d | p = %d | K = %d\n", N_regions, n_times, p, K))

# ==============================================================================
# 2. REGIÕES PARA MONITORAR
# ==============================================================================
cluster_ids <- apply(h_mat, 1, sum)
set.seed(42)
REGIONS_INTEREST <- unlist(lapply(1:K, function(cl) {
  regs <- which(cluster_ids == cl)
  if (length(regs) >= 3) sample(regs, 3) else regs
}))
cat("Regiões monitoradas:\n"); print(sort(REGIONS_INTEREST))

MU_MONITORS <- unlist(lapply(
  REGIONS_INTEREST,
  function(r) paste0("mu[", r, ", ", seq_len(n_times), "]")
))

# ==============================================================================
# 3. CONSTANTES E DADOS NIMBLE
# ==============================================================================
constants_spatial <- list(
  n_regions = N_regions, n_times = n_times, p = p, K = K,
  h = h_mat, mu_beta = rep(0, p),
  a_unif = 0.0, b_unif = 0.1,
  adj = adj_vec, num = num_vec, weights = weights_vec, n_adj = n_adj_val
)
constants_nonspatial <- constants_spatial[
  setdiff(names(constants_spatial), c("adj", "num", "weights", "n_adj"))
]

data_nimble <- list(Y = Y_mat, E = E_norm, x = x)

# ==============================================================================
# 4. INICIALIZAÇÕES (gamma agora é matriz K × T)
# ==============================================================================
set.seed(123)

# Função auxiliar para gerar uma coluna de gamma (para um dado t)
gerar_gamma_coluna <- function(K) {
  g <- numeric(K)
  g[1] <- runif(1, 0.02, 0.08)
  for (k in 2:K) {
    max_val <- 1 - sum(g[1:(k-1)])
    g[k] <- runif(1, 0, min(0.2, max_val))
  }
  g
}

gamma_mat1 <- replicate(n_times, gerar_gamma_coluna(K))  # K × T
gamma_mat2 <- replicate(n_times, gerar_gamma_coluna(K))

inits_spatial_1 <- list(
  beta0t  = rep(0, n_times),
  beta    = rep(0, p),
  gamma   = gamma_mat1,
  sigma_s = 0.5,
  s       = rep(0, N_regions)
)
inits_spatial_2 <- list(
  beta0t  = rnorm(n_times, 0, 0.3),
  beta    = rnorm(p, 0, 0.3),
  gamma   = gamma_mat2,
  sigma_s = 1.0,
  s       = rep(0, N_regions)
)
inits_nonspatial_1 <- list(
  beta0t = rep(0, n_times),
  beta   = rep(0, p),
  gamma  = gamma_mat1
)
inits_nonspatial_2 <- list(
  beta0t = rnorm(n_times, 0, 0.3),
  beta   = rnorm(p, 0, 0.3),
  gamma  = gamma_mat2
)

inits_list_spatial    <- list(inits_spatial_1,    inits_spatial_2)
inits_list_nonspatial <- list(inits_nonspatial_1, inits_nonspatial_2)

# ==============================================================================
# 5. FUNÇÃO WORKER
# ==============================================================================
run_model <- function(model_type, output_dir) {
  
  library(nimble); library(coda); library(dplyr)
  library(ggplot2); library(readr); library(stringr); library(tibble)
  
  # ── Código NIMBLE (gamma indexado como gamma[k, t]) ────────────────────────
  code_spatial <- nimbleCode({
    for (t in 1:n_times) {
      beta0t[t] ~ dnorm(0, sd = 1)
      gamma[1, t] ~ dunif(a_unif, b_unif)
      for (k in 2:K) {
        gamma[k, t] ~ dunif(0, 1 - sum(gamma[1:(k-1), t]))
      }
    }
    for (j in 1:p) {
      beta[j] ~ dnorm(0, sd = 1)
    }
    sigma_s ~ T(dt(0, 1, 1), 0, )
    tau_s   <- 1 / (sigma_s^2)
    s[1:n_regions] ~ dcar_normal(
      adj[1:n_adj], weights[1:n_adj], num[1:n_regions], tau_s, zero_mean = 1
    )
    for (i in 1:n_regions) {
      for (t in 1:n_times) {
        epsilon[i, t] <- 1 - inprod(h[i, 1:K], gamma[1:K, t])
        log(mu[i, t]) <- beta0t[t] + log(E[i, t]) + log(epsilon[i, t]) +
          inprod(beta[1:p], x[i, t, 1:p]) + s[i]
        Y[i, t]          ~ dpois(mu[i, t])
        logLik_Y[i, t]   <- dpois(Y[i, t], mu[i, t], log = TRUE)
      }
    }
  })
  
  code_nonspatial <- nimbleCode({
    for (t in 1:n_times) {
      beta0t[t] ~ dnorm(0, sd = 1)
      gamma[1, t] ~ dunif(a_unif, b_unif)
      for (k in 2:K) {
        gamma[k, t] ~ dunif(0, 1 - sum(gamma[1:(k-1), t]))
      }
    }
    for (j in 1:p) {
      beta[j] ~ dnorm(0, sd = 1)
    }
    for (i in 1:n_regions) {
      for (t in 1:n_times) {
        epsilon[i, t] <- 1 - inprod(h[i, 1:K], gamma[1:K, t])
        log(mu[i, t]) <- beta0t[t] + log(E[i, t]) + log(epsilon[i, t]) +
          inprod(beta[1:p], x[i, t, 1:p])
        Y[i, t]          ~ dpois(mu[i, t])
        logLik_Y[i, t]   <- dpois(Y[i, t], mu[i, t], log = TRUE)
      }
    }
  })
  
  is_spatial <- (model_type == "spatial")
  model_code <- if (is_spatial) code_spatial     else code_nonspatial
  constants  <- if (is_spatial) constants_spatial else constants_nonspatial
  inits_list <- if (is_spatial) inits_list_spatial else inits_list_nonspatial
  
  cat("\n=== Modelo:", model_type, "===\n")
  scenario_dir <- file.path(output_dir, model_type)
  dir.create(scenario_dir, recursive = TRUE, showWarnings = FALSE)
  
  # ── Compilar ──────────────────────────────────────────────────────────────
  model  <- nimbleModel(
    code = model_code, constants = constants,
    data = data_nimble, inits = inits_list[[1]], check = FALSE
  )
  Cmodel <- compileNimble(model)
  
  conf <- configureMCMC(model)
  # Amostradores para gamma: um AF_slice por coluna t (bloco de K parâmetros)
  conf$removeSamplers("gamma")
  for (t_idx in seq_len(n_times)) {
    conf$addSampler(
      target = paste0("gamma[", seq_len(K), ", ", t_idx, "]"),
      type = "AF_slice"
    )
  }
  # beta0t já é amostrado individualmente por padrão (RW)
  
  monitors_base <- c("beta0t", "beta", "gamma", "logLik_Y")
  if (is_spatial) monitors_base <- c(monitors_base, "s", "sigma_s", "tau_s")
  conf$addMonitors(monitors_base)
  conf$addMonitors(MU_MONITORS)
  
  Rmcmc <- buildMCMC(conf)
  Cmcmc <- compileNimble(Rmcmc, project = model)
  
  # ── MCMC ──────────────────────────────────────────────────────────────────
  niter   <- 50000; nburnin <- 10000; nchains <- 2; thin <- 10
  cat(sprintf("[%s] niter=%d | nburnin=%d | thin=%d\n",
              model_type, niter, nburnin, thin))
  
  samples <- runMCMC(
    Cmcmc, niter = niter, nburnin = nburnin, nchains = nchains, thin = thin,
    inits = inits_list, samplesAsCodaMCMC = TRUE, summary = FALSE, WAIC = FALSE
  )
  saveRDS(samples, file.path(scenario_dir, "samples.rds"))
  
  samples_mat    <- as.matrix(samples)
  mcmc_list_full <- mcmc.list(lapply(1:nchains, function(ch) as.mcmc(samples[[ch]])))
  rm(samples); gc()
  
  # ── Auxiliares ────────────────────────────────────────────────────────────
  safe_hpd <- function(sv) {
    if (var(sv) < 1e-12) return(c(NA_real_, NA_real_))
    as.numeric(HPDinterval(as.mcmc(sv), prob = 0.95))
  }
  safe_gelman <- function(obj) {
    tryCatch(
      gelman.diag(obj, autoburnin = FALSE)$psrf[, 1],
      error = function(e) rep(NA_real_, nvar(obj))
    )
  }
  
  beta0t_names <- paste0("beta0t[", 1:n_times, "]")
  beta_names   <- paste0("beta[",   1:p,        "]")
  # Nomes de gamma: gamma[1, 1], gamma[2, 1], ..., gamma[K, T]
  gamma_names_mat <- outer(
    1:K, 1:n_times,
    function(k, t) paste0("gamma[", k, ", ", t, "]")
  )  # K × T
  
  # ── Epsilon posterior (N × T) ─────────────────────────────────────────────
  n_draw        <- nrow(samples_mat)
  epsilon_draws <- array(NA_real_, dim = c(n_draw, N_regions, n_times))
  for (t in seq_len(n_times)) {
    g_t <- samples_mat[, gamma_names_mat[, t], drop = FALSE]  # n_draw × K
    epsilon_draws[, , t] <- 1 - g_t %*% t(h_mat)              # n_draw × N
  }
  
  # ── Sumário de gamma[k, t] (estilo do exemplo) ───────────────────────────
  gamma_summary <- do.call(rbind, lapply(1:K, function(k) {
    do.call(rbind, lapply(1:n_times, function(t) {
      nm  <- gamma_names_mat[k, t]
      sv  <- samples_mat[, nm]; hpd <- safe_hpd(sv)
      tibble(k = k, Time = t, Parameter = nm,
             Mean = mean(sv), SD = sd(sv), Lower = hpd[1], Upper = hpd[2],
             ESS  = as.numeric(effectiveSize(mcmc_list_full[, nm])),
             Rhat = safe_gelman(mcmc_list_full[, nm]))
    }))
  }))
  write_csv(gamma_summary, file.path(scenario_dir, "gamma_summary.csv"))
  
  # Painel temporal de gamma (facet por cluster k)
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
      labs(title = paste("Trajetória de gamma[k, t] (", model_type, ")"),
           x = "Tempo", y = expression(gamma[kt])),
    width = 10, height = 7
  )
  
  # ── Sumário de epsilon para regiões selecionadas ─────────────────────────
  # ── Sumário de epsilon por cluster ───────────────────────────────────────
  eps_full <- do.call(rbind, lapply(1:N_regions, function(i) {
    do.call(rbind, lapply(1:n_times, function(t) {
      sv  <- epsilon_draws[, i, t]; hpd <- safe_hpd(sv)
      tibble(Region = i, Cluster = cluster_ids[i], Time = t,
             Mean = mean(sv), Lower = hpd[1], Upper = hpd[2])
    }))
  }))
  write_csv(eps_full, file.path(scenario_dir, "epsilon_summary.csv"))
  
  # Agrega por cluster × tempo (pool de draws das regiões do cluster)
  eps_by_cluster <- do.call(rbind, lapply(1:K, function(k) {
    regs_k <- which(cluster_ids == k)
    do.call(rbind, lapply(1:n_times, function(t) {
      sv  <- as.vector(epsilon_draws[, regs_k, t])
      hpd <- safe_hpd(sv)
      tibble(k = k, Time = t, Mean = mean(sv), Lower = hpd[1], Upper = hpd[2])
    }))
  }))
  write_csv(eps_by_cluster, file.path(scenario_dir, "epsilon_cluster_summary.csv"))
  
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
      labs(title = paste("Trajetória de epsilon por cluster (", model_type, ")"),
           x = "Tempo", y = expression(epsilon[kt])),
    width = 10, height = 7
  )
  
  # ── Sumário de beta0t ────────────────────────────────────────────────────
  beta0t_summary <- do.call(rbind, lapply(1:n_times, function(t) {
    nm <- beta0t_names[t]; sv <- samples_mat[, nm]; hpd <- safe_hpd(sv)
    tibble(Time = t, Parameter = nm, Mean = mean(sv), SD = sd(sv),
           Lower = hpd[1], Upper = hpd[2],
           ESS  = as.numeric(effectiveSize(mcmc_list_full[, nm])),
           Rhat = safe_gelman(mcmc_list_full[, nm]))
  }))
  write_csv(beta0t_summary, file.path(scenario_dir, "beta0t_summary.csv"))
  
  anos_label <- colnames(data_nimble$Y)
  if (is.null(anos_label)) anos_label <- seq_len(n_times)
  
  ggsave(
    file.path(scenario_dir, "painel_beta0t.png"),
    ggplot(beta0t_summary, aes(x = Time)) +
      geom_ribbon(aes(ymin = Lower, ymax = Upper), fill = "grey60", alpha = 0.4) +
      geom_line(aes(y = Mean), color = "black", linewidth = 0.9) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
      scale_x_continuous(breaks = 1:n_times, labels = anos_label) +
      theme_bw(base_size = 12) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      labs(title = paste("Intercepto temporal beta0t (", model_type, ")"),
           x = "Ano", y = expression(beta[0][t])),
    width = 10, height = 5
  )
  
  # ── Sumário de beta ──────────────────────────────────────────────────────
  beta_summary <- do.call(rbind, lapply(beta_names, function(nm) {
    sv  <- samples_mat[, nm]; hpd <- safe_hpd(sv)
    tibble(Parameter = nm, Mean = mean(sv), SD = sd(sv),
           Lower = hpd[1], Upper = hpd[2],
           ESS  = as.numeric(effectiveSize(mcmc_list_full[, nm])),
           Rhat = safe_gelman(mcmc_list_full[, nm]))
  }))
  write_csv(beta_summary, file.path(scenario_dir, "beta_summary.csv"))
  
  # ── tau_s e s (espacial) ─────────────────────────────────────────────────
  ESS_tau <- NA_real_
  if (is_spatial) {
    tau_sv <- samples_mat[, "tau_s"]; hpd_t <- safe_hpd(tau_sv)
    tau_sum <- tibble(
      Parameter = "tau_s", Mean = mean(tau_sv), SD = sd(tau_sv),
      Lower = hpd_t[1], Upper = hpd_t[2],
      ESS  = as.numeric(effectiveSize(mcmc_list_full[, "tau_s"])),
      Rhat = safe_gelman(mcmc_list_full[, "tau_s"])
    )
    write_csv(tau_sum, file.path(scenario_dir, "tau_summary.csv"))
    ESS_tau <- tau_sum$ESS
    
    s_names   <- paste0("s[", 1:N_regions, "]")
    s_summary <- do.call(rbind, lapply(1:N_regions, function(i) {
      sv <- samples_mat[, s_names[i]]; hpd <- safe_hpd(sv)
      tibble(Region = i, Mean = mean(sv), SD = sd(sv),
             Lower = hpd[1], Upper = hpd[2],
             ESS = as.numeric(effectiveSize(mcmc_list_full[, s_names[i]])))
    }))
    write_csv(s_summary, file.path(scenario_dir, "s_summary.csv"))
    
    ggsave(
      file.path(scenario_dir, "s_posterior.png"),
      ggplot(s_summary, aes(x = Region, y = Mean)) +
        geom_point(size = 0.9) +
        geom_errorbar(aes(ymin = Lower, ymax = Upper), width = 0.4, linewidth = 0.3) +
        geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
        theme_bw() +
        labs(title = paste("Efeito espacial s[i] (", model_type, ")"),
             y = "s[i]", x = "Região"),
      width = 10, height = 5
    )
  }
  
  # ── mu para regiões selecionadas ─────────────────────────────────────────
  mu_summary <- do.call(rbind, lapply(REGIONS_INTEREST, function(i) {
    do.call(rbind, lapply(1:n_times, function(t) {
      nm  <- paste0("mu[", i, ", ", t, "]")
      sv  <- samples_mat[, nm]; hpd <- safe_hpd(sv)
      tibble(Region = i, Time = t, Mean = mean(sv),
             Lower = hpd[1], Upper = hpd[2], model = model_type)
    }))
  }))
  write_csv(mu_summary, file.path(scenario_dir, "mu_selected.csv"))
  
  mu_regs <- sort(unique(mu_summary$Region))
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
  ggsave(
    file.path(scenario_dir, "painel_mu.png"),
    ggplot(mu_summary, aes(x = Time)) +
      geom_ribbon(aes(ymin = Lower, ymax = Upper), fill = "darkorange", alpha = 0.3) +
      geom_line(aes(y = Mean), color = "darkorange", linewidth = 0.8) +
      facet_wrap(
        ~ Region, ncol = 3, scales = "free_y",
        labeller = labeller(Region = make_label(mu_regs))
      ) +
      theme_bw(base_size = 10) +
      labs(title = paste("Mu estimado (", model_type, ")"),
           x = "Tempo", y = expression(mu[it])),
    width = 12, height = 10
  )
  
  # ── WAIC e LPML ───────────────────────────────────────────────────────────
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
  
  # ── Diagnósticos ACF (organizados como no exemplo) ───────────────────────
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
  
  # ACF de beta0t
  acf_beta0t_df <- do.call(rbind, lapply(beta0t_names, function(nm) {
    ac <- acf(samples_mat[, nm], lag.max = 100, plot = FALSE)
    t_idx <- as.integer(str_extract(nm, "\\d+"))
    tibble(Time = t_idx, Lag = as.vector(ac$lag[-1]), ACF = as.vector(ac$acf[-1]))
  }))
  ggsave(
    file.path(scenario_dir, "acf_beta0t.png"),
    ggplot(acf_beta0t_df, aes(x = Lag, y = ACF)) +
      geom_col(width = 0.6, fill = "steelblue") +
      geom_hline(yintercept = c(-0.10, 0.10), linetype = "dashed",
                 color = "blue", linewidth = 0.5) +
      geom_hline(yintercept = c(-0.05, 0.05), linetype = "dotted",
                 color = "red",  linewidth = 0.5) +
      facet_wrap(~ Time, scales = "free_y", ncol = 6,
                 labeller = label_bquote(beta[0*.(Time)])) +
      theme_bw(base_size = 9) +
      labs(title    = paste("ACF de beta0t (", model_type, ")"),
           subtitle = "Azul: |0.10| | Vermelho: |0.05|"),
    width = 14, height = 10
  )
  
  # ACF de gamma[k, t] por cluster k
  invisible(lapply(1:K, function(k_val) {
    acf_g_df <- do.call(rbind, lapply(1:n_times, function(t) {
      nm <- gamma_names_mat[k_val, t]
      ac <- acf(samples_mat[, nm], lag.max = 100, plot = FALSE)
      tibble(Time = t, Lag = as.vector(ac$lag[-1]), ACF = as.vector(ac$acf[-1]))
    }))
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
  df_trace <- do.call(rbind, lapply(1:nchains, function(ch) {
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
           subtitle = "Linha grossa = média ergódica | Linha fina = cadeia",
           x = "Iteração (pós-burnin)", y = "Valor"),
    width = 10, height = max(5, 3 * ceiling(length(params_acf_beta) / 3))
  )
  
  # ── Retorno resumido ──────────────────────────────────────────────────────
  tibble(
    model = model_type, niter = niter, nburnin = nburnin, thin = thin,
    WAIC = waic, LPML = LPML,
    ESS_beta_min    = min(beta_summary$ESS,   na.rm = TRUE),
    ESS_gamma_min   = min(gamma_summary$ESS,  na.rm = TRUE),
    ESS_beta0t_min  = min(beta0t_summary$ESS, na.rm = TRUE),
    ESS_tau         = ESS_tau,
    Rhat_beta_max   = max(beta_summary$Rhat,   na.rm = TRUE),
    Rhat_gamma_max  = max(gamma_summary$Rhat,  na.rm = TRUE),
    Rhat_beta0t_max = max(beta0t_summary$Rhat, na.rm = TRUE)
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
  "resultados_intercepto_temporal_gamma_t"
)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

cl <- makeCluster(n_cores)
clusterExport(cl, c(
  "constants_spatial", "constants_nonspatial", "data_nimble",
  "inits_list_spatial", "inits_list_nonspatial",
  "N_regions", "n_times", "p", "K", "h_mat", "cluster_ids",
  "run_model", "output_dir", "REGIONS_INTEREST", "MU_MONITORS"
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

# Comparativo beta0t
beta0t_all <- lapply(model_types, function(m) {
  path <- file.path(output_dir, m, "beta0t_summary.csv")
  if (!file.exists(path)) return(NULL)
  read_csv(path, show_col_types = FALSE) %>% mutate(model = m)
}) |> bind_rows()

if (nrow(beta0t_all) > 0) {
  anos_label <- colnames(data_nimble$Y)
  if (is.null(anos_label)) anos_label <- seq_len(n_times)
  beta0t_all <- beta0t_all %>% mutate(Ano = rep(anos_label, length(model_types)))
  
  ggsave(
    file.path(output_dir, "beta0t_comparativo.png"),
    ggplot(beta0t_all, aes(x = Time, y = Mean, color = model, fill = model)) +
      geom_ribbon(aes(ymin = Lower, ymax = Upper), alpha = 0.15, color = NA) +
      geom_line(linewidth = 0.9) +
      scale_x_continuous(breaks = 1:n_times, labels = anos_label) +
      theme_bw(base_size = 12) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            legend.position = "bottom") +
      labs(title = "Comparação beta0t: espacial vs. não-espacial",
           x = "Ano", y = expression(beta[0][t]),
           color = "Modelo", fill = "Modelo"),
    width = 10, height = 5, dpi = 300
  )
}

# Comparativo gamma[k,t]
gamma_all <- lapply(model_types, function(m) {
  path <- file.path(output_dir, m, "gamma_summary.csv")
  if (!file.exists(path)) return(NULL)
  read_csv(path, show_col_types = FALSE) %>% mutate(model = m)
}) |> bind_rows()

if (nrow(gamma_all) > 0) {
  # Banda da priori estática [0, 0.1]
  prior_ribbon <- tibble(
    Time    = seq_len(n_times),
    k       = 1L,
    a_lower = constants_spatial$a_unif,
    b_upper = constants_spatial$b_unif
  )
  ggsave(
    file.path(output_dir, "gamma_comparativo.png"),
    ggplot(gamma_all, aes(x = Time, y = Mean, color = model, fill = model)) +
      geom_ribbon(aes(ymin = Lower, ymax = Upper), alpha = 0.12, color = NA) +
      geom_line(linewidth = 0.8) +
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
      labs(title = "Comparação gamma[k,t]: espacial vs. não-espacial (K=4)",
           subtitle = "Banda cinza em k=1: suporte da priori estática [0, 0.1]",
           x = "Tempo", y = expression(gamma[kt]),
           color = "Modelo", fill = "Modelo"),
    width = 10, height = 8, dpi = 300
  )
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
      labs(title = "Comparação epsilon por cluster: espacial vs. não-espacial (K=4)",
           subtitle = "Banda = IC 95% HPD",
           x = "Tempo", y = expression(epsilon[kt]),
           color = "Modelo", fill = "Modelo"),
    width = 10, height = 8, dpi = 300
  )
}

cat("\nTempo total de execução:\n")
print(Sys.time() - inicio_global)