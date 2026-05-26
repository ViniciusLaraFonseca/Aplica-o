# ==============================================================================
# ajuste_intercepto_temporal.R
# Modelo Poisson log-linear com intercepto temporal i.i.d.
#
# log(mu[i,t]) = beta0t[t] + log(E[i,t]) + log(epsilon[i])
#                + beta[1:p] %*% x[i,t,1:p]  (+ s[i]  se espacial)
#
# Priors:
#   beta0t[t]  ~ N(0, 1)  i.i.d. para cada t = 1,...,T
#   beta[j]    ~ N(0, 1)
#   gamma[1]   ~ Unif(0, 0.1)
#   gamma[j]   ~ Unif(0, 1 - sum(gamma[1:(j-1)]))
#   sigma_s    ~ Half-t(0, 1, 1)   [apenas espacial]
#   tau_s      = 1/sigma_s^2
#   s          ~ ICAR(tau_s)        [apenas espacial]
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
# 1. DADOS  (gerados por covariaveis.R)
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
# 2. REGIÕES PARA MONITORAR mu
# ==============================================================================
cluster_ids <- apply(h_mat, 1, sum)
set.seed(42)
REGIONS_INTEREST <- unlist(lapply(1:K, function(cl) {
  regs <- which(cluster_ids == cl)
  if (length(regs) >= 3) sample(regs, 3) else regs
}))
cat("Regioes monitoradas:\n"); print(sort(REGIONS_INTEREST))

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
# 4. INICIALIZAÇÕES
# ==============================================================================
set.seed(123)

inits_spatial_1 <- list(
  beta0t  = rep(0, n_times),
  beta    = rep(0, p),
  gamma   = c(0.05, 0.10, 0.10, 0.15),
  sigma_s = 0.5,
  s       = rep(0, N_regions)
)
inits_spatial_2 <- list(
  beta0t  = rnorm(n_times, 0, 0.3),
  beta    = rnorm(p, 0, 0.3),
  gamma   = c(0.04, 0.09, 0.09, 0.14),
  sigma_s = 1.0,
  s       = rep(0, N_regions)
)
inits_nonspatial_1 <- list(
  beta0t = rep(0, n_times),
  beta   = rep(0, p),
  gamma  = c(0.05, 0.10, 0.10, 0.15)
)
inits_nonspatial_2 <- list(
  beta0t = rnorm(n_times, 0, 0.3),
  beta   = rnorm(p, 0, 0.3),
  gamma  = c(0.04, 0.09, 0.09, 0.14)
)

inits_list_spatial    <- list(inits_spatial_1,    inits_spatial_2)
inits_list_nonspatial <- list(inits_nonspatial_1, inits_nonspatial_2)

# ==============================================================================
# 5. FUNÇÃO WORKER
# ==============================================================================
run_model <- function(model_type, output_dir) {

  library(nimble); library(coda); library(dplyr)
  library(ggplot2); library(readr); library(stringr); library(tibble)

  # ── Código NIMBLE ─────────────────────────────────────────────────────────
  code_spatial <- nimbleCode({
    for (t in 1:n_times) beta0t[t] ~ dnorm(0, sd = 1)
    for (j in 1:p)       beta[j]   ~ dnorm(0, sd = 1)
    gamma[1] ~ dunif(a_unif, b_unif)
    for (j in 2:K) gamma[j] ~ dunif(0, 1 - sum(gamma[1:(j-1)]))
    sigma_s ~ T(dt(0, 1, 1), 0, )
    tau_s   <- 1 / (sigma_s^2)
    s[1:n_regions] ~ dcar_normal(
      adj[1:n_adj], weights[1:n_adj], num[1:n_regions], tau_s, zero_mean = 1
    )
    for (i in 1:n_regions) {
      epsilon[i] <- 1 - inprod(h[i, 1:K], gamma[1:K])
      for (t in 1:n_times) {
        log(mu[i, t]) <- beta0t[t] + log(E[i, t]) + log(epsilon[i]) +
                         inprod(beta[1:p], x[i, t, 1:p]) + s[i]
        Y[i, t]          ~ dpois(mu[i, t])
        logLik_Y[i, t]   <- dpois(Y[i, t], mu[i, t], log = TRUE)
      }
    }
  })

  code_nonspatial <- nimbleCode({
    for (t in 1:n_times) beta0t[t] ~ dnorm(0, sd = 1)
    for (j in 1:p)       beta[j]   ~ dnorm(0, sd = 1)
    gamma[1] ~ dunif(a_unif, b_unif)
    for (j in 2:K) gamma[j] ~ dunif(0, 1 - sum(gamma[1:(j-1)]))
    for (i in 1:n_regions) {
      epsilon[i] <- 1 - inprod(h[i, 1:K], gamma[1:K])
      for (t in 1:n_times) {
        log(mu[i, t]) <- beta0t[t] + log(E[i, t]) + log(epsilon[i]) +
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
  conf$removeSamplers("gamma")
  conf$addSampler(target = "gamma", type = "AF_slice")
  # beta0t[1:T]: amostrados individualmente pelo RW padrão do NIMBLE
  # (bloco opcional — descomentar se convergência for lenta)
  # conf$removeSamplers("beta0t")
  # conf$addSampler(target = paste0("beta0t[", 1:n_times, "]"),
  #                 type = "RW_block")

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
  resumo_param <- function(nm_vec) {
    do.call(rbind, lapply(nm_vec, function(nm) {
      sv  <- samples_mat[, nm]
      hpd <- safe_hpd(sv)
      tibble(Parameter = nm, Mean = mean(sv), SD = sd(sv),
             HPD_Lower = hpd[1], HPD_Upper = hpd[2],
             ESS  = as.numeric(effectiveSize(mcmc_list_full[, nm])),
             Rhat = safe_gelman(mcmc_list_full[, nm]))
    }))
  }

  beta0t_names <- paste0("beta0t[", 1:n_times, "]")
  beta_names   <- paste0("beta[",   1:p,        "]")
  gamma_names  <- paste0("gamma[",  1:K,        "]")

  # ── Epsilon ────────────────────────────────────────────────────────────────
  gamma_draws   <- samples_mat[, gamma_names, drop = FALSE]
  epsilon_draws <- 1 - gamma_draws %*% t(h_mat)

  epsilon_summary <- do.call(rbind, lapply(1:N_regions, function(i) {
    hpd <- safe_hpd(epsilon_draws[, i])
    tibble(Region = i, Cluster = cluster_ids[i],
           Eps_Mean = mean(epsilon_draws[, i]),
           Eps_Lower = hpd[1], Eps_Upper = hpd[2])
  }))
  write_csv(epsilon_summary, file.path(scenario_dir, "epsilon_summary.csv"))

  # ── Parâmetros ─────────────────────────────────────────────────────────────
  beta0t_summary <- resumo_param(beta0t_names)
  beta_summary   <- resumo_param(beta_names)
  gamma_summary  <- resumo_param(gamma_names)
  write_csv(beta0t_summary, file.path(scenario_dir, "beta0t_summary.csv"))
  write_csv(beta_summary,   file.path(scenario_dir, "beta_summary.csv"))
  write_csv(gamma_summary,  file.path(scenario_dir, "gamma_summary.csv"))

  # ── Trajetória de beta0t (gráfico central deste modelo) ────────────────────
  anos_label <- colnames(data_nimble$Y)
  if (is.null(anos_label)) anos_label <- seq_len(n_times)
  
  beta0t_plot <- beta0t_summary %>%
    mutate(Tempo = seq_along(beta0t_names),
           Ano   = anos_label)
  
  ggsave(file.path(scenario_dir, "painel_beta0t.png"),
         ggplot(beta0t_plot, aes(x = Tempo)) +
           geom_ribbon(aes(ymin = HPD_Lower, ymax = HPD_Upper),
                       fill = "grey60", alpha = 0.4) +
           geom_line(aes(y = Mean), color = "black", linewidth = 0.9) +
           geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
           scale_x_continuous(breaks = seq_len(n_times), labels = anos_label) +
           theme_bw(base_size = 12) +
           theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
           labs(title = paste("Intercepto temporal beta0t (", model_type, ")"),
                subtitle = "Média posterior e HPD 95%",
                x = "Ano", y = expression(beta[0][t])),
         width = 10, height = 5)

  # ── tau_s e s (espacial) ───────────────────────────────────────────────────
  ESS_tau <- NA_real_
  if (is_spatial) {
    tau_sum <- resumo_param("tau_s")
    write_csv(tau_sum, file.path(scenario_dir, "tau_summary.csv"))
    ESS_tau <- tau_sum$ESS

    s_names   <- paste0("s[", 1:N_regions, "]")
    s_summary <- do.call(rbind, lapply(1:N_regions, function(i) {
      sv <- samples_mat[, s_names[i]]; hpd <- safe_hpd(sv)
      tibble(Region = i, Mean = mean(sv), SD = sd(sv),
             HPD_Lower = hpd[1], HPD_Upper = hpd[2],
             ESS = as.numeric(effectiveSize(mcmc_list_full[, s_names[i]])))
    }))
    write_csv(s_summary, file.path(scenario_dir, "s_summary.csv"))

    ggsave(file.path(scenario_dir, "s_posterior.png"),
      ggplot(s_summary, aes(x = Region, y = Mean)) +
        geom_point(size = 0.9) +
        geom_errorbar(aes(ymin = HPD_Lower, ymax = HPD_Upper),
                      width = 0.4, linewidth = 0.3) +
        geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
        theme_bw() +
        labs(title = paste("Efeito espacial s[i] —", model_type),
             y = "s[i]", x = "Regiao"),
      width = 10, height = 5)
  }

  # ── mu para regiões selecionadas ──────────────────────────────────────────
  mu_summary <- do.call(rbind, lapply(MU_MONITORS, function(nm) {
    idx <- str_match(nm, "mu\\[(\\d+),\\s*(\\d+)\\]")
    i <- as.integer(idx[2]); t <- as.integer(idx[3])
    sv  <- samples_mat[, nm]; hpd <- safe_hpd(sv)
    tibble(Region = i, Time = t, Mean = mean(sv),
           Lower = hpd[1], Upper = hpd[2], model = model_type)
  }))
  write_csv(mu_summary, file.path(scenario_dir, "mu_selected.csv"))

  mu_regs <- sort(unique(mu_summary$Region))
  reg_labels <- setNames(
    sprintf("Reg %d (C%d)\ne=%.3f",
            mu_regs, cluster_ids[mu_regs],
            epsilon_summary$Eps_Mean[mu_regs]),
    as.character(mu_regs)
  )
  ggsave(file.path(scenario_dir, "painel_mu.png"),
    ggplot(mu_summary, aes(x = Time)) +
      geom_ribbon(aes(ymin = Lower, ymax = Upper),
                  fill = "darkorange", alpha = 0.3) +
      geom_line(aes(y = Mean), color = "darkorange") +
      facet_wrap(~Region, scales = "free_y", ncol = 3,
                 labeller = labeller(Region = reg_labels)) +
      theme_bw(base_size = 10) +
      labs(title = paste("mu estimado —", model_type),
           x = "Tempo", y = expression(mu[i * t])),
    width = 12, height = 10)

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

  # ── Diagnósticos (Padrão Unificado) ───────────────────────────────────────
  params_beta   <- c(beta_names, if (is_spatial) "tau_s" else character(0))
  params_gamma  <- gamma_names
  params_struct <- c(params_beta, params_gamma)
  
  # 1. Numérico: beta0t salvo com o resto da estrutura
  ESS_b0t  <- effectiveSize(mcmc_list_full[, beta0t_names])
  Rhat_b0t <- safe_gelman(mcmc_list_full[, beta0t_names])
  acf_b0t  <- do.call(rbind, lapply(beta0t_names, function(nm) {
    ac   <- acf(samples_mat[, nm], lag.max = 200, plot = FALSE)
    lags <- as.vector(ac$lag[-1]); acfs <- as.vector(ac$acf[-1])
    tibble(Parameter = nm, ESS = ESS_b0t[nm], Rhat = Rhat_b0t[nm],
           lag_0.10 = { v <- lags[which(abs(acfs) < 0.10)[1]]; ifelse(is.na(v), Inf, v) },
           lag_0.05 = { v <- lags[which(abs(acfs) < 0.05)[1]]; ifelse(is.na(v), Inf, v) },
           acf_lag1 = acfs[1])
  }))
  
  ESS_struct  <- effectiveSize(mcmc_list_full[, params_struct])
  Rhat_struct <- safe_gelman(mcmc_list_full[, params_struct])
  acf_struct <- do.call(rbind, lapply(params_struct, function(nm) {
    ac   <- acf(samples_mat[, nm], lag.max = 200, plot = FALSE)
    lags <- as.vector(ac$lag[-1]); acfs <- as.vector(ac$acf[-1])
    tibble(Parameter = nm, ESS = ESS_struct[nm], Rhat = Rhat_struct[nm],
           lag_0.10 = { v <- lags[which(abs(acfs) < 0.10)[1]]; ifelse(is.na(v), Inf, v) },
           lag_0.05 = { v <- lags[which(abs(acfs) < 0.05)[1]]; ifelse(is.na(v), Inf, v) },
           acf_lag1 = acfs[1])
  }))
  
  write_csv(bind_rows(acf_b0t, acf_struct),
            file.path(scenario_dir, "acf_diagnostics.csv"))
  
  # 2. Gráficos ACF
  acf_beta_df <- do.call(rbind, lapply(params_beta, function(nm) {
    ac <- acf(samples_mat[, nm], lag.max = 100, plot = FALSE)
    tibble(Parameter = nm, Lag = as.vector(ac$lag[-1]), ACF = as.vector(ac$acf[-1]))
  }))
  ggsave(file.path(scenario_dir, "acf_beta.png"),
         ggplot(acf_beta_df, aes(x = Lag, y = ACF)) +
           geom_col(width = 0.6, fill = "grey50") +
           geom_hline(yintercept = c(-0.10, 0.10), linetype = "dashed", color = "blue", linewidth = 0.5) +
           geom_hline(yintercept = c(-0.05, 0.05), linetype = "dotted", color = "red", linewidth = 0.5) +
           facet_wrap(~Parameter, scales = "free_y") + theme_bw(base_size = 11) +
           labs(title = paste("ACF de beta (", model_type, ")"), subtitle = "Azul: |0.10| | Vermelho: |0.05|"),
         width = 10, height = 5)
  
  acf_gamma_df <- do.call(rbind, lapply(params_gamma, function(nm) {
    ac <- acf(samples_mat[, nm], lag.max = 100, plot = FALSE)
    tibble(Parameter = nm, Lag = as.vector(ac$lag[-1]), ACF = as.vector(ac$acf[-1]))
  }))
  ggsave(file.path(scenario_dir, "acf_gamma.png"),
         ggplot(acf_gamma_df, aes(x = Lag, y = ACF)) +
           geom_col(width = 0.6, fill = "darkorange") +
           geom_hline(yintercept = c(-0.10, 0.10), linetype = "dashed", color = "blue", linewidth = 0.5) +
           geom_hline(yintercept = c(-0.05, 0.05), linetype = "dotted", color = "red", linewidth = 0.5) +
           facet_wrap(~Parameter, scales = "free_y") + theme_bw(base_size = 11) +
           labs(title = paste("ACF de gamma (", model_type, ")"), subtitle = "Azul: |0.10| | Vermelho: |0.05|"),
         width = 10, height = 5)
  
  acf_b0t_df_plot <- do.call(rbind, lapply(beta0t_names, function(nm) {
    ac <- acf(samples_mat[, nm], lag.max = 100, plot = FALSE)
    t_idx <- as.integer(str_extract(nm, "\\d+"))
    tibble(Time = t_idx, Lag = as.vector(ac$lag[-1]), ACF = as.vector(ac$acf[-1]))
  }))
  ggsave(file.path(scenario_dir, "acf_beta0t.png"),
         ggplot(acf_b0t_df_plot, aes(x = Lag, y = ACF)) +
           geom_col(width = 0.6, fill = "steelblue") +
           geom_hline(yintercept = c(-0.10, 0.10), linetype = "dashed", color = "blue", linewidth = 0.5) +
           geom_hline(yintercept = c(-0.05, 0.05), linetype = "dotted", color = "red", linewidth = 0.5) +
           facet_wrap(~Time, scales = "free_y", ncol = 6, labeller = label_bquote(beta[0*.(Time)])) +
           theme_bw(base_size = 9) +
           labs(title = paste("ACF de beta0t (", model_type, ")"), subtitle = "Azul: |0.10| | Vermelho: |0.05|"),
         width = 14, height = 10)
  
  # 3. Traceplots (Apenas Beta e Gamma para manter limpo, igual ao M7)
  cores_cadeia <- c("Cadeia 1" = "#2166AC", "Cadeia 2" = "#D6604D")
  
  df_trace_beta <- do.call(rbind, lapply(1:nchains, function(ch) {
    cm <- as.matrix(mcmc_list_full[[ch]])
    do.call(rbind, lapply(params_beta, function(nm) {
      vals <- cm[, nm]
      tibble(Iter = seq_along(vals), Value = vals, ErgMedia = cumsum(vals) / seq_along(vals),
             Parameter = nm, Cadeia = paste0("Cadeia ", ch))
    }))
  }))
  ggsave(file.path(scenario_dir, "traceplots_beta.png"),
         ggplot(df_trace_beta, aes(x = Iter, color = Cadeia)) +
           geom_line(aes(y = Value), alpha = 0.25, linewidth = 0.20) +
           geom_line(aes(y = ErgMedia), alpha = 0.90, linewidth = 0.75) +
           scale_color_manual(values = cores_cadeia) + facet_wrap(~Parameter, scales = "free_y") +
           theme_bw(base_size = 11) + theme(legend.position = "bottom") +
           labs(title = paste("Traceplots beta (", model_type, ")"), x = "Iteracao (pos-burnin)", y = "Valor"),
         width = 10, height = max(5, 3 * ceiling(length(params_beta) / 3)))
  
  df_trace_gamma <- do.call(rbind, lapply(1:nchains, function(ch) {
    cm <- as.matrix(mcmc_list_full[[ch]])
    do.call(rbind, lapply(params_gamma, function(nm) {
      vals <- cm[, nm]
      tibble(Iter = seq_along(vals), Value = vals, ErgMedia = cumsum(vals) / seq_along(vals),
             Parameter = nm, Cadeia = paste0("Cadeia ", ch))
    }))
  }))
  ggsave(file.path(scenario_dir, "traceplots_gamma.png"),
         ggplot(df_trace_gamma, aes(x = Iter, color = Cadeia)) +
           geom_line(aes(y = Value), alpha = 0.25, linewidth = 0.20) +
           geom_line(aes(y = ErgMedia), alpha = 0.90, linewidth = 0.75) +
           scale_color_manual(values = cores_cadeia) + facet_wrap(~Parameter, scales = "free_y") +
           theme_bw(base_size = 11) + theme(legend.position = "bottom") +
           labs(title = paste("Traceplots gamma (", model_type, ")"), x = "Iteracao (pos-burnin)", y = "Valor"),
         width = 10, height = max(5, 3 * ceiling(length(params_gamma) / 3)))

  # ── Retorno ────────────────────────────────────────────────────────────────
  tibble(
    model             = model_type,
    niter             = niter, nburnin = nburnin, thin = thin,
    WAIC              = waic,  LPML    = LPML,
    ESS_beta0t_min    = min(acf_b0t$ESS,   na.rm = TRUE),
    ESS_beta0t_max    = max(acf_b0t$ESS,   na.rm = TRUE),
    Rhat_beta0t_max   = max(acf_b0t$Rhat,  na.rm = TRUE),
    ESS_beta_min      = min(acf_struct$ESS[acf_struct$Parameter %in% beta_names],  na.rm = TRUE),
    ESS_gamma_min     = min(acf_struct$ESS[acf_struct$Parameter %in% gamma_names], na.rm = TRUE),
    ESS_tau           = ESS_tau,
    ESS_global_min    = min(c(acf_b0t$ESS, ESS_struct), na.rm = TRUE),
    Rhat_max          = max(c(acf_b0t$Rhat, Rhat_struct), na.rm = TRUE),
    lag_max_0.10      = max(c(acf_b0t$lag_0.10, acf_struct$lag_0.10), na.rm = TRUE),
    lag_max_0.05      = max(c(acf_b0t$lag_0.05, acf_struct$lag_0.05), na.rm = TRUE)
  )
}

# ==============================================================================
# 6. PARALELO
# ==============================================================================
model_types <- c("spatial", "non_spatial")
n_cores     <- min(length(model_types), parallel::detectCores() - 1)
if (n_cores < 1) n_cores <- 1

output_dir <- file.path(
  "C:/Users/vlara/OneDrive/Estatistica UFMG/Mestrado/Pesquisa/Aplicação",
  "resultados_intercepto_temporal"
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
# 7. CONSOLIDAÇÃO
# ==============================================================================
resumo <- bind_rows(resultados)
write_csv(resumo, file.path(output_dir, "resumo_comparativo.csv"))
cat("\n=== RESUMO ===\n"); print(resumo)

# Comparativo de mu
all_mu <- lapply(model_types, function(m) {
  path <- file.path(output_dir, m, "mu_selected.csv")
  if (file.exists(path)) read_csv(path, show_col_types = FALSE) else NULL
}) |> bind_rows()

if (nrow(all_mu) > 0) {
  ggsave(file.path(output_dir, "mu_comparativo.png"),
    ggplot(all_mu, aes(x = Time, y = Mean, color = model, fill = model)) +
      geom_ribbon(aes(ymin = Lower, ymax = Upper), alpha = 0.15, color = NA) +
      geom_line(linewidth = 0.8) +
      facet_wrap(~Region, scales = "free_y", ncol = 3) +
      theme_bw(base_size = 12) + theme(legend.position = "bottom") +
      labs(title = "Comparacao mu: espacial vs. nao-espacial",
           x = "Tempo", color = "Modelo", fill = "Modelo"),
    width = 14, height = 10, dpi = 300)
}

# Comparativo de beta0t entre modelos
all_b0t <- lapply(model_types, function(m) {
  path <- file.path(output_dir, m, "beta0t_summary.csv")
  if (!file.exists(path)) return(NULL)
  read_csv(path, show_col_types = FALSE) %>%
    mutate(Tempo = seq_len(n()), model = m)
}) |> bind_rows()

if (nrow(all_b0t) > 0) {
  anos_label <- colnames(data_nimble$Y)
  if (!is.null(anos_label)) {
    all_b0t <- all_b0t %>% mutate(Ano = rep(anos_label, length(model_types)))
  }
  
  ggsave(file.path(output_dir, "beta0t_comparativo.png"),
         ggplot(all_b0t, aes(x = Tempo, y = Mean, color = model, fill = model)) +
           geom_ribbon(aes(ymin = HPD_Lower, ymax = HPD_Upper), alpha = 0.15, color = NA) +
           geom_line(linewidth = 0.9) +
           geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
           scale_x_continuous(breaks = seq_len(length(anos_label)), labels = anos_label) +
           theme_bw(base_size = 12) +
           theme(axis.text.x = element_text(angle = 45, hjust = 1),
                 legend.position = "bottom") +
           labs(title = "Comparacao de beta0t: espacial vs. nao-espacial",
                subtitle = "Media posterior e HPD 95%",
                x = "Ano", y = expression(beta[0][t]),
                color = "Modelo", fill = "Modelo"),
         width = 10, height = 5, dpi = 300)
}

cat("\nTempo total:\n"); print(Sys.time() - inicio_global)
