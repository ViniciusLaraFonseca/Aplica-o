# ==============================================================================
# ajuste_dinamico_BH_ARACUAI.R
# Modelo Bayesiano Dinâmico (espaço-estados) — séries univariadas
# Uma microrregião por ajuste: BH e Araçuaí em paralelo
#
# Diferenças em relação a ajuste_gamma_temporal.R:
#   - N_regions = 1 (série temporal única por ajuste)
#   - Sem efeito espacial (sem s[i], sigma_s, ICAR)
#   - Sem estrutura de clusters: epsilon[t] ~ dunif(a_eps, b_eps)
#       BH     : epsilon[t] ~ dunif(0.95, 1.00)
#       Araçuaí: epsilon[t] ~ dunif(0.75, 0.85)
#   - Preditor: log(mu[t]) = log(lambda[t]) + log(E[t]) + log(epsilon[t])
#                            + beta' x[t]
#   - FFBS simplificado: g_t = E[t] * epsilon[t] * exp(beta' x[t])
#   - slice em epsilon (23 escalares independentes)
# ==============================================================================

inicio_global <- Sys.time()

setwd("C:/Users/vlara/OneDrive/Estatistica UFMG/Mestrado/Pesquisa/Aplicação/main")

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
# 1. DADOS — via script de covariáveis
# ==============================================================================

source("Covariaveis.R")

# ── BH ────────────────────────────────────────────────────────────────────────
grouped_consultas_BH    <- grouped_consultas_MG   %>% filter(MICRO_ == "BELO_HORIZONTE_NOVA_LIMA_CAETE")
grouped_instrucao_BH    <- grouped_instrucao_MG   %>% filter(MICRO_ == "BELO_HORIZONTE_NOVA_LIMA_CAETE")
grouped_subnutridos_BH  <- grouped_subnutridos_MG %>% filter(MICRO_ == "BELO_HORIZONTE_NOVA_LIMA_CAETE")
grouped_E_BH            <- Grouped_total_MG       %>% filter(MICRO_ == "BELO_HORIZONTE_NOVA_LIMA_CAETE")
grouped_Y_BH            <- Grouped_Y              %>% filter(MICRO_ == "BELO_HORIZONTE_NOVA_LIMA_CAETE")

# ── Araçuaí ───────────────────────────────────────────────────────────────────
grouped_consultas_ARACUAI   <- grouped_consultas_MG   %>% filter(MICRO_ == "ARACUAI")
grouped_instrucao_ARACUAI   <- grouped_instrucao_MG   %>% filter(MICRO_ == "ARACUAI")
grouped_subnutridos_ARACUAI <- grouped_subnutridos_MG %>% filter(MICRO_ == "ARACUAI")
grouped_E_ARACUAI           <- Grouped_total_MG       %>% filter(MICRO_ == "ARACUAI")
grouped_Y_ARACUAI           <- Grouped_Y              %>% filter(MICRO_ == "ARACUAI")

# ── Monta arrays x[1, T, 3], vetores E[T] e Y[T] por região ──────────────────
make_data_region <- function(cons_df, inst_df, subn_df, E_df, Y_df) {
  anos  <- colnames(cons_df)[-1]           # colunas de ano (excl. MICRO_)
  n_t   <- length(anos)
  x_mat <- array(NA_real_, dim = c(1L, n_t, 2L),
                 dimnames = list(NULL, anos,
                                 c("prenatal","baixo_peso")))
  
  x_mat[1, , 1] <- as.numeric(cons_df[1, -1])
  x_mat[1, , 2] <- as.numeric(subn_df[1, -1])
  E_vec <- as.numeric(E_df[1, -1])
  Y_vec <- as.numeric(Y_df[1, -1])
  list(x = x_mat, E = matrix(E_vec, nrow = 1), Y = matrix(Y_vec, nrow = 1),
       n_times = n_t)
}

data_BH     <- make_data_region(grouped_consultas_BH,    grouped_instrucao_BH,
                                grouped_subnutridos_BH,  grouped_E_BH,    grouped_Y_BH)
data_ARACUAI <- make_data_region(grouped_consultas_ARACUAI, grouped_instrucao_ARACUAI,
                                 grouped_subnutridos_ARACUAI, grouped_E_ARACUAI, grouped_Y_ARACUAI)

p       <- 2L
n_times <- data_BH$n_times        # 23
stopifnot(data_ARACUAI$n_times == n_times)

# Normaliza E pelo próprio valor médio (análogo ao ajuste original)
norm_E <- function(E_mat) E_mat / mean(E_mat)
data_BH$E     <- norm_E(data_BH$E)
data_ARACUAI$E <- norm_E(data_ARACUAI$E)

# ==============================================================================
# 2. CONFIGURAÇÃO POR REGIÃO
# ==============================================================================

# Limites de epsilon por região
eps_bounds <- list(
  BH      = c(a_eps = 0.95, b_eps = 1.00),
  ARACUAI = c(a_eps = 0.75, b_eps = 0.85)
)

make_constants <- function(bounds) {
  list(
    n_times  = n_times,
    p        = p,
    mu_beta  = rep(0, p),
    a0       = 1.0, b0 = 1.0, w = 0.9,
    a_eps    = bounds["a_eps"],
    b_eps    = bounds["b_eps"]
  )
}

constants_BH      <- make_constants(eps_bounds$BH)
constants_ARACUAI <- make_constants(eps_bounds$ARACUAI)

# Dados nimble (listas separadas por região)
nimble_data_BH      <- list(Y = data_BH$Y,      E = data_BH$E,      x = data_BH$x)
nimble_data_ARACUAI <- list(Y = data_ARACUAI$Y, E = data_ARACUAI$E, x = data_ARACUAI$x)

# ==============================================================================
# 3. INICIALIZAÇÕES
# ==============================================================================

make_inits <- function(bounds, seed1 = 123, seed2 = 456) {
  set.seed(seed1)
  eps1 <- runif(n_times, bounds["a_eps"], bounds["b_eps"])
  set.seed(seed2)
  eps2 <- runif(n_times, bounds["a_eps"], bounds["b_eps"])
  list(
    list(beta = rep(0, p),        epsilon = eps1, lambda = rep(1.0, n_times)),
    list(beta = rnorm(p, 0, 0.3), epsilon = eps2, lambda = rgamma(n_times, 1, 1))
  )
}

inits_BH      <- make_inits(eps_bounds$BH)
inits_ARACUAI <- make_inits(eps_bounds$ARACUAI)

# ==============================================================================
# 4. FUNÇÃO WORKER
# ==============================================================================

run_model_univariate <- function(region_id, output_dir_base) {
  
  library(nimble); library(coda); library(dplyr)
  library(ggplot2); library(readr); library(stringr); library(tibble)
  
  # ── 4a. Recupera objetos exportados para o worker ────────────────────────────
  constants  <- get(paste0("constants_",  region_id))
  data_nim   <- get(paste0("nimble_data_", region_id))
  inits_list <- get(paste0("inits_",      region_id))
  
  n_times <- constants$n_times
  p       <- constants$p
  
  # ── 4b. Código NIMBLE ────────────────────────────────────────────────────────
  model_code <- nimbleCode({
    
    # Priori beta
    for (j in 1:p) beta[j] ~ dnorm(mu_beta[j], sd = 1)
    
    # epsilon[t]: probabilidade de notificação — priori uniforme por região
    for (t in 1:n_times) {
      epsilon[t] ~ dunif(a_eps, b_eps)
    }
    
    # lambda[t]: intercepto temporal — priori vaga (amostrado via FFBS)
    for (t in 1:n_times) lambda[t] ~ dgamma(1, 1)
    
    # Verossimilhança (N_regions = 1, índice i omitido)
    for (t in 1:n_times) {
      log(mu[t]) <- log(lambda[t]) + log(E[1, t]) + log(epsilon[t]) +
        inprod(beta[1:p], x[1, t, 1:p])
      Y[1, t]        ~ dpois(mu[t])
      logLik_Y[1, t] <- dpois(Y[1, t], mu[t], log = TRUE)
    }
  })
  
  # ── 4c. FFBS univariado ──────────────────────────────────────────────────────
  # g_t = E[t] * epsilon[t] * exp(beta' x[t])
  ffbs_univariate <- nimbleFunction(
    contains = sampler_BASE,
    setup = function(model, mvSaved, target, control) {
      n_times <- control$n_times
      p       <- control$p
      a0      <- control$a0
      b0      <- control$b0
      w       <- control$w
      at_buf  <- nimNumeric(n_times + 1, 0)
      bt_buf  <- nimNumeric(n_times + 1, 0)
      calcNodes   <- model$getDependencies(target, self = FALSE)
      targetNodes <- model$expandNodeNames(target)
      setupOutputs(at_buf, bt_buf)
    },
    run = function() {
      declare(t,      integer())
      declare(t_idx,  integer())
      declare(t_back, integer())
      declare(k,      integer())
      declare(prod_val, double())
      declare(g_t,    double())
      declare(nu,     double())
      
      at_buf[1] <<- a0
      bt_buf[1] <<- b0
      
      # Passagem para frente (forward)
      for (t in 1:n_times) {
        prod_val <- 0
        for (k in 1:p) prod_val <- prod_val + model$x[1, t, k] * model$beta[k]
        g_t <- model$E[1, t] * model$epsilon[t] * exp(prod_val)
        
        at_buf[t + 1] <<- w * at_buf[t] + model$Y[1, t]
        bt_buf[t + 1] <<- w * bt_buf[t] + g_t
      }
      
      # Amostragem para trás (backward sampling)
      model$lambda[n_times] <<- rgamma(1,
                                       shape = at_buf[n_times + 1],
                                       rate  = bt_buf[n_times + 1])
      for (t_idx in 1:(n_times - 1)) {
        t_back <- n_times - t_idx
        nu <- rgamma(1,
                     shape = (1 - w) * at_buf[t_back + 1],
                     rate  = bt_buf[t_back + 1])
        model$lambda[t_back] <<- nu + w * model$lambda[t_back + 1]
      }
      
      model$calculate(calcNodes)
      copy(from = model, to = mvSaved, row = 1,
           nodes = targetNodes, logProb = TRUE)
    },
    methods = list(reset = function() {})
  )
  
  # ── 4d. Construir e compilar ─────────────────────────────────────────────────
  cat("\n=== Iniciando região:", region_id, "===\n")
  scenario_dir <- file.path(output_dir_base, region_id)
  dir.create(scenario_dir, recursive = TRUE, showWarnings = FALSE)
  
  model  <- nimbleModel(code = model_code, constants = constants,
                        data = data_nim, inits = inits_list[[1]], check = FALSE)
  Cmodel <- compileNimble(model)
  
  conf <- configureMCMC(model)
  
  # FFBS para lambda
  conf$removeSamplers("lambda")
  conf$addSampler(
    target  = "lambda",
    type    = ffbs_univariate,
    control = list(n_times = n_times, p = p,
                   a0 = constants$a0, b0 = constants$b0, w = constants$w)
  )
  
  # slice escalar por epsilon[t]
  conf$removeSamplers("epsilon")
  for (t_idx in seq_len(n_times)) {
    conf$addSampler(
      target = paste0("epsilon[", t_idx, "]"),
      type   = "slice"
    )
  }
  
  conf$addMonitors(c("beta", "epsilon", "lambda", "logLik_Y"))
  conf$printSamplers()
  
  Rmcmc <- buildMCMC(conf)
  Cmcmc <- compileNimble(Rmcmc, project = model)
  
  # ── 4e. MCMC ─────────────────────────────────────────────────────────────────
  niter   <- 50000; nburnin <- 10000; nchains <- 2; thin <- 10
  
  cat(sprintf("[%s] niter=%d | nburnin=%d | thin=%d | cadeias=%d\n",
              region_id, niter, nburnin, thin, nchains))
  
  samples <- runMCMC(
    Cmcmc, niter = niter, nburnin = nburnin, nchains = nchains, thin = thin,
    inits = inits_list, samplesAsCodaMCMC = TRUE, summary = FALSE, WAIC = FALSE
  )
  saveRDS(samples, file.path(scenario_dir, "samples.rds"))
  cat("[", region_id, "] Amostras salvas.\n")
  
  samples_mat    <- as.matrix(samples)
  mcmc_list_full <- mcmc.list(lapply(seq_len(nchains),
                                     function(ch) as.mcmc(samples[[ch]])))
  rm(samples); gc()
  
  # ── 4f. Funções auxiliares ───────────────────────────────────────────────────
  safe_hpd <- function(sv) {
    if (var(sv) < 1e-12) return(c(NA_real_, NA_real_))
    as.numeric(HPDinterval(as.mcmc(sv), prob = 0.95))
  }
  safe_gelman <- function(obj) {
    tryCatch(gelman.diag(obj, autoburnin = FALSE)$psrf[, 1],
             error = function(e) rep(NA_real_, nvar(obj)))
  }
  
  beta_names    <- paste0("beta[",    seq_len(p),       "]")
  lambda_names  <- paste0("lambda[",  seq_len(n_times), "]")
  epsilon_names <- paste0("epsilon[", seq_len(n_times), "]")
  
  # ── 4g. Sumário de epsilon[t] ────────────────────────────────────────────────
  epsilon_summary <- do.call(rbind, lapply(seq_len(n_times), function(t) {
    nm  <- epsilon_names[t]
    sv  <- samples_mat[, nm]; hpd <- safe_hpd(sv)
    tibble(Time = t, Parameter = nm,
           Mean  = mean(sv), SD = sd(sv),
           Lower = hpd[1],  Upper = hpd[2],
           ESS   = as.numeric(effectiveSize(mcmc_list_full[, nm])),
           Rhat  = safe_gelman(mcmc_list_full[, nm]))
  }))
  write_csv(epsilon_summary, file.path(scenario_dir, "epsilon_summary.csv"))
  
  ggsave(
    file.path(scenario_dir, "painel_epsilon.png"),
    ggplot(epsilon_summary, aes(x = Time)) +
      geom_ribbon(aes(ymin = Lower, ymax = Upper), fill = "darkorange", alpha = 0.3) +
      geom_line(aes(y = Mean), color = "darkorange", linewidth = 0.8) +
      theme_bw(base_size = 12) +
      labs(title = paste("Trajetória de epsilon[t] —", region_id),
           x = "Tempo", y = expression(epsilon[t])),
    width = 9, height = 5
  )
  
  # ── 4h. Sumário de beta ───────────────────────────────────────────────────────
  beta_summary <- do.call(rbind, lapply(beta_names, function(nm) {
    sv  <- samples_mat[, nm]; hpd <- safe_hpd(sv)
    tibble(Parameter = nm, Mean = mean(sv), SD = sd(sv),
           HPD_Lower = hpd[1], HPD_Upper = hpd[2],
           ESS  = as.numeric(effectiveSize(mcmc_list_full[, nm])),
           Rhat = safe_gelman(mcmc_list_full[, nm]))
  }))
  write_csv(beta_summary, file.path(scenario_dir, "beta_summary.csv"))
  
  # ── 4i. Sumário de lambda[t] ──────────────────────────────────────────────────
  lambda_summary <- do.call(rbind, lapply(seq_len(n_times), function(t) {
    nm <- lambda_names[t]; sv <- samples_mat[, nm]; hpd <- safe_hpd(sv)
    tibble(Time = t, Mean = mean(sv), SD = sd(sv),
           Lower = hpd[1], Upper = hpd[2],
           ESS  = as.numeric(effectiveSize(mcmc_list_full[, nm])),
           Rhat = safe_gelman(mcmc_list_full[, nm]),
           region = region_id)
  }))
  write_csv(lambda_summary, file.path(scenario_dir, "lambda_summary.csv"))
  
  ggsave(
    file.path(scenario_dir, "painel_lambda.png"),
    ggplot(lambda_summary, aes(x = Time)) +
      geom_ribbon(aes(ymin = Lower, ymax = Upper), fill = "grey70", alpha = 0.5) +
      geom_line(aes(y = Mean), color = "black", linewidth = 0.8) +
      theme_bw(base_size = 12) +
      labs(title = paste("Intercepto temporal lambda[t] —", region_id),
           x = "Tempo", y = expression(lambda[t])),
    width = 9, height = 5
  )
  
  # ── 4j. theta[t] e mu[t] ─────────────────────────────────────────────────────
  beta_cols <- grep("^beta\\[", colnames(samples_mat), value = TRUE)
  n_draw    <- nrow(samples_mat)
  
  # x[1, t, ] para a única região
  x_reg <- matrix(NA_real_, nrow = n_times, ncol = p)
  for (t in seq_len(n_times))
    x_reg[t, ] <- data_nim$x[1, t, ]
  
  theta_summary <- do.call(rbind, lapply(seq_len(n_times), function(t) {
    ldraws <- samples_mat[, lambda_names[t]]
    bdraws <- samples_mat[, beta_cols, drop = FALSE]
    theta  <- ldraws * exp(as.vector(bdraws %*% x_reg[t, ]))
    hpd    <- safe_hpd(theta)
    tibble(Time = t, Mean = mean(theta), Lower = hpd[1], Upper = hpd[2],
           region = region_id)
  }))
  write_csv(theta_summary, file.path(scenario_dir, "theta_summary.csv"))
  
  ggsave(
    file.path(scenario_dir, "painel_theta.png"),
    ggplot(theta_summary, aes(x = Time)) +
      geom_ribbon(aes(ymin = Lower, ymax = Upper), fill = "steelblue", alpha = 0.3) +
      geom_line(aes(y = Mean), color = "steelblue", linewidth = 0.8) +
      theme_bw(base_size = 12) +
      labs(title = paste("Risco relativo theta[t] —", region_id),
           x = "Tempo", y = expression(theta[t])),
    width = 9, height = 5
  )
  
  mu_summary <- do.call(rbind, lapply(seq_len(n_times), function(t) {
    ldraws   <- samples_mat[, lambda_names[t]]
    bdraws   <- samples_mat[, beta_cols, drop = FALSE]
    eps_t    <- samples_mat[, epsilon_names[t]]
    mu_draws <- ldraws * exp(as.vector(bdraws %*% x_reg[t, ])) *
      data_nim$E[1, t] * eps_t
    hpd <- safe_hpd(mu_draws)
    tibble(Time = t, Mean = mean(mu_draws), Lower = hpd[1], Upper = hpd[2],
           Y_obs = data_nim$Y[1, t], region = region_id)
  }))
  write_csv(mu_summary, file.path(scenario_dir, "mu_summary.csv"))
  
  ggsave(
    file.path(scenario_dir, "painel_mu.png"),
    ggplot(mu_summary, aes(x = Time)) +
      geom_ribbon(aes(ymin = Lower, ymax = Upper), fill = "darkorange", alpha = 0.3) +
      geom_line(aes(y = Mean), color = "darkorange", linewidth = 0.8) +
      geom_point(aes(y = Y_obs), shape = 21, fill = "white",
                 color = "black", size = 2) +
      theme_bw(base_size = 12) +
      labs(title    = paste("Contagem esperada mu[t] —", region_id),
           subtitle = "Pontos = Y observado",
           x = "Tempo", y = expression(mu[t])),
    width = 9, height = 5
  )
  
  # ── 4k. WAIC e LPML ───────────────────────────────────────────────────────────
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
    cat(sprintf("[%s] WAIC = %.2f | LPML = %.2f\n", region_id, waic, LPML))
  }
  
  # ── 4l. Diagnósticos ACF ──────────────────────────────────────────────────────
  params_diag <- c(beta_names, lambda_names, epsilon_names)
  
  acf_results <- do.call(rbind, lapply(params_diag, function(nm) {
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
  
  # ACF de beta
  acf_beta_df <- do.call(rbind, lapply(beta_names, function(nm) {
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
      labs(title    = paste("ACF de beta —", region_id),
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
      labs(title    = paste("ACF de lambda[t] —", region_id),
           subtitle = "Azul: |0.10| | Vermelho: |0.05|"),
    width = 14, height = 10
  )
  
  # ACF de epsilon[t]
  acf_eps_df <- do.call(rbind, lapply(epsilon_names, function(nm) {
    ac <- acf(samples_mat[, nm], lag.max = 100, plot = FALSE)
    tibble(Time = as.integer(str_extract(nm, "\\d+")),
           Lag  = as.vector(ac$lag[-1]), ACF = as.vector(ac$acf[-1]))
  }))
  ggsave(
    file.path(scenario_dir, "acf_epsilon.png"),
    ggplot(acf_eps_df, aes(x = Lag, y = ACF)) +
      geom_col(width = 0.6, fill = "darkorange") +
      geom_hline(yintercept = c(-0.10, 0.10), linetype = "dashed",
                 color = "blue", linewidth = 0.5) +
      geom_hline(yintercept = c(-0.05, 0.05), linetype = "dotted",
                 color = "red",  linewidth = 0.5) +
      facet_wrap(~ Time, scales = "free_y", ncol = 6,
                 labeller = label_bquote(epsilon[.(Time)])) +
      theme_bw(base_size = 9) +
      labs(title    = paste("ACF de epsilon[t] —", region_id),
           subtitle = "Azul: |0.10| | Vermelho: |0.05|"),
    width = 14, height = 10
  )
  
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
      geom_line(aes(y = Value),    alpha = 0.25, linewidth = 0.20) +
      geom_line(aes(y = ErgMedia), alpha = 0.90, linewidth = 0.75) +
      scale_color_manual(values = cores_cadeia) +
      facet_wrap(~ Parameter, scales = "free_y") +
      theme_bw(base_size = 11) + theme(legend.position = "bottom") +
      labs(title    = paste("Traceplots beta —", region_id),
           subtitle = "Linha grossa = média ergódica | Linha fina = cadeia",
           x = "Iteração (pós-burnin)", y = "Valor"),
    width = 10, height = max(5, 3 * ceiling(length(beta_names) / 3))
  )
  
  # ── 4m. Retorno resumido ──────────────────────────────────────────────────────
  tibble(
    region  = region_id, niter = niter, nburnin = nburnin, thin = thin,
    a_eps   = constants$a_eps, b_eps = constants$b_eps,
    WAIC    = waic, LPML = LPML,
    ESS_beta_min    = min(beta_summary$ESS,    na.rm = TRUE),
    ESS_epsilon_min = min(epsilon_summary$ESS, na.rm = TRUE),
    ESS_lambda_min  = min(lambda_summary$ESS,  na.rm = TRUE),
    Rhat_beta_max    = max(beta_summary$Rhat,    na.rm = TRUE),
    Rhat_epsilon_max = max(epsilon_summary$Rhat, na.rm = TRUE),
    Rhat_lambda_max  = max(lambda_summary$Rhat,  na.rm = TRUE),
    lag_max_0.10 = max(acf_results$lag_0.10, na.rm = TRUE),
    lag_max_0.05 = max(acf_results$lag_0.05, na.rm = TRUE)
  )
}

# ==============================================================================
# 5. EXECUÇÃO PARALELA
# ==============================================================================
regions <- c("BH", "ARACUAI")
n_cores <- min(length(regions), parallel::detectCores() - 1)
if (n_cores < 1) n_cores <- 1

output_dir <- file.path(
  "C:/Users/vlara/OneDrive/Estatistica UFMG/Mestrado/Pesquisa/Aplicação",
  "resultados_dinamico_BH_ARACUAI_sem_instrucao"
)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

cl <- makeCluster(n_cores)
clusterExport(cl, c(
  "constants_BH", "constants_ARACUAI",
  "nimble_data_BH", "nimble_data_ARACUAI",
  "inits_BH", "inits_ARACUAI",
  "n_times", "p",
  "run_model_univariate", "output_dir"
))
clusterEvalQ(cl, {
  library(nimble); library(coda); library(dplyr)
  library(ggplot2); library(readr); library(stringr); library(tibble)
  Sys.setenv(OMP_NUM_THREADS = "1"); Sys.setenv(MKL_NUM_THREADS = "1")
  if (requireNamespace("RhpcBLASctl", quietly = TRUE))
    RhpcBLASctl::blas_set_num_threads(1)
})

resultados <- parLapply(cl, regions,
                        function(r) run_model_univariate(r, output_dir))
stopCluster(cl)

# ==============================================================================
# 6. CONSOLIDAÇÃO E GRÁFICOS COMPARATIVOS
# ==============================================================================
resumo <- bind_rows(resultados)
write_csv(resumo, file.path(output_dir, "resumo_comparativo.csv"))
cat("\n=== RESUMO COMPARATIVO ===\n"); print(resumo)

# ── Comparativo lambda[t] ──────────────────────────────────────────────────────
lambda_all <- lapply(regions, function(r) {
  path <- file.path(output_dir, r, "lambda_summary.csv")
  if (!file.exists(path)) return(NULL)
  read_csv(path, show_col_types = FALSE)
}) |> bind_rows()

if (nrow(lambda_all) > 0) {
  ggsave(
    file.path(output_dir, "lambda_comparativo.png"),
    ggplot(lambda_all, aes(x = Time, y = Mean, color = region, fill = region)) +
      geom_ribbon(aes(ymin = Lower, ymax = Upper), alpha = 0.15, color = NA) +
      geom_line(linewidth = 0.9) +
      theme_bw(base_size = 12) + theme(legend.position = "bottom") +
      labs(title = "Comparação lambda[t]: BH vs. Araçuaí",
           x = "Tempo", y = expression(lambda[t]),
           color = "Região", fill = "Região"),
    width = 10, height = 5, dpi = 300
  )
}

# ── Comparativo epsilon[t] ─────────────────────────────────────────────────────
epsilon_all <- lapply(regions, function(r) {
  path <- file.path(output_dir, r, "epsilon_summary.csv")
  if (!file.exists(path)) return(NULL)
  read_csv(path, show_col_types = FALSE) %>% mutate(region = r)
}) |> bind_rows()

if (nrow(epsilon_all) > 0) {
  ggsave(
    file.path(output_dir, "epsilon_comparativo.png"),
    ggplot(epsilon_all, aes(x = Time, y = Mean, color = region, fill = region)) +
      geom_ribbon(aes(ymin = Lower, ymax = Upper), alpha = 0.15, color = NA) +
      geom_line(linewidth = 0.9) +
      theme_bw(base_size = 12) + theme(legend.position = "bottom") +
      labs(title    = "Comparação epsilon[t]: BH vs. Araçuaí",
           subtitle = "Suportes distintos: BH ~ U(0.95,1) | Araçuaí ~ U(0.75,0.85)",
           x = "Tempo", y = expression(epsilon[t]),
           color = "Região", fill = "Região"),
    width = 10, height = 5, dpi = 300
  )
}

# ── Comparativo theta[t] e mu[t] ──────────────────────────────────────────────
for (tipo in c("theta", "mu")) {
  all_df <- lapply(regions, function(r) {
    path <- file.path(output_dir, r, paste0(tipo, "_summary.csv"))
    if (!file.exists(path)) return(NULL)
    read_csv(path, show_col_types = FALSE)
  }) |> bind_rows()
  
  if (nrow(all_df) > 0) {
    p_comp <- ggplot(all_df, aes(x = Time, y = Mean,
                                 color = region, fill = region)) +
      geom_ribbon(aes(ymin = Lower, ymax = Upper), alpha = 0.15, color = NA) +
      geom_line(linewidth = 0.9) +
      theme_bw(base_size = 12) + theme(legend.position = "bottom") +
      labs(title = sprintf("Comparação %s[t]: BH vs. Araçuaí", tipo),
           x = "Tempo", y = tipo,
           color = "Região", fill = "Região")
    
    # Adiciona pontos observados no gráfico de mu
    if (tipo == "mu" && "Y_obs" %in% colnames(all_df)) {
      p_comp <- p_comp +
        geom_point(aes(y = Y_obs), shape = 21, fill = "white",
                   color = "black", size = 1.5, show.legend = FALSE)
    }
    ggsave(file.path(output_dir, paste0(tipo, "_comparativo.png")),
           p_comp, width = 10, height = 5, dpi = 300)
  }
}

cat("\nTempo total de execução:\n")
print(Sys.time() - inicio_global)