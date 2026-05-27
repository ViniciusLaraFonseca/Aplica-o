# ==============================================================================
# ajuste_delta_temporal_gamma_t.R
# Modelo Poisson com intercepto temporal multiplicativo delta[t] ~ Gamma(1,1)
# e gamma variável no tempo.
#
# PREDITOR:
#   mu[i,t] = E[i,t] * epsilon[i,t] * delta[t] * exp(beta' x[i,t])  [* exp(s[i])]
#
#   log(mu[i,t]) = log(delta[t]) + log(E[i,t]) + log(epsilon[i,t])
#                  + beta' x[i,t]  [+ s[i]  se espacial]
#
# COMPONENTES:
#   delta[t]      : efeito temporal multiplicativo i.i.d. ~ Gamma(1, 1)
#   epsilon[i,t]  : 1 - sum_k h[i,k] * gamma[k,t]  (gamma VARIÁVEL no tempo)
#   gamma[k,t]    : matriz K × T
#
# PRIORIS:
#   delta[t]   ~ Gamma(1, 1)         i.i.d. para cada t = 1,...,T
#   beta[j]    ~ N(0, 1)
#   gamma[1,t] ~ Unif(0, 0.1)        para cada t
#   gamma[k,t] ~ Unif(0, 1 - sum_{j<k} gamma[j,t])  para k = 2,...,K, cada t
#   sigma_s    ~ Half-t(0, 1, 1)     [apenas espacial]
#   tau_s      = 1 / sigma_s^2
#   s          ~ ICAR(tau_s)          [apenas espacial]
#
# AMOSTRADOR DE delta[t] — GIBBS CONJUGADO EXATO:
#   Dado o restante, delta[t] é conjugado:
#     delta[t] | resto ~ Gamma(a_t, b_t)
#   onde:
#     a_t = 1 + sum_i Y[i,t]
#     b_t = 1 + sum_i g[i,t]
#     g[i,t] = E[i,t] * epsilon[i,t] * exp(beta' x[i,t])  [* exp(s[i])]
#
#   Um único amostrador (nimbleFunction) itera sobre todos os T tempos,
#   amostrando cada delta[t] diretamente da distribuição Gamma posterior —
#   sem rejeição, sem MH, sem estado latente temporal.
#
# INSPIRADO em M5 (ajuste_intercepto_temporal_gamma_t.R):
#   Mantém estrutura de código, diagnósticos, saídas e gráficos.
#   Substituição: beta0t[t] ~ N(0,1)  →  delta[t] ~ Gamma(1,1)
#   Predictor:    beta0t[t] + ...      →  log(delta[t]) + ...
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
# 4. INICIALIZAÇÕES  (delta[t] ~ Gamma(1,1): init = 1.0 é o valor esperado)
# ==============================================================================
set.seed(123)

# gamma[k,t]: mesmo esquema de M5 — coluna a coluna via função auxiliar
gerar_gamma_coluna <- function(K) {
  g    <- numeric(K)
  g[1] <- runif(1, 0.02, 0.08)
  for (k in 2:K) {
    max_val <- 1 - sum(g[1:(k - 1)])
    g[k]    <- runif(1, 0, min(0.2, max_val))
  }
  g
}

gamma_mat1 <- replicate(n_times, gerar_gamma_coluna(K))  # K × T
gamma_mat2 <- replicate(n_times, gerar_gamma_coluna(K))

inits_spatial_1 <- list(
  delta   = rep(1.0, n_times),          # E[delta[t]] sob Gamma(1,1)
  beta    = rep(0, p),
  gamma   = gamma_mat1,
  sigma_s = 0.5,
  s       = rep(0, N_regions)
)
inits_spatial_2 <- list(
  delta   = rgamma(n_times, 1, 1),      # amostras iniciais da priori
  beta    = rnorm(p, 0, 0.3),
  gamma   = gamma_mat2,
  sigma_s = 1.0,
  s       = rep(0, N_regions)
)
inits_nonspatial_1 <- list(
  delta = rep(1.0, n_times),
  beta  = rep(0, p),
  gamma = gamma_mat1
)
inits_nonspatial_2 <- list(
  delta = rgamma(n_times, 1, 1),
  beta  = rnorm(p, 0, 0.3),
  gamma = gamma_mat2
)

inits_list_spatial    <- list(inits_spatial_1,    inits_spatial_2)
inits_list_nonspatial <- list(inits_nonspatial_1, inits_nonspatial_2)

# ==============================================================================
# 5. FUNÇÃO WORKER
# ==============================================================================
run_model <- function(model_type, output_dir) {
  
  library(nimble); library(coda); library(dplyr)
  library(ggplot2); library(readr); library(stringr); library(tibble)
  
  # ── 5a. Código NIMBLE ────────────────────────────────────────────────────────
  #
  # delta[t] entra como log(delta[t]) no preditor — estrutura idêntica ao
  # ajuste_lambda_shared.R / ajuste_gamma_temporal.R, mas sem FFBS:
  # delta[t] é i.i.d. e será amostrado via Gibbs conjugado.
  
  code_spatial <- nimbleCode({
    # Covariáveis
    for (j in 1:p) beta[j] ~ dnorm(mu_beta[j], sd = 1)
    
    # Pesos de cluster (gamma variável no tempo — estrutura de M5)
    for (t in 1:n_times) {
      gamma[1, t] ~ dunif(a_unif, b_unif)
      for (k in 2:K)
        gamma[k, t] ~ dunif(0, 1 - sum(gamma[1:(k-1), t]))
    }
    
    # Intercepto temporal multiplicativo i.i.d.
    for (t in 1:n_times) delta[t] ~ dgamma(1, 1)
    
    # Efeito espacial ICAR
    sigma_s ~ T(dt(0, 1, 1), 0, )
    tau_s   <- 1 / (sigma_s^2)
    s[1:n_regions] ~ dcar_normal(
      adj[1:n_adj], weights[1:n_adj], num[1:n_regions], tau_s, zero_mean = 1
    )
    
    # Verossimilhança
    for (i in 1:n_regions) {
      for (t in 1:n_times) {
        epsilon[i, t] <- 1 - inprod(h[i, 1:K], gamma[1:K, t])
        log(mu[i, t]) <- log(delta[t]) + log(E[i, t]) + log(epsilon[i, t]) +
          inprod(beta[1:p], x[i, t, 1:p]) + s[i]
        Y[i, t]        ~ dpois(mu[i, t])
        logLik_Y[i, t] <- dpois(Y[i, t], mu[i, t], log = TRUE)
      }
    }
  })
  
  code_nonspatial <- nimbleCode({
    for (j in 1:p) beta[j] ~ dnorm(mu_beta[j], sd = 1)
    
    for (t in 1:n_times) {
      gamma[1, t] ~ dunif(a_unif, b_unif)
      for (k in 2:K)
        gamma[k, t] ~ dunif(0, 1 - sum(gamma[1:(k-1), t]))
    }
    
    for (t in 1:n_times) delta[t] ~ dgamma(1, 1)
    
    for (i in 1:n_regions) {
      for (t in 1:n_times) {
        epsilon[i, t] <- 1 - inprod(h[i, 1:K], gamma[1:K, t])
        log(mu[i, t]) <- log(delta[t]) + log(E[i, t]) + log(epsilon[i, t]) +
          inprod(beta[1:p], x[i, t, 1:p])
        Y[i, t]        ~ dpois(mu[i, t])
        logLik_Y[i, t] <- dpois(Y[i, t], mu[i, t], log = TRUE)
      }
    }
  })
  
  # ── 5b. Amostrador Gibbs conjugado para delta — versão ESPACIAL ──────────────
  #
  # Derivação:
  #   p(delta[t] | resto) ∝ prod_i Poisson(Y[i,t] | mu[i,t]) * Gamma(delta[t] | 1,1)
  #
  # Como mu[i,t] = E[i,t] * epsilon[i,t] * delta[t] * exp(beta'x[i,t] + s[i])
  # = g[i,t] * delta[t],  temos:
  #
  #   prod_i [delta[t]^Y[i,t] * exp(-g[i,t]*delta[t])] * delta[t]^0 * exp(-delta[t])
  #   = delta[t]^(sum_i Y[i,t]) * exp(-(1 + sum_i g[i,t]) * delta[t])
  #   ~ Gamma(a_t, b_t)
  #
  # com  a_t = 1 + sum_i Y[i,t]
  #      b_t = 1 + sum_i g[i,t]  (b_t é a TAXA — parametrização shape/rate)
  #
  # O amostrador percorre todos os T tempos em uma única chamada run(),
  # amostrando cada delta[t] diretamente — exato e sem rejeição.
  
  gibbs_delta_spatial <- nimbleFunction(
    contains = sampler_BASE,
    setup = function(model, mvSaved, target, control) {
      n_regions  <- control$n_regions
      n_times    <- control$n_times
      p          <- control$p
      calcNodes  <- model$getDependencies(target, self = FALSE)
      targetNodes <- model$expandNodeNames(target)
      setupOutputs()
    },
    run = function() {
      declare(i,        integer())
      declare(t,        integer())
      declare(k,        integer())
      declare(prod_val, double())
      declare(g_it,     double())
      declare(sum_Y_t,  double())
      declare(sum_g_t,  double())
      declare(a_t,      double())
      declare(b_t,      double())
      
      for (t in 1:n_times) {
        sum_Y_t <- 0
        sum_g_t <- 0
        for (i in 1:n_regions) {
          sum_Y_t  <- sum_Y_t + model$Y[i, t]
          # g[i,t] = E[i,t] * epsilon[i,t] * exp(beta'x[i,t] + s[i])
          prod_val <- 0
          for (k in 1:p) prod_val <- prod_val + model$x[i, t, k] * model$beta[k]
          g_it    <- model$E[i, t] * model$epsilon[i, t] * exp(prod_val + model$s[i])
          sum_g_t <- sum_g_t + g_it
        }
        # Gibbs exato: delta[t] ~ Gamma(a_t, b_t)   (parametrização shape/rate)
        a_t <- 1 + sum_Y_t
        b_t <- 1 + sum_g_t
        model$delta[t] <<- rgamma(1, shape = a_t, rate = b_t)
      }
      model$calculate(calcNodes)
      copy(from = model, to = mvSaved, row = 1, nodes = targetNodes, logProb = TRUE)
    },
    methods = list(reset = function() {})
  )
  
  # ── 5c. Amostrador Gibbs conjugado para delta — versão NÃO-ESPACIAL ──────────
  # Idêntica ao espacial, sem s[i] no cálculo de g[i,t].
  
  gibbs_delta_nonspatial <- nimbleFunction(
    contains = sampler_BASE,
    setup = function(model, mvSaved, target, control) {
      n_regions   <- control$n_regions
      n_times     <- control$n_times
      p           <- control$p
      calcNodes   <- model$getDependencies(target, self = FALSE)
      targetNodes <- model$expandNodeNames(target)
      setupOutputs()
    },
    run = function() {
      declare(i,        integer())
      declare(t,        integer())
      declare(k,        integer())
      declare(prod_val, double())
      declare(g_it,     double())
      declare(sum_Y_t,  double())
      declare(sum_g_t,  double())
      declare(a_t,      double())
      declare(b_t,      double())
      
      for (t in 1:n_times) {
        sum_Y_t <- 0
        sum_g_t <- 0
        for (i in 1:n_regions) {
          sum_Y_t  <- sum_Y_t + model$Y[i, t]
          prod_val <- 0
          for (k in 1:p) prod_val <- prod_val + model$x[i, t, k] * model$beta[k]
          g_it    <- model$E[i, t] * model$epsilon[i, t] * exp(prod_val)  # sem s[i]
          sum_g_t <- sum_g_t + g_it
        }
        a_t <- 1 + sum_Y_t
        b_t <- 1 + sum_g_t
        model$delta[t] <<- rgamma(1, shape = a_t, rate = b_t)
      }
      model$calculate(calcNodes)
      copy(from = model, to = mvSaved, row = 1, nodes = targetNodes, logProb = TRUE)
    },
    methods = list(reset = function() {})
  )
  
  # ── 5d. Seleção de objetos ───────────────────────────────────────────────────
  is_spatial <- (model_type == "spatial")
  model_code <- if (is_spatial) code_spatial      else code_nonspatial
  constants  <- if (is_spatial) constants_spatial  else constants_nonspatial
  inits_list <- if (is_spatial) inits_list_spatial else inits_list_nonspatial
  gibbs_fn   <- if (is_spatial) gibbs_delta_spatial else gibbs_delta_nonspatial
  
  cat("\n=== Modelo:", model_type, "===\n")
  scenario_dir <- file.path(output_dir, model_type)
  dir.create(scenario_dir, recursive = TRUE, showWarnings = FALSE)
  
  # ── 5e. Construir e compilar ─────────────────────────────────────────────────
  model  <- nimbleModel(
    code = model_code, constants = constants,
    data = data_nimble, inits = inits_list[[1]], check = FALSE
  )
  Cmodel <- compileNimble(model)
  
  conf <- configureMCMC(model)
  
  # Substitui amostrador padrão de delta pelo Gibbs conjugado (T tempos de uma vez)
  conf$removeSamplers("delta")
  conf$addSampler(
    target  = "delta",
    type    = gibbs_fn,
    control = list(n_regions = N_regions, n_times = n_times, p = p)
  )
  
  # gamma: um AF_slice por coluna-t (K parâmetros por bloco) — idêntico a M5
  conf$removeSamplers("gamma")
  for (t_idx in seq_len(n_times)) {
    conf$addSampler(
      target = paste0("gamma[", seq_len(K), ", ", t_idx, "]"),
      type   = "AF_slice"
    )
  }
  
  monitors_base <- c("delta", "beta", "gamma", "logLik_Y")
  if (is_spatial) monitors_base <- c(monitors_base, "s", "sigma_s", "tau_s")
  conf$addMonitors(monitors_base)
  conf$addMonitors(MU_MONITORS)
  conf$printSamplers()
  
  Rmcmc <- buildMCMC(conf)
  Cmcmc <- compileNimble(Rmcmc, project = model)
  
  # ── 5f. MCMC ─────────────────────────────────────────────────────────────────
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
  
  # ── 5g. Funções auxiliares ───────────────────────────────────────────────────
  safe_hpd <- function(sv) {
    if (var(sv) < 1e-12) return(c(NA_real_, NA_real_))
    as.numeric(HPDinterval(as.mcmc(sv), prob = 0.95))
  }
  safe_gelman <- function(obj) {
    tryCatch(gelman.diag(obj, autoburnin = FALSE)$psrf[, 1],
             error = function(e) rep(NA_real_, nvar(obj)))
  }
  
  beta_names  <- paste0("beta[",  seq_len(p),       "]")
  delta_names <- paste0("delta[", seq_len(n_times), "]")
  
  # Nomes de gamma: gamma[k, t]
  gamma_names_mat <- outer(
    seq_len(K), seq_len(n_times),
    function(k, t) paste0("gamma[", k, ", ", t, "]")
  )  # K × T
  
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
  
  # ── 5j. Epsilon por cluster ───────────────────────────────────────────────────
  eps_full <- do.call(rbind, lapply(seq_len(N_regions), function(i) {
    do.call(rbind, lapply(seq_len(n_times), function(t) {
      sv  <- epsilon_draws[, i, t]; hpd <- safe_hpd(sv)
      tibble(Region = i, Cluster = cluster_ids[i], Time = t,
             Mean = mean(sv), Lower = hpd[1], Upper = hpd[2])
    }))
  }))
  write_csv(eps_full, file.path(scenario_dir, "epsilon_summary.csv"))
  
  eps_by_cluster <- do.call(rbind, lapply(seq_len(K), function(k) {
    regs_k <- which(cluster_ids == k)
    do.call(rbind, lapply(seq_len(n_times), function(t) {
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
  
  # ── 5k. Sumário e painel de delta[t] ─────────────────────────────────────────
  delta_summary <- do.call(rbind, lapply(seq_len(n_times), function(t) {
    nm <- delta_names[t]; sv <- samples_mat[, nm]; hpd <- safe_hpd(sv)
    tibble(Time = t, Parameter = nm,
           Mean = mean(sv), SD = sd(sv), Lower = hpd[1], Upper = hpd[2],
           ESS  = as.numeric(effectiveSize(mcmc_list_full[, nm])),
           Rhat = safe_gelman(mcmc_list_full[, nm]),
           model = model_type)
  }))
  write_csv(delta_summary, file.path(scenario_dir, "delta_summary.csv"))
  
  anos_label <- colnames(data_nimble$Y)
  if (is.null(anos_label)) anos_label <- seq_len(n_times)
  
  ggsave(
    file.path(scenario_dir, "painel_delta.png"),
    ggplot(delta_summary, aes(x = Time)) +
      geom_ribbon(aes(ymin = Lower, ymax = Upper), fill = "grey70", alpha = 0.5) +
      geom_line(aes(y = Mean), color = "black", linewidth = 0.8) +
      geom_hline(yintercept = 1, linetype = "dashed", color = "red",
                 linewidth = 0.5) +                         # E[delta] = 1 sob Gamma(1,1)
      scale_x_continuous(breaks = seq_len(n_times), labels = anos_label) +
      theme_bw(base_size = 12) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      labs(title    = paste("Intercepto temporal delta[t] ~ Gamma(1,1) (", model_type, ")"),
           subtitle = "Média posterior e HPD 95% | Linha vermelha: E[delta] = 1",
           x = "Ano", y = expression(delta[t])),
    width = 10, height = 5
  )
  
  # ── 5l. Sumário de beta ───────────────────────────────────────────────────────
  beta_summary <- do.call(rbind, lapply(beta_names, function(nm) {
    sv  <- samples_mat[, nm]; hpd <- safe_hpd(sv)
    tibble(Parameter = nm, Mean = mean(sv), SD = sd(sv),
           HPD_Lower = hpd[1], HPD_Upper = hpd[2],
           ESS  = as.numeric(effectiveSize(mcmc_list_full[, nm])),
           Rhat = safe_gelman(mcmc_list_full[, nm]))
  }))
  write_csv(beta_summary, file.path(scenario_dir, "beta_summary.csv"))
  
  # ── 5m. tau_s e s[i] (apenas espacial) ───────────────────────────────────────
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
      sv <- samples_mat[, s_names[i]]; hpd <- safe_hpd(sv)
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
             y = "s[i]", x = "Região"),
      width = 10, height = 5
    )
  }
  
  # ── 5n. mu para regiões selecionadas ─────────────────────────────────────────
  eps_mean_region <- eps_full |>
    group_by(Region) |>
    summarise(MeanEps = mean(Mean), .groups = "drop")
  
  make_label <- function(regs) {
    setNames(
      sprintf("Reg %d (C%d)\ne=%.3f", regs, cluster_ids[regs],
              eps_mean_region$MeanEps[match(regs, eps_mean_region$Region)]),
      as.character(regs)
    )
  }
  
  mu_summary <- do.call(rbind, lapply(REGIONS_INTEREST, function(i) {
    do.call(rbind, lapply(seq_len(n_times), function(t) {
      nm  <- paste0("mu[", i, ", ", t, "]")
      sv  <- samples_mat[, nm]; hpd <- safe_hpd(sv)
      tibble(Region = i, Time = t, Mean = mean(sv),
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
  params_struct <- c(beta_names, if (is_spatial) "tau_s" else character(0))
  all_diag      <- c(delta_names, params_struct, as.vector(gamma_names_mat))
  
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
           lag_0.10  = { v <- lags[which(abs(acfs) < 0.10)[1]]; ifelse(is.na(v), Inf, v) },
           lag_0.05  = { v <- lags[which(abs(acfs) < 0.05)[1]]; ifelse(is.na(v), Inf, v) },
           acf_lag1  = acfs[1])
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
  
  # ACF de delta[t] — painel por tempo (análogo ao painel de beta0t em M5
  #                   e ao painel de lambda em M3/M4)
  acf_delta_df <- do.call(rbind, lapply(delta_names, function(nm) {
    ac <- acf(samples_mat[, nm], lag.max = 100, plot = FALSE)
    tibble(Time = as.integer(str_extract(nm, "\\d+")),
           Lag  = as.vector(ac$lag[-1]), ACF = as.vector(ac$acf[-1]))
  }))
  ggsave(
    file.path(scenario_dir, "acf_delta.png"),
    ggplot(acf_delta_df, aes(x = Lag, y = ACF)) +
      geom_col(width = 0.6, fill = "steelblue") +
      geom_hline(yintercept = c(-0.10, 0.10), linetype = "dashed",
                 color = "blue", linewidth = 0.5) +
      geom_hline(yintercept = c(-0.05, 0.05), linetype = "dotted",
                 color = "red",  linewidth = 0.5) +
      facet_wrap(~ Time, scales = "free_y", ncol = 6,
                 labeller = label_bquote(delta[.(Time)])) +
      theme_bw(base_size = 9) +
      labs(title    = paste("ACF de delta[t] (", model_type, ")"),
           subtitle = "Gibbs conjugado — autocorrelação esperada baixa"),
    width = 14, height = 10
  )
  
  # ACF de gamma[k, t] — um arquivo por cluster k
  invisible(lapply(seq_len(K), function(k_val) {
    acf_g_df <- do.call(rbind, lapply(seq_len(n_times), function(t) {
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
           subtitle = "Linha grossa = média ergódica | Linha fina = cadeia",
           x = "Iteração (pós-burnin)", y = "Valor"),
    width = 10, height = max(5, 3 * ceiling(length(params_acf_beta) / 3))
  )
  
  # ── 5q. Retorno resumido ──────────────────────────────────────────────────────
  tibble(
    model = model_type, niter = niter, nburnin = nburnin, thin = thin,
    WAIC = waic, LPML = LPML,
    ESS_beta_min    = min(beta_summary$ESS,   na.rm = TRUE),
    ESS_gamma_min   = min(gamma_summary$ESS,  na.rm = TRUE),
    ESS_delta_min   = min(delta_summary$ESS,  na.rm = TRUE),
    ESS_delta_max   = max(delta_summary$ESS,  na.rm = TRUE),
    ESS_tau         = ESS_tau,
    Rhat_beta_max   = max(Rhat_struct[beta_names], na.rm = TRUE),
    Rhat_gamma_max  = max(gamma_summary$Rhat,  na.rm = TRUE),
    Rhat_delta_max  = max(delta_summary$Rhat,  na.rm = TRUE),
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
  "resultados_delta_temporal_gamma_t"
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

# Comparativo delta[t]
delta_all <- lapply(model_types, function(m) {
  path <- file.path(output_dir, m, "delta_summary.csv")
  if (!file.exists(path)) return(NULL)
  read_csv(path, show_col_types = FALSE)
}) |> bind_rows()

if (nrow(delta_all) > 0) {
  ggsave(
    file.path(output_dir, "delta_comparativo.png"),
    ggplot(delta_all, aes(x = Time, y = Mean, color = model, fill = model)) +
      geom_ribbon(aes(ymin = Lower, ymax = Upper), alpha = 0.15, color = NA) +
      geom_line(linewidth = 0.9) +
      geom_hline(yintercept = 1, linetype = "dashed", color = "grey40") +
      theme_bw(base_size = 12) + theme(legend.position = "bottom") +
      labs(title    = "Comparação delta[t]: espacial vs. não-espacial",
           subtitle = "Linha tracejada: E[delta] = 1 (média da priori Gamma(1,1))",
           x = "Tempo", y = expression(delta[t]),
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
      labs(title    = "Comparação gamma[k,t]: espacial vs. não-espacial (K=4)",
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
      labs(title    = "Comparação epsilon por cluster: espacial vs. não-espacial (K=4)",
           subtitle = "Banda = IC 95% HPD",
           x = "Tempo", y = expression(epsilon[kt]),
           color = "Modelo", fill = "Modelo"),
    width = 10, height = 8, dpi = 300
  )
}

# Comparativo mu
mu_all <- lapply(model_types, function(m) {
  path <- file.path(output_dir, m, "mu_selected.csv")
  if (!file.exists(path)) return(NULL)
  read_csv(path, show_col_types = FALSE)
}) |> bind_rows()

if (nrow(mu_all) > 0) {
  ggsave(
    file.path(output_dir, "mu_comparativo.png"),
    ggplot(mu_all, aes(x = Time, y = Mean, color = model, fill = model)) +
      geom_ribbon(aes(ymin = Lower, ymax = Upper), alpha = 0.15, color = NA) +
      geom_line(linewidth = 0.8) +
      facet_wrap(~ Region, scales = "free_y", ncol = 3) +
      theme_bw(base_size = 12) + theme(legend.position = "bottom") +
      labs(title = "Comparação mu: espacial vs. não-espacial",
           x = "Tempo", color = "Modelo", fill = "Modelo"),
    width = 14, height = 10, dpi = 300
  )
}

cat("\nTempo total de execução:\n")
print(Sys.time() - inicio_global)
