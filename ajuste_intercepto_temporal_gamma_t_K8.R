# ==============================================================================
# ajuste_intercepto_temporal_gamma_t_K8.R
# Modelo Poisson log-linear com intercepto temporal i.i.d. e gamma temporal
# K = 8 clusters
#
# PREDITOR:
#   log(mu[i,t]) = beta0t[t] + log(E[i,t]) + log(epsilon[i,t])
#                  + beta' x[i,t]  [+ s[i]  se espacial]
#
# COMPONENTES:
#   beta0t[t]     : intercepto temporal i.i.d. ~ N(0, 1)
#   epsilon[i,t]  : 1 - sum_k h[i,k] * gamma[k,t]   (gamma VARIÁVEL no tempo)
#   gamma[k,t]    : matriz K × T
#
# PRIORIS TEMPO-VARIANTES PARA gamma[1,t]  — HERDADO DO PADRÃO OURO:
#   t = 1..10  → Unif(0.05, 0.15)  : âncora restrita  (início da série)
#   t = 11..17 → Unif(0.00, 0.10)  : âncora moderada  (período intermediário)
#   t = 18..23 → Unif(0.00, 0.05)  : âncora difusa    (período mais recente)
#   gamma[k,t] ~ Unif(0, 1 - sum_{j<k} gamma[j,t])   k = 2..K
#
# INICIALIZAÇÃO ÓTIMA:
#   gamma[1,t] no ponto médio de [a_unif[t], b_unif[t]] para cada t;
#   gamma[2..K,t] com pesos decrescentes no simplex residual;
#   restrições a_unif[t] ≤ gamma[1,t] ≤ b_unif[t] e sum < 1 verificadas.
#
# DEMAIS PRIORIS:
#   beta0t[t]  ~ N(0, 1)
#   beta[j]    ~ N(0, 1)
#   sigma_s    ~ Half-t(0, 1, 1)   [apenas espacial]
#   tau_s      = 1/sigma_s^2
#   s          ~ ICAR(tau_s)        [apenas espacial]
#
# K = 8  (lido automaticamente de _dataCaseStudy_K8.R)
# Estrutura, diagnósticos e saídas seguem ajuste_gamma_temporal_K8.R.
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
# 1. DADOS
# ==============================================================================

source("_dataCaseStudy_K8.R")   # Y_mat, E, x, data$adj/num/sumNumNeigh/hAI (K=8)

stopifnot(
  "Y_mat deve ter 75 linhas"       = nrow(Y_mat) == 75,
  "E deve ter 75 linhas"           = nrow(E)     == 75,
  "x deve ter dimensão correta"    = all(dim(x) == c(75, 23, 3)),
  "Y_mat e E devem ter mesma dim"  = identical(dim(Y_mat), dim(E)),
  "linhas de Y_mat e E alinhadas"  = identical(row.names(Y_mat), row.names(E)),
  "linhas de Y_mat e x alinhadas"  = identical(row.names(Y_mat), dimnames(x)[[1]]),
  "estrutura espacial consistente" = sum(data$num) == data$sumNumNeigh,
  "K deve ser 8"                   = ncol(data$hAI) == 8
)
cat("Verificações de consistência: OK\n")

E_norm    <- E / mean(E)
N_regions <- nrow(Y_mat)          # 75
n_times   <- ncol(Y_mat)          # 23
p         <- dim(x)[3]            # 3
K         <- ncol(data$hAI)       # 8  (lido automaticamente)
cat(sprintf("N = %d | T = %d | p = %d | K = %d\n", N_regions, n_times, p, K))

adj_vec     <- as.integer(data$adj)
num_vec     <- as.integer(data$num)
n_adj_val   <- as.integer(data$sumNumNeigh)
weights_vec <- rep(1.0, n_adj_val)
h_mat       <- data$hAI   # 75 × 8

# ==============================================================================
# 2. PRIORI TEMPO-VARIANTE PARA gamma[1,t]  — IDÊNTICA AO PADRÃO OURO
# ==============================================================================
#
#   t = 1..10  → Unif(0.05, 0.15)  : âncora restrita
#   t = 11..17 → Unif(0.00, 0.10)  : âncora moderada
#   t = 18..23 → Unif(0.00, 0.05)  : âncora difusa

a_unif_vec <- c(rep(0.05, 10), rep(0.00, 13))
b_unif_vec <- c(rep(0.15, 10), rep(0.10,  7), rep(0.05, 6))

stopifnot(
  "a_unif_vec deve ter comprimento T" = length(a_unif_vec) == n_times,
  "b_unif_vec deve ter comprimento T" = length(b_unif_vec) == n_times,
  "b_unif > a_unif em todo t"         = all(b_unif_vec > a_unif_vec)
)
cat("Priori de gamma[1,t] verificada: OK\n")
cat(sprintf("  a_unif: %.2f (t=1..10) | %.2f (t=11..17) | %.2f (t=18..23)\n",
            a_unif_vec[1], a_unif_vec[11], a_unif_vec[18]))
cat(sprintf("  b_unif: %.2f (t=1..10) | %.2f (t=11..17) | %.2f (t=18..23)\n",
            b_unif_vec[1], b_unif_vec[11], b_unif_vec[18]))

# ==============================================================================
# 3. REGIÕES DE INTERESSE  (3 por cluster → até 24 regiões)
# ==============================================================================

cluster_ids <- apply(h_mat, 1, sum)   # cluster de cada região (1..K)

set.seed(42)
REGIONS_INTEREST <- unlist(lapply(1:K, function(cl) {
  regs <- which(cluster_ids == cl)
  if (length(regs) >= 3) sample(regs, 3) else regs
}))
cat("Regiões selecionadas:\n"); print(sort(REGIONS_INTEREST))

MU_MONITORS <- unlist(lapply(
  REGIONS_INTEREST,
  function(r) paste0("mu[", r, ", ", seq_len(n_times), "]")
))

# ==============================================================================
# 4. CONSTANTES E DADOS NIMBLE
# ==============================================================================

constants_spatial <- list(
  n_regions = N_regions, n_times = n_times, p = p, K = K,
  h         = h_mat,
  mu_beta   = rep(0, p),
  # Vetores de comprimento T (priori tempo-variante, igual ao padrão ouro)
  a_unif    = a_unif_vec,
  b_unif    = b_unif_vec,
  adj     = adj_vec,    num     = num_vec,
  weights = weights_vec, n_adj  = n_adj_val
)
constants_nonspatial <- constants_spatial[
  setdiff(names(constants_spatial), c("adj", "num", "weights", "n_adj"))
]
data_nimble <- list(Y = Y_mat, E = E_norm, x = x)

# ==============================================================================
# 5. INICIALIZAÇÕES — K = 8, gamma K × T
# ==============================================================================
#
# Para cada t:
#   gamma[1,t] = ponto médio de [a_unif[t], b_unif[t]]
#   Cadeia 2:   0.9 × ponto médio  (diversidade)
#   gamma[2..K,t]: pesos decrescentes em 50 % do simplex residual
#
# Verificação formal para t ∈ {1, 11, 18, 23}:
#   a_unif[t] ≤ gamma[1,t] ≤ b_unif[t]  e  sum(gamma[,t]) < 1

set.seed(123)

gamma_init1 <- matrix(0, K, n_times)
gamma_init2 <- matrix(0, K, n_times)

for (t in seq_len(n_times)) {
  g1_mid <- (a_unif_vec[t] + b_unif_vec[t]) / 2
  
  gamma_init1[1, t] <- g1_mid
  gamma_init2[1, t] <- 0.9 * g1_mid
  
  pesos <- seq(0.08, 0.02, length.out = K - 1)
  pesos <- pesos / sum(pesos)
  
  R1 <- 1 - gamma_init1[1, t]
  R2 <- 1 - gamma_init2[1, t]
  
  gamma_init1[2:K, t] <- R1 * 0.5 * pesos
  gamma_init2[2:K, t] <- R2 * 0.5 * pesos
}

# Verificação formal das restrições para t representativos
for (t_check in c(1, 11, 18, n_times)) {
  s1 <- sum(gamma_init1[, t_check])
  s2 <- sum(gamma_init2[, t_check])
  stopifnot(
    "gamma_init1 soma < 1"          = s1 < 1,
    "gamma_init2 soma < 1"          = s2 < 1,
    "gamma_init1[1,t] >= a_unif[t]" = gamma_init1[1, t_check] >= a_unif_vec[t_check],
    "gamma_init1[1,t] <= b_unif[t]" = gamma_init1[1, t_check] <= b_unif_vec[t_check],
    "gamma_init2[1,t] >= a_unif[t]" = gamma_init2[1, t_check] >= a_unif_vec[t_check],
    "gamma_init2[1,t] <= b_unif[t]" = gamma_init2[1, t_check] <= b_unif_vec[t_check]
  )
}
cat("Inicializações de gamma[K×T] verificadas: OK\n")

inits_list_spatial <- list(
  list(beta0t  = rep(0, n_times),
       beta    = rep(0, p),
       gamma   = gamma_init1,
       sigma_s = 0.5,
       s       = rep(0, N_regions)),
  list(beta0t  = rnorm(n_times, 0, 0.3),
       beta    = rnorm(p, 0, 0.3),
       gamma   = gamma_init2,
       sigma_s = 1.0,
       s       = rep(0, N_regions))
)
inits_list_nonspatial <- list(
  list(beta0t = rep(0, n_times),
       beta   = rep(0, p),
       gamma  = gamma_init1),
  list(beta0t = rnorm(n_times, 0, 0.3),
       beta   = rnorm(p, 0, 0.3),
       gamma  = gamma_init2)
)

# ==============================================================================
# 6. FUNÇÃO WORKER
# ==============================================================================

run_model <- function(model_type, output_dir) {
  
  library(nimble); library(coda); library(dplyr)
  library(ggplot2); library(readr); library(stringr); library(tibble)
  
  # ── 6a. Código NIMBLE ────────────────────────────────────────────────────────
  #
  # gamma[1,t] ~ dunif(a_unif[t], b_unif[t])  — vetores indexados por t
  # (idêntico ao padrão ouro; apenas lambda[t] é removido e substituído por beta0t[t])
  
  code_spatial <- nimbleCode({
    # Intercepto temporal i.i.d.
    for (t in 1:n_times) beta0t[t] ~ dnorm(0, sd = 1)
    # Covariáveis
    for (j in 1:p) beta[j] ~ dnorm(mu_beta[j], sd = 1)
    # Pesos de cluster — priori tempo-variante (âncora k=1 com três regimes)
    for (t in 1:n_times) {
      gamma[1, t] ~ dunif(min = a_unif[t], max = b_unif[t])
      for (k in 2:K)
        gamma[k, t] ~ dunif(min = 0, max = (1 - sum(gamma[1:(k-1), t])))
    }
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
        log(mu[i, t]) <- beta0t[t] + log(E[i, t]) + log(epsilon[i, t]) +
          inprod(beta[1:p], x[i, t, 1:p]) + s[i]
        Y[i, t]        ~ dpois(mu[i, t])
        logLik_Y[i, t] <- dpois(Y[i, t], mu[i, t], log = TRUE)
      }
    }
  })
  
  code_nonspatial <- nimbleCode({
    for (t in 1:n_times) beta0t[t] ~ dnorm(0, sd = 1)
    for (j in 1:p) beta[j] ~ dnorm(mu_beta[j], sd = 1)
    for (t in 1:n_times) {
      gamma[1, t] ~ dunif(min = a_unif[t], max = b_unif[t])
      for (k in 2:K)
        gamma[k, t] ~ dunif(min = 0, max = (1 - sum(gamma[1:(k-1), t])))
    }
    for (i in 1:n_regions) {
      for (t in 1:n_times) {
        epsilon[i, t] <- 1 - inprod(h[i, 1:K], gamma[1:K, t])
        log(mu[i, t]) <- beta0t[t] + log(E[i, t]) + log(epsilon[i, t]) +
          inprod(beta[1:p], x[i, t, 1:p])
        Y[i, t]        ~ dpois(mu[i, t])
        logLik_Y[i, t] <- dpois(Y[i, t], mu[i, t], log = TRUE)
      }
    }
  })
  
  # ── 6b. Seleção de objetos ───────────────────────────────────────────────────
  
  is_spatial <- (model_type == "spatial")
  model_code <- if (is_spatial) code_spatial      else code_nonspatial
  constants  <- if (is_spatial) constants_spatial  else constants_nonspatial
  inits_list <- if (is_spatial) inits_list_spatial else inits_list_nonspatial
  
  cat("\n=== Iniciando modelo:", model_type, "===\n")
  scenario_dir <- file.path(output_dir, model_type)
  dir.create(scenario_dir, recursive = TRUE, showWarnings = FALSE)
  
  # ── 6c. Construir e compilar ─────────────────────────────────────────────────
  
  model  <- nimbleModel(code = model_code, constants = constants,
                        data = data_nimble, inits = inits_list[[1]], check = FALSE)
  Cmodel <- compileNimble(model)
  
  conf <- configureMCMC(model)
  # gamma: um AF_slice por coluna-t (23 blocos de K=8 parâmetros — padrão ouro)
  conf$removeSamplers("gamma")
  for (t_idx in seq_len(n_times)) {
    conf$addSampler(
      target = paste0("gamma[", seq_len(K), ", ", t_idx, "]"),
      type   = "AF_slice"
    )
  }
  # beta0t[t]: amostrador RW escalar padrão do NIMBLE (um por tempo)
  
  monitors_base <- c("beta0t", "beta", "gamma", "logLik_Y")
  if (is_spatial) monitors_base <- c(monitors_base, "s", "sigma_s", "tau_s")
  conf$addMonitors(monitors_base)
  conf$addMonitors(MU_MONITORS)
  conf$printSamplers()
  
  Rmcmc <- buildMCMC(conf)
  Cmcmc <- compileNimble(Rmcmc, project = model)
  
  # ── 6d. MCMC ─────────────────────────────────────────────────────────────────
  
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
  
  # ── 6e. Funções auxiliares ───────────────────────────────────────────────────
  
  safe_hpd <- function(sv) {
    if (var(sv) < 1e-12) return(c(NA_real_, NA_real_))
    as.numeric(HPDinterval(as.mcmc(sv), prob = 0.95))
  }
  safe_gelman <- function(obj) {
    tryCatch(gelman.diag(obj, autoburnin = FALSE)$psrf[, 1],
             error = function(e) rep(NA_real_, nvar(obj)))
  }
  
  beta0t_names <- paste0("beta0t[", seq_len(n_times), "]")
  beta_names   <- paste0("beta[",   seq_len(p),       "]")
  
  # Nomes "gamma[k, t]" — matriz K × T
  gamma_names_mat <- outer(
    seq_len(K), seq_len(n_times),
    function(k, t) paste0("gamma[", k, ", ", t, "]")
  )
  
  # ── 6f. Epsilon posterior — array (n_draw × N × T) ──────────────────────────
  
  n_draw        <- nrow(samples_mat)
  epsilon_draws <- array(NA_real_, dim = c(n_draw, N_regions, n_times))
  for (t in seq_len(n_times)) {
    g_t <- samples_mat[, gamma_names_mat[, t], drop = FALSE]  # n_draw × K
    epsilon_draws[, , t] <- 1 - g_t %*% t(h_mat)              # n_draw × N
  }
  
  # ── 6g. Sumário de gamma[k, t] ───────────────────────────────────────────────
  
  gamma_summary <- do.call(rbind, lapply(seq_len(K), function(k) {
    do.call(rbind, lapply(seq_len(n_times), function(t) {
      nm  <- gamma_names_mat[k, t]
      sv  <- samples_mat[, nm]; hpd <- safe_hpd(sv)
      tibble(k = k, Time = t, Parameter = nm,
             Mean = mean(sv), SD = sd(sv), Lower = hpd[1], Upper = hpd[2],
             # Limites da priori registrados para diagnóstico (só k=1)
             a_unif_t = if (k == 1) a_unif_vec[t] else 0,
             b_unif_t = if (k == 1) b_unif_vec[t] else NA_real_,
             ESS  = as.numeric(effectiveSize(mcmc_list_full[, nm])),
             Rhat = safe_gelman(mcmc_list_full[, nm]))
    }))
  }))
  write_csv(gamma_summary, file.path(scenario_dir, "gamma_summary.csv"))
  
  # Painel temporal de gamma — linhas verticais marcam mudanças de regime
  ggsave(
    file.path(scenario_dir, "painel_gamma.png"),
    ggplot(gamma_summary, aes(x = Time)) +
      geom_ribbon(aes(ymin = Lower, ymax = Upper), fill = "steelblue", alpha = 0.3) +
      geom_line(aes(y = Mean), color = "steelblue", linewidth = 0.8) +
      geom_vline(xintercept = c(10.5, 17.5), linetype = "dashed",
                 color = "firebrick", linewidth = 0.5) +
      facet_wrap(
        ~ k, ncol = 2, scales = "free_y",
        labeller = labeller(k = function(x) paste0("gamma[", x, ", t]"))
      ) +
      theme_bw(base_size = 11) +
      labs(title    = paste("Trajetória de gamma[k, t] (", model_type, ") — K=8"),
           subtitle = "Linhas vermelhas: mudança de regime da priori de gamma[1,t]",
           x = "Tempo", y = expression(gamma[kt])),
    width = 12, height = 10
  )
  
  # ── 6h. Epsilon por cluster (agrega regiões do mesmo cluster) ────────────────
  
  eps_full <- do.call(rbind, lapply(seq_len(N_regions), function(i) {
    do.call(rbind, lapply(seq_len(n_times), function(t) {
      sv  <- epsilon_draws[, i, t]; hpd <- safe_hpd(sv)
      tibble(Region = i, Cluster = cluster_ids[i], Time = t,
             Mean = mean(sv), Lower = hpd[1], Upper = hpd[2])
    }))
  }))
  write_csv(eps_full, file.path(scenario_dir, "epsilon_summary.csv"))
  
  # Média por cluster × tempo (pool de draws × regiões do cluster)
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
      geom_vline(xintercept = c(10.5, 17.5), linetype = "dashed",
                 color = "firebrick", linewidth = 0.5) +
      facet_wrap(
        ~ k, ncol = 2, scales = "free_y",
        labeller = labeller(k = function(x) paste0("epsilon[", x, ", t]"))
      ) +
      theme_bw(base_size = 11) +
      labs(title    = paste("Trajetória de epsilon por cluster (", model_type, ") — K=8"),
           subtitle = "Linhas vermelhas: mudança de regime da priori de gamma[1,t]",
           x = "Tempo", y = expression(epsilon[kt])),
    width = 12, height = 10
  )
  
  # ── 6i. Sumário e painel de beta0t ───────────────────────────────────────────
  
  beta0t_summary <- do.call(rbind, lapply(seq_len(n_times), function(t) {
    nm <- beta0t_names[t]; sv <- samples_mat[, nm]; hpd <- safe_hpd(sv)
    tibble(Time = t, Parameter = nm,
           Mean = mean(sv), SD = sd(sv), Lower = hpd[1], Upper = hpd[2],
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
      scale_x_continuous(breaks = seq_len(n_times), labels = anos_label) +
      theme_bw(base_size = 12) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      labs(title    = paste("Intercepto temporal beta0t (", model_type, ") — K=8"),
           subtitle = "Média posterior e HPD 95%",
           x = "Ano", y = expression(beta[0][t])),
    width = 10, height = 5
  )
  
  # ── 6j. Sumário de beta ───────────────────────────────────────────────────────
  
  beta_summary <- do.call(rbind, lapply(beta_names, function(nm) {
    sv  <- samples_mat[, nm]; hpd <- safe_hpd(sv)
    tibble(Parameter = nm, Mean = mean(sv), SD = sd(sv),
           HPD_Lower = hpd[1], HPD_Upper = hpd[2],
           ESS  = as.numeric(effectiveSize(mcmc_list_full[, nm])),
           Rhat = safe_gelman(mcmc_list_full[, nm]))
  }))
  write_csv(beta_summary, file.path(scenario_dir, "beta_summary.csv"))
  
  # ── 6k. tau_s e s[i] (apenas espacial) ───────────────────────────────────────
  
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
        labs(title = paste("Efeito espacial s[i] (", model_type, ") — K=8"),
             y = "s[i]", x = "Região"),
      width = 10, height = 5
    )
  }
  
  # ── 6l. mu para regiões selecionadas ─────────────────────────────────────────
  
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
      labs(title = paste("Mu estimado (", model_type, ") — K=8"),
           x = "Tempo", y = expression(mu[it])),
    width = 14, height = 16
  )
  
  # ── 6m. WAIC e LPML ───────────────────────────────────────────────────────────
  
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
  
  # ── 6n. Diagnósticos ACF ──────────────────────────────────────────────────────
  
  # Parâmetros estruturais (não temporais): beta + tau_s
  params_struct <- c(beta_names, if (is_spatial) "tau_s" else character(0))
  # Diagnóstico completo: beta0t + beta + gamma + tau_s
  all_diag <- c(beta0t_names, params_struct, as.vector(gamma_names_mat))
  
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
  
  # ACF de beta (e tau_s) — padrão ouro
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
      labs(title    = paste("ACF de beta (", model_type, ") — K=8"),
           subtitle = "Azul: |0.10| | Vermelho: |0.05|"),
    width = 10, height = 5
  )
  
  # ACF de beta0t — painel por tempo (padrão ouro)
  acf_beta0t_df <- do.call(rbind, lapply(beta0t_names, function(nm) {
    ac <- acf(samples_mat[, nm], lag.max = 100, plot = FALSE)
    tibble(Time = as.integer(str_extract(nm, "\\d+")),
           Lag  = as.vector(ac$lag[-1]), ACF = as.vector(ac$acf[-1]))
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
      labs(title    = paste("ACF de beta0t (", model_type, ") — K=8"),
           subtitle = "Azul: |0.10| | Vermelho: |0.05|"),
    width = 14, height = 10
  )
  
  # ACF de gamma[k,t] — um arquivo por cluster k (padrão ouro)
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
      labs(title    = paste0("ACF de gamma[", k_local, ", t] (", model_type, ") — K=8"),
           subtitle = "Azul: |0.10| | Vermelho: |0.05|",
           x = "Lag", y = "ACF")
    ggsave(file.path(scenario_dir, sprintf("acf_gamma_k%d.png", k_local)),
           p_acf, width = 14, height = 10)
  }))
  
  # Traceplots de beta (e tau_s) — padrão ouro
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
      labs(title    = paste("Traceplots beta (", model_type, ") — K=8"),
           subtitle = "Linha grossa = média ergódica | Linha fina = cadeia",
           x = "Iteração (pós-burnin)", y = "Valor"),
    width = 10, height = max(5, 3 * ceiling(length(params_acf_beta) / 3))
  )
  
  # ── 6o. Retorno resumido ──────────────────────────────────────────────────────
  
  tibble(
    model = model_type, K = K, niter = niter, nburnin = nburnin, thin = thin,
    WAIC = waic, LPML = LPML,
    ESS_beta_min    = min(beta_summary$ESS,    na.rm = TRUE),
    ESS_gamma_min   = min(gamma_summary$ESS,   na.rm = TRUE),
    ESS_beta0t_min  = min(beta0t_summary$ESS,  na.rm = TRUE),
    ESS_beta0t_max  = max(beta0t_summary$ESS,  na.rm = TRUE),
    ESS_tau         = ESS_tau,
    Rhat_beta_max   = max(Rhat_struct[beta_names], na.rm = TRUE),
    Rhat_gamma_max  = max(gamma_summary$Rhat,  na.rm = TRUE),
    Rhat_beta0t_max = max(beta0t_summary$Rhat, na.rm = TRUE),
    lag_max_0.10    = max(acf_results$lag_0.10, na.rm = TRUE),
    lag_max_0.05    = max(acf_results$lag_0.05, na.rm = TRUE)
  )
}

# ==============================================================================
# 7. EXECUÇÃO PARALELA
# ==============================================================================

model_types <- c("spatial", "non_spatial")
n_cores     <- min(length(model_types), parallel::detectCores() - 1)
if (n_cores < 1) n_cores <- 1

output_dir <- file.path(
  "C:/Users/vlara/OneDrive/Estatistica UFMG/Mestrado/Pesquisa/Aplicação",
  "resultados_intercepto_temporal_gamma_t_K8"
)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

cl <- makeCluster(n_cores)
clusterExport(cl, c(
  "constants_spatial", "constants_nonspatial", "data_nimble",
  "inits_list_spatial", "inits_list_nonspatial",
  "N_regions", "n_times", "p", "K",
  "h_mat", "cluster_ids",
  "a_unif_vec", "b_unif_vec",        # vetores da priori tempo-variante
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
# 8. CONSOLIDAÇÃO E GRÁFICOS COMPARATIVOS
# ==============================================================================

resumo <- bind_rows(resultados)
write_csv(resumo, file.path(output_dir, "resumo_comparativo.csv"))
cat("\n=== RESUMO COMPARATIVO (M6 — K=8) ===\n"); print(resumo)

# Comparativo beta0t
beta0t_all <- lapply(model_types, function(m) {
  path <- file.path(output_dir, m, "beta0t_summary.csv")
  if (!file.exists(path)) return(NULL)
  read_csv(path, show_col_types = FALSE) |> mutate(model = m)
}) |> bind_rows()

if (nrow(beta0t_all) > 0) {
  anos_label <- colnames(data_nimble$Y)
  if (is.null(anos_label)) anos_label <- seq_len(n_times)
  ggsave(
    file.path(output_dir, "beta0t_comparativo.png"),
    ggplot(beta0t_all, aes(x = Time, y = Mean, color = model, fill = model)) +
      geom_ribbon(aes(ymin = Lower, ymax = Upper), alpha = 0.15, color = NA) +
      geom_line(linewidth = 0.9) +
      scale_x_continuous(breaks = seq_len(n_times), labels = anos_label) +
      theme_bw(base_size = 12) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            legend.position = "bottom") +
      labs(title = "Comparação beta0t: espacial vs. não-espacial — K=8",
           x = "Ano", y = expression(beta[0][t]),
           color = "Modelo", fill = "Modelo"),
    width = 10, height = 5, dpi = 300
  )
}

# Comparativo gamma[k,t] com bandas da priori (padrão ouro)
gamma_all <- lapply(model_types, function(m) {
  path <- file.path(output_dir, m, "gamma_summary.csv")
  if (!file.exists(path)) return(NULL)
  read_csv(path, show_col_types = FALSE) |> mutate(model = m)
}) |> bind_rows()

if (nrow(gamma_all) > 0) {
  prior_ribbon <- tibble(
    Time    = seq_len(n_times),
    k       = 1L,
    a_lower = a_unif_vec,
    b_upper = b_unif_vec
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
      geom_vline(xintercept = c(10.5, 17.5), linetype = "dashed",
                 color = "firebrick", linewidth = 0.4) +
      facet_wrap(
        ~ k, ncol = 2, scales = "free_y",
        labeller = labeller(k = function(x) paste0("gamma[", x, ", t]"))
      ) +
      theme_bw(base_size = 12) + theme(legend.position = "bottom") +
      labs(title    = "Comparação gamma[k,t]: espacial vs. não-espacial — K=8",
           subtitle = "Banda cinza em k=1: suporte da priori | Linhas vermelhas: mudança de regime",
           x = "Tempo", y = expression(gamma[kt]),
           color = "Modelo", fill = "Modelo"),
    width = 12, height = 14, dpi = 300
  )
}

# Comparativo epsilon por cluster
eps_cluster_all <- lapply(model_types, function(m) {
  path <- file.path(output_dir, m, "epsilon_cluster_summary.csv")
  if (!file.exists(path)) return(NULL)
  read_csv(path, show_col_types = FALSE) |> mutate(model = m)
}) |> bind_rows()

if (nrow(eps_cluster_all) > 0) {
  ggsave(
    file.path(output_dir, "epsilon_comparativo.png"),
    ggplot(eps_cluster_all,
           aes(x = Time, y = Mean, color = model, fill = model)) +
      geom_ribbon(aes(ymin = Lower, ymax = Upper), alpha = 0.12, color = NA) +
      geom_line(linewidth = 0.8) +
      geom_vline(xintercept = c(10.5, 17.5), linetype = "dashed",
                 color = "firebrick", linewidth = 0.4) +
      facet_wrap(
        ~ k, ncol = 2, scales = "free_y",
        labeller = labeller(k = function(x) paste0("epsilon[", x, ", t]"))
      ) +
      theme_bw(base_size = 12) + theme(legend.position = "bottom") +
      labs(title    = "Comparação epsilon por cluster: espacial vs. não-espacial — K=8",
           subtitle = "Banda = IC 95% HPD | Linhas vermelhas: mudança de regime da priori",
           x = "Tempo", y = expression(epsilon[kt]),
           color = "Modelo", fill = "Modelo"),
    width = 12, height = 14, dpi = 300
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
      labs(title = "Comparação mu: espacial vs. não-espacial — K=8",
           x = "Tempo", color = "Modelo", fill = "Modelo"),
    width = 14, height = 16, dpi = 300
  )
}

cat("\nTempo total de execução:\n")
print(Sys.time() - inicio_global)
