# ==============================================================================
# ajuste_intercepto_temporal_gamma_proporcional.R
#
# Modelo Poisson log-linear com epsilon proporcional ao melhor cluster
#
# log(mu[i,t]) = beta0t[t] + log(E[i,t]) + log(epsilon[i,t])
#                + beta[1:p] %*% x[i,t,1:p]   (+ s[i]  se espacial)
#
# epsilon[i,t] = 1 - gamma_1[t] * inprod(h[i, 1:K], kappa[1:K])
#
# ── Processo temporal gamma_1[t] (melhor cluster) ──────────────────────────
#   gamma_1[t] ~ U(L[t] - delta_g, L[t] + delta_g)
#   L[t]: sequência decrescente fixada (ex.: seq(0.12, 0.02, T))
#   delta_g: meia-largura fixa; requer L[T] > delta_g
#   Interpretação: probabilidade de notificação melhora ao longo do tempo
#
# ── Multiplicadores estáticos kappa[k] ────────────────────────────────────
#   kappa[1] = 1  (melhor cluster, fixo)
#   kappa[k] = kappa[k-1] + incr[k-1],  k = 2,...,K
#   incr[k]  ~ dunif(0, M_incr)
#   => garante 1 = kappa[1] < kappa[2] < kappa[3] < kappa[4] por construção
#
# ── Restrição de viabilidade (epsilon > 0) ────────────────────────────────
#   constraint_kappa ~ dconstraint(kappa[K] * (L_max + delta_g) < 1)
#   Qualquer proposta que torne epsilon <= 0 é imediatamente rejeitada pelo
#   dconstraint, sem risco de log(epsilon) = -Inf travar as cadeias.
#
# ── Priors restantes ──────────────────────────────────────────────────────
#   beta0t[t] ~ N(0,1) i.i.d.
#   beta[j]   ~ N(0,1)
#   sigma_s   ~ Half-t(0,1,1)  [apenas espacial]
#   s         ~ ICAR(tau_s)    [apenas espacial]
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
# 2. HIPERPARÂMETROS DO PROCESSO TEMPORAL gamma_1[t]
# ==============================================================================
#
# L[t]: centros do prior uniforme; sequência linear decrescente.
# delta_g: meia-largura da uniforme.  Deve satisfazer L[T] > delta_g.
#
# Ajuste esses valores conforme o seu conhecimento de domínio:
#   gamma_1t,max ≈ 0.12  (cluster melhor, tempo inicial)
#   gamma_1t,min ≈ 0.02  (cluster melhor, tempo final)
#   delta_g      ≈ 0.015 (intervalo [L[t]-0.015, L[t]+0.015] em cada t)
# ------------------------------------------------------------------------------
gamma_1_max <- 0.12    # L[1]: centro do prior em t = 1
gamma_1_min <- 0.02    # L[T]: centro do prior em t = T
delta_g     <- 0.015   # meia-largura fixa; altere aqui se necessário

# Verificação: L[T] - delta_g > 0 é necessário para que a uniforme seja válida
stopifnot("delta_g deve ser menor que gamma_1_min (L[T] > delta_g)" =
            delta_g < gamma_1_min)

L_seq <- seq(gamma_1_max, gamma_1_min, length.out = n_times)
L_max <- L_seq[1]   # = gamma_1_max  (pior caso para a restrição de kappa)

# Limite individual dos incrementos (cada incr[k] individualmente).
# A restrição conjunta é tratada pelo dconstraint, mas M_incr evita
# que o amostrador explore regiões absurdas de forma desnecessária.
# kappa[K] < 1/(L_max + delta_g) => máximo total dos incrementos < 1/(L_max+delta_g) - 1
kappa_max_hard <- 1 / (L_max + delta_g)          # e.g. ≈ 7.41 para L_max=0.12, delta_g=0.015
M_incr         <- ceiling(kappa_max_hard - 1)    # e.g. 7 (limite por incremento)

cat(sprintf(
  "L[1]=%.3f | L[T]=%.3f | delta_g=%.3f\nkappa[K] max: %.3f | M_incr=%d\n",
  L_max, tail(L_seq, 1), delta_g, kappa_max_hard, M_incr
))

# ==============================================================================
# 3. REGIÕES PARA MONITORAR
# ==============================================================================
# cluster_ids[i]: cluster dominante da região i (argmax de h_mat[i,])
cluster_ids <- apply(h_mat, 1, which.max)

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
# 4. CONSTANTES E DADOS NIMBLE
# ==============================================================================
constants_base <- list(
  n_regions = N_regions, n_times = n_times, p = p, K = K,
  h = h_mat,
  L      = L_seq,    # vetor de centros decrescentes para gamma_1[t]
  delta_g = delta_g, # meia-largura (escalar)
  L_max  = L_max,    # = L[1], usado na restrição de kappa
  M_incr = M_incr    # limite superior individual de cada incr[k]
)

constants_spatial <- c(constants_base, list(
  adj = adj_vec, num = num_vec, weights = weights_vec, n_adj = n_adj_val
))
constants_nonspatial <- constants_base

# constraint_kappa = 1L indica ao NIMBLE que a condição deve ser VERDADEIRA.
# Qualquer proposta que a viole recebe log-verossimilhança = -Inf e é rejeitada.
data_nimble <- list(
  Y = Y_mat, E = E_norm, x = x,
  constraint_kappa = 1L
)

# ==============================================================================
# 5. INICIALIZAÇÕES
# ==============================================================================
set.seed(123)

# gamma_1 inicializado nos centros L[t] (interior válido do prior)
gamma_1_init1 <- L_seq
gamma_1_init2 <- pmax(
  L_seq - delta_g * 0.8,
  pmin(L_seq + delta_g * 0.8, L_seq + rnorm(n_times, 0, delta_g / 3))
)

# incr inicializado para kappa = (1, 2, 3, 4)  [incr = (1, 1, 1)]
# Verificação: kappa[K] * (L_max + delta_g) = 4 * 0.135 = 0.54 < 1  ✓
incr_init1 <- rep(1,   K - 1)
incr_init2 <- rep(0.5, K - 1)

inits_spatial_1 <- list(
  beta0t = rep(0, n_times), beta = rep(0, p),
  gamma_1 = gamma_1_init1, incr = incr_init1,
  sigma_s = 0.5, s = rep(0, N_regions)
)
inits_spatial_2 <- list(
  beta0t = rnorm(n_times, 0, 0.3), beta = rnorm(p, 0, 0.3),
  gamma_1 = gamma_1_init2, incr = incr_init2,
  sigma_s = 1.0, s = rep(0, N_regions)
)
inits_nonspatial_1 <- list(
  beta0t = rep(0, n_times), beta = rep(0, p),
  gamma_1 = gamma_1_init1, incr = incr_init1
)
inits_nonspatial_2 <- list(
  beta0t = rnorm(n_times, 0, 0.3), beta = rnorm(p, 0, 0.3),
  gamma_1 = gamma_1_init2, incr = incr_init2
)

inits_list_spatial    <- list(inits_spatial_1,    inits_spatial_2)
inits_list_nonspatial <- list(inits_nonspatial_1, inits_nonspatial_2)

# ==============================================================================
# 6. FUNÇÃO WORKER
# ==============================================================================
run_model <- function(model_type, output_dir) {

  library(nimble); library(coda); library(dplyr)
  library(ggplot2); library(readr); library(stringr); library(tibble)

  # ── Código NIMBLE ────────────────────────────────────────────────────────
  #
  #  Estrutura do modelo:
  #
  #  1. gamma_1[t] ~ U(L[t] - delta_g, L[t] + delta_g)
  #     Processo do melhor cluster; prior decrescente garantido pelos centros L[t].
  #
  #  2. incr[k] ~ dunif(0, M_incr)   (k = 1,...,K-1)
  #     kappa[1] = 1  (fixo, melhor cluster)
  #     kappa[k] = kappa[k-1] + incr[k-1]
  #     => 1 = kappa[1] < kappa[2] < ... < kappa[K]  por construção.
  #
  #  3. constraint_kappa ~ dconstraint(kappa[K] * (L_max + delta_g) < 1)
  #     Rejeição imediata de qualquer proposta que tornaria epsilon <= 0.
  #
  #  4. epsilon[i,t] = 1 - gamma_1[t] * inprod(h[i,], kappa)
  # -------------------------------------------------------------------------

  code_spatial <- nimbleCode({
    # Interceptos temporais i.i.d.
    for (t in 1:n_times) { beta0t[t] ~ dnorm(0, sd = 1) }

    # Covariáveis
    for (j in 1:p) { beta[j] ~ dnorm(0, sd = 1) }

    # Processo temporal: gamma_1[t] ~ U(L[t]-delta_g, L[t]+delta_g)
    # Os centros L[t] são decrescentes, codificando a melhora de notificação
    for (t in 1:n_times) {
      gamma_1[t] ~ dunif(L[t] - delta_g, L[t] + delta_g)
    }

    # Multiplicadores: reparametrização por incrementos (garante ordem sem rejeição)
    for (k in 1:(K - 1)) { incr[k] ~ dunif(0, M_incr) }
    kappa[1] <- 1                                   # melhor cluster, fixo
    for (k in 2:K) { kappa[k] <- kappa[k - 1] + incr[k - 1] }

    # Restrição dura de viabilidade: epsilon[i,t] > 0 para todo i,t.
    # O pior caso é kappa[K] × max(gamma_1) = kappa[K] × (L_max + delta_g).
    # dconstraint retorna densidade 0 (log = -Inf) se a condição for falsa,
    # rejeitando automaticamente a proposta no passo MH.
    constraint_kappa ~ dconstraint(kappa[K] * (L_max + delta_g) < 1)

    # Efeito espacial ICAR
    sigma_s ~ T(dt(0, 1, 1), 0, )
    tau_s   <- 1 / (sigma_s ^ 2)
    s[1:n_regions] ~ dcar_normal(
      adj[1:n_adj], weights[1:n_adj], num[1:n_regions], tau_s, zero_mean = 1
    )

    # Verossimilhança
    for (i in 1:n_regions) {
      for (t in 1:n_times) {
        epsilon[i, t]    <- 1 - gamma_1[t] * inprod(h[i, 1:K], kappa[1:K])
        log(mu[i, t])    <- beta0t[t] + log(E[i, t]) + log(epsilon[i, t]) +
                            inprod(beta[1:p], x[i, t, 1:p]) + s[i]
        Y[i, t]          ~ dpois(mu[i, t])
        logLik_Y[i, t]   <- dpois(Y[i, t], mu[i, t], log = TRUE)
      }
    }
  })

  code_nonspatial <- nimbleCode({
    for (t in 1:n_times) { beta0t[t] ~ dnorm(0, sd = 1) }
    for (j in 1:p)        { beta[j]   ~ dnorm(0, sd = 1) }

    for (t in 1:n_times) {
      gamma_1[t] ~ dunif(L[t] - delta_g, L[t] + delta_g)
    }

    for (k in 1:(K - 1)) { incr[k] ~ dunif(0, M_incr) }
    kappa[1] <- 1
    for (k in 2:K) { kappa[k] <- kappa[k - 1] + incr[k - 1] }

    constraint_kappa ~ dconstraint(kappa[K] * (L_max + delta_g) < 1)

    for (i in 1:n_regions) {
      for (t in 1:n_times) {
        epsilon[i, t]    <- 1 - gamma_1[t] * inprod(h[i, 1:K], kappa[1:K])
        log(mu[i, t])    <- beta0t[t] + log(E[i, t]) + log(epsilon[i, t]) +
                            inprod(beta[1:p], x[i, t, 1:p])
        Y[i, t]          ~ dpois(mu[i, t])
        logLik_Y[i, t]   <- dpois(Y[i, t], mu[i, t], log = TRUE)
      }
    }
  })

  is_spatial <- (model_type == "spatial")
  model_code <- if (is_spatial) code_spatial      else code_nonspatial
  constants  <- if (is_spatial) constants_spatial  else constants_nonspatial
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

  # Amostrador AF_slice em BLOCO para incr[1:(K-1)]:
  #   - Propõe todos os incrementos simultaneamente, respeitando a restrição
  #     conjunta via dconstraint de forma muito mais eficiente que RW individuais.
  #   - Evita alta rejeição ("travar" das cadeias) que ocorreria com truncamentos
  #     dinâmicos individuais.
  conf$removeSamplers("incr")
  conf$addSampler(
    target = paste0("incr[", seq_len(K - 1), "]"),
    type   = "AF_slice"
  )

  # gamma_1[t]: RW individual com suporte restrito pelo prior uniforme.
  # O NIMBLE ajusta automaticamente a escala do RW para cada t.
  conf$removeSamplers("gamma_1")
  for (t_idx in seq_len(n_times)) {
    conf$addSampler(target = paste0("gamma_1[", t_idx, "]"), type = "RW")
  }

  monitors_base <- c("beta0t", "beta", "gamma_1", "kappa", "incr", "logLik_Y")
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
  gamma1_names <- paste0("gamma_1[", 1:n_times,  "]")
  kappa_names  <- paste0("kappa[",   1:K,         "]")
  incr_names   <- paste0("incr[",    1:(K - 1),   "]")

  # ── Epsilon posterior (N × T) ─────────────────────────────────────────────
  # epsilon[i,t] = 1 - gamma_1[t] * inprod(h[i,], kappa)
  #
  # Estratégia eficiente:
  #   h_kappa[i, d] = inprod(h[i,], kappa[d,])  →  N × n_draw
  #   epsilon[d, i, t] = 1 - h_kappa[i, d] * gamma_1[d, t]
  n_draw       <- nrow(samples_mat)
  kappa_draws  <- samples_mat[, kappa_names, drop = FALSE]  # n_draw × K
  h_kappa      <- h_mat %*% t(kappa_draws)                  # N × n_draw

  epsilon_draws <- array(NA_real_, dim = c(n_draw, N_regions, n_times))
  for (t in seq_len(n_times)) {
    g1_t <- samples_mat[, gamma1_names[t]]           # n_draw
    # sweep multiplica coluna d de h_kappa pelo draw g1_t[d]
    epsilon_draws[, , t] <- t(1 - sweep(h_kappa, 2, g1_t, `*`))   # n_draw × N
  }

  # ── Sumário de gamma_1[t] ─────────────────────────────────────────────────
  gamma1_summary <- do.call(rbind, lapply(1:n_times, function(t) {
    nm  <- gamma1_names[t]; sv <- samples_mat[, nm]; hpd <- safe_hpd(sv)
    tibble(Time = t, Parameter = nm,
           Mean = mean(sv), SD = sd(sv), Lower = hpd[1], Upper = hpd[2],
           L_center = L_seq[t],
           ESS  = as.numeric(effectiveSize(mcmc_list_full[, nm])),
           Rhat = safe_gelman(mcmc_list_full[, nm]))
  }))
  write_csv(gamma1_summary, file.path(scenario_dir, "gamma1_summary.csv"))

  ggsave(
    file.path(scenario_dir, "painel_gamma1.png"),
    ggplot(gamma1_summary, aes(x = Time)) +
      geom_ribbon(aes(ymin = Lower, ymax = Upper), fill = "steelblue", alpha = 0.3) +
      geom_line(aes(y = Mean),     color = "steelblue", linewidth = 0.9) +
      geom_line(aes(y = L_center), color = "black",     linewidth = 0.7, linetype = "dashed") +
      theme_bw(base_size = 12) +
      labs(title    = paste("Processo temporal gamma_1[t] (", model_type, ")"),
           subtitle = "Linha tracejada = centro do prior L[t]",
           x = "Tempo", y = expression(gamma[1][t])),
    width = 8, height = 5
  )

  # ── Sumário de kappa[k] (multiplicadores) ────────────────────────────────
  kappa_summary <- do.call(rbind, lapply(1:K, function(k) {
    nm  <- kappa_names[k]; sv <- samples_mat[, nm]; hpd <- safe_hpd(sv)
    tibble(k = k, Parameter = nm,
           Mean = mean(sv), SD = sd(sv), Lower = hpd[1], Upper = hpd[2],
           ESS  = as.numeric(effectiveSize(mcmc_list_full[, nm])),
           Rhat = safe_gelman(mcmc_list_full[, nm]))
  }))
  write_csv(kappa_summary, file.path(scenario_dir, "kappa_summary.csv"))

  ggsave(
    file.path(scenario_dir, "kappa_posterior.png"),
    ggplot(kappa_summary, aes(x = factor(k), y = Mean)) +
      geom_point(size = 3, color = "steelblue") +
      geom_errorbar(aes(ymin = Lower, ymax = Upper),
                    width = 0.3, linewidth = 0.7, color = "steelblue") +
      geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
      theme_bw(base_size = 12) +
      labs(title    = paste("Multiplicadores kappa[k] (", model_type, ")"),
           subtitle = "kappa[1] = 1 fixo | kappa[k] > kappa[k-1] por construção",
           x = "Cluster k", y = expression(kappa[k])),
    width = 6, height = 5
  )

  # ── Sumário de epsilon para regiões selecionadas ─────────────────────────
  eps_full <- do.call(rbind, lapply(1:N_regions, function(i) {
    do.call(rbind, lapply(1:n_times, function(t) {
      sv  <- epsilon_draws[, i, t]; hpd <- safe_hpd(sv)
      tibble(Region = i, Cluster = cluster_ids[i], Time = t,
             Mean = mean(sv), Lower = hpd[1], Upper = hpd[2])
    }))
  }))
  write_csv(eps_full, file.path(scenario_dir, "epsilon_summary.csv"))

  eps_ri     <- sort(REGIONS_INTEREST)
  eps_labels <- setNames(
    sprintf("Reg %d (C%d)", eps_ri, cluster_ids[eps_ri]),
    as.character(eps_ri)
  )
  ggsave(
    file.path(scenario_dir, "painel_epsilon.png"),
    ggplot(eps_full %>% filter(Region %in% eps_ri), aes(x = Time)) +
      geom_ribbon(aes(ymin = Lower, ymax = Upper), fill = "darkorange", alpha = 0.3) +
      geom_line(aes(y = Mean), color = "darkorange", linewidth = 0.8) +
      facet_wrap(~ Region, ncol = 3, scales = "free_y",
                 labeller = labeller(Region = eps_labels)) +
      theme_bw(base_size = 10) +
      labs(title = paste("Epsilon[i,t] estimado (", model_type, ")"),
           x = "Tempo", y = expression(epsilon[it])),
    width = 12, height = 10
  )

  # ── Sumário de beta0t ────────────────────────────────────────────────────
  anos_label <- colnames(data_nimble$Y)
  if (is.null(anos_label)) anos_label <- seq_len(n_times)

  beta0t_summary <- do.call(rbind, lapply(1:n_times, function(t) {
    nm <- beta0t_names[t]; sv <- samples_mat[, nm]; hpd <- safe_hpd(sv)
    tibble(Time = t, Parameter = nm, Mean = mean(sv), SD = sd(sv),
           Lower = hpd[1], Upper = hpd[2],
           ESS  = as.numeric(effectiveSize(mcmc_list_full[, nm])),
           Rhat = safe_gelman(mcmc_list_full[, nm]))
  }))
  write_csv(beta0t_summary, file.path(scenario_dir, "beta0t_summary.csv"))

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
    tau_sv  <- samples_mat[, "tau_s"]; hpd_t <- safe_hpd(tau_sv)
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
      facet_wrap(~ Region, ncol = 3, scales = "free_y",
                 labeller = labeller(Region = make_label(mu_regs))) +
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

  # ── Diagnósticos ACF ─────────────────────────────────────────────────────
  params_acf_beta <- c(beta_names, if (is_spatial) "tau_s" else character(0))

  acf_beta_df <- do.call(rbind, lapply(params_acf_beta, function(nm) {
    ac <- acf(samples_mat[, nm], lag.max = 100, plot = FALSE)
    tibble(Parameter = nm, Lag = as.vector(ac$lag[-1]), ACF = as.vector(ac$acf[-1]))
  }))
  ggsave(
    file.path(scenario_dir, "acf_beta.png"),
    ggplot(acf_beta_df, aes(x = Lag, y = ACF)) +
      geom_col(width = 0.6, fill = "grey50") +
      geom_hline(yintercept = c(-0.10,  0.10), linetype = "dashed",
                 color = "blue", linewidth = 0.5) +
      geom_hline(yintercept = c(-0.05,  0.05), linetype = "dotted",
                 color = "red",  linewidth = 0.5) +
      facet_wrap(~ Parameter, scales = "free_y") +
      theme_bw(base_size = 11) +
      labs(title    = paste("ACF de beta (", model_type, ")"),
           subtitle = "Azul: |0.10| | Vermelho: |0.05|"),
    width = 10, height = 5
  )

  # ACF de gamma_1[t]
  acf_g1_df <- do.call(rbind, lapply(1:n_times, function(t) {
    nm <- gamma1_names[t]
    ac <- acf(samples_mat[, nm], lag.max = 100, plot = FALSE)
    tibble(Time = t, Lag = as.vector(ac$lag[-1]), ACF = as.vector(ac$acf[-1]))
  }))
  ggsave(
    file.path(scenario_dir, "acf_gamma1.png"),
    ggplot(acf_g1_df, aes(x = Lag, y = ACF)) +
      geom_col(width = 0.6, fill = "steelblue") +
      geom_hline(yintercept = c(-0.10,  0.10), linetype = "dashed",
                 color = "blue", linewidth = 0.5) +
      geom_hline(yintercept = c(-0.05,  0.05), linetype = "dotted",
                 color = "red",  linewidth = 0.5) +
      facet_wrap(~ Time, scales = "free_y", ncol = 6,
                 labeller = label_bquote(gamma[1][.(Time)])) +
      theme_bw(base_size = 9) +
      labs(title    = paste("ACF de gamma_1[t] (", model_type, ")"),
           subtitle = "Azul: |0.10| | Vermelho: |0.05|"),
    width = 14, height = 10
  )

  # ACF dos incrementos incr[k] (e kappa derivados)
  acf_incr_df <- do.call(rbind, lapply(seq_len(K - 1), function(k) {
    nm <- incr_names[k]
    ac <- acf(samples_mat[, nm], lag.max = 100, plot = FALSE)
    tibble(Parameter = nm, Lag = as.vector(ac$lag[-1]), ACF = as.vector(ac$acf[-1]))
  }))
  ggsave(
    file.path(scenario_dir, "acf_incr_kappa.png"),
    ggplot(acf_incr_df, aes(x = Lag, y = ACF)) +
      geom_col(width = 0.6, fill = "purple") +
      geom_hline(yintercept = c(-0.10,  0.10), linetype = "dashed",
                 color = "blue", linewidth = 0.5) +
      geom_hline(yintercept = c(-0.05,  0.05), linetype = "dotted",
                 color = "red",  linewidth = 0.5) +
      facet_wrap(~ Parameter, scales = "free_y") +
      theme_bw(base_size = 11) +
      labs(title    = paste("ACF dos incrementos incr[k] (", model_type, ")"),
           subtitle = "Azul: |0.10| | Vermelho: |0.05|"),
    width = 8, height = 5
  )

  # ACF de beta0t
  acf_beta0t_df <- do.call(rbind, lapply(beta0t_names, function(nm) {
    ac    <- acf(samples_mat[, nm], lag.max = 100, plot = FALSE)
    t_idx <- as.integer(str_extract(nm, "\\d+"))
    tibble(Time = t_idx, Lag = as.vector(ac$lag[-1]), ACF = as.vector(ac$acf[-1]))
  }))
  ggsave(
    file.path(scenario_dir, "acf_beta0t.png"),
    ggplot(acf_beta0t_df, aes(x = Lag, y = ACF)) +
      geom_col(width = 0.6, fill = "steelblue") +
      geom_hline(yintercept = c(-0.10,  0.10), linetype = "dashed",
                 color = "blue", linewidth = 0.5) +
      geom_hline(yintercept = c(-0.05,  0.05), linetype = "dotted",
                 color = "red",  linewidth = 0.5) +
      facet_wrap(~ Time, scales = "free_y", ncol = 6,
                 labeller = label_bquote(beta[0][.(Time)])) +
      theme_bw(base_size = 9) +
      labs(title    = paste("ACF de beta0t (", model_type, ")"),
           subtitle = "Azul: |0.10| | Vermelho: |0.05|"),
    width = 14, height = 10
  )

  # Traceplots de kappa (nós determinísticos derivados de incr)
  cores_cadeia <- c("Cadeia 1" = "#2166AC", "Cadeia 2" = "#D6604D")
  df_trace_kappa <- do.call(rbind, lapply(1:nchains, function(ch) {
    cm <- as.matrix(mcmc_list_full[[ch]])
    do.call(rbind, lapply(kappa_names, function(nm) {
      vals <- cm[, nm]
      tibble(Iter = seq_along(vals), Value = vals,
             ErgMedia = cumsum(vals) / seq_along(vals),
             Parameter = nm, Cadeia = paste0("Cadeia ", ch))
    }))
  }))
  ggsave(
    file.path(scenario_dir, "traceplots_kappa.png"),
    ggplot(df_trace_kappa, aes(x = Iter, color = Cadeia)) +
      geom_line(aes(y = Value),    alpha = 0.25, linewidth = 0.20) +
      geom_line(aes(y = ErgMedia), alpha = 0.90, linewidth = 0.75) +
      scale_color_manual(values = cores_cadeia) +
      facet_wrap(~ Parameter, scales = "free_y") +
      theme_bw(base_size = 11) + theme(legend.position = "bottom") +
      labs(title    = paste("Traceplots kappa[k] (", model_type, ")"),
           subtitle = "Linha grossa = média ergódica | Linha fina = cadeia",
           x = "Iteração (pós-burnin)", y = "Valor"),
    width = 10, height = 6
  )

  # ── Retorno resumido ──────────────────────────────────────────────────────
  tibble(
    model = model_type, niter = niter, nburnin = nburnin, thin = thin,
    WAIC = waic, LPML = LPML,
    ESS_beta_min    = min(beta_summary$ESS,    na.rm = TRUE),
    ESS_gamma1_min  = min(gamma1_summary$ESS,  na.rm = TRUE),
    ESS_kappa_min   = min(kappa_summary$ESS,   na.rm = TRUE),
    ESS_beta0t_min  = min(beta0t_summary$ESS,  na.rm = TRUE),
    ESS_tau         = ESS_tau,
    Rhat_beta_max   = max(beta_summary$Rhat,   na.rm = TRUE),
    Rhat_gamma1_max = max(gamma1_summary$Rhat, na.rm = TRUE),
    Rhat_kappa_max  = max(kappa_summary$Rhat,  na.rm = TRUE),
    Rhat_beta0t_max = max(beta0t_summary$Rhat, na.rm = TRUE)
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
  "resultados_gamma_proporcional"
)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

cl <- makeCluster(n_cores)
clusterExport(cl, c(
  "constants_spatial", "constants_nonspatial", "data_nimble",
  "inits_list_spatial", "inits_list_nonspatial",
  "N_regions", "n_times", "p", "K", "h_mat", "cluster_ids",
  "run_model", "output_dir", "REGIONS_INTEREST", "MU_MONITORS",
  "L_seq", "L_max", "delta_g", "gamma_1_max", "gamma_1_min", "M_incr"
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
cat("\n=== RESUMO COMPARATIVO ===\n"); print(resumo)

# ── Comparativo gamma_1[t] ────────────────────────────────────────────────
gamma1_all <- lapply(model_types, function(m) {
  path <- file.path(output_dir, m, "gamma1_summary.csv")
  if (!file.exists(path)) return(NULL)
  read_csv(path, show_col_types = FALSE) %>% mutate(model = m)
}) |> bind_rows()

if (nrow(gamma1_all) > 0) {
  ggsave(
    file.path(output_dir, "gamma1_comparativo.png"),
    ggplot(gamma1_all, aes(x = Time, y = Mean, color = model, fill = model)) +
      geom_ribbon(aes(ymin = Lower, ymax = Upper), alpha = 0.15, color = NA) +
      geom_line(linewidth = 0.9) +
      geom_line(
        data = gamma1_all %>% filter(model == model_types[1]),
        aes(y = L_center), color = "black", linetype = "dashed",
        linewidth = 0.6, inherit.aes = FALSE
      ) +
      theme_bw(base_size = 12) + theme(legend.position = "bottom") +
      labs(title    = "Comparação gamma_1[t]: espacial vs. não-espacial",
           subtitle = "Tracejado = centro do prior L[t]",
           x = "Tempo", y = expression(gamma[1][t]),
           color = "Modelo", fill = "Modelo"),
    width = 8, height = 5, dpi = 300
  )
}

# ── Comparativo kappa[k] ─────────────────────────────────────────────────
kappa_all <- lapply(model_types, function(m) {
  path <- file.path(output_dir, m, "kappa_summary.csv")
  if (!file.exists(path)) return(NULL)
  read_csv(path, show_col_types = FALSE) %>% mutate(model = m)
}) |> bind_rows()

if (nrow(kappa_all) > 0) {
  ggsave(
    file.path(output_dir, "kappa_comparativo.png"),
    ggplot(kappa_all, aes(x = factor(k), y = Mean, color = model)) +
      geom_point(size = 3, position = position_dodge(0.4)) +
      geom_errorbar(aes(ymin = Lower, ymax = Upper), width = 0.3,
                    linewidth = 0.7, position = position_dodge(0.4)) +
      geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
      theme_bw(base_size = 12) + theme(legend.position = "bottom") +
      labs(title = "Comparação kappa[k]: espacial vs. não-espacial",
           x = "Cluster k", y = expression(kappa[k]),
           color = "Modelo"),
    width = 7, height = 5, dpi = 300
  )
}

# ── Comparativo beta0t ────────────────────────────────────────────────────
beta0t_all <- lapply(model_types, function(m) {
  path <- file.path(output_dir, m, "beta0t_summary.csv")
  if (!file.exists(path)) return(NULL)
  read_csv(path, show_col_types = FALSE) %>% mutate(model = m)
}) |> bind_rows()

if (nrow(beta0t_all) > 0) {
  anos_label <- colnames(data_nimble$Y)
  if (is.null(anos_label)) anos_label <- seq_len(n_times)
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

# ── Comparativo mu ────────────────────────────────────────────────────────
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

# ── Comparativo epsilon ───────────────────────────────────────────────────
eps_all <- lapply(model_types, function(m) {
  path <- file.path(output_dir, m, "epsilon_summary.csv")
  if (!file.exists(path)) return(NULL)
  read_csv(path, show_col_types = FALSE) %>%
    filter(Region %in% REGIONS_INTEREST) %>%
    mutate(model = m)
}) |> bind_rows()

if (nrow(eps_all) > 0) {
  ggsave(
    file.path(output_dir, "epsilon_comparativo.png"),
    ggplot(eps_all, aes(x = Time, y = Mean, color = model, fill = model)) +
      geom_ribbon(aes(ymin = Lower, ymax = Upper), alpha = 0.15, color = NA) +
      geom_line(linewidth = 0.8) +
      facet_wrap(~ Region, scales = "free_y", ncol = 1) +
      theme_bw(base_size = 12) + theme(legend.position = "bottom") +
      labs(title = "Comparação epsilon: espacial vs. não-espacial",
           x = "Tempo", y = expression(epsilon[it]),
           color = "Modelo", fill = "Modelo"),
    width = 10, height = 3 * length(REGIONS_INTEREST), dpi = 300
  )
}

cat("\nTempo total de execução:\n")
print(Sys.time() - inicio_global)
