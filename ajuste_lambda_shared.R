# ==============================================================================
# ajuste_lambda_shared.R
# Modelo Bayesiano Espaço-Estados — dados reais MG
# A = 75 microrregiões | T = 23 anos | K = 4 clusters | p = 3 covariáveis
#
# DIFERENÇA-CHAVE em relação a ajuste_dados_reais.R:
#   Em vez de lambda[i,t] (efeito latente por região × tempo),
#   usamos lambda[t] (efeito latente COMPARTILHADO entre regiões, variando no tempo).
#
#   Equação: log(mu[i,t]) = log(lambda[t]) + log(E[i,t]) + log(epsilon[i])
#                           + beta' x[i,t] + s[i]    (s[i] apenas no modelo espacial)
#
# FFBS adaptado: no passo de filtragem, a suficiência de lambda[t] é
#   obtida somando Y[i,t] e g[i,t] sobre todas as regiões i a cada tempo t.
#   Forward:  a[t+1] = w * a[t] + sum_i Y[i,t]
#             b[t+1] = w * b[t] + sum_i g[i,t]
#   Backward: idem ao caso original (Gamma conjugado).
#
# Modelos: espacial (ICAR) vs. não-espacial — execução paralela em 2 núcleos.
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

# ── Verificações de consistência ───────────────────────────────────────────────
stopifnot(
  "Y_mat deve ter 75 linhas"      = nrow(Y_mat) == 75,
  "E deve ter 75 linhas"          = nrow(E)     == 75,
  "x deve ter dimensão correta"   = all(dim(x) == c(75, 23, 3)),
  "Y_mat e E devem ter mesma dim" = identical(dim(Y_mat), dim(E)),
  "linhas de Y_mat e E alinhadas" = identical(row.names(Y_mat), row.names(E)),
  "linhas de Y_mat e x alinhadas" = identical(row.names(Y_mat), dimnames(x)[[1]]),
  "estrutura espacial consistente"= sum(data$num) == data$sumNumNeigh
)
cat("Verificacoes de consistencia: OK\n")

# ── Normalizar E (estabilidade numérica) ──────────────────────────────────────
E_norm <- E / mean(E)

# ── Dimensões ──────────────────────────────────────────────────────────────────
N_regions <- nrow(Y_mat)    # 75
n_times   <- ncol(Y_mat)    # 23
p         <- dim(x)[3]      # 3
K         <- ncol(data$hAI) # 4

cat(sprintf("N = %d | T = %d | p = %d | K = %d\n", N_regions, n_times, p, K))

# ── Estrutura espacial ────────────────────────────────────────────────────────
adj_vec     <- as.integer(data$adj)
num_vec     <- as.integer(data$num)
n_adj_val   <- as.integer(data$sumNumNeigh)
weights_vec <- rep(1.0, n_adj_val)
h_mat       <- data$hAI   # 75 × 4

# ==============================================================================
# 2. SELEÇÃO DE REGIÕES PARA MONITORAMENTO DE THETA E MU
#    (lambda[t] já é monitorado para todos os t — apenas 23 nós)
# ==============================================================================
cluster_ids <- apply(h_mat, 1, sum)   # cluster dominante de cada região
set.seed(42)
REGIONS_INTEREST <- unlist(lapply(1:K, function(cl) {
  regs <- which(cluster_ids == cl)
  if (length(regs) >= 3) sample(regs, 3) else regs
}))
cat("Regioes selecionadas para monitoramento de theta/mu:\n")
print(sort(REGIONS_INTEREST))

# Monitores de lambda: apenas os T tempos (sem dimensão de região)
LAMBDA_MONITORS <- paste0("lambda[", seq_len(n_times), "]")

# ==============================================================================
# 3. CONSTANTES E DADOS NIMBLE
# ==============================================================================
constants_spatial <- list(
  n_regions = N_regions,
  n_times   = n_times,
  p         = p,
  K         = K,
  h         = h_mat,
  mu_beta   = rep(0, p),
  a_unif    = 0.0,
  b_unif    = 0.1,
  a0        = 1.0,
  b0        = 1.0,
  w         = 0.9,
  adj       = adj_vec,
  num       = num_vec,
  weights   = weights_vec,
  n_adj     = n_adj_val
)

constants_nonspatial <- constants_spatial[
  setdiff(names(constants_spatial), c("adj", "num", "weights", "n_adj"))
]

data_nimble <- list(
  Y = Y_mat,
  E = E_norm,
  x = x
)

# ==============================================================================
# 4. INICIALIZAÇÕES
#    lambda agora é um vetor de comprimento n_times (não mais uma matriz)
# ==============================================================================
set.seed(123)

inits_list_spatial <- list(
  list(
    beta    = rep(0, p),
    gamma   = c(0.05, 0.10, 0.10, 0.15),
    lambda  = rep(1.0, n_times),        # <-- vetor T
    sigma_s = 0.5,
    s       = rep(0, N_regions)
  ),
  list(
    beta    = rnorm(p, 0, 0.3),
    gamma   = c(0.04, 0.09, 0.09, 0.14),
    lambda  = rgamma(n_times, 1, 1),    # <-- vetor T
    sigma_s = 1.0,
    s       = rep(0, N_regions)
  )
)

inits_list_nonspatial <- list(
  list(
    beta   = rep(0, p),
    gamma  = c(0.05, 0.10, 0.10, 0.15),
    lambda = rep(1.0, n_times)          # <-- vetor T
  ),
  list(
    beta   = rnorm(p, 0, 0.3),
    gamma  = c(0.04, 0.09, 0.09, 0.14),
    lambda = rgamma(n_times, 1, 1)      # <-- vetor T
  )
)

# ==============================================================================
# 5. FUNÇÃO WORKER — roda um modelo em processo separado
# ==============================================================================
run_model <- function(model_type, output_dir) {
  
  library(nimble); library(coda); library(dplyr)
  library(ggplot2); library(readr); library(stringr); library(tibble)
  
  # ── 5a. Código NIMBLE ────────────────────────────────────────────────────────
  #
  # Mudança principal: lambda é indexado apenas em t, não em (i,t).
  # O loop externo é sobre t (lambda[t] é compartilhado entre regiões).
  # ─────────────────────────────────────────────────────────────────────────────
  
  code_spatial <- nimbleCode({
    # Priori dos parâmetros estruturais
    for (j in 1:p) beta[j] ~ dnorm(mu_beta[j], sd = 1)
    gamma[1] ~ dunif(min = a_unif, max = b_unif)
    for (j in 2:K) gamma[j] ~ dunif(min = 0, max = (1 - sum(gamma[1:(j-1)])))
    
    # Efeito espacial ICAR
    sigma_s ~ T(dt(0, 1, 1), 0, )
    tau_s   <- 1 / (sigma_s^2)
    s[1:n_regions] ~ dcar_normal(
      adj[1:n_adj], weights[1:n_adj], num[1:n_regions], tau_s, zero_mean = 1
    )
    
    # lambda[t]: efeito latente temporal compartilhado entre TODAS as regiões
    for (t in 1:n_times) {
      lambda[t] ~ dgamma(1, 1)
    }
    
    # Verossimilhança
    for (i in 1:n_regions) {
      epsilon[i] <- 1 - inprod(h[i, 1:K], gamma[1:K])
      for (t in 1:n_times) {
        log(mu[i, t]) <- log(lambda[t]) + log(E[i, t]) + log(epsilon[i]) +
          inprod(beta[1:p], x[i, t, 1:p]) + s[i]
        Y[i, t]        ~ dpois(mu[i, t])
        logLik_Y[i, t] <- dpois(Y[i, t], mu[i, t], log = TRUE)
      }
    }
  })
  
  code_nonspatial <- nimbleCode({
    for (j in 1:p) beta[j] ~ dnorm(mu_beta[j], sd = 1)
    gamma[1] ~ dunif(min = a_unif, max = b_unif)
    for (j in 2:K) gamma[j] ~ dunif(min = 0, max = (1 - sum(gamma[1:(j-1)])))
    
    for (t in 1:n_times) {
      lambda[t] ~ dgamma(1, 1)
    }
    
    for (i in 1:n_regions) {
      epsilon[i] <- 1 - inprod(h[i, 1:K], gamma[1:K])
      for (t in 1:n_times) {
        log(mu[i, t]) <- log(lambda[t]) + log(E[i, t]) + log(epsilon[i]) +
          inprod(beta[1:p], x[i, t, 1:p])
        Y[i, t]        ~ dpois(mu[i, t])
        logLik_Y[i, t] <- dpois(Y[i, t], mu[i, t], log = TRUE)
      }
    }
  })
  
  # ── 5b. FFBS espacial — lambda[t] compartilhado ──────────────────────────────
  #
  # Forward: para cada tempo t, acumula suficiência somando sobre todas as
  #          regiões i: sum_i Y[i,t]  e  sum_i g[i,t].
  #   g[i,t] = E[i,t] * epsilon[i] * exp(beta' x[i,t] + s[i])
  #
  # Modelo dinâmico conjugado:
  #   lambda[t] | D_{t-1} ~ Gamma(a_{t-1}^pred, b_{t-1}^pred)
  #     com  a^pred = w * a,  b^pred = w * b   (desconto)
  #   Atualização:
  #     a[t] = w * a[t-1] + sum_i Y[i,t]
  #     b[t] = w * b[t-1] + sum_i g[i,t]
  #
  # Backward:
  #   lambda[T] ~ Gamma(a[T], b[T])
  #   nu ~ Gamma((1-w)*a[t], b[t])  → lambda[t] = nu + w * lambda[t+1]
  # ─────────────────────────────────────────────────────────────────────────────
  
  ffbs_spatial <- nimbleFunction(
    contains = sampler_BASE,
    setup = function(model, mvSaved, target, control) {
      n_regions <- control$n_regions
      n_times   <- control$n_times
      p         <- control$p
      a0        <- control$a0
      b0        <- control$b0
      w         <- control$w
      # Buffers de tamanho T+1
      at_buf <- nimNumeric(n_times + 1, 0)
      bt_buf <- nimNumeric(n_times + 1, 0)
      calcNodes   <- model$getDependencies(target, self = FALSE)
      targetNodes <- model$expandNodeNames(target)
      setupOutputs(at_buf, bt_buf)
    },
    run = function() {
      declare(i,        integer())
      declare(t,        integer())
      declare(t_idx,    integer())
      declare(t_back,   integer())
      declare(k,        integer())
      declare(prod_val, double())
      declare(g_it,     double())
      declare(sum_Y_t,  double())
      declare(sum_g_t,  double())
      declare(nu,       double())
      
      # ── Passo forward ─────────────────────────────────────────────────────
      at_buf[1] <<- a0
      bt_buf[1] <<- b0
      
      for (t in 1:n_times) {
        sum_Y_t <- 0
        sum_g_t <- 0
        for (i in 1:n_regions) {
          sum_Y_t <- sum_Y_t + model$Y[i, t]
          prod_val <- 0
          for (k in 1:p) prod_val <- prod_val + model$x[i, t, k] * model$beta[k]
          # g[i,t] inclui efeito espacial s[i]
          g_it    <- model$E[i, t] * model$epsilon[i] * exp(prod_val + model$s[i])
          sum_g_t <- sum_g_t + g_it
        }
        at_buf[t + 1] <<- w * at_buf[t] + sum_Y_t
        bt_buf[t + 1] <<- w * bt_buf[t] + sum_g_t
      }
      
      # ── Passo backward ────────────────────────────────────────────────────
      model$lambda[n_times] <<- rgamma(1, shape = at_buf[n_times + 1],
                                       rate  = bt_buf[n_times + 1])
      for (t_idx in 1:(n_times - 1)) {
        t_back <- n_times - t_idx
        nu <- rgamma(1, shape = (1 - w) * at_buf[t_back + 1],
                     rate  = bt_buf[t_back + 1])
        model$lambda[t_back] <<- nu + w * model$lambda[t_back + 1]
      }
      
      model$calculate(calcNodes)
      copy(from = model, to = mvSaved, row = 1,
           nodes = targetNodes, logProb = TRUE)
    },
    methods = list(reset = function() {})
  )
  
  # ── 5c. FFBS não-espacial — lambda[t] compartilhado ──────────────────────────
  # Idêntico ao espacial, mas sem s[i] no cálculo de g[i,t].
  
  ffbs_nonspatial <- nimbleFunction(
    contains = sampler_BASE,
    setup = function(model, mvSaved, target, control) {
      n_regions <- control$n_regions
      n_times   <- control$n_times
      p         <- control$p
      a0        <- control$a0
      b0        <- control$b0
      w         <- control$w
      at_buf <- nimNumeric(n_times + 1, 0)
      bt_buf <- nimNumeric(n_times + 1, 0)
      calcNodes   <- model$getDependencies(target, self = FALSE)
      targetNodes <- model$expandNodeNames(target)
      setupOutputs(at_buf, bt_buf)
    },
    run = function() {
      declare(i,        integer())
      declare(t,        integer())
      declare(t_idx,    integer())
      declare(t_back,   integer())
      declare(k,        integer())
      declare(prod_val, double())
      declare(g_it,     double())
      declare(sum_Y_t,  double())
      declare(sum_g_t,  double())
      declare(nu,       double())
      
      at_buf[1] <<- a0
      bt_buf[1] <<- b0
      
      for (t in 1:n_times) {
        sum_Y_t <- 0
        sum_g_t <- 0
        for (i in 1:n_regions) {
          sum_Y_t <- sum_Y_t + model$Y[i, t]
          prod_val <- 0
          for (k in 1:p) prod_val <- prod_val + model$x[i, t, k] * model$beta[k]
          g_it    <- model$E[i, t] * model$epsilon[i] * exp(prod_val)  # sem s[i]
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
      copy(from = model, to = mvSaved, row = 1,
           nodes = targetNodes, logProb = TRUE)
    },
    methods = list(reset = function() {})
  )
  
  # ── 5d. Seleção de objetos conforme model_type ───────────────────────────────
  is_spatial <- (model_type == "spatial")
  model_code <- if (is_spatial) code_spatial      else code_nonspatial
  constants  <- if (is_spatial) constants_spatial  else constants_nonspatial
  inits_list <- if (is_spatial) inits_list_spatial else inits_list_nonspatial
  ffbs_fn    <- if (is_spatial) ffbs_spatial       else ffbs_nonspatial
  
  cat("\n=== Iniciando modelo:", model_type, "===\n")
  scenario_dir <- file.path(output_dir, model_type)
  dir.create(scenario_dir, recursive = TRUE, showWarnings = FALSE)
  
  # ── 5e. Construir e compilar ─────────────────────────────────────────────────
  model  <- nimbleModel(
    code      = model_code,
    constants = constants,
    data      = data_nimble,
    inits     = inits_list[[1]],
    check     = FALSE
  )
  Cmodel <- compileNimble(model)
  
  conf <- configureMCMC(model)
  conf$removeSamplers("lambda")
  conf$addSampler(
    target  = "lambda",
    type    = ffbs_fn,
    control = list(
      n_regions = N_regions, n_times = n_times, p = p,
      a0 = constants$a0, b0 = constants$b0, w = constants$w
    )
  )
  conf$removeSamplers("gamma")
  conf$addSampler(target = "gamma", type = "AF_slice")
  
  monitors_base <- c("beta", "gamma", "logLik_Y", "lambda")   # lambda[t]: 23 nós
  if (is_spatial) monitors_base <- c(monitors_base, "s", "sigma_s", "tau_s")
  conf$addMonitors(monitors_base)
  conf$printSamplers()
  
  Rmcmc <- buildMCMC(conf)
  Cmcmc <- compileNimble(Rmcmc, project = model)
  
  # ── 5f. MCMC ─────────────────────────────────────────────────────────────────
  niter   <- 50000
  nburnin <- 10000
  nchains <- 2
  thin    <- 10
  
  cat(sprintf("[%s] niter=%d | nburnin=%d | thin=%d | cadeias=%d\n",
              model_type, niter, nburnin, thin, nchains))
  
  samples <- runMCMC(
    Cmcmc,
    niter             = niter,
    nburnin           = nburnin,
    nchains           = nchains,
    thin              = thin,
    inits             = inits_list,
    samplesAsCodaMCMC = TRUE,
    summary           = FALSE,
    WAIC              = FALSE
  )
  
  saveRDS(samples, file.path(scenario_dir, "samples.rds"))
  cat("[", model_type, "] Amostras salvas.\n")
  
  samples_mat    <- as.matrix(samples)
  mcmc_list_full <- mcmc.list(lapply(1:nchains, function(ch) as.mcmc(samples[[ch]])))
  rm(samples); gc()
  
  # ── 5g. Funções auxiliares ───────────────────────────────────────────────────
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
  
  beta_names   <- paste0("beta[",   1:p, "]")
  gamma_names  <- paste0("gamma[",  1:K, "]")
  lambda_names <- paste0("lambda[", 1:n_times, "]")   # vetor T, não matriz
  
  # ── 5h. Epsilon posterior ────────────────────────────────────────────────────
  gamma_draws   <- samples_mat[, gamma_names, drop = FALSE]
  epsilon_draws <- 1 - gamma_draws %*% t(h_mat)   # n_draw × N
  
  epsilon_summary <- do.call(rbind, lapply(1:N_regions, function(i) {
    hpd <- safe_hpd(epsilon_draws[, i])
    tibble(Region    = i,
           Cluster   = cluster_ids[i],
           Eps_Mean  = mean(epsilon_draws[, i]),
           Eps_Lower = hpd[1],
           Eps_Upper = hpd[2])
  }))
  write_csv(epsilon_summary, file.path(scenario_dir, "epsilon_summary.csv"))
  
  make_region_label <- function(regions_vec) {
    setNames(
      sprintf("Reg %d (C%d)\ne=%.3f",
              regions_vec,
              cluster_ids[regions_vec],
              epsilon_summary$Eps_Mean[regions_vec]),
      as.character(regions_vec)
    )
  }
  
  # ── 5i. Sumário de beta, gamma ───────────────────────────────────────────────
  resumo_param <- function(nm_vec) {
    do.call(rbind, lapply(nm_vec, function(nm) {
      sv  <- samples_mat[, nm]
      hpd <- safe_hpd(sv)
      tibble(Parameter = nm,
             Mean      = mean(sv),
             SD        = sd(sv),
             HPD_Lower = hpd[1],
             HPD_Upper = hpd[2],
             ESS       = as.numeric(effectiveSize(mcmc_list_full[, nm])),
             Rhat      = safe_gelman(mcmc_list_full[, nm]))
    }))
  }
  
  beta_summary  <- resumo_param(beta_names)
  gamma_summary <- resumo_param(gamma_names)
  write_csv(beta_summary,  file.path(scenario_dir, "beta_summary.csv"))
  write_csv(gamma_summary, file.path(scenario_dir, "gamma_summary.csv"))
  
  # ── 5j. tau_s (apenas espacial) ──────────────────────────────────────────────
  ESS_tau <- NA_real_
  if (is_spatial) {
    tau_sv  <- samples_mat[, "tau_s"]
    hpd_t   <- safe_hpd(tau_sv)
    tau_sum <- tibble(Parameter = "tau_s",
                      Mean      = mean(tau_sv),
                      SD        = sd(tau_sv),
                      HPD_Lower = hpd_t[1],
                      HPD_Upper = hpd_t[2],
                      ESS       = as.numeric(effectiveSize(mcmc_list_full[, "tau_s"])),
                      Rhat      = safe_gelman(mcmc_list_full[, "tau_s"]))
    write_csv(tau_sum, file.path(scenario_dir, "tau_summary.csv"))
    ESS_tau <- tau_sum$ESS
    
    s_names   <- paste0("s[", 1:N_regions, "]")
    s_summary <- do.call(rbind, lapply(1:N_regions, function(i) {
      sv  <- samples_mat[, s_names[i]]
      hpd <- safe_hpd(sv)
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
  
  # ── 5k. lambda[t] — sumário completo para todos os tempos ────────────────────
  #
  # Como lambda[t] é COMPARTILHADO, seu sumário tem apenas T linhas
  # (uma por tempo), sem dimensão de região.
  # ─────────────────────────────────────────────────────────────────────────────
  lambda_summary <- do.call(rbind, lapply(seq_len(n_times), function(t) {
    nm  <- lambda_names[t]
    sv  <- samples_mat[, nm]
    hpd <- safe_hpd(sv)
    tibble(Time  = t,
           Mean  = mean(sv),
           SD    = sd(sv),
           Lower = hpd[1],
           Upper = hpd[2],
           ESS   = as.numeric(effectiveSize(mcmc_list_full[, nm])),
           Rhat  = safe_gelman(mcmc_list_full[, nm]),
           model = model_type)
  }))
  write_csv(lambda_summary, file.path(scenario_dir, "lambda_summary.csv"))
  
  # Gráfico da trajetória temporal de lambda[t]
  ggsave(
    file.path(scenario_dir, "lambda_temporal.png"),
    ggplot(lambda_summary, aes(x = Time)) +
      geom_ribbon(aes(ymin = Lower, ymax = Upper), fill = "grey70", alpha = 0.5) +
      geom_line(aes(y = Mean), color = "black", linewidth = 0.8) +
      theme_bw(base_size = 12) +
      labs(title = paste("Lambda compartilhado lambda[t] (", model_type, ")"),
           x = "Tempo", y = expression(lambda[t])),
    width = 9, height = 5
  )
  
  # ── 5l. theta[i,t] e mu[i,t] para regiões selecionadas ──────────────────────
  #
  # theta[i,t] = lambda[t] * exp(beta' x[i,t])
  # mu[i,t]    = lambda[t] * exp(beta' x[i,t]) * E[i,t] * epsilon[i]
  #
  # lambda[t] é o mesmo para todas as regiões; a variação regional vem de
  # x[i,t], E[i,t] e epsilon[i].
  # ─────────────────────────────────────────────────────────────────────────────
  beta_cols <- grep("^beta\\[", colnames(samples_mat), value = TRUE)
  
  theta_summary <- do.call(rbind, lapply(REGIONS_INTEREST, function(i) {
    do.call(rbind, lapply(seq_len(n_times), function(t) {
      ldraws <- samples_mat[, lambda_names[t]]
      bdraws <- samples_mat[, beta_cols, drop = FALSE]
      x_it   <- data_nimble$x[i, t, ]
      theta  <- ldraws * exp(as.vector(bdraws %*% x_it))
      hpd    <- safe_hpd(theta)
      tibble(Region = i, Time = t, Mean = mean(theta),
             Lower = hpd[1], Upper = hpd[2], model = model_type)
    }))
  }))
  write_csv(theta_summary, file.path(scenario_dir, "theta_selected.csv"))
  
  mu_summary <- do.call(rbind, lapply(REGIONS_INTEREST, function(i) {
    do.call(rbind, lapply(seq_len(n_times), function(t) {
      ldraws         <- samples_mat[, lambda_names[t]]
      bdraws         <- samples_mat[, beta_cols, drop = FALSE]
      x_it           <- data_nimble$x[i, t, ]
      epsilon_i_draw <- epsilon_draws[, i]
      mu_draws <- ldraws * exp(as.vector(bdraws %*% x_it)) *
        data_nimble$E[i, t] * epsilon_i_draw
      hpd <- safe_hpd(mu_draws)
      tibble(Region = i, Time = t, Mean = mean(mu_draws),
             Lower = hpd[1], Upper = hpd[2], model = model_type)
    }))
  }))
  write_csv(mu_summary, file.path(scenario_dir, "mu_selected.csv"))
  
  # ── 5m. WAIC e LPML ───────────────────────────────────────────────────────────
  loglik_names <- grep("logLik_Y", colnames(samples_mat), value = TRUE)
  waic <- NA_real_; LPML <- NA_real_
  if (length(loglik_names) > 0) {
    lm     <- samples_mat[, loglik_names, drop = FALSE]
    lppd   <- sum(apply(lm, 2, function(x) { mx <- max(x); mx + log(mean(exp(x - mx))) }))
    p_waic <- sum(apply(lm, 2, var))
    waic   <- -2 * (lppd - p_waic)
    LPML   <- sum(log(1 / apply(lm, 2, function(x) mean(exp(-x)))))
    write_csv(
      tibble(WAIC = waic, LPML = LPML, lppd = lppd, pWAIC = p_waic),
      file.path(scenario_dir, "criteria.csv")
    )
    cat(sprintf("[%s] WAIC = %.2f | LPML = %.2f\n", model_type, waic, LPML))
  }
  
  # ── 5n. Diagnósticos ACF e ESS (Padrão Unificado) ────────────────────────────
  params_beta   <- c(beta_names, if (is_spatial) "tau_s" else character(0))
  params_gamma  <- gamma_names
  params_struct <- c(params_beta, params_gamma, lambda_names)
  
  ESS_struct  <- effectiveSize(mcmc_list_full[, params_struct])
  Rhat_struct <- safe_gelman(mcmc_list_full[, params_struct])
  
  acf_results <- do.call(rbind, lapply(params_struct, function(nm) {
    ac   <- acf(samples_mat[, nm], lag.max = 200, plot = FALSE)
    lags <- as.vector(ac$lag[-1]); acfs <- as.vector(ac$acf[-1])
    tibble(
      Parameter = nm, ESS = ESS_struct[nm], Rhat = Rhat_struct[nm],
      lag_0.10  = { v <- lags[which(abs(acfs) < 0.10)[1]]; ifelse(is.na(v), Inf, v) },
      lag_0.05  = { v <- lags[which(abs(acfs) < 0.05)[1]]; ifelse(is.na(v), Inf, v) },
      acf_lag1  = acfs[1]
    )
  }))
  write_csv(acf_results, file.path(scenario_dir, "acf_diagnostics.csv"))
  
  # ── 5o. Gráficos ACF Separados ───────────────────────────────────────────────
  
  # 1. ACF de Beta e Tau_s
  acf_beta_df <- do.call(rbind, lapply(params_beta, function(nm) {
    ac <- acf(samples_mat[, nm], lag.max = 100, plot = FALSE)
    tibble(Parameter = nm, Lag = as.vector(ac$lag[-1]), ACF = as.vector(ac$acf[-1]))
  }))
  ggsave(file.path(scenario_dir, "acf_beta.png"),
         ggplot(acf_beta_df, aes(x = Lag, y = ACF)) +
           geom_col(width = 0.6, fill = "grey50") +
           geom_hline(yintercept = c(-0.10, 0.10), linetype = "dashed", color = "blue",  linewidth = 0.5) +
           geom_hline(yintercept = c(-0.05, 0.05), linetype = "dotted", color = "red",   linewidth = 0.5) +
           facet_wrap(~Parameter, scales = "free_y") +
           theme_bw(base_size = 11) +
           labs(title = paste("ACF de beta (", model_type, ")"), subtitle = "Azul: |0.10| | Vermelho: |0.05|"),
         width = 10, height = 5)
  
  # 2. ACF de Gamma
  acf_gamma_df <- do.call(rbind, lapply(params_gamma, function(nm) {
    ac <- acf(samples_mat[, nm], lag.max = 100, plot = FALSE)
    tibble(Parameter = nm, Lag = as.vector(ac$lag[-1]), ACF = as.vector(ac$acf[-1]))
  }))
  ggsave(file.path(scenario_dir, "acf_gamma.png"),
         ggplot(acf_gamma_df, aes(x = Lag, y = ACF)) +
           geom_col(width = 0.6, fill = "darkorange") +
           geom_hline(yintercept = c(-0.10, 0.10), linetype = "dashed", color = "blue",  linewidth = 0.5) +
           geom_hline(yintercept = c(-0.05, 0.05), linetype = "dotted", color = "red",   linewidth = 0.5) +
           facet_wrap(~Parameter, scales = "free_y") +
           theme_bw(base_size = 11) +
           labs(title = paste("ACF de gamma (", model_type, ")"), subtitle = "Azul: |0.10| | Vermelho: |0.05|"),
         width = 10, height = 5)
  
  # 3. ACF de Lambda[t] (Mantido inalterado com cor steelblue)
  acf_lambda_df <- do.call(rbind, lapply(lambda_names, function(nm) {
    ac <- acf(samples_mat[, nm], lag.max = 100, plot = FALSE)
    tibble(Time = as.integer(str_extract(nm, "\\d+")),
           Lag  = as.vector(ac$lag[-1]), ACF  = as.vector(ac$acf[-1]))
  }))
  ggsave(file.path(scenario_dir, "acf_lambda.png"),
         ggplot(acf_lambda_df, aes(x = Lag, y = ACF)) +
           geom_col(width = 0.6, fill = "steelblue") +
           geom_hline(yintercept = c(-0.10, 0.10), linetype = "dashed", color = "blue",  linewidth = 0.5) +
           geom_hline(yintercept = c(-0.05, 0.05), linetype = "dotted", color = "red",   linewidth = 0.5) +
           facet_wrap(~Time, scales = "free_y", ncol = 6, labeller = label_bquote(lambda[.(Time)])) +
           theme_bw(base_size = 9) +
           labs(title = paste("ACF de lambda[t] (", model_type, ")"), subtitle = "Azul: |0.10| | Vermelho: |0.05|"),
         width = 14, height = 10)
  
  # ── 5p. Traceplots Separados ─────────────────────────────────────────────────
  cores_cadeia <- c("Cadeia 1" = "#2166AC", "Cadeia 2" = "#D6604D")
  
  # Traceplots Beta
  df_trace_beta <- do.call(rbind, lapply(seq_len(nchains), function(ch) {
    cm <- as.matrix(mcmc_list_full[[ch]])
    do.call(rbind, lapply(params_beta, function(nm) {
      vals <- cm[, nm]; erg_mean <- cumsum(vals) / seq_along(vals)
      tibble(Iter = seq_along(vals), Value = vals, ErgMedia = erg_mean,
             Parameter = nm, Cadeia = paste0("Cadeia ", ch))
    }))
  }))
  ggsave(file.path(scenario_dir, "traceplots_beta.png"),
         ggplot(df_trace_beta, aes(x = Iter, color = Cadeia, fill = Cadeia)) +
           geom_line(aes(y = Value),    alpha = 0.25, linewidth = 0.20) +
           geom_line(aes(y = ErgMedia), alpha = 0.90, linewidth = 0.75) +
           scale_color_manual(values = cores_cadeia) + scale_fill_manual(values  = cores_cadeia) +
           facet_wrap(~Parameter, scales = "free_y") +
           theme_bw(base_size = 11) + theme(legend.position = "bottom") +
           labs(title = paste("Traceplots beta (", model_type, ")"),
                subtitle = "Linha grossa = media ergodica | Linha fina = cadeia",
                x = "Iteracao (pos-burnin)", y = "Valor"),
         width  = 10, height = max(5, 3 * ceiling(length(params_beta) / 3)))
  
  # Traceplots Gamma
  df_trace_gamma <- do.call(rbind, lapply(seq_len(nchains), function(ch) {
    cm <- as.matrix(mcmc_list_full[[ch]])
    do.call(rbind, lapply(params_gamma, function(nm) {
      vals <- cm[, nm]; erg_mean <- cumsum(vals) / seq_along(vals)
      tibble(Iter = seq_along(vals), Value = vals, ErgMedia = erg_mean,
             Parameter = nm, Cadeia = paste0("Cadeia ", ch))
    }))
  }))
  ggsave(file.path(scenario_dir, "traceplots_gamma.png"),
         ggplot(df_trace_gamma, aes(x = Iter, color = Cadeia, fill = Cadeia)) +
           geom_line(aes(y = Value),    alpha = 0.25, linewidth = 0.20) +
           geom_line(aes(y = ErgMedia), alpha = 0.90, linewidth = 0.75) +
           scale_color_manual(values = cores_cadeia) + scale_fill_manual(values  = cores_cadeia) +
           facet_wrap(~Parameter, scales = "free_y") +
           theme_bw(base_size = 11) + theme(legend.position = "bottom") +
           labs(title = paste("Traceplots gamma (", model_type, ")"),
                subtitle = "Linha grossa = media ergodica | Linha fina = cadeia",
                x = "Iteracao (pos-burnin)", y = "Valor"),
         width  = 10, height = max(5, 3 * ceiling(length(params_gamma) / 3)))
  
  # ── 5q. Painel theta e mu por região ──────────────────────────────────────────
  theta_regs <- sort(unique(theta_summary$Region))
  ggsave(
    file.path(scenario_dir, "painel_theta.png"),
    ggplot(theta_summary, aes(x = Time)) +
      geom_ribbon(aes(ymin = Lower, ymax = Upper), fill = "steelblue", alpha = 0.3) +
      geom_line(aes(y = Mean), color = "steelblue") +
      facet_wrap(~Region, scales = "free_y", ncol = 3,
                 labeller = labeller(Region = make_region_label(theta_regs))) +
      theme_bw(base_size = 10) +
      labs(title = paste("Theta estimado (", model_type, ")"),
           x = "Tempo", y = expression(theta[i * t])),
    width = 12, height = 10
  )
  
  mu_regs <- sort(unique(mu_summary$Region))
  ggsave(
    file.path(scenario_dir, "painel_mu.png"),
    ggplot(mu_summary, aes(x = Time)) +
      geom_ribbon(aes(ymin = Lower, ymax = Upper), fill = "darkorange", alpha = 0.3) +
      geom_line(aes(y = Mean), color = "darkorange") +
      facet_wrap(~Region, scales = "free_y", ncol = 3,
                 labeller = labeller(Region = make_region_label(mu_regs))) +
      theme_bw(base_size = 10) +
      labs(title = paste("Mu estimado (", model_type, ")"),
           x = "Tempo", y = expression(mu[i * t])),
    width = 12, height = 10
  )
  
  # ── 5r. Epsilon posterior ─────────────────────────────────────────────────────
  ggsave(
    file.path(scenario_dir, "epsilon_posterior.png"),
    ggplot(epsilon_summary, aes(x = Region, y = Eps_Mean,
                                color = factor(Cluster))) +
      geom_point(size = 0.9) +
      geom_errorbar(aes(ymin = Eps_Lower, ymax = Eps_Upper),
                    width = 0.3, linewidth = 0.3) +
      scale_color_discrete(name = "Cluster") +
      theme_bw() +
      labs(title = paste("Epsilon posterior (", model_type, ")"),
           y = expression(epsilon[i]), x = "Regiao"),
    width = 10, height = 4
  )
  
  # ── 5s. Retorno resumido ──────────────────────────────────────────────────────
  params_diag <- c(beta_names, gamma_names, if (is_spatial) "tau_s" else character(0))
  tibble(
    model           = model_type,
    niter           = niter,
    nburnin         = nburnin,
    thin            = thin,
    WAIC            = waic,
    LPML            = LPML,
    ESS_beta_min    = min(beta_summary$ESS,    na.rm = TRUE),
    ESS_gamma_min   = min(gamma_summary$ESS,   na.rm = TRUE),
    ESS_lambda_min  = min(lambda_summary$ESS,  na.rm = TRUE),
    ESS_tau         = ESS_tau,
    Rhat_max        = max(Rhat_struct[params_diag], na.rm = TRUE),
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
  "resultados_lambda_shared"
)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

cl <- makeCluster(n_cores)
clusterExport(cl, c(
  "constants_spatial", "constants_nonspatial", "data_nimble",
  "inits_list_spatial", "inits_list_nonspatial",
  "N_regions", "n_times", "p", "K",
  "h_mat", "cluster_ids",
  "run_model", "output_dir",
  "REGIONS_INTEREST", "LAMBDA_MONITORS"
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

cat("\n=== RESUMO COMPARATIVO ===\n")
print(resumo)

# ── Gráfico comparativo de lambda[t] entre modelos ────────────────────────────
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
      theme_bw(base_size = 12) +
      theme(legend.position = "bottom") +
      labs(title = "Comparacao de lambda[t]: espacial vs. nao-espacial",
           x = "Tempo (ano)", y = expression(lambda[t]),
           color = "Modelo", fill = "Modelo"),
    width = 10, height = 5, dpi = 300
  )
}

# ── Gráficos comparativos de theta e mu ───────────────────────────────────────
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
        facet_wrap(~Region, scales = "free_y", ncol = 3) +
        theme_bw(base_size = 12) +
        theme(legend.position = "bottom") +
        labs(title = sprintf("Comparacao de %s: espacial vs. nao-espacial", tipo),
             x = "Tempo", color = "Modelo", fill = "Modelo"),
      width = 14, height = 10, dpi = 300
    )
  }
}

cat("\nTempo total de execucao:\n")
print(Sys.time() - inicio_global)
