# ==============================================================================
# ajuste_lambda_ffbs_gamma_proporcional.R
#
# Modelo Poisson log-linear com:
#
#   mu[i,t] = lambda[t] * E[i,t] * epsilon[i,t] * exp(x[i,t,]' beta + s[i])
#
# ── lambda[t]: intercepto temporal DINÂMICO (Gama-DLM) ──────────────────────
#   Evolução multiplicativa: lambda[t] = w^{-1} * lambda[t-1] * xi[t]
#     xi[t] ~ Beta(w*a_{t-1}/lambda[t-1], ...)  (perturbação estocástica)
#   Na prática amostrado por FFBS Gama-conjugado (West & Harrison):
#     forward:  a[t+1] = w*a[t] + sum_i Y[i,t]
#               b[t+1] = w*b[t] + sum_i g[i,t]
#               g[i,t] = E[i,t] * epsilon[i,t] * exp(x'beta + s[i])
#     backward: lambda[T] ~ Gamma(a[T+1], b[T+1])
#               nu ~ Gamma((1-w)*a[t], b[t])
#               lambda[t] = nu + w * lambda[t+1]
#   Prior inicial: lambda[t] ~ Gamma(a0, b0)  (nó declarativo p/ o grafo NIMBLE;
#   o amostrador padrão é removido e substituído pelo FFBS)
#
# ── epsilon[i,t]: falha de notificação proporcional ─────────────────────────
#   epsilon[i,t] = 1 - gamma_1[t] * inprod(h[i,], kappa)
#   gamma_1[t] ~ Unif(L[t]-delta_g, L[t]+delta_g)
#     L[t]: sequência decrescente fixa  (ex.: 0.12 → 0.02)
#   kappa[1] = 1 (melhor cluster, fixo)
#   kappa[k] = kappa[k-1] + incr[k-1],  incr[k] ~ Unif(0, M_incr)
#   dconstraint(kappa[K]*(L_max+delta_g) < 1)  garante epsilon > 0
#
# ── Demais parâmetros ────────────────────────────────────────────────────────
#   beta[j] ~ N(0,1)
#   sigma_s ~ Half-t(0,1,1)   [apenas espacial]
#   s ~ ICAR(tau_s)            [apenas espacial]
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
# 2. HIPERPARÂMETROS
# ==============================================================================

# ── Gama-DLM (FFBS) ──────────────────────────────────────────────────────────
# a0, b0: prior inicial de lambda (lambda[0] ~ Gamma(a0, b0), media = a0/b0)
# w:      fator de desconto (0 < w < 1). Quanto maior w, mais suave a trajetória.
#         w = 0.85 é um ponto de partida comum para séries epidemiológicas.
a0 <- 1.0   # shape inicial
b0 <- 1.0   # rate inicial   →  E[lambda_0] = 1
w  <- 0.85  # fator de desconto

# ── Processo gamma_1[t] ───────────────────────────────────────────────────────
gamma_1_max <- 0.12   # L[1]: centro do prior em t = 1
gamma_1_min <- 0.02   # L[T]: centro do prior em t = T
delta_g     <- 0.015  # meia-largura da uniforme; exige L[T] > delta_g

stopifnot("delta_g deve ser < gamma_1_min" = delta_g < gamma_1_min)

L_seq <- seq(gamma_1_max, gamma_1_min, length.out = n_times)
L_max <- L_seq[1]

# Limite por incremento: kappa[K] < 1/(L_max + delta_g)
kappa_max_hard <- 1 / (L_max + delta_g)
M_incr         <- ceiling(kappa_max_hard - 1)

cat(sprintf(
  "a0=%.1f | b0=%.1f | w=%.2f\nL[1]=%.3f | L[T]=%.3f | delta_g=%.3f | M_incr=%d\n",
  a0, b0, w, L_max, tail(L_seq, 1), delta_g, M_incr
))

# ==============================================================================
# 3. REGIÕES PARA MONITORAR
# ==============================================================================
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
  a0 = a0, b0 = b0, w = w,
  L       = L_seq,
  delta_g = delta_g,
  L_max   = L_max,
  M_incr  = M_incr
)

constants_spatial <- c(constants_base, list(
  adj = adj_vec, num = num_vec, weights = weights_vec, n_adj = n_adj_val
))
constants_nonspatial <- constants_base

data_nimble <- list(
  Y = Y_mat, E = E_norm, x = x,
  constraint_kappa = 1L     # ativa o dconstraint de viabilidade de epsilon
)

# ==============================================================================
# 5. INICIALIZAÇÕES
# ==============================================================================
set.seed(123)

# lambda: inicializado em 1 (interior da Gama prior com a0=b0=1)
lambda_init <- rep(1.0, n_times)

# gamma_1: nos centros L[t]
gamma_1_init1 <- L_seq
gamma_1_init2 <- pmax(L_seq - delta_g * 0.5, L_seq + rnorm(n_times, 0, delta_g / 4))
gamma_1_init2 <- pmin(pmax(gamma_1_init2, L_seq - delta_g + 1e-4),
                      L_seq + delta_g - 1e-4)

# incr: kappa = (1, 2, 3, 4) -> incr = (1, 1, 1)
incr_init1 <- rep(1.0, K - 1)
incr_init2 <- rep(0.5, K - 1)

inits_spatial_1 <- list(
  lambda  = lambda_init,
  beta    = rep(0, p),
  gamma_1 = gamma_1_init1, incr = incr_init1,
  sigma_s = 0.5, s = rep(0, N_regions)
)
inits_spatial_2 <- list(
  lambda  = lambda_init * runif(n_times, 0.8, 1.2),
  beta    = rnorm(p, 0, 0.3),
  gamma_1 = gamma_1_init2, incr = incr_init2,
  sigma_s = 1.0, s = rep(0, N_regions)
)
inits_nonspatial_1 <- list(
  lambda  = lambda_init,
  beta    = rep(0, p),
  gamma_1 = gamma_1_init1, incr = incr_init1
)
inits_nonspatial_2 <- list(
  lambda  = lambda_init * runif(n_times, 0.8, 1.2),
  beta    = rnorm(p, 0, 0.3),
  gamma_1 = gamma_1_init2, incr = incr_init2
)

inits_list_spatial    <- list(inits_spatial_1, inits_spatial_2)
inits_list_nonspatial <- list(inits_nonspatial_1, inits_nonspatial_2)

# ==============================================================================
# 6. FUNÇÃO WORKER
# ==============================================================================
run_model <- function(model_type, output_dir) {
  
  library(nimble); library(coda); library(dplyr)
  library(ggplot2); library(readr); library(stringr); library(tibble)
  
  # ── 6a. Código NIMBLE ──────────────────────────────────────────────────────
  #
  #  mu[i,t] = lambda[t] * E[i,t] * epsilon[i,t] * exp(x'beta + s[i])
  #
  #  lambda[t] ~ Gamma(a0, b0)  ← nó declarativo para o grafo;
  #    o amostrador padrão é removido e substituído pelo FFBS abaixo.
  #    O FFBS Gama-conjugado garante que lambda[t] > 0 e captura a dinâmica
  #    multiplicativa w^{-1} * lambda[t-1] * xi[t], xi[t] ~ Beta.
  #
  #  epsilon[i,t] = 1 - gamma_1[t] * inprod(h[i,], kappa)
  #    gamma_1[t] ~ Unif(L[t]-delta_g, L[t]+delta_g)   (decrescente em t)
  #    kappa[1] = 1;  kappa[k] = kappa[k-1] + incr[k-1]
  #    dconstraint: kappa[K]*(L_max+delta_g) < 1  =>  epsilon > 0
  # ---------------------------------------------------------------------------
  
  code_spatial <- nimbleCode({
    
    # ── lambda[t]: nó declarativo (prior Gama para o grafo NIMBLE)
    # O amostrador FFBS substitui o RW padrão e amostra diretamente
    # toda a trajetória {lambda[1],...,lambda[T]} em um único passo.
    for (t in 1:n_times) {
      lambda[t] ~ dgamma(shape = a0, rate = b0)
    }
    
    # ── Covariáveis fixas
    for (j in 1:p) { beta[j] ~ dnorm(0, sd = 1) }
    
    # ── Processo gamma_1[t]: prior uniforme centrado em L[t] (decrescente)
    for (t in 1:n_times) {
      gamma_1[t] ~ dunif(L[t] - delta_g, L[t] + delta_g)
    }
    
    # ── Multiplicadores kappa via incrementos (ordem garantida por construção)
    for (k in 1:(K - 1)) { incr[k] ~ dunif(0, M_incr) }
    kappa[1] <- 1
    for (k in 2:K) { kappa[k] <- kappa[k - 1] + incr[k - 1] }
    
    # ── Restrição de viabilidade: epsilon[i,t] > 0 para todo (i,t).
    # Pior caso: kappa[K] * max(gamma_1[t]) = kappa[K] * (L_max + delta_g).
    # dconstraint = 0 (log = -Inf) se violada → proposta rejeitada no MH.
    constraint_kappa ~ dconstraint(kappa[K] * (L_max + delta_g) < 1)
    
    # ── Efeito espacial ICAR
    sigma_s ~ T(dt(0, 1, 1), 0, )
    tau_s   <- 1 / (sigma_s ^ 2)
    s[1:n_regions] ~ dcar_normal(
      adj[1:n_adj], weights[1:n_adj], num[1:n_regions], tau_s, zero_mean = 1
    )
    
    # ── Verossimilhança e quantidade derivada epsilon
    for (i in 1:n_regions) {
      for (t in 1:n_times) {
        epsilon[i, t]  <- 1 - gamma_1[t] * inprod(h[i, 1:K], kappa[1:K])
        mu[i, t]       <- lambda[t] * E[i, t] * epsilon[i, t] *
          exp(inprod(beta[1:p], x[i, t, 1:p]) + s[i])
        Y[i, t]        ~ dpois(mu[i, t])
        logLik_Y[i, t] <- dpois(Y[i, t], mu[i, t], log = TRUE)
      }
    }
  })
  
  code_nonspatial <- nimbleCode({
    
    for (t in 1:n_times) {
      lambda[t] ~ dgamma(shape = a0, rate = b0)
    }
    
    for (j in 1:p) { beta[j] ~ dnorm(0, sd = 1) }
    
    for (t in 1:n_times) {
      gamma_1[t] ~ dunif(L[t] - delta_g, L[t] + delta_g)
    }
    
    for (k in 1:(K - 1)) { incr[k] ~ dunif(0, M_incr) }
    kappa[1] <- 1
    for (k in 2:K) { kappa[k] <- kappa[k - 1] + incr[k - 1] }
    
    constraint_kappa ~ dconstraint(kappa[K] * (L_max + delta_g) < 1)
    
    for (i in 1:n_regions) {
      for (t in 1:n_times) {
        epsilon[i, t]  <- 1 - gamma_1[t] * inprod(h[i, 1:K], kappa[1:K])
        mu[i, t]       <- lambda[t] * E[i, t] * epsilon[i, t] *
          exp(inprod(beta[1:p], x[i, t, 1:p]))
        Y[i, t]        ~ dpois(mu[i, t])
        logLik_Y[i, t] <- dpois(Y[i, t], mu[i, t], log = TRUE)
      }
    }
  })
  
  # ── 6b. FFBS espacial ──────────────────────────────────────────────────────
  #
  # Amostrador Gama-conjugado para {lambda[1],...,lambda[T]}.
  #
  # Filtro forward (acumulação de estatísticas suficientes):
  #   a[t+1] = w * a[t] + sum_i Y[i,t]          (soma de contagens)
  #   b[t+1] = w * b[t] + sum_i g[i,t]          (soma de exposições ajustadas)
  #   g[i,t] = E[i,t] * epsilon[i,t] * exp(x'beta + s[i])
  #
  # Amostrador backward:
  #   lambda[T] ~ Gamma(a[T+1], b[T+1])
  #   nu        ~ Gamma((1-w)*a[t], b[t])
  #   lambda[t] = nu + w * lambda[t+1]     (garante lambda[t] > 0)
  # ---------------------------------------------------------------------------
  
  ffbs_spatial <- nimbleFunction(
    contains = sampler_BASE,
    setup = function(model, mvSaved, target, control) {
      n_regions <- control$n_regions; n_times <- control$n_times
      p  <- control$p
      a0 <- control$a0; b0 <- control$b0; w <- control$w
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
      
      # Forward: acumula contagens e exposições ajustadas por epsilon e beta
      # CORREÇÃO 2: clip de epsilon para evitar g_it <= 0 em draws transitórios
      # antes do dconstraint rejeitar a proposta. bt_buf <= 0 quebra rgamma.
      for (t in 1:n_times) {
        sum_Y_t <- 0; sum_g_t <- 0
        for (i in 1:n_regions) {
          sum_Y_t <- sum_Y_t + model$Y[i, t]
          prod_val <- 0
          for (k in 1:p) prod_val <- prod_val + model$x[i, t, k] * model$beta[k]
          g_it <- model$E[i, t] *
            max(model$epsilon[i, t], 1e-10) *
            exp(prod_val + model$s[i])
          sum_g_t <- sum_g_t + g_it
        }
        at_buf[t + 1] <<- w * at_buf[t] + sum_Y_t
        bt_buf[t + 1] <<- w * bt_buf[t] + sum_g_t
      }
      
      # Backward: amostra a trajetória completa de lambda (toda de uma vez)
      # CORREÇÃO 3: bt_buf deve ser > 0 para rgamma; garante numericamente.
      model$lambda[n_times] <<- rgamma(1, shape = at_buf[n_times + 1],
                                       rate  = max(bt_buf[n_times + 1], 1e-10))
      for (t_idx in 1:(n_times - 1)) {
        t_back <- n_times - t_idx
        # nu captura a inovação local; w*lambda[t+1] carrega a memória do futuro
        nu <- rgamma(1, shape = (1 - w) * at_buf[t_back + 1],
                     rate  = max(bt_buf[t_back + 1], 1e-10))
        model$lambda[t_back] <<- nu + w * model$lambda[t_back + 1]
      }
      
      model$calculate(calcNodes)
      copy(from = model, to = mvSaved, row = 1, nodes = targetNodes, logProb = TRUE)
    },
    methods = list(reset = function() {})
  )
  
  # ── 6c. FFBS não-espacial ──────────────────────────────────────────────────
  # Idêntico ao espacial, exceto que g[i,t] não inclui s[i].
  
  ffbs_nonspatial <- nimbleFunction(
    contains = sampler_BASE,
    setup = function(model, mvSaved, target, control) {
      n_regions <- control$n_regions; n_times <- control$n_times
      p  <- control$p
      a0 <- control$a0; b0 <- control$b0; w <- control$w
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
          sum_Y_t <- sum_Y_t + model$Y[i, t]
          prod_val <- 0
          for (k in 1:p) prod_val <- prod_val + model$x[i, t, k] * model$beta[k]
          # CORREÇÃO 2: clip de epsilon (mesma razão do FFBS espacial)
          g_it    <- model$E[i, t] * max(model$epsilon[i, t], 1e-10) * exp(prod_val)
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
  
  # ── 6d. Selecionar objetos por modelo ──────────────────────────────────────
  is_spatial <- (model_type == "spatial")
  model_code <- if (is_spatial) code_spatial      else code_nonspatial
  constants  <- if (is_spatial) constants_spatial  else constants_nonspatial
  inits_list <- if (is_spatial) inits_list_spatial else inits_list_nonspatial
  ffbs_fn    <- if (is_spatial) ffbs_spatial       else ffbs_nonspatial
  
  cat("\n=== Iniciando modelo:", model_type, "===\n")
  scenario_dir <- file.path(output_dir, model_type)
  dir.create(scenario_dir, recursive = TRUE, showWarnings = FALSE)
  
  # ── 6e. Construir e compilar ───────────────────────────────────────────────
  model  <- nimbleModel(code = model_code, constants = constants,
                        data = data_nimble, inits = inits_list[[1]], check = FALSE)
  
  # CORREÇÃO 1: forçar o cálculo de todos os nós determinísticos
  # (epsilon, kappa, mu) antes da compilação.
  # Sem isso, model$epsilon[i,t] fica NA na primeira chamada do FFBS,
  # propagando NA para bt_buf e quebrando rgamma com "TRUE/FALSE necessário".
  model$calculate()
  
  Cmodel <- compileNimble(model)
  
  conf <- configureMCMC(model)
  
  # lambda: remove o amostrador padrão (RW) e insere o FFBS para toda a trajetória
  conf$removeSamplers("lambda")
  conf$addSampler(
    target  = "lambda",
    type    = ffbs_fn,
    control = list(
      n_regions = N_regions, n_times = n_times, p = p,
      a0 = constants$a0, b0 = constants$b0, w = constants$w
    )
  )
  
  # incr: AF_slice em bloco — propõe todos os incrementos juntos,
  # respeitando a restrição conjunta via dconstraint sem alta rejeição
  conf$removeSamplers("incr")
  conf$addSampler(
    target = paste0("incr[", seq_len(K - 1), "]"),
    type   = "AF_slice"
  )
  
  # gamma_1[t]: RW individual por t (suporte restrito pelo prior uniforme)
  conf$removeSamplers("gamma_1")
  for (t_idx in seq_len(n_times)) {
    conf$addSampler(target = paste0("gamma_1[", t_idx, "]"), type = "RW")
  }
  
  monitors_base <- c("lambda", "beta", "gamma_1", "kappa", "incr", "logLik_Y")
  if (is_spatial) monitors_base <- c(monitors_base, "s", "sigma_s", "tau_s")
  conf$addMonitors(monitors_base)
  conf$addMonitors(MU_MONITORS)
  
  Rmcmc <- buildMCMC(conf)
  Cmcmc <- compileNimble(Rmcmc, project = model)
  
  # ── 6f. MCMC ───────────────────────────────────────────────────────────────
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
  
  # ── 6g. Auxiliares ─────────────────────────────────────────────────────────
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
  
  lambda_names <- paste0("lambda[", 1:n_times, "]")
  beta_names   <- paste0("beta[",   1:p,        "]")
  gamma1_names <- paste0("gamma_1[", 1:n_times, "]")
  kappa_names  <- paste0("kappa[",   1:K,        "]")
  incr_names   <- paste0("incr[",    1:(K - 1),  "]")
  
  # ── 6h. Epsilon posterior (N × T) ──────────────────────────────────────────
  # epsilon[i,t] = 1 - gamma_1[t] * inprod(h[i,], kappa)
  n_draw      <- nrow(samples_mat)
  kappa_draws <- samples_mat[, kappa_names, drop = FALSE]  # n_draw × K
  h_kappa     <- h_mat %*% t(kappa_draws)                  # N × n_draw
  
  epsilon_draws <- array(NA_real_, dim = c(n_draw, N_regions, n_times))
  for (t in seq_len(n_times)) {
    g1_t <- samples_mat[, gamma1_names[t]]                         # n_draw
    epsilon_draws[, , t] <- t(1 - sweep(h_kappa, 2, g1_t, `*`))   # n_draw × N
  }
  
  # ── 6i. Sumário de lambda[t] ───────────────────────────────────────────────
  anos_label <- colnames(data_nimble$Y)
  if (is.null(anos_label)) anos_label <- seq_len(n_times)
  
  lambda_summary <- do.call(rbind, lapply(1:n_times, function(t) {
    nm  <- lambda_names[t]; sv <- samples_mat[, nm]; hpd <- safe_hpd(sv)
    tibble(Time = t, Parameter = nm,
           Mean = mean(sv), SD = sd(sv), Lower = hpd[1], Upper = hpd[2],
           ESS  = as.numeric(effectiveSize(mcmc_list_full[, nm])),
           Rhat = safe_gelman(mcmc_list_full[, nm]))
  }))
  write_csv(lambda_summary, file.path(scenario_dir, "lambda_summary.csv"))
  
  ggsave(
    file.path(scenario_dir, "painel_lambda.png"),
    ggplot(lambda_summary, aes(x = Time)) +
      geom_ribbon(aes(ymin = Lower, ymax = Upper), fill = "grey70", alpha = 0.5) +
      geom_line(aes(y = Mean), color = "black", linewidth = 0.9) +
      geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
      scale_x_continuous(breaks = 1:n_times, labels = anos_label) +
      theme_bw(base_size = 12) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      labs(title    = paste("Intensidade dinâmica lambda[t] (", model_type, ")"),
           subtitle = paste0("Gama-DLM  |  w = ", w,
                             "  |  E[lambda_0] = ", round(a0 / b0, 2)),
           x = "Tempo", y = expression(lambda[t])),
    width = 9, height = 5
  )
  
  # ── 6j. Sumário de gamma_1[t] ──────────────────────────────────────────────
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
      geom_line(aes(y = L_center), color = "black",       linewidth = 0.7,
                linetype = "dashed") +
      scale_x_continuous(breaks = 1:n_times, labels = anos_label) +
      theme_bw(base_size = 12) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      labs(title    = paste("Processo gamma_1[t] — melhor cluster (", model_type, ")"),
           subtitle = "Tracejado = centro do prior L[t]",
           x = "Tempo", y = expression(gamma[1][t])),
    width = 9, height = 5
  )
  
  # ── 6k. Sumário de kappa[k] ────────────────────────────────────────────────
  # kappa[1] = 1 é determinístico: ESS e Rhat não se aplicam.
  # kappa[k≥2] são derivados de incr[], portanto têm variância positiva.
  kappa_summary <- do.call(rbind, lapply(1:K, function(k) {
    nm  <- kappa_names[k]
    sv  <- samples_mat[, nm]
    hpd <- safe_hpd(sv)
    is_const <- (var(sv) < 1e-12)   # kappa[1] = 1 fixo
    tibble(k = k, Parameter = nm,
           Mean = mean(sv), SD = sd(sv), Lower = hpd[1], Upper = hpd[2],
           ESS  = if (is_const) NA_real_ else
             as.numeric(effectiveSize(mcmc_list_full[, nm])),
           Rhat = if (is_const) NA_real_ else
             safe_gelman(mcmc_list_full[, nm]))
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
           subtitle = "kappa[1] = 1 fixo | ordem garantida por incrementos",
           x = "Cluster k", y = expression(kappa[k])),
    width = 6, height = 5
  )
  
  # ── 6l. Sumário de epsilon ─────────────────────────────────────────────────
  # ── 6l. Sumário de epsilon por cluster ─────────────────────────────────────
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
      scale_x_continuous(breaks = 1:n_times, labels = anos_label) +
      theme_bw(base_size = 11) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      labs(title    = paste("Trajetória de epsilon por cluster (", model_type, ")"),
           x = "Tempo", y = expression(epsilon[kt])),
    width = 10, height = ceiling(K / 2) * 4
  )
  
  # ── 6m. Sumário de beta ────────────────────────────────────────────────────
  beta_summary <- do.call(rbind, lapply(beta_names, function(nm) {
    sv  <- samples_mat[, nm]; hpd <- safe_hpd(sv)
    tibble(Parameter = nm, Mean = mean(sv), SD = sd(sv),
           Lower = hpd[1], Upper = hpd[2],
           ESS  = as.numeric(effectiveSize(mcmc_list_full[, nm])),
           Rhat = safe_gelman(mcmc_list_full[, nm]))
  }))
  write_csv(beta_summary, file.path(scenario_dir, "beta_summary.csv"))
  
  # ── 6n. tau_s e s (espacial) ───────────────────────────────────────────────
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
  
  # ── 6o. mu para regiões selecionadas ───────────────────────────────────────
  # ── 6o. mu para regiões selecionadas ───────────────────────────────────────
  mu_summary <- do.call(rbind, lapply(REGIONS_INTEREST, function(i) {
    do.call(rbind, lapply(1:n_times, function(t) {
      nm  <- paste0("mu[", i, ", ", t, "]")
      sv  <- samples_mat[, nm]; hpd <- safe_hpd(sv)
      tibble(Region = i, Time = t, Mean = mean(sv),
             Lower = hpd[1], Upper = hpd[2], model = model_type)
    }))
  }))
  write_csv(mu_summary, file.path(scenario_dir, "mu_selected.csv"))
  
  ggsave(
    file.path(scenario_dir, "painel_mu.png"),
    ggplot(mu_summary, aes(x = Time)) +
      geom_ribbon(aes(ymin = Lower, ymax = Upper), fill = "darkorange", alpha = 0.3) +
      geom_line(aes(y = Mean), color = "darkorange", linewidth = 0.8) +
      facet_wrap(~ Region, ncol = 3, scales = "free_y",
                 labeller = labeller(Region = function(x) paste("Região", x))) +
      scale_x_continuous(breaks = 1:n_times, labels = anos_label) +
      theme_bw(base_size = 11) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      labs(title = paste("Mu estimado (", model_type, ")"),
           x = "Tempo", y = expression(mu[it])),
    width = 12, height = 10
  )
  
  # ── 6p. WAIC e LPML ────────────────────────────────────────────────────────
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
  
  # ── 6q. Diagnósticos ───────────────────────────────────────────────────────
  cores_cadeia <- c("Cadeia 1" = "#2166AC", "Cadeia 2" = "#D6604D")
  
  # ACF de lambda[t]
  acf_lambda_df <- do.call(rbind, lapply(1:n_times, function(t) {
    nm <- lambda_names[t]
    ac <- acf(samples_mat[, nm], lag.max = 100, plot = FALSE)
    tibble(Time = t, Lag = as.vector(ac$lag[-1]), ACF = as.vector(ac$acf[-1]))
  }))
  ggsave(
    file.path(scenario_dir, "acf_lambda.png"),
    ggplot(acf_lambda_df, aes(x = Lag, y = ACF)) +
      geom_col(width = 0.6, fill = "darkorange") +
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
  
  # ACF de gamma_1[t]
  acf_g1_df <- do.call(rbind, lapply(1:n_times, function(t) {
    nm <- gamma1_names[t]
    ac <- acf(samples_mat[, nm], lag.max = 100, plot = FALSE)
    tibble(Time = t, Lag = as.vector(ac$lag[-1]), ACF = as.vector(ac$acf[-1]))
  }))
  ggsave(
    file.path(scenario_dir, "acf_gamma1.png"),
    ggplot(acf_g1_df, aes(x = Lag, y = ACF)) +
      geom_col(width = 0.6, fill = "darkorange") +
      geom_hline(yintercept = c(-0.10, 0.10), linetype = "dashed",
                 color = "blue", linewidth = 0.5) +
      geom_hline(yintercept = c(-0.05, 0.05), linetype = "dotted",
                 color = "red",  linewidth = 0.5) +
      facet_wrap(~ Time, scales = "free_y", ncol = 6,
                 labeller = label_bquote(gamma[1][.(Time)])) +
      theme_bw(base_size = 9) +
      labs(title    = paste("ACF de gamma_1[t] (", model_type, ")"),
           subtitle = "Azul: |0.10| | Vermelho: |0.05|"),
    width = 14, height = 10
  )
  
  # ACF de incr[k] (incrementos dos multiplicadores)
  acf_incr_df <- do.call(rbind, lapply(seq_len(K - 1), function(k) {
    nm <- incr_names[k]
    ac <- acf(samples_mat[, nm], lag.max = 100, plot = FALSE)
    tibble(Parameter = nm, Lag = as.vector(ac$lag[-1]), ACF = as.vector(ac$acf[-1]))
  }))
  ggsave(
    file.path(scenario_dir, "acf_incr_kappa.png"),
    ggplot(acf_incr_df, aes(x = Lag, y = ACF)) +
      geom_col(width = 0.6, fill = "darkorange") +
      geom_hline(yintercept = c(-0.10, 0.10), linetype = "dashed",
                 color = "blue", linewidth = 0.5) +
      geom_hline(yintercept = c(-0.05, 0.05), linetype = "dotted",
                 color = "red",  linewidth = 0.5) +
      facet_wrap(~ Parameter, scales = "free_y") +
      theme_bw(base_size = 11) +
      labs(title    = paste("ACF dos incrementos incr[k] (", model_type, ")"),
           subtitle = "Azul: |0.10| | Vermelho: |0.05|"),
    width = 8, height = 5
  )
  
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
  
  # Traceplots de lambda[t] (trajetória completa)
  df_trace_lambda <- do.call(rbind, lapply(1:nchains, function(ch) {
    cm <- as.matrix(mcmc_list_full[[ch]])
    do.call(rbind, lapply(lambda_names, function(nm) {
      vals <- cm[, nm]
      tibble(Iter = seq_along(vals), Value = vals,
             ErgMedia = cumsum(vals) / seq_along(vals),
             Parameter = nm, Cadeia = paste0("Cadeia ", ch))
    }))
  }))
  ggsave(
    file.path(scenario_dir, "traceplots_lambda.png"),
    ggplot(df_trace_lambda, aes(x = Iter, color = Cadeia)) +
      geom_line(aes(y = Value),    alpha = 0.25, linewidth = 0.20) +
      geom_line(aes(y = ErgMedia), alpha = 0.90, linewidth = 0.75) +
      scale_color_manual(values = cores_cadeia) +
      facet_wrap(~ Parameter, scales = "free_y", ncol = 6) +
      theme_bw(base_size = 9) + theme(legend.position = "bottom") +
      labs(title    = paste("Traceplots lambda[t] (", model_type, ")"),
           subtitle = "Linha grossa = média ergódica | Linha fina = cadeia",
           x = "Iteração (pós-burnin)", y = expression(lambda[t])),
    width = 14, height = 10
  )
  
  # Traceplots de kappa[k]
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
           x = "Iteração (pós-burnin)", y = expression(kappa[k])),
    width = 10, height = 6
  )
  
  # ── 6r. Retorno resumido ───────────────────────────────────────────────────
  # ESS de kappa via incr (nós estocásticos reais; kappa[1]=1 é constante)
  incr_summary_ess <- do.call(rbind, lapply(incr_names, function(nm) {
    sv <- samples_mat[, nm]
    tibble(Parameter = nm,
           ESS  = as.numeric(effectiveSize(mcmc_list_full[, nm])),
           Rhat = safe_gelman(mcmc_list_full[, nm]))
  }))
  
  tibble(
    model = model_type, niter = niter, nburnin = nburnin, thin = thin,
    WAIC = waic, LPML = LPML,
    ESS_beta_min    = min(beta_summary$ESS,       na.rm = TRUE),
    ESS_lambda_min  = min(lambda_summary$ESS,     na.rm = TRUE),
    ESS_gamma1_min  = min(gamma1_summary$ESS,     na.rm = TRUE),
    ESS_incr_min    = min(incr_summary_ess$ESS,   na.rm = TRUE),
    ESS_tau         = ESS_tau,
    Rhat_beta_max   = max(beta_summary$Rhat,      na.rm = TRUE),
    Rhat_lambda_max = max(lambda_summary$Rhat,    na.rm = TRUE),
    Rhat_gamma1_max = max(gamma1_summary$Rhat,    na.rm = TRUE),
    Rhat_incr_max   = max(incr_summary_ess$Rhat,  na.rm = TRUE)
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
  "resultados_lambda_ffbs_gamma_proporcional"
)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

cl <- makeCluster(n_cores)
clusterExport(cl, c(
  "constants_spatial", "constants_nonspatial", "data_nimble",
  "inits_list_spatial", "inits_list_nonspatial",
  "N_regions", "n_times", "p", "K", "h_mat", "cluster_ids",
  "run_model", "output_dir", "REGIONS_INTEREST", "MU_MONITORS",
  "L_seq", "L_max", "delta_g", "M_incr", "a0", "b0", "w"
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

read_model_csv <- function(filename) {
  lapply(model_types, function(m) {
    path <- file.path(output_dir, m, filename)
    if (!file.exists(path)) return(NULL)
    read_csv(path, show_col_types = FALSE) %>% mutate(model = m)
  }) |> bind_rows()
}

anos_label <- colnames(data_nimble$Y)
if (is.null(anos_label)) anos_label <- seq_len(n_times)

# Comparativo lambda[t]
lambda_all <- read_model_csv("lambda_summary.csv")
if (nrow(lambda_all) > 0) {
  ggsave(
    file.path(output_dir, "lambda_comparativo.png"),
    ggplot(lambda_all, aes(x = Time, y = Mean, color = model, fill = model)) +
      geom_ribbon(aes(ymin = Lower, ymax = Upper), alpha = 0.15, color = NA) +
      geom_line(linewidth = 0.9) +
      geom_hline(yintercept = 1, linetype = "dashed", color = "grey40") +
      scale_x_continuous(breaks = 1:n_times, labels = anos_label) +
      theme_bw(base_size = 12) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            legend.position = "bottom") +
      labs(title    = "Comparação lambda[t]: espacial vs. não-espacial",
           subtitle = paste0("Gama-DLM  |  w = ", w),
           x = "Tempo", y = expression(lambda[t]),
           color = "Modelo", fill = "Modelo"),
    width = 9, height = 5, dpi = 300
  )
}

# Comparativo gamma_1[t]
gamma1_all <- read_model_csv("gamma1_summary.csv")
if (nrow(gamma1_all) > 0) {
  ggsave(
    file.path(output_dir, "gamma1_comparativo.png"),
    ggplot(gamma1_all, aes(x = Time, y = Mean, color = model, fill = model)) +
      geom_ribbon(aes(ymin = Lower, ymax = Upper), alpha = 0.15, color = NA) +
      geom_line(linewidth = 0.9) +
      geom_line(data = gamma1_all %>% filter(model == model_types[1]),
                aes(x = Time, y = L_center), color = "black", linetype = "dashed",
                linewidth = 0.6, inherit.aes = FALSE) +
      scale_x_continuous(breaks = 1:n_times, labels = anos_label) +
      theme_bw(base_size = 12) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            legend.position = "bottom") +
      labs(title    = "Comparação gamma_1[t]: espacial vs. não-espacial",
           subtitle = "Tracejado = centro do prior L[t]",
           x = "Tempo", y = expression(gamma[1][t]),
           color = "Modelo", fill = "Modelo"),
    width = 9, height = 5, dpi = 300
  )
}

# Comparativo kappa[k]
kappa_all <- read_model_csv("kappa_summary.csv")
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
           x = "Cluster k", y = expression(kappa[k]), color = "Modelo"),
    width = 7, height = 5, dpi = 300
  )
}

# Comparativo mu
mu_all <- read_model_csv("mu_selected.csv")
if (nrow(mu_all) > 0) {
  anos_label_cmp <- colnames(data_nimble$Y)
  if (is.null(anos_label_cmp)) anos_label_cmp <- seq_len(n_times)
  
  ggsave(
    file.path(output_dir, "mu_comparativo.png"),
    ggplot(mu_all, aes(x = Time, y = Mean, color = model, fill = model)) +
      geom_ribbon(aes(ymin = Lower, ymax = Upper), alpha = 0.15, color = NA) +
      geom_line(linewidth = 0.8) +
      facet_wrap(~ Region, scales = "free_y", ncol = 3,
                 labeller = labeller(Region = function(x) paste("Região", x))) +
      scale_x_continuous(breaks = 1:n_times, labels = anos_label_cmp) +
      theme_bw(base_size = 12) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            legend.position = "bottom") +
      labs(title = "Comparação mu[i,t]: espacial vs. não-espacial",
           x = "Tempo", y = expression(mu[it]),
           color = "Modelo", fill = "Modelo"),
    width = 14, height = 10, dpi = 300
  )
}

# Comparativo epsilon por cluster
eps_cluster_all <- lapply(model_types, function(m) {
  path <- file.path(output_dir, m, "epsilon_cluster_summary.csv")
  if (!file.exists(path)) return(NULL)
  read_csv(path, show_col_types = FALSE) %>% mutate(model = m)
}) |> bind_rows()

if (nrow(eps_cluster_all) > 0) {
  anos_label_cmp <- colnames(data_nimble$Y)
  if (is.null(anos_label_cmp)) anos_label_cmp <- seq_len(n_times)
  
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
      scale_x_continuous(breaks = 1:n_times, labels = anos_label_cmp) +
      theme_bw(base_size = 12) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            legend.position = "bottom") +
      labs(title    = "Comparação epsilon por cluster: espacial vs. não-espacial",
           subtitle = "Banda = IC 95% HPD",
           x = "Tempo", y = expression(epsilon[kt]),
           color = "Modelo", fill = "Modelo"),
    width = 10, height = ceiling(K / 2) * 4, dpi = 300
  )
}

cat("\nTempo total de execução:\n")
print(Sys.time() - inicio_global)
