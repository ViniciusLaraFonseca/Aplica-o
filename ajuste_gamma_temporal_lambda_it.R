# ==============================================================================
# ajuste_gamma_temporal_lambda_it.R
# Modelo Bayesiano Espaço-Estados — dados reais MG
# A = 75 microrregiões | T = 23 anos | K = 4 clusters | p = 3 covariáveis
#
# EXTENSÃO em relação ao original:
#   gamma: VETOR K  →  MATRIZ K × T   (gamma[k, t])
#   epsilon: VETOR N  →  MATRIZ N × T  (epsilon[i, t])
#   lambda: VETOR T → MATRIZ N × T (lambda[i, t])
#   epsilon[i, t] = 1 - sum_k h[i, k] * gamma[k, t]
#
#   Preditor:
#     log(mu[i,t]) = log(lambda[i,t]) + log(E[i,t]) + log(epsilon[i,t])
#                    + beta' x[i,t] + s[i]   (s[i] apenas no modelo espacial)
#
#   Priori de gamma por coluna-t (restrição de soma independente em cada t):
#     gamma[1, t] ~ dunif(a_unif, b_unif)
#     gamma[j, t] ~ dunif(0, 1 - sum(gamma[1:(j-1), t]))   j = 2..K
#
#   Amostrador: AF_slice separado por coluna-t de gamma (23 blocos de tamanho K).
#   lambda[i,t]: FFBS por região
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

# Nomes das regiões e anos
region_names <- rownames(Y_mat)
anos_label   <- colnames(Y_mat)
if (is.null(anos_label)) anos_label <- as.character(seq_len(n_times))

# ==============================================================================
# 2. REGIÕES DE INTERESSE (incluindo BH e Araçuaí)
# ==============================================================================
cluster_ids <- apply(h_mat, 1, sum)

# Encontrar índices de BH e Araçuaí
bh_idx      <- grep("BELO_HORIZONTE", region_names, ignore.case = TRUE)
aracuai_idx <- grep("ARACUAI", region_names, ignore.case = TRUE)

cat("BH encontrado na posição:", bh_idx, "-", region_names[bh_idx], "\n")
cat("Araçuaí encontrado na posição:", aracuai_idx, "-", region_names[aracuai_idx], "\n")

# Selecionar regiões de interesse: BH, Araçuaí + 3 aleatórias por cluster
set.seed(42)
REGIONS_INTEREST <- c(bh_idx, aracuai_idx)
REGIONS_INTEREST <- c(REGIONS_INTEREST,
                      unlist(lapply(1:K, function(cl) {
                        regs <- which(cluster_ids == cl)
                        regs <- setdiff(regs, REGIONS_INTEREST)  # remove as já incluídas
                        if (length(regs) >= 3) sample(regs, 3) else regs
                      }))
)
REGIONS_INTEREST <- sort(unique(REGIONS_INTEREST))
cat("Regiões selecionadas:\n")
print(data.frame(Indice = REGIONS_INTEREST, Nome = region_names[REGIONS_INTEREST]))

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

lambda_init1 <- matrix(1.0, nrow = N_regions, ncol = n_times)
lambda_init2 <- matrix(rgamma(N_regions * n_times, 1, 1), 
                       nrow = N_regions, ncol = n_times)

inits_list_spatial <- list(
  list(beta = rep(0, p),        gamma = gamma_init1,
       lambda = lambda_init1, sigma_s = 0.5, s = rep(0, N_regions)),
  list(beta = rnorm(p, 0, 0.3), gamma = gamma_init2,
       lambda = lambda_init2, sigma_s = 1.0, s = rep(0, N_regions))
)
inits_list_nonspatial <- list(
  list(beta = rep(0, p),        gamma = gamma_init1, lambda = lambda_init1),
  list(beta = rnorm(p, 0, 0.3), gamma = gamma_init2, lambda = lambda_init2)
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
    
    # Prior para lambda[i,t]
    for (i in 1:n_regions) {
      for (t in 1:n_times) {
        lambda[i, t] ~ dgamma(a0, b0)
      }
    }
    
    for (i in 1:n_regions) {
      for (t in 1:n_times) {
        epsilon[i, t] <- 1 - inprod(h[i, 1:K], gamma[1:K, t])
        log(mu[i, t]) <- log(lambda[i, t]) + log(E[i, t]) + log(epsilon[i, t]) +
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
    
    for (i in 1:n_regions) {
      for (t in 1:n_times) {
        lambda[i, t] ~ dgamma(a0, b0)
      }
    }
    
    for (i in 1:n_regions) {
      for (t in 1:n_times) {
        epsilon[i, t] <- 1 - inprod(h[i, 1:K], gamma[1:K, t])
        log(mu[i, t]) <- log(lambda[i, t]) + log(E[i, t]) + log(epsilon[i, t]) +
          inprod(beta[1:p], x[i, t, 1:p])
        Y[i, t]        ~ dpois(mu[i, t])
        logLik_Y[i, t] <- dpois(Y[i, t], mu[i, t], log = TRUE)
      }
    }
  })
  
  # ── 5b. FFBS para lambda[i,t] por região ────────────────────────────────────
  ffbs_spatial <- nimbleFunction(
    contains = sampler_BASE,
    setup = function(model, mvSaved, target, control) {
      n_times   <- control$n_times
      p         <- control$p
      a0        <- control$a0
      b0        <- control$b0
      w         <- control$w
      region_i  <- control$region_i
      
      at_buf <- nimNumeric(n_times + 1, 0)
      bt_buf <- nimNumeric(n_times + 1, 0)
      
      calcNodes   <- model$getDependencies(target, self = FALSE)
      targetNodes <- model$expandNodeNames(target)
      setupOutputs(at_buf, bt_buf)
    },
    run = function() {
      declare(t, integer()); declare(t_idx, integer())
      declare(t_back, integer()); declare(k, integer())
      declare(prod_val, double()); declare(g_it, double())
      declare(nu, double())
      
      at_buf[1] <<- a0; bt_buf[1] <<- b0
      
      for (t in 1:n_times) {
        prod_val <- 0.0
        for (k in 1:p) prod_val <- prod_val + model$x[region_i, t, k] * model$beta[k]
        g_it <- model$E[region_i, t] * model$epsilon[region_i, t] * 
          exp(prod_val + model$s[region_i])
        
        at_buf[t + 1] <<- w * at_buf[t] + model$Y[region_i, t]
        bt_buf[t + 1] <<- w * bt_buf[t] + g_it
      }
      
      model$lambda[region_i, n_times] <<- rgamma(1, shape = at_buf[n_times + 1],
                                                 rate  = bt_buf[n_times + 1])
      for (t_idx in 1:(n_times - 1)) {
        t_back <- n_times - t_idx
        nu <- rgamma(1, shape = (1 - w) * at_buf[t_back + 1],
                     rate  = bt_buf[t_back + 1])
        model$lambda[region_i, t_back] <<- nu + w * model$lambda[region_i, t_back + 1]
      }
      
      model$calculate(calcNodes)
      copy(from = model, to = mvSaved, row = 1, nodes = targetNodes, logProb = TRUE)
    },
    methods = list(reset = function() {})
  )
  
  ffbs_nonspatial <- nimbleFunction(
    contains = sampler_BASE,
    setup = function(model, mvSaved, target, control) {
      n_times   <- control$n_times
      p         <- control$p
      a0        <- control$a0
      b0        <- control$b0
      w         <- control$w
      region_i  <- control$region_i
      
      at_buf <- nimNumeric(n_times + 1, 0)
      bt_buf <- nimNumeric(n_times + 1, 0)
      
      calcNodes   <- model$getDependencies(target, self = FALSE)
      targetNodes <- model$expandNodeNames(target)
      setupOutputs(at_buf, bt_buf)
    },
    run = function() {
      declare(t, integer()); declare(t_idx, integer())
      declare(t_back, integer()); declare(k, integer())
      declare(prod_val, double()); declare(g_it, double())
      declare(nu, double())
      
      at_buf[1] <<- a0; bt_buf[1] <<- b0
      
      for (t in 1:n_times) {
        prod_val <- 0.0
        for (k in 1:p) prod_val <- prod_val + model$x[region_i, t, k] * model$beta[k]
        g_it <- model$E[region_i, t] * model$epsilon[region_i, t] * exp(prod_val)
        
        at_buf[t + 1] <<- w * at_buf[t] + model$Y[region_i, t]
        bt_buf[t + 1] <<- w * bt_buf[t] + g_it
      }
      
      model$lambda[region_i, n_times] <<- rgamma(1, shape = at_buf[n_times + 1],
                                                 rate  = bt_buf[n_times + 1])
      for (t_idx in 1:(n_times - 1)) {
        t_back <- n_times - t_idx
        nu <- rgamma(1, shape = (1 - w) * at_buf[t_back + 1],
                     rate  = bt_buf[t_back + 1])
        model$lambda[region_i, t_back] <<- nu + w * model$lambda[region_i, t_back + 1]
      }
      
      model$calculate(calcNodes)
      copy(from = model, to = mvSaved, row = 1, nodes = targetNodes, logProb = TRUE)
    },
    methods = list(reset = function() {})
  )
  
  # ── 5c. Selecionar objetos ───────────────────────────────────────────────────
  is_spatial <- (model_type == "spatial")
  model_code <- if (is_spatial) code_spatial      else code_nonspatial
  constants  <- if (is_spatial) constants_spatial  else constants_nonspatial
  inits_list <- if (is_spatial) inits_list_spatial else inits_list_nonspatial
  ffbs_fn    <- if (is_spatial) ffbs_spatial       else ffbs_nonspatial
  
  cat("\n=== Iniciando modelo:", model_type, "===\n")
  scenario_dir <- file.path(output_dir, model_type)
  dir.create(scenario_dir, recursive = TRUE, showWarnings = FALSE)
  
  # ── 5d. Construir e compilar ─────────────────────────────────────────────────
  model  <- nimbleModel(code = model_code, constants = constants,
                        data = data_nimble, inits = inits_list[[1]], check = FALSE)
  Cmodel <- compileNimble(model)
  
  conf <- configureMCMC(model)
  
  # FFBS para lambda[i,t] por região
  conf$removeSamplers("lambda")
  for (i in seq_len(N_regions)) {
    conf$addSampler(
      target = paste0("lambda[", i, ", 1:", n_times, "]"),
      type   = ffbs_fn,
      control = list(
        n_times   = n_times,
        p         = p,
        a0        = constants$a0,
        b0        = constants$b0,
        w         = constants$w,
        region_i  = i
      )
    )
  }
  
  # AF_slice por coluna-t de gamma
  conf$removeSamplers("gamma")
  for (t_idx in seq_len(n_times)) {
    conf$addSampler(
      target = paste0("gamma[", seq_len(K), ", ", t_idx, "]"),
      type   = "AF_slice"
    )
  }
  
  monitors_base <- c("beta", "gamma", "logLik_Y", "lambda")
  if (is_spatial) monitors_base <- c(monitors_base, "s", "sigma_s", "tau_s")
  
  # Adicionar mu para regiões selecionadas
  mu_monitors <- unlist(lapply(REGIONS_INTEREST, 
                               function(r) paste0("mu[", r, ", ", seq_len(n_times), "]")))
  
  conf$addMonitors(monitors_base)
  conf$addMonitors(mu_monitors)
  
  conf$printSamplers()
  
  Rmcmc <- buildMCMC(conf)
  Cmcmc <- compileNimble(Rmcmc, project = model)
  
  # ── 5e. MCMC ─────────────────────────────────────────────────────────────────
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
  
  # ── 5f. Funções auxiliares ───────────────────────────────────────────────────
  safe_hpd <- function(sv) {
    if (var(sv) < 1e-12) return(c(NA_real_, NA_real_))
    as.numeric(HPDinterval(as.mcmc(sv), prob = 0.95))
  }
  safe_gelman <- function(obj) {
    tryCatch(gelman.diag(obj, autoburnin = FALSE)$psrf[, 1],
             error = function(e) rep(NA_real_, nvar(obj)))
  }
  
  beta_names   <- paste0("beta[",   seq_len(p),       "]")
  
  # Nomes "gamma[k, t]" na matriz de amostras
  gamma_names_mat <- outer(
    seq_len(K), seq_len(n_times),
    function(k, t) paste0("gamma[", k, ", ", t, "]")
  )
  
  # Nomes "lambda[i, t]" na matriz de amostras
  lambda_names_mat <- outer(
    seq_len(N_regions), seq_len(n_times),
    function(i, t) paste0("lambda[", i, ", ", t, "]")
  )
  
  # ── 5g. Epsilon posterior — array (n_draw × N × T) ──────────────────────────
  n_draw        <- nrow(samples_mat)
  epsilon_draws <- array(NA_real_, dim = c(n_draw, N_regions, n_times))
  for (t in seq_len(n_times)) {
    g_t <- samples_mat[, gamma_names_mat[, t], drop = FALSE]
    epsilon_draws[, , t] <- 1 - g_t %*% t(h_mat)
  }
  
  # ── 5h. Sumário de gamma[k, t] ───────────────────────────────────────────────
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
  
  # Painel temporal de gamma[k, t]
  ggsave(
    file.path(scenario_dir, "painel_gamma.png"),
    ggplot(gamma_summary, aes(x = Time)) +
      geom_ribbon(aes(ymin = Lower, ymax = Upper), fill = "steelblue", alpha = 0.3) +
      geom_line(aes(y = Mean), color = "steelblue", linewidth = 0.8) +
      facet_wrap(
        ~ k, ncol = 2, scales = "free_y",
        labeller = labeller(k = function(x) paste0("Cluster ", x))
      ) +
      scale_x_continuous(breaks = 1:n_times, labels = anos_label) +
      theme_bw(base_size = 11) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      labs(title = paste("Trajetória de gamma[k, t] (", model_type, ")"),
           subtitle = "Banda = IC 95% HPD",
           x = "Ano", y = expression(gamma[kt])),
    width = 12, height = 8
  )
  
  # ── 5i. Sumário de epsilon[i, t] por região e por cluster ────────────────────
  
  # Sumário completo (N × T)
  eps_full <- do.call(rbind, lapply(seq_len(N_regions), function(i) {
    do.call(rbind, lapply(seq_len(n_times), function(t) {
      sv  <- epsilon_draws[, i, t]; hpd <- safe_hpd(sv)
      tibble(Region = i, Region_Name = region_names[i], 
             Cluster = cluster_ids[i], Time = t,
             Mean = mean(sv), Lower = hpd[1], Upper = hpd[2])
    }))
  }))
  write_csv(eps_full, file.path(scenario_dir, "epsilon_summary.csv"))
  
  # Agrega por cluster × tempo
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
        labeller = labeller(k = function(x) paste0("Cluster ", x))
      ) +
      scale_x_continuous(breaks = 1:n_times, labels = anos_label) +
      theme_bw(base_size = 11) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      labs(title = paste("Trajetória de epsilon por cluster (", model_type, ")"),
           subtitle = "Banda = IC 95% HPD",
           x = "Ano", y = expression(epsilon[kt])),
    width = 12, height = 8
  )
  
  # ── 5j. Sumário de beta ───────────────────────────────────────────────────────
  beta_summary <- do.call(rbind, lapply(beta_names, function(nm) {
    sv  <- samples_mat[, nm]; hpd <- safe_hpd(sv)
    tibble(Parameter = nm, Mean = mean(sv), SD = sd(sv),
           HPD_Lower = hpd[1], HPD_Upper = hpd[2],
           ESS  = as.numeric(effectiveSize(mcmc_list_full[, nm])),
           Rhat = safe_gelman(mcmc_list_full[, nm]))
  }))
  write_csv(beta_summary, file.path(scenario_dir, "beta_summary.csv"))
  
  # ── 5k. lambda[i,t] — sumários e painéis ─────────────────────────────────────
  
  # Sumário completo (N × T)
  lambda_full <- do.call(rbind, lapply(seq_len(N_regions), function(i) {
    do.call(rbind, lapply(seq_len(n_times), function(t) {
      nm <- lambda_names_mat[i, t]
      sv <- samples_mat[, nm]; hpd <- safe_hpd(sv)
      tibble(Region = i, Region_Name = region_names[i], Time = t,
             Mean = mean(sv), Lower = hpd[1], Upper = hpd[2],
             ESS = as.numeric(effectiveSize(mcmc_list_full[, nm])),
             Rhat = safe_gelman(mcmc_list_full[, nm]))
    }))
  }))
  write_csv(lambda_full, file.path(scenario_dir, "lambda_summary.csv"))
  
  # Painel lambda para regiões selecionadas (com nomes)
  lambda_sel <- lambda_full %>% filter(Region %in% REGIONS_INTEREST)
  
  # Criar labels com nome da região (abreviado para caber no gráfico)
  make_region_label <- function(names) {
    sapply(names, function(n) {
      n_clean <- gsub("_", " ", n)
      if (nchar(n_clean) > 25) {
        paste0(substr(n_clean, 1, 22), "...")
      } else {
        n_clean
      }
    })
  }
  
  lambda_sel$Region_Label <- make_region_label(lambda_sel$Region_Name)
  
  ggsave(
    file.path(scenario_dir, "painel_lambda.png"),
    ggplot(lambda_sel, aes(x = Time)) +
      geom_ribbon(aes(ymin = Lower, ymax = Upper), fill = "grey70", alpha = 0.5) +
      geom_line(aes(y = Mean), color = "black", linewidth = 0.8) +
      facet_wrap(~ Region_Label, scales = "free_y", ncol = 4) +
      scale_x_continuous(breaks = 1:n_times, labels = anos_label) +
      theme_bw(base_size = 11) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      labs(title = paste("lambda[i,t] - Regiões selecionadas (", model_type, ")"),
           subtitle = "Banda = IC 95% HPD",
           x = "Ano", y = expression(lambda[it])),
    width = 16, height = 12
  )
  
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
      tibble(Region = i, Region_Name = region_names[i], 
             Mean = mean(sv), SD = sd(sv),
             HPD_Lower = hpd[1], HPD_Upper = hpd[2],
             ESS = as.numeric(effectiveSize(mcmc_list_full[, s_names[i]])))
    }))
    write_csv(s_summary, file.path(scenario_dir, "s_summary.csv"))
    
    # Destacar BH e Araçuaí no gráfico de efeitos espaciais
    s_summary$highlight <- ifelse(s_summary$Region %in% c(bh_idx, aracuai_idx), 
                                  "Destaque", "Normal")
    
    ggsave(
      file.path(scenario_dir, "s_posterior.png"),
      ggplot(s_summary, aes(x = reorder(Region_Name, Mean), y = Mean, color = highlight)) +
        geom_point(size = 1.5) +
        geom_errorbar(aes(ymin = HPD_Lower, ymax = HPD_Upper),
                      width = 0.4, linewidth = 0.3) +
        geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
        scale_color_manual(values = c("Destaque" = "red", "Normal" = "black")) +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6),
              legend.position = "none") +
        labs(title = paste("Efeito espacial s[i] (", model_type, ")"),
             subtitle = "Vermelho: BH e Araçuaí",
             y = "s[i]", x = "Microrregião"),
      width = 14, height = 6
    )
  }
  
  # ── 5m. mu[i, t] para regiões selecionadas ───────────────────────────────────
  beta_cols <- grep("^beta\\[", colnames(samples_mat), value = TRUE)
  
  # Calcular mu_summary a partir dos monitores diretos
  mu_summary <- do.call(rbind, lapply(REGIONS_INTEREST, function(i) {
    do.call(rbind, lapply(seq_len(n_times), function(t) {
      nm <- paste0("mu[", i, ", ", t, "]")
      if (nm %in% colnames(samples_mat)) {
        sv <- samples_mat[, nm]; hpd <- safe_hpd(sv)
        tibble(Region = i, Region_Name = region_names[i], Time = t, 
               Mean = mean(sv), Lower = hpd[1], Upper = hpd[2], model = model_type)
      } else {
        # Calcular manualmente se não estiver nos monitores
        ldraws   <- samples_mat[, lambda_names_mat[i, t]]
        bdraws   <- samples_mat[, beta_cols, drop = FALSE]
        eps_it   <- epsilon_draws[, i, t]
        s_draws  <- if (is_spatial) samples_mat[, paste0("s[", i, "]")] else 0
        mu_draws <- ldraws * exp(as.vector(bdraws %*% data_nimble$x[i, t, ])) *
          data_nimble$E[i, t] * eps_it * exp(s_draws)
        hpd <- safe_hpd(mu_draws)
        tibble(Region = i, Region_Name = region_names[i], Time = t,
               Mean = mean(mu_draws), Lower = hpd[1], Upper = hpd[2], model = model_type)
      }
    }))
  }))
  write_csv(mu_summary, file.path(scenario_dir, "mu_selected.csv"))
  
  mu_summary$Region_Label <- make_region_label(mu_summary$Region_Name)
  
  ggsave(
    file.path(scenario_dir, "painel_mu.png"),
    ggplot(mu_summary, aes(x = Time)) +
      geom_ribbon(aes(ymin = Lower, ymax = Upper), fill = "darkorange", alpha = 0.3) +
      geom_line(aes(y = Mean), color = "darkorange", linewidth = 0.8) +
      facet_wrap(~ Region_Label, scales = "free_y", ncol = 4) +
      scale_x_continuous(breaks = 1:n_times, labels = anos_label) +
      theme_bw(base_size = 10) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      labs(title = paste("Mu estimado (", model_type, ")"),
           subtitle = "Banda = IC 95% HPD",
           x = "Ano", y = expression(mu[it])),
    width = 16, height = 12
  )
  
  # ── 5n. WAIC e LPML ───────────────────────────────────────────────────────────
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
  
  # ── 5o. Diagnósticos ACF ──────────────────────────────────────────────────────
  cores_cadeia <- c("Cadeia 1" = "#2166AC", "Cadeia 2" = "#D6604D")
  
  params_struct <- c(beta_names,
                     if (is_spatial) "tau_s" else character(0))
  
  # ACF de beta
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
      labs(title    = paste("ACF de beta e tau_s (", model_type, ")"),
           subtitle = "Azul: |0.10| | Vermelho: |0.05|"),
    width = 10, height = 5
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
  
  # Traceplots de beta
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
      labs(title    = paste("Traceplots beta e tau_s (", model_type, ")"),
           subtitle = "Linha grossa = média ergódica | Linha fina = cadeia",
           x = "Iteração (pós-burnin)", y = "Valor"),
    width = 10, height = max(5, 3 * ceiling(length(params_acf_beta) / 3))
  )
  
  # ── 5p. Retorno resumido ──────────────────────────────────────────────────────
  tibble(
    model = model_type, niter = niter, nburnin = nburnin, thin = thin,
    WAIC = waic, LPML = LPML,
    ESS_beta_min    = min(beta_summary$ESS,   na.rm = TRUE),
    ESS_gamma_min   = min(gamma_summary$ESS,  na.rm = TRUE),
    ESS_lambda_min  = min(lambda_full$ESS, na.rm = TRUE),
    ESS_tau         = ESS_tau,
    Rhat_beta_max   = max(beta_summary$Rhat,   na.rm = TRUE),
    Rhat_gamma_max  = max(gamma_summary$Rhat, na.rm = TRUE),
    Rhat_lambda_max = max(lambda_full$Rhat, na.rm = TRUE)
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
  "h_mat", "cluster_ids", "region_names", "anos_label",
  "run_model", "output_dir", "REGIONS_INTEREST", 
  "bh_idx", "aracuai_idx"
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

# Função auxiliar para ler CSVs
read_model_csv <- function(filename) {
  lapply(model_types, function(m) {
    path <- file.path(output_dir, m, filename)
    if (!file.exists(path)) return(NULL)
    read_csv(path, show_col_types = FALSE) %>% mutate(model = m)
  }) |> bind_rows()
}

# Comparativo gamma[k, t]
gamma_all <- read_model_csv("gamma_summary.csv")
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
      facet_wrap(~ k, ncol = 2, scales = "free_y",
                 labeller = labeller(k = function(x) paste0("Cluster ", x))) +
      scale_x_continuous(breaks = 1:n_times, labels = anos_label) +
      theme_bw(base_size = 12) + 
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            legend.position = "bottom") +
      labs(title = "Comparação gamma[k,t]: espacial vs. não-espacial",
           subtitle = "Banda cinza em k=1: suporte da priori estática [0, 0.1]",
           x = "Ano", y = expression(gamma[kt]),
           color = "Modelo", fill = "Modelo"),
    width = 12, height = 8, dpi = 300
  )
}

# Comparativo lambda
lambda_all <- read_model_csv("lambda_summary.csv")
if (nrow(lambda_all) > 0) {
  # Filtrar regiões selecionadas e criar labels
  lambda_sel <- lambda_all %>% 
    filter(Region %in% REGIONS_INTEREST)
  
  make_region_label <- function(names) {
    sapply(names, function(n) {
      n_clean <- gsub("_", " ", n)
      if (nchar(n_clean) > 25) {
        paste0(substr(n_clean, 1, 22), "...")
      } else {
        n_clean
      }
    })
  }
  
  lambda_sel$Region_Label <- make_region_label(lambda_sel$Region_Name)
  
  ggsave(
    file.path(output_dir, "lambda_comparativo.png"),
    ggplot(lambda_sel, aes(x = Time, y = Mean, color = model, fill = model)) +
      geom_ribbon(aes(ymin = Lower, ymax = Upper), alpha = 0.15, color = NA) +
      geom_line(linewidth = 0.8) +
      facet_wrap(~ Region_Label, scales = "free_y", ncol = 4) +
      scale_x_continuous(breaks = 1:n_times, labels = anos_label) +
      theme_bw(base_size = 11) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
            legend.position = "bottom") +
      labs(title = "Comparação lambda[i,t]: espacial vs. não-espacial",
           subtitle = "Microrregiões selecionadas (incluindo BH e Araçuaí)",
           x = "Ano", y = expression(lambda[it]),
           color = "Modelo", fill = "Modelo"),
    width = 16, height = 12, dpi = 300
  )
}

# Comparativo mu
mu_all <- read_model_csv("mu_selected.csv")
if (nrow(mu_all) > 0) {
  mu_all$Region_Label <- make_region_label(mu_all$Region_Name)
  
  ggsave(
    file.path(output_dir, "mu_comparativo.png"),
    ggplot(mu_all, aes(x = Time, y = Mean, color = model, fill = model)) +
      geom_ribbon(aes(ymin = Lower, ymax = Upper), alpha = 0.15, color = NA) +
      geom_line(linewidth = 0.8) +
      facet_wrap(~ Region_Label, scales = "free_y", ncol = 4) +
      scale_x_continuous(breaks = 1:n_times, labels = anos_label) +
      theme_bw(base_size = 11) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
            legend.position = "bottom") +
      labs(title = "Comparação mu[i,t]: espacial vs. não-espacial",
           subtitle = "Microrregiões selecionadas (incluindo BH e Araçuaí)",
           x = "Ano", y = expression(mu[it]),
           color = "Modelo", fill = "Modelo"),
    width = 16, height = 12, dpi = 300
  )
}

# Comparativo epsilon por cluster
eps_cluster_all <- read_model_csv("epsilon_cluster_summary.csv")
if (nrow(eps_cluster_all) > 0) {
  ggsave(
    file.path(output_dir, "epsilon_comparativo.png"),
    ggplot(eps_cluster_all, aes(x = Time, y = Mean, color = model, fill = model)) +
      geom_ribbon(aes(ymin = Lower, ymax = Upper), alpha = 0.12, color = NA) +
      geom_line(linewidth = 0.8) +
      facet_wrap(~ k, ncol = 2, scales = "free_y",
                 labeller = labeller(k = function(x) paste0("Cluster ", x))) +
      scale_x_continuous(breaks = 1:n_times, labels = anos_label) +
      theme_bw(base_size = 12) + 
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            legend.position = "bottom") +
      labs(title = "Comparação epsilon por cluster: espacial vs. não-espacial",
           subtitle = "Banda = IC 95% HPD",
           x = "Ano", y = expression(epsilon[kt]),
           color = "Modelo", fill = "Modelo"),
    width = 12, height = 8, dpi = 300
  )
}

# Comparativo beta
beta_all <- read_model_csv("beta_summary.csv")
if (nrow(beta_all) > 0) {
  ggsave(
    file.path(output_dir, "beta_comparativo.png"),
    ggplot(beta_all, aes(x = Parameter, y = Mean, color = model)) +
      geom_point(position = position_dodge(width = 0.5), size = 3) +
      geom_errorbar(aes(ymin = HPD_Lower, ymax = HPD_Upper),
                    position = position_dodge(width = 0.5), width = 0.2) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
      theme_bw(base_size = 12) + theme(legend.position = "bottom") +
      labs(title = "Comparação dos coeficientes beta: espacial vs. não-espacial",
           subtitle = "Barras = IC 95% HPD",
           x = "Covariável", y = expression(beta),
           color = "Modelo"),
    width = 8, height = 5, dpi = 300
  )
}

cat("\n========================================\n")
cat("Tempo total de execução:\n")
print(Sys.time() - inicio_global)
cat("\nResultados salvos em:", output_dir, "\n")
cat("========================================\n")