# ==============================================================================
# ajuste_cluster1_lambda_global_vs_local.R
# Modelo Bayesiano para Cluster 1: lambda[t] vs lambda[i,t]
#
# COMPARAÇÃO:
#   Modelo A (lambda_global): lambda[t] compartilhado entre todas as regiões
#   Modelo B (lambda_local):  lambda[i,t] específico por região
#
# AMBOS OS MODELOS:
#   - Apenas regiões do Cluster 1 (melhor reportação)
#   - epsilon[t] ~ U(0.95, 1) i.i.d.
#   - Efeito espacial ICAR: s[i]
#   - Covariáveis fixas: beta
#
# PREDITOR:
#   Modelo A: mu[i,t] = lambda[t]   * E[i,t] * epsilon[t] * exp(x'beta + s[i])
#   Modelo B: mu[i,t] = lambda[i,t] * E[i,t] * epsilon[t] * exp(x'beta + s[i])
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

source("_dataCaseStudy.r")

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

# ==============================================================================
# 2. SELECIONAR REGIÕES DO CLUSTER 1
# ==============================================================================

h_mat       <- data$hAI
cluster_ids <- apply(h_mat, 1, sum)

# Apenas cluster 1
cluster1_idx <- which(cluster_ids == 1)
cat(sprintf("\nRegiões no Cluster 1: %d de %d\n", length(cluster1_idx), N_regions))
cat("Regiões selecionadas:\n")
print(data.frame(
  Indice_Original = cluster1_idx,
  Nome = rownames(Y_mat)[cluster1_idx]
))

# Filtrar dados
Y_c1 <- Y_mat[cluster1_idx, , drop = FALSE]
E_c1 <- E_norm[cluster1_idx, , drop = FALSE]
x_c1 <- x[cluster1_idx, , , drop = FALSE]

N_c1 <- nrow(Y_c1)
cat(sprintf("N_cluster1 = %d | T = %d | p = %d\n", N_c1, n_times, p))

# ==============================================================================
# 3. MATRIZ DE ADJACÊNCIA PARA CLUSTER 1
# ==============================================================================

adj_original <- data$adj
num_original <- data$num

old_to_new <- setNames(seq_along(cluster1_idx), cluster1_idx)

adj_c1 <- c()
num_c1 <- integer(N_c1)
weights_c1 <- c()

for (i_new in seq_len(N_c1)) {
  i_old <- cluster1_idx[i_new]
  
  if (num_original[i_old] > 0) {
    start_idx <- if (i_old == 1) 1 else sum(num_original[1:(i_old-1)]) + 1
    end_idx   <- start_idx + num_original[i_old] - 1
    neighbors_old <- adj_original[start_idx:end_idx]
    
    neighbors_new <- old_to_new[as.character(neighbors_old)]
    neighbors_new <- neighbors_new[!is.na(neighbors_new)]
    
    num_c1[i_new] <- length(neighbors_new)
    if (length(neighbors_new) > 0) {
      adj_c1 <- c(adj_c1, neighbors_new)
      weights_c1 <- c(weights_c1, rep(1, length(neighbors_new)))
    }
  }
}

n_adj_c1 <- length(adj_c1)
cat(sprintf("Subgrafo Cluster 1: %d regiões, %d arestas\n", N_c1, n_adj_c1))

# Verificar regiões isoladas
isoladas <- which(num_c1 == 0)
if (length(isoladas) > 0) {
  cat("ATENÇÃO: Regiões isoladas no cluster 1:\n")
  print(data.frame(Indice = isoladas, Nome = rownames(Y_c1)[isoladas]))
}

# ==============================================================================
# 4. REGIÕES DE INTERESSE
# ==============================================================================

region_names_c1 <- rownames(Y_c1)
anos_label <- colnames(Y_mat)
if (is.null(anos_label)) anos_label <- as.character(seq_len(n_times))

# Encontrar BH
bh_idx_c1 <- grep("BELO_HORIZONTE", region_names_c1, ignore.case = TRUE)
if (length(bh_idx_c1) > 0) {
  cat(sprintf("BH: índice %d - %s\n", bh_idx_c1[1], region_names_c1[bh_idx_c1[1]]))
}

# Selecionar regiões para monitorar
set.seed(42)
n_monitor <- min(12, N_c1)
REGIONS_INTEREST <- if (length(bh_idx_c1) > 0) {
  c(bh_idx_c1[1], sample(setdiff(seq_len(N_c1), bh_idx_c1), 
                         min(n_monitor - 1, N_c1 - length(bh_idx_c1))))
} else {
  sample(seq_len(N_c1), n_monitor)
}
REGIONS_INTEREST <- sort(REGIONS_INTEREST)

cat("\nRegiões monitoradas:\n")
print(data.frame(
  Indice = REGIONS_INTEREST,
  Nome = region_names_c1[REGIONS_INTEREST]
))

# ==============================================================================
# 5. HIPERPARÂMETROS E CONSTANTES
# ==============================================================================

# Gama-DLM
a0 <- 1.0
b0 <- 1.0
w  <- 0.85

# Epsilon
eps_lower <- 0.95
eps_upper <- 1.00

# Constantes comuns
constants_base <- list(
  n_regions = N_c1,
  n_times   = n_times,
  p         = p,
  a0 = a0, b0 = b0, w = w,
  eps_lower = eps_lower,
  eps_upper = eps_upper,
  adj     = adj_c1,
  num     = num_c1,
  weights = weights_c1,
  n_adj   = n_adj_c1
)

data_nimble <- list(Y = Y_c1, E = E_c1, x = x_c1)

# ==============================================================================
# 6. INICIALIZAÇÕES
# ==============================================================================

set.seed(123)
epsilon_init <- rep(0.975, n_times)

# ── Modelo lambda global ─────────────────────────────────────────────────
lambda_global_init1 <- rep(1.0, n_times)
lambda_global_init2 <- rgamma(n_times, 1, 1)

inits_global_1 <- list(
  lambda  = lambda_global_init1,
  beta    = rep(0, p),
  epsilon = epsilon_init,
  sigma_s = 0.5,
  s       = rep(0, N_c1)
)
inits_global_2 <- list(
  lambda  = lambda_global_init2,
  beta    = rnorm(p, 0, 0.3),
  epsilon = runif(n_times, 0.96, 0.99),
  sigma_s = 1.0,
  s       = rep(0, N_c1)
)

# ── Modelo lambda local ──────────────────────────────────────────────────
lambda_local_init1 <- matrix(1.0, nrow = N_c1, ncol = n_times)
lambda_local_init2 <- matrix(rgamma(N_c1 * n_times, 1, 1), 
                             nrow = N_c1, ncol = n_times)

inits_local_1 <- list(
  lambda  = lambda_local_init1,
  beta    = rep(0, p),
  epsilon = epsilon_init,
  sigma_s = 0.5,
  s       = rep(0, N_c1)
)
inits_local_2 <- list(
  lambda  = lambda_local_init2,
  beta    = rnorm(p, 0, 0.3),
  epsilon = runif(n_times, 0.96, 0.99),
  sigma_s = 1.0,
  s       = rep(0, N_c1)
)

inits_list_global <- list(inits_global_1, inits_global_2)
inits_list_local  <- list(inits_local_1, inits_local_2)

# ==============================================================================
# 7. FUNÇÃO WORKER
# ==============================================================================

run_model <- function(model_type, output_dir) {
  
  library(nimble); library(coda); library(dplyr)
  library(ggplot2); library(readr); library(stringr); library(tibble)
  
  is_global <- (model_type == "lambda_global")
  model_label <- if (is_global) "Lambda Global" else "Lambda Local"
  
  # ── 7a. Código NIMBLE ──────────────────────────────────────────────────
  
  code_global <- nimbleCode({
    
    for (j in 1:p) {
      beta[j] ~ dnorm(0, sd = 1)
    }
    
    for (t in 1:n_times) {
      epsilon[t] ~ dunif(eps_lower, eps_upper)
    }
    
    sigma_s ~ T(dt(0, 1, 1), 0, )
    tau_s   <- 1 / (sigma_s^2)
    s[1:n_regions] ~ dcar_normal(
      adj[1:n_adj], weights[1:n_adj], num[1:n_regions], tau_s, zero_mean = 1
    )
    
    # Lambda[t] global - prior declarativo (substituído por FFBS)
    for (t in 1:n_times) {
      lambda[t] ~ dgamma(a0, b0)
    }
    
    for (i in 1:n_regions) {
      for (t in 1:n_times) {
        mu[i, t] <- lambda[t] * E[i, t] * epsilon[t] *
          exp(inprod(beta[1:p], x[i, t, 1:p]) + s[i])
        Y[i, t]        ~ dpois(mu[i, t])
        logLik_Y[i, t] <- dpois(Y[i, t], mu[i, t], log = TRUE)
      }
    }
  })
  
  code_local <- nimbleCode({
    
    for (j in 1:p) {
      beta[j] ~ dnorm(0, sd = 1)
    }
    
    for (t in 1:n_times) {
      epsilon[t] ~ dunif(eps_lower, eps_upper)
    }
    
    sigma_s ~ T(dt(0, 1, 1), 0, )
    tau_s   <- 1 / (sigma_s^2)
    s[1:n_regions] ~ dcar_normal(
      adj[1:n_adj], weights[1:n_adj], num[1:n_regions], tau_s, zero_mean = 1
    )
    
    # Lambda[i,t] local - prior declarativo (substituído por FFBS)
    for (i in 1:n_regions) {
      for (t in 1:n_times) {
        lambda[i, t] ~ dgamma(a0, b0)
      }
    }
    
    for (i in 1:n_regions) {
      for (t in 1:n_times) {
        mu[i, t] <- lambda[i, t] * E[i, t] * epsilon[t] *
          exp(inprod(beta[1:p], x[i, t, 1:p]) + s[i])
        Y[i, t]        ~ dpois(mu[i, t])
        logLik_Y[i, t] <- dpois(Y[i, t], mu[i, t], log = TRUE)
      }
    }
  })
  
  # ── 7b. FFBS para lambda global ────────────────────────────────────────
  
  ffbs_global <- nimbleFunction(
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
      declare(i, integer())
      declare(t, integer())
      declare(t_idx, integer())
      declare(t_back, integer())
      declare(k, integer())
      declare(prod_val, double())
      declare(g_it, double())
      declare(sum_Y, double())
      declare(sum_g, double())
      declare(nu, double())
      
      at_buf[1] <<- a0
      bt_buf[1] <<- b0
      
      # Forward: agrega sobre todas as regiões
      for (t in 1:n_times) {
        sum_Y <- 0.0
        sum_g <- 0.0
        for (i in 1:n_regions) {
          sum_Y <- sum_Y + model$Y[i, t]
          prod_val <- 0.0
          for (k in 1:p) {
            prod_val <- prod_val + model$x[i, t, k] * model$beta[k]
          }
          g_it <- model$E[i, t] * model$epsilon[t] * 
            exp(prod_val + model$s[i])
          sum_g <- sum_g + g_it
        }
        at_buf[t + 1] <<- w * at_buf[t] + sum_Y
        bt_buf[t + 1] <<- w * bt_buf[t] + sum_g
      }
      
      # Backward
      model$lambda[n_times] <<- rgamma(1, shape = at_buf[n_times + 1],
                                       rate  = max(bt_buf[n_times + 1], 1e-10))
      for (t_idx in 1:(n_times - 1)) {
        t_back <- n_times - t_idx
        nu <- rgamma(1, shape = (1 - w) * at_buf[t_back + 1],
                     rate  = max(bt_buf[t_back + 1], 1e-10))
        model$lambda[t_back] <<- nu + w * model$lambda[t_back + 1]
      }
      
      model$calculate(calcNodes)
      copy(from = model, to = mvSaved, row = 1, nodes = targetNodes, logProb = TRUE)
    },
    methods = list(reset = function() {})
  )
  
  # ── 7c. FFBS para lambda local ─────────────────────────────────────────
  
  ffbs_local <- nimbleFunction(
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
      declare(t, integer())
      declare(t_idx, integer())
      declare(t_back, integer())
      declare(k, integer())
      declare(prod_val, double())
      declare(g_it, double())
      declare(nu, double())
      
      at_buf[1] <<- a0
      bt_buf[1] <<- b0
      
      # Forward: apenas para a região i
      for (t in 1:n_times) {
        prod_val <- 0.0
        for (k in 1:p) {
          prod_val <- prod_val + model$x[region_i, t, k] * model$beta[k]
        }
        g_it <- model$E[region_i, t] * model$epsilon[t] * 
          exp(prod_val + model$s[region_i])
        
        at_buf[t + 1] <<- w * at_buf[t] + model$Y[region_i, t]
        bt_buf[t + 1] <<- w * bt_buf[t] + g_it
      }
      
      # Backward
      model$lambda[region_i, n_times] <<- rgamma(
        1, shape = at_buf[n_times + 1],
        rate = max(bt_buf[n_times + 1], 1e-10)
      )
      for (t_idx in 1:(n_times - 1)) {
        t_back <- n_times - t_idx
        nu <- rgamma(1, shape = (1 - w) * at_buf[t_back + 1],
                     rate  = max(bt_buf[t_back + 1], 1e-10))
        model$lambda[region_i, t_back] <<- nu + w * model$lambda[region_i, t_back + 1]
      }
      
      model$calculate(calcNodes)
      copy(from = model, to = mvSaved, row = 1, nodes = targetNodes, logProb = TRUE)
    },
    methods = list(reset = function() {})
  )
  
  # ── 7d. Configurar modelo ──────────────────────────────────────────────
  
  model_code <- if (is_global) code_global else code_local
  inits_list <- if (is_global) inits_list_global else inits_list_local
  ffbs_fn    <- if (is_global) ffbs_global else ffbs_local
  
  cat(sprintf("\n=== Modelo: %s ===\n", model_label))
  scenario_dir <- file.path(output_dir, model_type)
  dir.create(scenario_dir, recursive = TRUE, showWarnings = FALSE)
  
  model  <- nimbleModel(code = model_code, constants = constants_base,
                        data = data_nimble, inits = inits_list[[1]], check = FALSE)
  Cmodel <- compileNimble(model)
  
  conf <- configureMCMC(model)
  
  # FFBS para lambda
  conf$removeSamplers("lambda")
  if (is_global) {
    # Um único FFBS para lambda[1:T]
    conf$addSampler(
      target  = paste0("lambda[1:", n_times, "]"),
      type    = ffbs_fn,
      control = list(
        n_regions = N_c1,
        n_times   = n_times,
        p         = p,
        a0        = a0,
        b0        = b0,
        w         = w
      )
    )
  } else {
    # Um FFBS por região
    for (i in seq_len(N_c1)) {
      conf$addSampler(
        target  = paste0("lambda[", i, ", 1:", n_times, "]"),
        type    = ffbs_fn,
        control = list(
          n_times   = n_times,
          p         = p,
          a0        = a0,
          b0        = b0,
          w         = w,
          region_i  = i
        )
      )
    }
  }
  
  # Epsilon: slice por t
  conf$removeSamplers("epsilon")
  for (t in seq_len(n_times)) {
    conf$addSampler(target = paste0("epsilon[", t, "]"), type = "slice")
  }
  
  # Monitores
  monitors_base <- c("beta", "epsilon", "logLik_Y", "lambda", "s", "sigma_s", "tau_s")
  
  # mu para regiões de interesse
  mu_monitors <- unlist(lapply(
    REGIONS_INTEREST,
    function(r) paste0("mu[", r, ", ", seq_len(n_times), "]")
  ))
  
  conf$addMonitors(monitors_base)
  conf$addMonitors(mu_monitors)
  
  Rmcmc <- buildMCMC(conf)
  Cmcmc <- compileNimble(Rmcmc, project = model)
  
  # ── 7e. Executar MCMC ──────────────────────────────────────────────────
  
  niter   <- 50000
  nburnin <- 10000
  nchains <- 2
  thin    <- 10
  
  cat(sprintf("[%s] niter=%d | nburnin=%d | thin=%d\n",
              model_type, niter, nburnin, thin))
  
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
  
  # ── 7f. Funções auxiliares ─────────────────────────────────────────────
  
  safe_hpd <- function(sv) {
    if (var(sv) < 1e-12) return(c(NA_real_, NA_real_))
    as.numeric(HPDinterval(as.mcmc(sv), prob = 0.95))
  }
  
  safe_gelman <- function(obj) {
    tryCatch(gelman.diag(obj, autoburnin = FALSE)$psrf[, 1],
             error = function(e) rep(NA_real_, nvar(obj)))
  }
  
  make_region_label <- function(names) {
    sapply(names, function(n) {
      n_clean <- gsub("_", " ", n)
      if (nchar(n_clean) > 25) paste0(substr(n_clean, 1, 22), "...") else n_clean
    })
  }
  
  beta_names    <- paste0("beta[",    seq_len(p),       "]")
  epsilon_names <- paste0("epsilon[", seq_len(n_times), "]")
  
  if (is_global) {
    lambda_names <- paste0("lambda[", seq_len(n_times), "]")
  } else {
    lambda_names <- as.vector(outer(
      seq_len(N_c1), seq_len(n_times),
      function(i, t) paste0("lambda[", i, ", ", t, "]")
    ))
  }
  
  # ── 7g. Sumário de epsilon[t] ─────────────────────────────────────────
  
  epsilon_summary <- do.call(rbind, lapply(seq_len(n_times), function(t) {
    nm  <- epsilon_names[t]
    sv  <- samples_mat[, nm]
    hpd <- safe_hpd(sv)
    tibble(
      Time = t,
      Mean = mean(sv), SD = sd(sv),
      Lower = hpd[1], Upper = hpd[2],
      ESS  = as.numeric(effectiveSize(mcmc_list_full[, nm])),
      Rhat = safe_gelman(mcmc_list_full[, nm]),
      model = model_type
    )
  }))
  write_csv(epsilon_summary, file.path(scenario_dir, "epsilon_summary.csv"))
  
  # ── 7h. Sumário de beta ───────────────────────────────────────────────
  
  beta_summary <- do.call(rbind, lapply(beta_names, function(nm) {
    sv  <- samples_mat[, nm]
    hpd <- safe_hpd(sv)
    tibble(
      Parameter = nm,
      Mean = mean(sv), SD = sd(sv),
      HPD_Lower = hpd[1], HPD_Upper = hpd[2],
      ESS  = as.numeric(effectiveSize(mcmc_list_full[, nm])),
      Rhat = safe_gelman(mcmc_list_full[, nm]),
      model = model_type
    )
  }))
  write_csv(beta_summary, file.path(scenario_dir, "beta_summary.csv"))
  
  # ── 7i. Sumário de lambda ─────────────────────────────────────────────
  
  if (is_global) {
    lambda_summary <- do.call(rbind, lapply(seq_len(n_times), function(t) {
      nm <- lambda_names[t]
      sv <- samples_mat[, nm]
      hpd <- safe_hpd(sv)
      tibble(
        Region = NA_integer_,
        Region_Name = "Global",
        Time = t,
        Mean = mean(sv), Lower = hpd[1], Upper = hpd[2],
        ESS  = as.numeric(effectiveSize(mcmc_list_full[, nm])),
        Rhat = safe_gelman(mcmc_list_full[, nm]),
        model = model_type
      )
    }))
  } else {
    lambda_summary <- do.call(rbind, lapply(seq_len(N_c1), function(i) {
      do.call(rbind, lapply(seq_len(n_times), function(t) {
        nm <- paste0("lambda[", i, ", ", t, "]")
        sv <- samples_mat[, nm]
        hpd <- safe_hpd(sv)
        tibble(
          Region = i,
          Region_Name = region_names_c1[i],
          Time = t,
          Mean = mean(sv), Lower = hpd[1], Upper = hpd[2],
          ESS  = as.numeric(effectiveSize(mcmc_list_full[, nm])),
          Rhat = safe_gelman(mcmc_list_full[, nm]),
          model = model_type
        )
      }))
    }))
  }
  write_csv(lambda_summary, file.path(scenario_dir, "lambda_summary.csv"))
  
  # ── 7j. Efeitos espaciais ──────────────────────────────────────────────
  
  tau_sv  <- samples_mat[, "tau_s"]
  hpd_t   <- safe_hpd(tau_sv)
  tau_sum <- tibble(
    Parameter = "tau_s",
    Mean = mean(tau_sv), SD = sd(tau_sv),
    HPD_Lower = hpd_t[1], HPD_Upper = hpd_t[2],
    ESS  = as.numeric(effectiveSize(mcmc_list_full[, "tau_s"])),
    Rhat = safe_gelman(mcmc_list_full[, "tau_s"]),
    model = model_type
  )
  write_csv(tau_sum, file.path(scenario_dir, "tau_summary.csv"))
  
  s_names   <- paste0("s[", seq_len(N_c1), "]")
  s_summary <- do.call(rbind, lapply(seq_len(N_c1), function(i) {
    sv  <- samples_mat[, s_names[i]]
    hpd <- safe_hpd(sv)
    tibble(
      Region = i, Region_Name = region_names_c1[i],
      Mean = mean(sv), SD = sd(sv),
      HPD_Lower = hpd[1], HPD_Upper = hpd[2],
      ESS = as.numeric(effectiveSize(mcmc_list_full[, s_names[i]])),
      model = model_type
    )
  }))
  write_csv(s_summary, file.path(scenario_dir, "s_summary.csv"))
  
  # ── 7k. mu para regiões monitoradas ───────────────────────────────────
  
  mu_summary <- do.call(rbind, lapply(REGIONS_INTEREST, function(i) {
    do.call(rbind, lapply(seq_len(n_times), function(t) {
      nm <- paste0("mu[", i, ", ", t, "]")
      sv <- samples_mat[, nm]
      hpd <- safe_hpd(sv)
      tibble(
        Region = i, Region_Name = region_names_c1[i],
        Time = t,
        Mean = mean(sv), Lower = hpd[1], Upper = hpd[2],
        model = model_type
      )
    }))
  }))
  write_csv(mu_summary, file.path(scenario_dir, "mu_selected.csv"))
  
  # ── 7l. Painéis ───────────────────────────────────────────────────────
  
  # Epsilon
  ggsave(
    file.path(scenario_dir, "painel_epsilon.png"),
    ggplot(epsilon_summary, aes(x = Time)) +
      geom_ribbon(aes(ymin = Lower, ymax = Upper), fill = "darkorange", alpha = 0.3) +
      geom_line(aes(y = Mean), color = "darkorange", linewidth = 0.8) +
      geom_hline(yintercept = eps_lower, linetype = "dashed", color = "grey50") +
      geom_hline(yintercept = eps_upper, linetype = "dashed", color = "grey50") +
      scale_x_continuous(breaks = 1:n_times, labels = anos_label) +
      theme_bw(base_size = 12) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      labs(title = paste("epsilon[t] ~ U(0.95,1) -", model_label),
           subtitle = "Linhas tracejadas: limites da priori",
           x = "Ano", y = expression(epsilon[t])),
    width = 10, height = 5
  )
  
  # Lambda
  if (is_global) {
    ggsave(
      file.path(scenario_dir, "painel_lambda.png"),
      ggplot(lambda_summary, aes(x = Time)) +
        geom_ribbon(aes(ymin = Lower, ymax = Upper), fill = "grey70", alpha = 0.5) +
        geom_line(aes(y = Mean), color = "black", linewidth = 0.9) +
        scale_x_continuous(breaks = 1:n_times, labels = anos_label) +
        theme_bw(base_size = 12) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        labs(title = paste("lambda[t] Global -", model_label),
             subtitle = "Banda = IC 95% HPD",
             x = "Ano", y = expression(lambda[t])),
      width = 10, height = 5
    )
  } else {
    lambda_sel <- lambda_summary %>% filter(Region %in% REGIONS_INTEREST)
    lambda_sel$Region_Label <- make_region_label(lambda_sel$Region_Name)
    
    ggsave(
      file.path(scenario_dir, "painel_lambda.png"),
      ggplot(lambda_sel, aes(x = Time)) +
        geom_ribbon(aes(ymin = Lower, ymax = Upper), fill = "grey70", alpha = 0.5) +
        geom_line(aes(y = Mean), color = "black", linewidth = 0.8) +
        facet_wrap(~ Region_Label, scales = "free_y", ncol = 4) +
        scale_x_continuous(breaks = 1:n_times, labels = anos_label) +
        theme_bw(base_size = 10) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        labs(title = paste("lambda[i,t] Local -", model_label),
             subtitle = "Banda = IC 95% HPD",
             x = "Ano", y = expression(lambda[it])),
      width = 16, height = 12
    )
  }
  
  # Efeitos espaciais
  s_summary$highlight <- ifelse(s_summary$Region %in% bh_idx_c1, "BH", "Outras")
  
  ggsave(
    file.path(scenario_dir, "s_posterior.png"),
    ggplot(s_summary, aes(x = reorder(Region_Name, Mean), y = Mean, color = highlight)) +
      geom_point(size = 1.5) +
      geom_errorbar(aes(ymin = HPD_Lower, ymax = HPD_Upper), width = 0.4, linewidth = 0.3) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
      scale_color_manual(values = c("BH" = "red", "Outras" = "black")) +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6),
            legend.position = "none") +
      labs(title = paste("Efeito espacial s[i] -", model_label),
           subtitle = if (length(bh_idx_c1) > 0) "Vermelho: BH" else "",
           y = "s[i]", x = "Microrregião"),
    width = 14, height = 6
  )
  
  # Mu
  mu_summary$Region_Label <- make_region_label(mu_summary$Region_Name)
  
  ggsave(
    file.path(scenario_dir, "painel_mu.png"),
    ggplot(mu_summary, aes(x = Time)) +
      geom_ribbon(aes(ymin = Lower, ymax = Upper), fill = "steelblue", alpha = 0.3) +
      geom_line(aes(y = Mean), color = "steelblue", linewidth = 0.8) +
      facet_wrap(~ Region_Label, scales = "free_y", ncol = 4) +
      scale_x_continuous(breaks = 1:n_times, labels = anos_label) +
      theme_bw(base_size = 10) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      labs(title = paste("Mu estimado -", model_label),
           subtitle = "Banda = IC 95% HPD",
           x = "Ano", y = expression(mu[it])),
    width = 16, height = 12
  )
  
  # ── 7m. WAIC e LPML ───────────────────────────────────────────────────
  
  loglik_names <- grep("logLik_Y", colnames(samples_mat), value = TRUE)
  waic <- NA_real_; LPML <- NA_real_
  if (length(loglik_names) > 0) {
    lm     <- samples_mat[, loglik_names, drop = FALSE]
    lppd   <- sum(apply(lm, 2, function(x) { mx <- max(x); mx + log(mean(exp(x - mx))) }))
    p_waic <- sum(apply(lm, 2, var))
    waic   <- -2 * (lppd - p_waic)
    LPML   <- sum(log(1 / apply(lm, 2, function(x) mean(exp(-x)))))
    write_csv(
      tibble(WAIC = waic, LPML = LPML, lppd = lppd, pWAIC = p_waic, model = model_type),
      file.path(scenario_dir, "criteria.csv")
    )
    cat(sprintf("[%s] WAIC = %.2f | LPML = %.2f\n", model_type, waic, LPML))
  }
  
  # ── 7n. Retorno resumido ───────────────────────────────────────────────
  
  tibble(
    model = model_type,
    N_regions = N_c1,
    niter = niter, nburnin = nburnin, thin = thin,
    WAIC = waic, LPML = LPML,
    ESS_beta_min   = min(beta_summary$ESS, na.rm = TRUE),
    ESS_eps_min    = min(epsilon_summary$ESS, na.rm = TRUE),
    ESS_lambda_min = min(lambda_summary$ESS, na.rm = TRUE),
    ESS_tau        = tau_sum$ESS[1],
    Rhat_beta_max   = max(beta_summary$Rhat, na.rm = TRUE),
    Rhat_eps_max    = max(epsilon_summary$Rhat, na.rm = TRUE),
    Rhat_lambda_max = max(lambda_summary$Rhat, na.rm = TRUE)
  )
}

# ==============================================================================
# 8. EXECUÇÃO PARALELA: lambda_global vs lambda_local
# ==============================================================================

model_types <- c("lambda_global", "lambda_local")
n_cores     <- min(length(model_types), parallel::detectCores() - 1)
if (n_cores < 1) n_cores <- 1

output_dir <- file.path(
  "C:/Users/vlara/OneDrive/Estatistica UFMG/Mestrado/Pesquisa/Aplicação",
  "resultados_cluster1_lambda_comparacao"
)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

cl <- makeCluster(n_cores)
clusterExport(cl, c(
  "constants_base", "data_nimble",
  "inits_list_global", "inits_list_local",
  "N_c1", "n_times", "p",
  "run_model", "output_dir", "REGIONS_INTEREST",
  "region_names_c1", "anos_label", "bh_idx_c1",
  "eps_lower", "eps_upper", "a0", "b0", "w"
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
# 9. CONSOLIDAÇÃO E GRÁFICOS COMPARATIVOS
# ==============================================================================

resumo <- bind_rows(resultados)
write_csv(resumo, file.path(output_dir, "resumo_comparativo.csv"))
cat("\n=== RESUMO COMPARATIVO: lambda[t] vs lambda[i,t] ===\n")
print(resumo)

make_region_label <- function(names) {
  sapply(names, function(n) {
    n_clean <- gsub("_", " ", n)
    if (nchar(n_clean) > 25) paste0(substr(n_clean, 1, 22), "...") else n_clean
  })
}

read_model_csv <- function(filename) {
  lapply(model_types, function(m) {
    path <- file.path(output_dir, m, filename)
    if (!file.exists(path)) return(NULL)
    read_csv(path, show_col_types = FALSE) %>% mutate(model = m)
  }) |> bind_rows()
}

# ── Comparativo epsilon ───────────────────────────────────────────────────
eps_all <- read_model_csv("epsilon_summary.csv")
if (nrow(eps_all) > 0) {
  ggsave(
    file.path(output_dir, "epsilon_comparativo.png"),
    ggplot(eps_all, aes(x = Time, y = Mean, color = model, fill = model)) +
      geom_ribbon(aes(ymin = Lower, ymax = Upper), alpha = 0.15, color = NA) +
      geom_line(linewidth = 0.9) +
      geom_hline(yintercept = eps_lower, linetype = "dashed", color = "grey50") +
      geom_hline(yintercept = eps_upper, linetype = "dashed", color = "grey50") +
      scale_x_continuous(breaks = 1:n_times, labels = anos_label) +
      theme_bw(base_size = 12) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            legend.position = "bottom") +
      labs(title = "Comparação epsilon[t]: lambda global vs local",
           subtitle = "Cluster 1 | Priori: U(0.95, 1)",
           x = "Ano", y = expression(epsilon[t]),
           color = "Modelo", fill = "Modelo"),
    width = 10, height = 5, dpi = 300
  )
}

# ── Comparativo lambda ────────────────────────────────────────────────────
lambda_all <- read_model_csv("lambda_summary.csv")
if (nrow(lambda_all) > 0) {
  # Para lambda global, replicar para todas as regiões monitoradas
  lambda_global_df <- lambda_all %>% filter(model == "lambda_global")
  lambda_local_df  <- lambda_all %>% filter(model == "lambda_local")
  
  lambda_local_sel <- lambda_local_df %>% filter(Region %in% REGIONS_INTEREST)
  lambda_local_sel$Region_Label <- make_region_label(lambda_local_sel$Region_Name)
  
  ggsave(
    file.path(output_dir, "lambda_comparativo.png"),
    ggplot() +
      # Lambda global (mesma curva em todos os painéis)
      geom_ribbon(data = lambda_global_df %>% 
                    slice(rep(1:n(), length(REGIONS_INTEREST))) %>%
                    mutate(Region_Label = rep(
                      make_region_label(region_names_c1[REGIONS_INTEREST]), 
                      each = n_times)),
                  aes(x = Time, ymin = Lower, ymax = Upper), 
                  fill = "blue", alpha = 0.15) +
      geom_line(data = lambda_global_df %>% 
                  slice(rep(1:n(), length(REGIONS_INTEREST))) %>%
                  mutate(Region_Label = rep(
                    make_region_label(region_names_c1[REGIONS_INTEREST]), 
                    each = n_times)),
                aes(x = Time, y = Mean, color = "Lambda Global"), 
                linewidth = 0.8) +
      # Lambda local
      geom_ribbon(data = lambda_local_sel,
                  aes(x = Time, ymin = Lower, ymax = Upper), 
                  fill = "red", alpha = 0.15) +
      geom_line(data = lambda_local_sel,
                aes(x = Time, y = Mean, color = "Lambda Local"), 
                linewidth = 0.8) +
      facet_wrap(~ Region_Label, scales = "free_y", ncol = 4) +
      scale_color_manual(
        values = c("Lambda Global" = "blue", "Lambda Local" = "red"),
        name = "Modelo"
      ) +
      scale_x_continuous(breaks = 1:n_times, labels = anos_label) +
      theme_bw(base_size = 10) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
            legend.position = "bottom") +
      labs(title = "Comparação lambda: Global vs Local",
           subtitle = paste("Cluster 1 -", N_c1, "microrregiões | Azul: lambda[t] | Vermelho: lambda[i,t]"),
           x = "Ano", y = expression(lambda)),
    width = 16, height = 12, dpi = 300
  )
}

# ── Comparativo mu ────────────────────────────────────────────────────────
mu_all <- read_model_csv("mu_selected.csv")
if (nrow(mu_all) > 0) {
  mu_all$Region_Label <- make_region_label(mu_all$Region_Name)
  mu_all$Modelo <- ifelse(mu_all$model == "lambda_global", 
                          "Lambda Global", "Lambda Local")
  
  ggsave(
    file.path(output_dir, "mu_comparativo.png"),
    ggplot(mu_all, aes(x = Time, y = Mean, color = Modelo, fill = Modelo)) +
      geom_ribbon(aes(ymin = Lower, ymax = Upper), alpha = 0.15, color = NA) +
      geom_line(linewidth = 0.8) +
      facet_wrap(~ Region_Label, scales = "free_y", ncol = 4) +
      scale_color_manual(values = c("Lambda Global" = "blue", "Lambda Local" = "red")) +
      scale_fill_manual(values = c("Lambda Global" = "blue", "Lambda Local" = "red")) +
      scale_x_continuous(breaks = 1:n_times, labels = anos_label) +
      theme_bw(base_size = 10) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
            legend.position = "bottom") +
      labs(title = "Comparação mu[i,t]: Lambda Global vs Local",
           subtitle = paste("Cluster 1 -", N_c1, "microrregiões"),
           x = "Ano", y = expression(mu[it])),
    width = 16, height = 12, dpi = 300
  )
}

# ── Comparativo beta ──────────────────────────────────────────────────────
beta_all <- read_model_csv("beta_summary.csv")
if (nrow(beta_all) > 0) {
  beta_all$Modelo <- ifelse(beta_all$model == "lambda_global", 
                            "Lambda Global", "Lambda Local")
  
  ggsave(
    file.path(output_dir, "beta_comparativo.png"),
    ggplot(beta_all, aes(x = Parameter, y = Mean, color = Modelo)) +
      geom_point(position = position_dodge(width = 0.5), size = 3) +
      geom_errorbar(aes(ymin = HPD_Lower, ymax = HPD_Upper),
                    position = position_dodge(width = 0.5), width = 0.2) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
      scale_color_manual(values = c("Lambda Global" = "blue", "Lambda Local" = "red")) +
      theme_bw(base_size = 12) +
      theme(legend.position = "bottom") +
      labs(title = "Comparação beta: Lambda Global vs Local",
           subtitle = paste("Cluster 1 -", N_c1, "microrregiões"),
           x = "Covariável", y = expression(beta)),
    width = 8, height = 5, dpi = 300
  )
}

# ── Comparativo efeitos espaciais ──────────────────────────────────────────
s_all <- read_model_csv("s_summary.csv")
if (nrow(s_all) > 0) {
  s_all$Modelo <- ifelse(s_all$model == "lambda_global", 
                         "Lambda Global", "Lambda Local")
  
  ggsave(
    file.path(output_dir, "s_comparativo.png"),
    ggplot(s_all, aes(x = Mean_lambda_global, y = Mean_lambda_local)) +
      geom_point() +
      geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
      theme_bw() +
      labs(title = "Efeitos espaciais: Lambda Global vs Local",
           x = "s[i] - Lambda Global", y = "s[i] - Lambda Local"),
    width = 7, height = 7
  )
  
  # Versão wide para scatter plot
  s_wide <- s_all %>%
    select(Region, Region_Name, model, Mean) %>%
    tidyr::pivot_wider(names_from = model, values_from = Mean, 
                       names_prefix = "Mean_")
  
  if (all(c("Mean_lambda_global", "Mean_lambda_local") %in% colnames(s_wide))) {
    ggsave(
      file.path(output_dir, "s_comparativo.png"),
      ggplot(s_wide, aes(x = Mean_lambda_global, y = Mean_lambda_local)) +
        geom_point(size = 2, alpha = 0.7) +
        geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
        theme_bw(base_size = 12) +
        labs(title = "Efeitos espaciais: Lambda Global vs Local",
             subtitle = "Cada ponto = uma microrregião do Cluster 1",
             x = "s[i] - Lambda Global", y = "s[i] - Lambda Local"),
      width = 7, height = 7, dpi = 300
    )
  }
}

cat("\n========================================\n")
cat("Tempo total de execução:\n")
print(Sys.time() - inicio_global)
cat("\nResultados salvos em:", output_dir, "\n")
cat("========================================\n")