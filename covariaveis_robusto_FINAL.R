# ==============================================================================
# covariaveis_robusto_FINAL.R
# Versão Final com Matching Topológico Robusto
# ==============================================================================

rm(list = ls())
setwd("C:/Users/vlara/OneDrive/Estatistica UFMG/Mestrado/Pesquisa/PesquisaMestrado")

# ==============================================================================
# 0. PACOTES
# ==============================================================================
pkgs <- c("dplyr", "tidyr", "readr", "stringr", "tibble", "sf", "spdep")
for (p in pkgs) {
  if (!require(p, character.only = TRUE)) install.packages(p)
  library(p, character.only = TRUE)
}

# ==============================================================================
# 1. FUNÇÕES AUXILIARES
# ==============================================================================

check_unique_key <- function(df, keys, nome_df) {
  dup <- df %>%
    dplyr::count(dplyr::across(dplyr::all_of(keys))) %>%
    dplyr::filter(n > 1)
  
  if (nrow(dup) > 0) {
    print(dup)
    stop(paste0("❌ Duplicatas encontradas em ", nome_df))
  }
}

safe_numeric <- function(df) {
  df %>%
    dplyr::mutate(dplyr::across(
      -COD_Municipio_Reduzido,
      ~ readr::parse_number(
        dplyr::na_if(as.character(.), "-")
      )
    ))
}

extrair_codigo <- function(df) {
  df %>%
    dplyr::mutate(
      COD_Municipio_Reduzido = stringr::str_extract(trimws(Município), "^\\d+")
    ) %>%
    dplyr::filter(!is.na(COD_Municipio_Reduzido)) %>%
    dplyr::select(-Município)
}

# ==============================================================================
# 2. BASE DE MICRORREGIÕES
# ==============================================================================

cat("\n📊 Carregando dados das microrregiões...\n")

cod_micro <- read.table("DATASUS/dados_codMICRO.txt", header = TRUE) %>%
  dplyr::rename(COD_Municipio_Reduzido = COD_MUN) %>%
  dplyr::mutate(
    COD_Municipio_Reduzido = as.character(COD_Municipio_Reduzido),
    MICRO_ = as.character(MICRO_)
  )

cat("✅ Microrregiões carregadas:", length(unique(cod_micro$MICRO_)), "\n")

# ==============================================================================
# 3. MATCHING TOPOLÓGICO
# ==============================================================================

cat("\n🔍 Iniciando matching topológico...\n")

# 3.1 Dados do modelo
adj_modelo <- list(
  c(2,3,5,12,17,18), c(1,4,5,6,17), c(1,5,7,12,19), c(2,6,17,31),
  c(1,2,3,6,7), c(2,4,5,7), c(3,5,6), c(9,12,15,19,21,30,58),
  c(8,15,19), c(11,18,27,28,39,58,74), c(10,17,18,20,28),
  c(1,3,8,18,19,58), c(14,17,20), c(13,20,23,42), c(8,9,21,36),
  c(22,26,45,68), c(1,2,4,11,13,18,20,31), c(1,10,11,12,17,58),
  c(3,8,9,12), c(11,13,14,17,23,28,67,73), c(8,15,27,30,36,58),
  c(16,26,42,45), c(14,20,32,35,42,67,73), c(57,64,71),
  c(36,41,49,51,63,66), c(16,22,33,45,47,62,65,68,70),
  c(10,21,36,39,58), c(10,11,20,39,55,67,74), c(41,48),
  c(8,21,58), c(4,17), c(23,35,55,67), c(26,44,65,68,70),
  c(40,54), c(23,32,42,45,47,55,56), c(15,21,25,27,37,39,51,66),
  c(36,51,71), c(40,43,52,54,72), c(10,27,28,36,53,60,61,66,74),
  c(34,38,46,52,54), c(25,29,48,49,51,54,64,72),
  c(14,22,23,35,45), c(38,52,59,72), c(33,68),
  c(16,22,26,35,42,47), c(40,52),
  c(26,35,45,49,50,56,61,62,74), c(29,41,54),
  c(25,41,47,50,59,61,63,69,72), c(47,49,52,59,62,70),
  c(25,36,37,41,64,71), c(38,40,43,46,50,59,70),
  c(39,61,63,66,69), c(34,38,40,41,48,72),
  c(28,32,35,56,67,74), c(35,47,55,61,74), c(24,71),
  c(8,10,12,18,21,27,30), c(43,49,50,52,72), c(39,61,74),
  c(39,47,49,53,56,60,69,74,75), c(26,47,50,65,70),
  c(25,49,53,66,69), c(24,41,51,71), c(26,33,62,70),
  c(25,36,39,53,63), c(20,23,28,32,55,73), c(16,26,33,44),
  c(49,53,61,63), c(26,33,50,52,62,65), c(24,37,51,57,64),
  c(38,41,43,49,54,59), c(20,23,67),
  c(10,28,39,47,55,56,61,60,75), c(61,74)
)

# 3.2 Tentar obter shapefile dos municípios
matriz_file <- "DATASUS/matriz_adj_municipal_MG.rda"

if(file.exists(matriz_file)) {
  cat("📂 Carregando matriz de adjacência salva...\n")
  load(matriz_file)
} else {
  cat("📥 Tentando obter shapefile dos municípios...\n")
  
  # Tentar múltiplos métodos
  mun_sf <- NULL
  
  # Método 1: geobr
  tryCatch({
    cat("  Tentando via geobr...\n")
    options(timeout = 300)
    mun_sf <- geobr::read_municipality(code_muni = "MG", year = 2010)
    if(inherits(mun_sf, "sf")) {
      cat("  ✅ Sucesso via geobr!\n")
    }
  }, error = function(e) {
    cat("  ❌ geobr falhou:", e$message, "\n")
  })
  
  # Método 2: Download direto do IBGE
  if(is.null(mun_sf)) {
    tryCatch({
      cat("  Tentando download direto do IBGE...\n")
      temp_dir <- tempdir()
      temp_zip <- tempfile(fileext = ".zip")
      
      # URL alternativa (mais estável)
      url <- "https://geoftp.ibge.gov.br/organizacao_do_territorio/malhas_territoriais/malhas_municipais/municipio_2015/MG/mg_municipios.zip"
      
      download.file(url, temp_zip, mode = "wb", timeout = 300)
      unzip(temp_zip, exdir = temp_dir)
      
      # Procurar pelo arquivo .shp
      shp_file <- list.files(temp_dir, pattern = "\\.shp$", full.names = TRUE, recursive = TRUE)
      if(length(shp_file) > 0) {
        mun_sf <- sf::read_sf(shp_file[1])
        cat("  ✅ Sucesso via download IBGE!\n")
      }
    }, error = function(e) {
      cat("  ❌ Download IBGE falhou:", e$message, "\n")
    })
  }
  
  # Método 3: Tentar ano diferente
  if(is.null(mun_sf)) {
    tryCatch({
      cat("  Tentando ano alternativo...\n")
      mun_sf <- geobr::read_municipality(code_muni = "MG", year = 2000)
      if(inherits(mun_sf, "sf")) {
        cat("  ✅ Sucesso com ano 2000!\n")
      }
    }, error = function(e) {
      cat("  ❌ Ano alternativo falhou:", e$message, "\n")
    })
  }
  
  if(!is.null(mun_sf) && inherits(mun_sf, "sf")) {
    cat("🔗 Gerando matriz de adjacência municipal...\n")
    
    # Garantir que é sf
    mun_sf <- sf::st_as_sf(mun_sf)
    
    # Extrair código do município (6 dígitos)
    if("code_muni" %in% names(mun_sf)) {
      mun_sf$cod_6dig <- substr(mun_sf$code_muni, 1, 6)
    } else if("CD_GEOCODM" %in% names(mun_sf)) {
      mun_sf$cod_6dig <- substr(mun_sf$CD_GEOCODM, 1, 6)
    } else {
      # Tentar encontrar coluna com código
      code_col <- grep("cod|CD_GEO", names(mun_sf), value = TRUE, ignore.case = TRUE)[1]
      if(!is.na(code_col)) {
        mun_sf$cod_6dig <- substr(mun_sf[[code_col]], 1, 6)
      }
    }
    
    # Criar adjacência
    vizinhos <- spdep::poly2nb(sf::st_geometry(mun_sf), queen = TRUE)
    
    # Salvar
    save(vizinhos, mun_sf, file = matriz_file)
    cat("💾 Matriz salva em", matriz_file, "\n")
  }
}

# 3.3 Construir adjacência entre microrregiões
if(exists("vizinhos") && exists("mun_sf")) {
  cat("🏗️  Construindo adjacência entre microrregiões...\n")
  
  municipios_por_micro <- split(cod_micro$COD_Municipio_Reduzido, cod_micro$MICRO_)
  micro_codigos <- names(municipios_por_micro)
  n_micro <- length(micro_codigos)
  
  # Verificar matching de códigos
  codigos_municipais <- mun_sf$cod_6dig
  seus_codigos <- unique(cod_micro$COD_Municipio_Reduzido)
  
  cat("  Códigos no shapefile:", length(codigos_municipais), "\n")
  cat("  Seus códigos:", length(seus_codigos), "\n")
  
  # Verificar interseção
  match_count <- sum(seus_codigos %in% codigos_municipais)
  cat("  Municípios com match:", match_count, "de", length(seus_codigos), "\n")
  
  if(match_count > 0) {
    # Construir adjacência
    adj_micro_r <- vector("list", n_micro)
    names(adj_micro_r) <- micro_codigos
    
    for(i in 1:n_micro) {
      micro_i <- micro_codigos[i]
      muns_i <- municipios_por_micro[[micro_i]]
      
      vizinhos_micro <- character(0)
      
      for(j in 1:n_micro) {
        if(i == j) next
        
        micro_j <- micro_codigos[j]
        muns_j <- municipios_por_micro[[micro_j]]
        
        # Verificar vizinhança
        viz_encontrado <- FALSE
        for(mun_i in muns_i) {
          idx_i <- which(codigos_municipais == mun_i)
          if(length(idx_i) == 1 && idx_i <= length(vizinhos)) {
            viz_i <- codigos_municipais[vizinhos[[idx_i]]]
            if(any(viz_i %in% muns_j)) {
              vizinhos_micro <- c(vizinhos_micro, micro_j)
              viz_encontrado <- TRUE
              break
            }
          }
        }
      }
      
      adj_micro_r[[micro_i]] <- vizinhos_micro
    }
    
    # Criar assinaturas
    grau_r <- lengths(adj_micro_r)
    grau_modelo <- lengths(adj_modelo)
    
    assinatura_r <- data.frame(
      MICRO_ = micro_codigos,
      grau = grau_r,
      soma_graus = sapply(adj_micro_r, function(v) {
        if(length(v) > 0) sum(grau_r[match(v, micro_codigos)]) else 0
      })
    )
    
    assinatura_modelo <- data.frame(
      ordem_modelo = 1:75,
      grau = grau_modelo,
      soma_graus = sapply(adj_modelo, function(v) {
        if(length(v) > 0) sum(grau_modelo[v]) else 0
      })
    )
    
    # Matching
    mapeamento <- merge(assinatura_modelo, assinatura_r,
                        by = c("grau", "soma_graus"))
    
    cat("\n📊 Resultado do matching:", nrow(mapeamento), "de 75 microrregiões\n")
    
    if(nrow(mapeamento) == 75) {
      cat("✅ MATCHING PERFEITO!\n")
      ordem_modelo <- mapeamento$MICRO_[order(mapeamento$ordem_modelo)]
    } else {
      cat("⚠️  Matching parcial. Usando ordem original como fallback.\n")
      ordem_modelo <- unique(cod_micro$MICRO_)
    }
  } else {
    cat("⚠️  Nenhum match de códigos. Usando ordem original.\n")
    ordem_modelo <- unique(cod_micro$MICRO_)
  }
} else {
  cat("⚠️  Shapefile não disponível. Usando ordem original.\n")
  ordem_modelo <- unique(cod_micro$MICRO_)
}

cat("\n✅ Ordem definida! Primeiras 5 microrregiões:",
    paste(ordem_modelo[1:5], collapse = ", "), "\n\n")

# ==============================================================================
# 4. PROCESSAMENTO DOS DADOS
# ==============================================================================

cat("📊 Processando dados do DATASUS...\n")

# 4.1 Total de nascidos vivos
total_raw <- read.table("DATASUS/nascidos_vivos_total.csv", sep = ";", header = TRUE,
                        fileEncoding = "Latin1", check.names = FALSE, stringsAsFactors = FALSE)

total <- total_raw[, -25] %>%
  dplyr::filter(!grepl("IGNORADO", Município, ignore.case = TRUE)) %>%
  extrair_codigo() %>%
  safe_numeric() %>%
  tidyr::pivot_longer(
    cols = -COD_Municipio_Reduzido,
    names_to = "Ano",
    values_to = "Nascidos_vivos"
  )

check_unique_key(total, c("COD_Municipio_Reduzido", "Ano"), "total")

# 4.2 Funções de agregação
agg_micro <- function(df, var_name) {
  check_unique_key(df, c("COD_Municipio_Reduzido", "Ano"), "covariavel")
  
  df %>%
    dplyr::inner_join(cod_micro, by = "COD_Municipio_Reduzido") %>%
    dplyr::group_by(MICRO_, Ano) %>%
    dplyr::summarise(
      valor = sum(.data[[var_name]], na.rm = TRUE),
      total = sum(Nascidos_vivos, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::mutate(taxa = valor / total) %>%
    dplyr::select(MICRO_, Ano, taxa) %>%
    tidyr::pivot_wider(names_from = Ano, values_from = taxa) %>%
    dplyr::mutate(MICRO_ = as.character(MICRO_)) %>%
    dplyr::slice(match(ordem_modelo, MICRO_))
}

process_cov <- function(path) {
  df_raw <- read.table(path, sep = ";", header = TRUE, fileEncoding = "Latin1",
                       check.names = FALSE, stringsAsFactors = FALSE)
  
  df <- df_raw[, -25] %>%
    dplyr::filter(!grepl("IGNORADO", Município, ignore.case = TRUE)) %>%
    extrair_codigo() %>%
    safe_numeric() %>%
    tidyr::pivot_longer(
      cols = -COD_Municipio_Reduzido,
      names_to = "Ano",
      values_to = "valor"
    )
  
  check_unique_key(df, c("COD_Municipio_Reduzido", "Ano"), "covariavel_pre_join")
  
  df <- df %>% dplyr::inner_join(total, by = c("COD_Municipio_Reduzido", "Ano"))
  return(df)
}

# 4.3 Processar covariáveis
cat("  Processando consultas pré-natal...\n")
consultas_long <- process_cov("DATASUS/nascidos_vivos_consultas.csv")

cat("  Processando instrução...\n")
instrucao_long <- process_cov("DATASUS/nascidos_vivos_instrução2.csv")

cat("  Processando baixo peso...\n")
subnutridos_long <- process_cov("DATASUS/nascidos_vivos_menos_2500g.csv")

grouped_consultas_MG <- agg_micro(consultas_long, "valor")
grouped_instrucao_MG <- agg_micro(instrucao_long, "valor")
grouped_subnutridos_MG <- agg_micro(subnutridos_long, "valor")

# ==============================================================================
# 5. CONSTRUÇÃO DAS MATRIZES FINAIS
# ==============================================================================

cat("\n🏗️  Construindo matrizes finais...\n")

anos <- colnames(grouped_consultas_MG)[-1]

# Array X
x <- array(NA, dim = c(75, 23, 3),
           dimnames = list(ordem_modelo, anos, c("prenatal", "instrucao", "baixo_peso")))
x[,,1] <- as.matrix(grouped_consultas_MG[, -1])
x[,,2] <- as.matrix(grouped_instrucao_MG[, -1])
x[,,3] <- as.matrix(grouped_subnutridos_MG[, -1])

# Matriz E (offset)
Grouped_total_MG <- total %>%
  dplyr::inner_join(cod_micro, by = "COD_Municipio_Reduzido") %>%
  dplyr::group_by(MICRO_, Ano) %>%
  dplyr::summarise(SumNascidos = sum(Nascidos_vivos), .groups = "drop") %>%
  tidyr::pivot_wider(names_from = Ano, values_from = SumNascidos) %>%
  dplyr::mutate(MICRO_ = as.character(MICRO_)) %>%
  dplyr::slice(match(ordem_modelo, MICRO_))

E <- as.matrix(Grouped_total_MG[, -1])
rownames(E) <- ordem_modelo

# Matriz Y (óbitos)
obitos_raw <- read.table("DATASUS/morte_neonatal_precoce.csv", sep = ";", header = TRUE,
                         fileEncoding = "Latin1", check.names = FALSE, stringsAsFactors = FALSE)

anos_cols <- as.character(2000:2022)

obitos_clean <- obitos_raw %>%
  dplyr::filter(!grepl("^\\s*Total", Município, ignore.case = TRUE),
                !grepl("IGNORADO", Município, ignore.case = TRUE)) %>%
  extrair_codigo() %>%
  dplyr::select(COD_Municipio_Reduzido, dplyr::all_of(anos_cols)) %>%
  safe_numeric() %>%
  tidyr::pivot_longer(cols = -COD_Municipio_Reduzido,
                      names_to = "Ano", values_to = "Obitos")

check_unique_key(obitos_clean, c("COD_Municipio_Reduzido", "Ano"), "obitos")

Grouped_Y <- obitos_clean %>%
  dplyr::inner_join(cod_micro, by = "COD_Municipio_Reduzido") %>%
  dplyr::group_by(MICRO_, Ano) %>%
  dplyr::summarise(Y = sum(Obitos, na.rm = TRUE), .groups = "drop") %>%
  tidyr::pivot_wider(names_from = Ano, values_from = Y) %>%
  dplyr::mutate(MICRO_ = as.character(MICRO_)) %>%
  dplyr::slice(match(ordem_modelo, MICRO_))

Y_mat <- as.matrix(Grouped_Y[, -1])
rownames(Y_mat) <- ordem_modelo

# ==============================================================================
# 6. INFORMAÇÕES DE CLUSTER
# ==============================================================================

cat("\n📊 Adicionando informações de cluster...\n")

clAI <- c(4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,3,3,3,3,3,3,3,3,3,3,3,3,3,3,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1)

AI <- c(-83.84,-75.69,-72.98,-40.01,-39.93,-31.28,-27.47,-23.89,3.71,5.04,5.26,5.72,8.03,9.01,13.17,13.25,18.89,24.08,27.06,28.38,33.62,35.44,38.70,39.14,39.92,44.24,44.49,48.59,50.00,50.00,50.00,50.21,54.52,55.80,56.14,57.16,57.61,59.11,56.86,60.07,60.88,60.95,61.74,63.74,65.77,65.86,69.77,71.45,71.95,72.95,73.83,75.64,76.61,77.17,77.83,78.04,78.17,78.19,78.30,79.86,82.15,82.79,84.47,86.83,86.94,88.11,88.23,89.65,89.75,90.69,91.94,92.20,92.73,98.61,100.00)

cluster_info <- data.frame(
  MICRO_ = ordem_modelo,
  ordem_modelo = 1:75,
  AI = AI,
  cluster = clAI,
  cluster_nome = dplyr::case_when(
    clAI == 1 ~ "Alta adequação",
    clAI == 2 ~ "Média-alta adequação",
    clAI == 3 ~ "Média-baixa adequação",
    clAI == 4 ~ "Baixa adequação"
  )
)

cat("Distribuição por cluster:\n")
print(table(cluster_info$cluster_nome))

# ==============================================================================
# 7. VERIFICAÇÕES FINAIS
# ==============================================================================

cat("\n🔍 Verificações finais...\n")

stopifnot(all(dim(Y_mat) == c(75, 23)))
stopifnot(all(dim(E) == c(75, 23)))
stopifnot(all(dim(x) == c(75, 23, 3)))

stopifnot(identical(rownames(Y_mat), ordem_modelo))
stopifnot(identical(rownames(E), ordem_modelo))
stopifnot(identical(dimnames(x)[[1]], ordem_modelo))

if(anyNA(Y_mat)) warning("⚠️  NAs em Y_mat")
if(anyNA(E)) warning("⚠️  NAs em E")
if(anyNA(x)) warning("⚠️  NAs em x")

cat("\n✅ DADOS PRONTOS PARA O MODELO!\n")
cat("📁 Dimensões finais:\n")
cat("  Y:", paste(dim(Y_mat), collapse = " × "), "\n")
cat("  E:", paste(dim(E), collapse = " × "), "\n")
cat("  x:", paste(dim(x), collapse = " × "), "\n")

# ==============================================================================
# 8. SALVAR
# ==============================================================================

setwd("C:/Users/vlara/OneDrive/Estatistica UFMG/Mestrado/Pesquisa/Aplicação/main")

save(Y_mat, E, x, ordem_modelo, cluster_info, 
     file = "dados_modelo_final.RData")

cat("\n💾 Dados salvos em 'dados_modelo_final.RData'\n")
cat("✅ Pronto para usar no modelo!\n")