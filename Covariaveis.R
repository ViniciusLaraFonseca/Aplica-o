# ==============================================================================
# covariaveis_robusto_FINAL.R
# ==============================================================================

rm(list = ls())
setwd("C:/Users/vlara/OneDrive/Estatistica UFMG/Mestrado/Pesquisa/PesquisaMestrado")
# ==============================================================================
# 0. PACOTES
# ==============================================================================
pkgs <- c("dplyr", "tidyr", "readr", "stringr", "tibble")
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

cod_micro <- read.table("DATASUS/dados_codMICRO.txt", header = TRUE) %>%
  dplyr::rename(COD_Municipio_Reduzido = COD_MUN) %>%
  dplyr::mutate(
    COD_Municipio_Reduzido = as.character(COD_Municipio_Reduzido),
    MICRO_ = as.character(MICRO_)
  )

ordem_micro <- cod_micro %>%
  dplyr::distinct(MICRO_) %>%
  dplyr::arrange(MICRO_) %>%
  dplyr::pull(MICRO_)

# ==============================================================================
# 3. OFFSET (TOTAL NASCIDOS VIVOS)
# ==============================================================================

total_raw <- read.table("DATASUS/nascidos_vivos_total.csv", sep = ";", header = TRUE,
                        fileEncoding = "Latin1",
                        check.names = FALSE,
                        stringsAsFactors = FALSE)

total <- total_raw[, -25] %>%  # remove coluna "Total"
  dplyr::filter(!grepl("IGNORADO", Município, ignore.case = TRUE)) %>%
  extrair_codigo() %>%
  safe_numeric() %>%
  tidyr::pivot_longer(
    cols = -COD_Municipio_Reduzido,
    names_to = "Ano",
    values_to = "Nascidos_vivos"
  )

check_unique_key(total, c("COD_Municipio_Reduzido", "Ano"), "total")

# ==============================================================================
# 4. FUNÇÃO DE AGREGAÇÃO
# ==============================================================================

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
    dplyr::slice(match(ordem_micro, MICRO_))
}

# ==============================================================================
# 5. FUNÇÃO PADRÃO PARA COVARIÁVEIS
# ==============================================================================

process_cov <- function(path) {
  
  df_raw <- df_raw <- read.table(
    path,
    sep = ";",
    header = TRUE,
    fileEncoding = "Latin1",
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  
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
  
  df <- df %>%
    dplyr::inner_join(total, by = c("COD_Municipio_Reduzido", "Ano"))
  
  return(df)
}

# ==============================================================================
# 6. COVARIÁVEIS
# ==============================================================================

consultas_long   <- process_cov("DATASUS/nascidos_vivos_consultas.csv")
instrucao_long   <- process_cov("DATASUS/nascidos_vivos_instrução2.csv")
subnutridos_long <- process_cov("DATASUS/nascidos_vivos_menos_2500g.csv")

grouped_consultas_MG   <- agg_micro(consultas_long, "valor")
grouped_instrucao_MG   <- agg_micro(instrucao_long, "valor")
grouped_subnutridos_MG <- agg_micro(subnutridos_long, "valor")

# ==============================================================================
# 7. ARRAY X
# ==============================================================================

anos <- colnames(grouped_consultas_MG)[-1]

x <- array(
  NA,
  dim      = c(75, 23, 3),
  dimnames = list(ordem_micro, anos, c("prenatal", "instrucao", "baixo_peso"))
)

x[,,1] <- as.matrix(grouped_consultas_MG[, -1])
x[,,2] <- as.matrix(grouped_instrucao_MG[, -1])
x[,,3] <- as.matrix(grouped_subnutridos_MG[, -1])

# ==============================================================================
# 8. OFFSET E
# ==============================================================================

Grouped_total_MG <- total %>%
  dplyr::inner_join(cod_micro, by = "COD_Municipio_Reduzido") %>%
  dplyr::group_by(MICRO_, Ano) %>%
  dplyr::summarise(SumNascidos = sum(Nascidos_vivos), .groups = "drop") %>%
  tidyr::pivot_wider(names_from = Ano, values_from = SumNascidos) %>%
  dplyr::mutate(MICRO_ = as.character(MICRO_)) %>%
  dplyr::slice(match(ordem_micro, MICRO_))

E <- as.matrix(Grouped_total_MG[, -1])
rownames(E) <- ordem_micro

# ==============================================================================
# 9. Y (ÓBITOS)
# ==============================================================================

obitos_raw <- read.table(
  "DATASUS/morte_neonatal_precoce.csv",
  sep = ";",
  header = TRUE,
  fileEncoding = "Latin1",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

anos_cols <- as.character(2000:2022)

obitos_clean <- obitos_raw %>%
  dplyr::filter(
    !grepl("^\\s*Total", Município, ignore.case = TRUE),
    !grepl("IGNORADO", Município, ignore.case = TRUE)
  ) %>%
  extrair_codigo() %>%
  dplyr::select(COD_Municipio_Reduzido, dplyr::all_of(anos_cols)) %>%
  safe_numeric() %>%
  tidyr::pivot_longer(
    cols = -COD_Municipio_Reduzido,
    names_to = "Ano",
    values_to = "Obitos"
  )

check_unique_key(obitos_clean, c("COD_Municipio_Reduzido", "Ano"), "obitos")

Grouped_Y <- obitos_clean %>%
  dplyr::inner_join(cod_micro, by = "COD_Municipio_Reduzido") %>%
  dplyr::group_by(MICRO_, Ano) %>%
  dplyr::summarise(Y = sum(Obitos, na.rm = TRUE), .groups = "drop") %>%
  tidyr::pivot_wider(names_from = Ano, values_from = Y) %>%
  dplyr::mutate(MICRO_ = as.character(MICRO_)) %>%
  dplyr::slice(match(ordem_micro, MICRO_))

Y_mat <- as.matrix(Grouped_Y[, -1])
rownames(Y_mat) <- ordem_micro

# ==============================================================================
# 10. CHECAGEM FINAL
# ==============================================================================

if (!all(dim(Y_mat) == c(75, 23))) stop("❌ Dimensão de Y incorreta")
if (!all(dim(E) == c(75, 23))) stop("❌ Dimensão de E incorreta")
if (!all(dim(x) == c(75, 23, 3))) stop("❌ Dimensão de x incorreta")

cat("✅ Dados prontos e consistentes para o modelo.\n")

setwd("C:/Users/vlara/OneDrive/Estatistica UFMG/Mestrado/Pesquisa/Aplicação/main")
