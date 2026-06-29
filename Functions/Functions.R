library(here)


#define a directory to upload the required data and save the figures
dir.create(here("SOD1_ddPCA_manuscript"), showWarnings = FALSE)


#FUNCTION: convert hgvs_pro to ID format####
aa_dict <- c(
  Ala="A", Arg="R", Asn="N", Asp="D", Cys="C",
  Glu="E", Gln="Q", Gly="G", His="H", Ile="I",
  Leu="L", Lys="K", Met="M", Phe="F", Pro="P",
  Ser="S", Thr="T", Trp="W", Tyr="Y", Val="V",
  Ter="*", X="*"  # Stop codons
)

aa_short <- c(
  A="A", R="R", N="N", D="D", C="C", E="E", Q="Q", G="G",
  H="H", I="I", L="L", K="K", M="M", F="F", P="P",
  S="S", T="T", W="W", Y="Y", V="V", "*"="*"
)

convert_hgvsp_to_id <- function(variant) {
  
  
  variant <- gsub("^p\\.", "", variant)
  
  # -------------------
  # 1) Deletions
  # -------------------
  if (grepl("del$", variant)) {
    matches <- regmatches(variant, gregexpr("[A-Za-z]+|[0-9]+", variant))[[1]]
    if (length(matches) >= 2) {
      ref <- matches[1]
      pos <- matches[2]
      ref <- ifelse(ref %in% names(aa_dict), aa_dict[ref], aa_short[ref])
      return(paste0(ref, pos, "d"))
    }
  }
  
  # -------------------
  # 2) Stop codon
  # -------------------
  if (grepl("X$", variant) || grepl("Ter$", variant)) {
    matches <- regmatches(variant, gregexpr("[A-Za-z]+|[0-9]+", variant))[[1]]
    if (length(matches) >= 2) {
      ref <- matches[1]
      pos <- matches[2]
      ref <- ifelse(ref %in% names(aa_dict), aa_dict[ref], aa_short[ref])
      return(paste0(ref, pos, "*"))
    }
  }
  
  # -------------------
  # 3) Synonymous
  # -------------------
  if (grepl("=$", variant)) {
    matches <- regmatches(variant, gregexpr("[A-Za-z]+|[0-9]+", variant))[[1]]
    if (length(matches) >= 2) {
      ref <- matches[1]
      pos <- matches[2]
      ref <- ifelse(ref %in% names(aa_dict), aa_dict[ref], aa_short[ref])
      return(paste0(ref, pos, ref))
    }
  }
  
  # -------------------
  # 4) Substitution
  # -------------------
  matches <- regmatches(variant, gregexpr("[A-Za-z]+|[0-9]+", variant))[[1]]
  
  if (length(matches) == 3) {
    ref <- matches[1]
    pos <- matches[2]
    alt <- matches[3]
    
    ref <- ifelse(ref %in% names(aa_dict), aa_dict[ref], aa_short[ref])
    alt <- ifelse(alt %in% names(aa_dict), aa_dict[alt], aa_short[alt])
    
    return(paste0(ref, pos, alt))
  }
  
  # -------------------
  # Unkown
  # -------------------
  warning(paste("Unknown format:", variant))
  return(NA)
}

#FUNCTION: to convert ID to hgvs_pro format####
id_to_hgvs_pro <- function(ID) {
  
  aa_map <- c(
    A = "Ala", R = "Arg", N = "Asn", D = "Asp", C = "Cys",
    Q = "Gln", E = "Glu", G = "Gly", H = "His", I = "Ile",
    L = "Leu", K = "Lys", M = "Met", F = "Phe", P = "Pro",
    S = "Ser", T = "Thr", W = "Trp", Y = "Tyr", V = "Val"
  )
  
  aa_wt  <- substr(ID, 1, 1)
  aa_mut <- substr(ID, nchar(ID), nchar(ID))
  pos    <- as.numeric(gsub("[A-Z]", "", ID))
  
  paste0(
    "p.",
    aa_map[aa_wt],
    pos,
    aa_map[aa_mut]
  )
}


#FUNCTION: substrate 1 position to ID####
decrement_id_pos <- function(id) {
  matches <- regmatches(id, regexec("^([A-Z])([0-9]+)([A-Z\\*=d])$", id))[[1]]
  
  if (length(matches) == 4) {
    aa <- matches[2]         
    pos <- as.numeric(matches[3]) - 1  
    alt <- matches[4]        
    return(paste0(aa, pos, alt))
  } else {
    warning(paste("Formato desconocido:", id))
    return(NA)
  }
}

#FUNCTION: sum 1 position to ID####
sum_id_pos <- function(id) {
  matches <- regmatches(id, regexec("^([A-Z])([0-9]+)([A-Z\\*=d])$", id))[[1]]
  
  if (length(matches) == 4) {
    aa <- matches[2]         
    pos <- as.numeric(matches[3]) + 1  
    alt <- matches[4]        
    return(paste0(aa, pos, alt))
  } else {
    warning(paste("Formato desconocido:", id))
    return(NA)
  }
}
