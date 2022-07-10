# create parameter table

# ----
alpha = c(0.1, 0.9) # values from which to calculate the width of the signal
compute_duration_assembly <- function(alpha, simul, ipath){
  if (paste0('ini1_p', ipath, '_r12') %in% colnames(simul)){
    timeval <- simul[ , paste0('kin2_p', ipath, '_r12')] 
  }
  else{
    timeval <- simul[ , paste0('kin2_p', ipath, '_r12')]
  }
  nH <- simul[ , paste0('kin1_p', ipath, '_r12')] 
  dT <- timeval*(alpha[2]/(1-alpha[2]))^(1/nH) - timeval*(alpha[1]/(1-alpha[1]))^(1/nH)
  return(dT)
}
# ------
current_path <- getActiveDocumentContext()$path
pathvar.workdir <- dirname(current_path)
setwd(pathvar.workdir)
pathvar.CIparfile <- '../../results/singlepore/20210417_bgsub_nuclei_dynamic_usestd1/20210417_single_npc_logn_hill_CI__par.txt'
pathvar.parfile <- '../../results/singlepore/20210417_bgsub_nuclei_dynamic_usestd1/20210417_single_npc_logn_hill_2min__par.txt'
pathvar.outfile_xlsx <- '../../results/singlepore/20210417_bgsub_nuclei_dynamic_usestd1/20210417_single_npc_logn_hill_CI__par.xlsx'
pathvar.outfile <- '../../results/singlepore/20210417_bgsub_nuclei_dynamic_usestd1/20210417_single_npc_logn_hill_CI__par_complete.txt'

CItab <- read.csv(pathvar.CIparfile, sep = '\t')
partab <- read.csv(pathvar.parfile, sep = '\t')

# ----
outframe <- data.frame()
par_names <- c('n_p', 'K_p_min', 'n_i', 'K_i_(min)', 'f_n', 'f_c', 'dT_p_min', 'dT_i_min')
df <- data.frame('POI', 'n_p', 'K_p_min', 'n_i', 'K_i_min', 'f_n', 'f_c', 'dT_p_min', 'dT_i_min')
for (poiv in unique(CItab$POI)){
  poidf_CI <- subset(CItab, POI == poiv)
  poidf <- subset(partab, POI == poiv)

  df <- data.frame('POI' = poiv)
  par_CI_indexes <- unique(poidf_CI$par_varied)
  for (ipar in seq_along(par_CI_indexes)){

    idx <- par_CI_indexes[ipar]
    poidf_CI_loc <- subset(poidf_CI, POI == poiv & par_varied == idx)
    col_name <- colnames(partab)[idx+1]
    par_mean <- poidf[idx+1]
    par_low <- poidf_CI_loc[1, idx+1]
    par_high <- poidf_CI_loc[2, idx+1]
    partab[partab$POI == poiv, paste0(col_name, "_lowCI")] <- par_low
    partab[partab$POI == poiv, paste0(col_name, "_highCI")] <- par_high
    str_out <-paste0(round(par_mean,2), " [", round(par_low,2), ", ", round(par_high,2), "]")
    df[par_names[ipar]] <- str_out
    
  }
  poidf_CI_loc <- subset(poidf_CI, POI == poiv & par_varied == 9)
  dT_p_CI <- c(compute_duration_assembly(alpha, poidf_CI_loc[1,], 1), compute_duration_assembly(alpha, poidf_CI_loc[2,], 1))
  poidf_CI_loc <- subset(poidf_CI, POI == poiv & par_varied == 10)
  dT_p_CI <- c(dT_p_CI, compute_duration_assembly(alpha, poidf_CI_loc[1,], 1), compute_duration_assembly(alpha, poidf_CI_loc[2,], 1))
  
  dT_p <- compute_duration_assembly(alpha, poidf, 1)
  str_out <- paste0(round(dT_p,2), " [", round(min(dT_p_CI),2), ", ", round(max(dT_p_CI),2), "]")
  df[par_names[7]] <- str_out
  poidf_CI_loc <- subset(poidf_CI, POI == poiv & par_varied == 11)
  dT_i_CI <- c(compute_duration_assembly(alpha, poidf_CI_loc[1,], 2), compute_duration_assembly(alpha, poidf_CI_loc[2,], 2))
  poidf_CI_loc <- subset(poidf_CI, POI == poiv & par_varied == 12)
  dT_i_CI <- c(dT_i_CI, compute_duration_assembly(alpha, poidf_CI_loc[1,], 2), compute_duration_assembly(alpha, poidf_CI_loc[2,], 2))
  
  dT_i <- compute_duration_assembly(alpha, poidf, 2)
  str_out <- paste0(round(dT_i,2), " [", round(min(dT_i_CI),2), ", ", round(max(dT_i_CI),2), "]")
  df[par_names[8]] <- str_out
  outframe <- rbind(outframe, df)
  partab[partab$POI == poiv, "dTp"] <- dT_p
  partab[partab$POI == poiv, "dTp_lowCI"] <-  min(dT_p_CI)
  partab[partab$POI == poiv, "dTp_highCI"] <- max(dT_p_CI)
  
  partab[partab$POI == poiv, "dTi"] <- dT_i
  partab[partab$POI == poiv, "dTi_lowCI"] <- min(dT_i_CI)
  partab[partab$POI == poiv, "dTi_highCI"] <- max(dT_i_CI)
}
write.table(x = partab, file = pathvar.outfile, row.names = FALSE, quote = FALSE, sep = '\t')

# ----
tab_test <- read.csv(pathvar.outfile)
xlsx::write.xlsx(x = outframe, file = pathvar.outfile_xlsx, row.names = FALSE)