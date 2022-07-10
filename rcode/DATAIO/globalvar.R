# ----- global variables for processing 
pathvar.datadir <- '../../expdata/livecellCalibrated/'
pathvar.xlswt <- '../../expdata/livecellCalibrated/HeLa4D-core-noncore-normalized-summary-190225_20210417.xlsx'
pathvar.tabwt <- '../../expdata/livecellCalibrated/20210417/HeLa4D-core-noncore-normalized-summary-190225.txt'
pathvar.tabwt_norm_bgsub <- '../../expdata/livecellCalibrated/20210417/HeLa4D-core-noncore-normalized-summary-190225_norm_bgsub.txt'
pathvar.tabwt_norm_bgsub_inter_meta <- '../../expdata/livecellCalibrated/20210417/HeLa4D-core-noncore-normalized-summary-190225_norm_bgsub_inter_meta.txt'
pathvar.tabwt_norm <- '../../expdata/livecellCalibrated/20210417/HeLa4D-core-noncore-normalized-summary-190225_norm.txt'
pathvar.tabwt_norm_bgsub_singlecell <- '../../expdata/livecellCalibrated/20210417/HeLa4D-core-noncore-normalized-summary-190225_norm_bgsub_singlecell.txt'
pathvar.xlsx_shotaro <- '../../expdata/livecellCalibrated/20210417/HeLa4D-core-noncore-normalized-summary-190225_norm_bgsub_inter_meta.xlsx'
fraction<- list('noncore' = 0.857, 'core' = 0.295) # fraction of postmitotic assembly



#' interphase_post_mitotic
#' Compute interphase and postmitotic assembly curves from data core and noncore regions and fraction of interphase in noncore and core
#' @param core an array with 1st column mean, 2nd column sd. Rows are time points. The data is normalized
#' @param noncore an array with 1st column mean, 2nd column sd.  Rows are time points. The data is normalized
#' @param df_id cell and poi id 
#' @param fraction fraction  of postmitotic assembly in noncore and core. Named list. example 
#'                  list('noncore' = 0.857, 'core' = 0.295) # fraction of postmitotic assembly
#'
#' @return A dataframe with separated interphase and postmitotic assembly
#' @export
#'
#' @examples
interphase_post_mitotic <- function(core, noncore, df_id, fraction){
  
  denominator <- fraction[['noncore']] - fraction[['core']]
  interphase <- (fraction[['noncore']]*core[,1] - fraction[['core']]*noncore[,1])/denominator
  
  postmitotic <- ((1-fraction[['core']])*noncore[,1]- (1-fraction[['noncore']])*core[,1])/denominator

  interphase.sd <- sqrt((fraction[['noncore']]^2*core[,2]^2  + fraction[['core']]^2*noncore[,2]^2)/denominator^2)
  postmitotic.sd <- sqrt(((1-fraction[['core']])^2*noncore[,2]^2 + (1-fraction[['noncore']])^2*core[,2]^2)/denominator^2)
    
 
  df <- data.frame('time_min' = df_id$time_min, 
                   'tot_int' = interphase,
                   'tot_int_sd' = interphase.sd,
                   'process' = 'interphase', 
                   'cellid' = df_id$cellid, 
                   'POI' = df_id$POI)
  df <- rbind(df, data.frame('time_min' = df_id$time_min, 
                             'tot_int' = postmitotic,
                             'tot_int_sd' = postmitotic.sd,
                             'process' = 'postmitotic', 
                             'cellid' = df_id$cellid, 
                             'POI' = df_id$POI))
  return(df)
}

# interphase_post_mitotic <- function(core, noncore, fraction){
#   denominator <- fraction[['noncore']] - fraction[['core']]
#   interphase <- (fraction[['noncore']]*core$tot_bgsub_norm - fraction[['core']]*noncore$tot_bgsub_norm)/denominator
#   
#   metaphase <- ((1-fraction[['core']])*noncore$tot_bgsub_norm - (1-fraction[['noncore']])*core$tot_bgsub_norm)/denominator
#   if ('tot_int_sd_norm' %in% colnames(core)){
#     interphase.sd <- sqrt((fraction[['noncore']]^2*core$tot_int_sd_norm^2  + fraction[['core']]^2*noncore$tot_int_sd_norm^2)/denominator^2)
#     metaphase.sd <- sqrt(((1-fraction[['core']])^2*noncore$tot_int_sd_norm^2 + (1-fraction[['noncore']])^2*noncore$tot_int_sd_norm^2)/denominator^2)
#     
#   } else {
#     interphase.sd <- 1
#     metaphase.sd <- 1
#   }
#   df <- data.frame('time_min' = core$time_min, 
#                    'tot_bgsub_norm' = interphase,
#                    'tot_sd_bgsub_norm' = interphase.sd,
#                    'process' = 'interphase', 
#                    'cellid' = core$cellid, 
#                    'POI' = core$POI)
#   df <- rbind(df, data.frame('time_min' = core$time_min, 
#                              'tot_bgsub_norm' = metaphase,
#                              'tot_sd_bgsub_norm' = metaphase.sd,
#                              'process' = 'metaphase', 
#                              'cellid' = core$cellid, 
#                              'POI' = core$POI))
#   return(df)
# }