# Flatten xls file from Shotaro Otsuka data into a txt file
# Recompute the a background subtraction and a normalization for the core region
# compute interphase and metaphase assembly according to model
# -----
library(rstudioapi)
current_path <- getActiveDocumentContext()$path
pathvar.workdir <- dirname(current_path)
setwd(pathvar.workdir)
source('./globalvar.R')
library('openxlsx')
library('ggplot2')
wb <- loadWorkbook(file = pathvar.xlswt)
names_sheets <- wb$sheet_names

# ----
# Read mean values only
names_sheets_split <- strsplit(names_sheets, split = '_')
outframe <- data.frame()
# pattern to recognize name of POI and number of cells
pattern1 <- '^(?<POI>\\w+)\\.\\(N=(?<n>\\d+)\\)'
maxprot <- 10 # number of proteins

for (name in names_sheets) {
  # Find the entry Sum that contains all proteins
  #name_split <- strsplit(name, split = '_')
  if (name == "Sum") {
    tabs <- list('noncore' = read.xlsx(wb, sheet = name, cols = seq(2, maxprot*2+2), startRow = 2), 
                 'core' = read.xlsx(wb, sheet = name, cols = c(2, seq(maxprot*2+4, maxprot*4+4)), startRow = 2))
    
    # For noncore and core region
    for (region in names(tabs)) {
      tab <- tabs[[region]]
      d <- dim(tab)
      
      # for each poi
      for (ipoi in seq(2, d[2], 2)){
        name_poi = colnames(tab)[ipoi]
        parsed <- regexpr(pattern = pattern1, name_poi , perl= T)
        
        startmat <- attr(parsed, "capture.start")
        endmat <- attr(parsed, "capture.start") + attr(parsed, "capture.length") - 1
        
        poi.name <- substr(name_poi, startmat[, 'POI'], endmat[, 'POI'])
        poi.n <- as.numeric(substr(name_poi, startmat[, 'n'], endmat[, 'n']))
        # start from 5th row (time starts at 3.5 min)
        x <- tab[seq(5, d[1]), c(1, ipoi, ipoi + 1)]
        x <- as.vector(x[complete.cases(x), ])
        
        cell_df <- data.frame('time_min' = as.double(x[,1]), 
                              'tot_int' = as.double(x[,2]), 
                              'tot_int_sd' = as.double(x[,3]))
        idx_to_use <- which((!is.na(cell_df$tot_int)) & (cell_df$time_min > 100 & cell_df$time_min < 120))
        norm_fac <- mean(cell_df$tot_int[idx_to_use])
        cell_df$tot_int_norm <- cell_df$tot_int/norm_fac
        cell_df$tot_int_sd_norm <- cell_df$tot_int_sd/norm_fac
        
        if (region == 'core') {
          # only perform a background subtraction for core region
          # find first 3 time points that are not NA
          idx_to_use <- which(!is.na(cell_df$tot_int))
          cell_df$tot_int_bgsub_norm <- cell_df$tot_int - mean(cell_df$tot_int[idx_to_use[c(1,2,3)]])
          
          #cell_df$tot_bgsub_norm <- cell_df$tot_int_mean
        } else {
          cell_df$tot_int_bgsub_norm <- cell_df$tot_int
        }
        
        idx_to_use <- which((!is.na(cell_df$tot_int)) & (cell_df$time_min > 100 & cell_df$time_min < 120))
       
        norm_fac <- mean(cell_df$tot_int_bgsub_norm[idx_to_use])
        cell_df$tot_int_bgsub_norm <- cell_df$tot_int_bgsub_norm/norm_fac
        cell_df$tot_int_sd_bgsub_norm <- cell_df$tot_int_sd/norm_fac
        
        
        
        
        cell_df$cellid = 1
        cell_df$region = region
        cell_df$POI =  poi.name 
        outframe <- rbind(outframe, 
                          cell_df)
        
      }
    }
  }
}






write.table(outframe, 
            file = pathvar.tabwt_norm_bgsub, 
            row.names = F, 
            quote = F, sep = '\t')



# ---- 
# initial background by subtracting the average of the first 3 points core region
outframe_inter_meta <- data.frame()
outframe_inter_meta_bgsub <- data.frame()
for (poiname in unique(outframe$POI)){
  poidf <- subset(outframe, POI == poiname)
  for (idcell in unique(poidf$cellid)){
    # reassign values for consistency
    core <- subset(poidf, region == 'core')
    noncore <- subset(poidf, region == 'noncore')
    
    cell_df <- interphase_post_mitotic(core = cbind(core$tot_int_norm, core$tot_int_sd_norm),
                                       noncore =  cbind(noncore$tot_int_norm, noncore$tot_int_sd_norm), 
                                       df_id = core,
                                       fraction = fraction)
    cell_df_bgsub <- interphase_post_mitotic(core = cbind(core$tot_int_bgsub_norm, core$tot_int_sd_bgsub_norm),
                                       noncore =  cbind(noncore$tot_int_bgsub_norm, noncore$tot_int_sd_bgsub_norm), 
                                       df_id = core,
                                       fraction = fraction)
    outframe_inter_meta <- rbind(outframe_inter_meta, cell_df)
    outframe_inter_meta_bgsub <- rbind(outframe_inter_meta_bgsub, cell_df_bgsub)
    
  }
}
write.table(outframe_inter_meta_bgsub, 
            file = pathvar.tabwt_norm_bgsub_inter_meta, 
            row.names = F, 
            quote = F, sep = '\t')

# ---- 
# create a table for shotaro with a sheet per protein 
options(java.parameters = c("-XX:+UseConcMarkSweepGC", "-Xmx8192m"))
gc()
file.remove( pathvar.xlsx_shotaro)

for (poiname in unique(outframe$POI)){
  poidf <- subset(outframe, POI == poiname)
  df_core <- subset(poidf, region == 'core')[, seq(2,7)]
  colnames(df_core) <- paste0('core_', colnames(df_core))
  df_noncore <- subset(poidf, region == 'noncore')[, seq(2,7)]
  colnames(df_noncore) <- paste0('noncore_', colnames(df_noncore))
  
  poidf <- subset(outframe_inter_meta, POI == poiname)
  df_interphase <- subset(poidf, process == 'interphase')[, seq(2,3)]
  colnames(df_interphase) <- paste0('interphase_norm', colnames(df_interphase))
  df_postmitotic <- subset(poidf, process == 'postmitotic')[, seq(2,3)]
  colnames(df_postmitotic) <- paste0('postmitotic_norm', colnames(df_postmitotic))
  
  poidf <- subset(outframe_inter_meta_bgsub, POI == poiname)
  df_interphase_bgsub <- subset(poidf, process == 'interphase')[, seq(2,3)]
  colnames(df_interphase_bgsub) <- paste0('interphase_bgsub_', colnames(df_interphase_bgsub))
  df_postmitotic_bgsub <- subset(poidf, process == 'postmitotic')[, seq(2,3)]
  colnames(df_postmitotic_bgsub) <- paste0('postmitotic_bgsub_', colnames(df_postmitotic_bgsub) )
  
  df <- cbind('time_min' = subset(poidf, process == 'interphase')$time_min, df_core, df_noncore, 
              df_interphase, df_postmitotic, 
              df_interphase_bgsub, df_postmitotic_bgsub)
  
  xlsx::write.xlsx2(file = pathvar.xlsx_shotaro, x = df, sheetName = poiname,row.names = FALSE, append = TRUE)
}

# ----
apoi = 'Nup358'
p <- ggplot(data = subset(outframe_inter_meta, POI == apoi))
p <- p + geom_ribbon(aes( x = time_min, ymin = tot_bgsub_norm - tot_sd_bgsub_norm, 
                 ymax = tot_bgsub_norm + tot_sd_bgsub_norm,  fill = process, alpha = 0.5))
p <- p + geom_line(aes(x = time_min, y =  tot_bgsub_norm, color = process))
p <- p + geom_line(data = subset(outframe, POI == apoi), aes(x = time_min, y =  tot_bgsub_norm, color = region))
#p <- p + geom_ribbon(aes(x = time_min, y = tot_bgsub_norm, color = process))

p