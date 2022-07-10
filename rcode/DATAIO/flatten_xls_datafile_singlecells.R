# Flatten xls file from Shotaro Otsuka data into a txt file
# Recompute the a background subtraction and a normalization for the core region
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

# ------
# Read single cell values
outframe <- data.frame()
# pattern to recognize sheet name
pattern_sheet <- '^total_(?<POI>\\w+)_\\w+'
maxprot <- 10 # number of proteins
pattern <-list('noncore' =  'nonCore\\d+AvgTot\\w+', 
               'core'  = 'innerCore\\dTotInt\\w+')

for (name in names_sheets) {
  # Find the entry Sum that contains all proteins
  #name_split <- strsplit(name, split = '_')
  parsed <- regexpr(pattern = pattern_sheet, name, perl= T)
  
  if (attr(parsed, "capture.start") == -1) {
    next
  }
  tab <- read.xlsx(wb, sheet = name)
  d <- dim(tab)
  x <- tab[seq(13, d[1]), ]
  #x <- as.vector(x[complete.cases(x), ])
  startmat <- attr(parsed, "capture.start")
  endmat <- attr(parsed, "capture.start") + attr(parsed, "capture.length") - 1
  
  poi.name <- substr(name, startmat[, 'POI'], endmat[, 'POI'])
  idx_cell <- list('noncore' = 1, 'core' = 1)
  
  for (col_name in colnames(tab)){
    for (region in names(pattern)){
      parsed_region <- regexpr(pattern = pattern[[region]], col_name, perl= T)
      if (parsed_region[1] != -1){
         # normalize and background subtract
         
         cell_df <- data.frame('time_min' = as.double(x[,1]), 
                               'tot' = as.double(x[,col_name]))
         if (region == 'core') {
           # only perform a background subtraction for core region
           # find first 3 time points that are not NA
           idx_to_use <- which(!is.na(cell_df$tot))
           cell_df$tot_bgsub_norm <- cell_df$tot - mean(cell_df$tot[idx_to_use[c(1,2,3)]])
         } else {
            cell_df$tot_bgsub_norm <- cell_df$tot
         }
         idx_to_use <- which((!is.na(cell_df$tot)) & (cell_df$time_min > 100 & cell_df$time_min < 120))
         norm_fac <- mean(cell_df$tot_bgsub_norm[idx_to_use])
         cell_df$tot_bgsub_norm <- cell_df$tot_bgsub_norm/norm_fac
         cell_df$cellid = idx_cell[[region]]
         cell_df$region = region
         cell_df$POI =  poi.name 
         cell_df$perturb = 'none'
         
         idx_cell[[region]] <- idx_cell[[region]] + 1
         
         outframe <- rbind(outframe, cell_df)

    }
  }
  }
}

# ---- 
# compute for each cell the interphase and postmitotic assemble
outframe_inter_meta <- data.frame()
for (poiname in unique(outframe$POI)){
  poidf <- subset(outframe, POI == poiname)
  for (idcell in unique(poidf$cellid)){
    cell_df <- interphase_post_mitotic(core = subset(poidf,  cellid == idcell & region == 'core'),
                                   noncore = subset(poidf, cellid == idcell & region == 'noncore'), 
                                   fraction = fraction)
    outframe_inter_meta <- rbind(outframe_inter_meta, cell_df)
  }
}

# ----

p <- ggplot(data = subset(outframe_inter_meta, POI == 'Nup107'), aes(x = time_min, y = tot_bgsub_norm, color = process))
p <- p + stat_summary(fun.min = function(z) { quantile(z,0.25) },
                      fun.max = function(z) { quantile(z,0.75) },
                      fun = mean, na.rm = TRUE)
#p <- p + geom_ribbon(aes(x = time_min, y = tot_bgsub_norm, color = process))

p
# ------
#write.table(outframe, 
#            file = pathvar.tabwt, 
#            row.names = F, 
#            quote = F, sep = '\t')


