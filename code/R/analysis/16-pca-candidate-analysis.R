source('code/R/functions/session_path.R')
## restore session
session::restore.session(session_path("15"))

## load parameters
rda_parameters <- "code/parameters/rda.toml" %>%
  RcppTOML::parseTOML() %>%
  `[[`(MODE)

## packages
library(gdsfmt)
library(SNPRelate)
library(ggplot2)

## load functions
source("code/R/functions/rsq.R")

## read data
  snp_admix_data <- readRDS(snp_pop_admix_path)
  x <- snp_metadata[[1]][grep('Valdemanco',snp_metadata[[1]]$locality_id),]$order_id
  snp_admix_data[[1]] <- snp_admix_data[[1]][!snp_admix_data[[1]]$ind.names %in% x, ]
  #snp_admix_data[[1]] <- snp_admix_data[[1]][snp_admix_data[[1]]$ind.names %in% snp_metadata[[1]]$order_id,]
  load('data/final/results.rda')
  #bayenv_results <- readRDS('data/final/bayenv_cand.rds')
  
## create directory to save Venn diagrams results
unlink("data/intermediate/pca-adapt", force = TRUE, recursive = TRUE)
dir.create("data/intermediate/pca-adapt", showWarnings = FALSE,
           recursive = TRUE)
## create directory to save LFMM results
curr_dir <- paste0("data/intermediate/pca-adapt/",
                   spp_parameters$species_atlas_names[i])
dir.create(curr_dir, showWarnings = FALSE, recursive = TRUE)

{## id of candidates
## id of OA candidates
oa_cand <- list(BS,SE,PA)
  oa_cand <- unlist(oa_cand)
  oa_con <- oa_cand[duplicated(oa_cand)]
  oa_all <- levels(as.factor(oa_cand))
## id of EAA candidates
ea_cand <- list(BE,SB,LF,RD)
  ea_cand <- unlist(ea_cand)
  ea_con <- ea_cand[duplicated(ea_cand)]
  ea_all <- levels(as.factor(ea_cand))
## id of MEM candidates
#mm_cand <- MM
## id of neutral candidates
nt_cand <- unique(c(ea_cand,
                    oa_cand
#                    mm_cand
                    ))
  nt_cand <- snp_admix_data[[1]][,!snp_admix_data[[1]]$loc.names %in% nt_cand]$loc.names
## Random SNP subset
#x <- sample.int(
#  length(snp_admix_data[[1]]$loc.names), # Among all SNPs
#  50, # Select 50 SNPs
#  replace = TRUE) # Do not duplicate
#rnd_cand <- snp_admix_data[[1]]$loc.names[x]

## Save candidates as lists
cand <- list(nt_cand, 
             #oa_con, oa_all, 
             BS, SE, PA, 
             #ea_con, ea_all, 
             BE, LF, RD 
#             mm_cand
#             rnd_cand
)
cand_name <- list("Neutral", 
                  #"OA_Concurrent", "OA_All", 
                  "Bayescan", "SelEstim", "PCAdapt",
                  #"EA_Concurrent", "EA_All", 
                  "Bayenv2", "LFMM", "RDA" 
#                  "MSOD",
#                  "Random"
)
cand_name <- unlist(cand_name)
}
{## define envi.cluster colours
snp_clust <- snp_metadata[[1]][snp_metadata[[1]]$order_id %in% snp_admix_data[[1]]$ind.names,]
  snp_clust$order_id <- as.numeric(snp_clust$order_id)
  snp_clust <- snp_clust[order(snp_clust$order_id),]
    snp_clust$cluster_id <- as.factor(snp_clust$cluster_id)
  levels(snp_clust$cluster_id) <- rda_parameters$env
  bg <- rda_parameters$col
  }
{## read matrix
snp_mx <- readRDS(snp_mx_path)
}
## limit input data to candidates alleles 
snp_cand <- lapply(seq_along(cand), function(i) {
  snp_cand <- snp_mx[,colnames(snp_mx) %in% cand[[i]] == TRUE] 
  ## RDA requires complete data frames (i.e., no missing data). 
  ## Using the most common genotype at each SNP across all individuals.
  apply(snp_cand, 2, function(x)
          replace(x, is.na(x), as.numeric(names(which.max(table(x))))))
})
## generate a single column data frame
snp_df <- lapply(seq_along(cand), function(i) {
  t(snp_cand[[i]])
})
## save and reload as geno object
lapply(seq_along(cand), function(i) {
  write.geno(snp_df[[i]], output.file = paste0(curr_dir,'/snp_',cand_name[i],'.geno'))
  })  
snp_geno <- lapply(seq_along(cand), function(i) {
    read.geno(paste0(curr_dir,'/snp_',cand_name[i],'.geno'))
  })
## generate .gds file
cand_gds <- lapply(seq_along(cand), function(i) {
  SNPRelate::snpgdsCreateGeno(gds.fn = paste0(curr_dir,'/snp_cand_',cand_name[i],'.gds'), genmat = snp_geno[[i]],
                              sample.id = rownames(snp_cand[[i]]), snp.id = colnames(snp_cand[[i]]),
                              snpfirstdim = TRUE)
  snpgdsOpen(paste0(curr_dir,'/snp_cand_',cand_name[i],'.gds'))
  })
  ## perform PCA
  cand_pca <- lapply(seq_along(cand), function(i) {
    snpgdsPCA(cand_gds[[i]], snp.id = colnames(snp_cand[[i]]),
              sample.id = rownames(snp_cand[[i]]), num.thread = 2)
})
  ## variance proportion (%)
  pc.percent <- lapply(seq_along(cand), function(i) {
    cand_pca[[i]]$varprop*100
  })
{## Plot OA results
    ## Get sample id
    sample.id <- read.gdsn(index.gdsn(cand_gds[[i]], "sample.id"))
    ## Get population information
    pop_code <- snp_admix_data[[1]]$pop
    #pop_code <- snp_oa$pop
      levels(pop_code) <- c("East","West")
      lineage <- c(21,22)
  }
## Make a data.frame
tab <- lapply(seq_along(cand), function(i) {
  data.frame(sample.id = cand_pca[[i]]$sample.id,
             pop = factor(pop_code)[match(cand_pca[[i]]$sample.id, sample.id)],
             clust = snp_clust$cluster_id,
             EV1 = cand_pca[[i]]$eigenvect[,1],    # the first eigenvector
             EV2 = cand_pca[[i]]$eigenvect[,2],    # the second eigenvector
             EV3 = cand_pca[[i]]$eigenvect[,3],    # the second eigenvector
             stringsAsFactors = FALSE)
})
  lapply(seq_along(cand), function(i) {
  write.csv(tab[[i]], paste0(curr_dir,'/cand_eigenvecter_',cand_name[i],'.csv'))
  })
## Print 1st PC value to spatially visualize described variation.
pc_plot <- lapply(seq_along(cand), function(i) {
    cbind(tab[[i]]$sample.id,sambada_env$lat,sambada_env$lon,tab[[i]]$EV1, tab[[i]]$EV2)
  })
  lapply(seq_along(cand), function(i) {
    colnames(pc_plot[[i]]) <- c("snp","lat","lon","pc","pc2")
    write.csv(pc_plot[[i]],paste0(curr_dir,'/PC_',cand_name[i],'.csv'))
  })
  

{### Procrustes distances 
procrat_name <- list()                              # name of comparisons
idx = 0                                           # start index to 0
for (i in seq_along(cand_name)){                    # for each list of candidates
    #not_i = seq_along(cand)[seq_along(cand_name)!=i]# every alternative list of candidates 
    #for (j in not_i){
    for (j in i:length(cand_name)){
      idx <- idx + 1                            # increase index
      procrat_name[idx] <- paste0(cand_name[i],".",cand_name[j])
    }
}
procrat_name <- unlist(procrat_name)

procrat <- list()                                 # list of Proc distances
idx = 0                                         # start index to 0
for (i in seq_along(cand)){                       # for each list of candidates
  #not_i = seq_along(cand)[seq_along(cand)!=i]     # every alternative list of candidates 
  #for (j in not_i){
  for (j in i:length(cand_name)){
    idx <- idx + 1                            # increase index
    # Calculate procrustes distance
    procrat[[idx]] <- procrustes(cand_pca[[i]]$eigenvect[,1:5],
                                 cand_pca[[j]]$eigenvect[,1:5],
                                 scores = 'sites', scale = T)#, symmetric = T)
  }
}
}
{### Estimate the significance of the Procrustes statistic.
prot <- list()                                 # list of Proc distances
  idx = 0                                         # start index to 0
for (i in seq_along(cand)){                       # for each list of candidates
  #not_i = seq_along(cand)[seq_along(cand)!=i]     # every alternative list of candidates 
  #for (j in not_i){
  for (j in i:length(cand_name)){
    idx <- idx + 1                            # increase index
    # Calculate significance
    prot[[idx]] <- protest(cand_pca[[i]]$eigenvect[,1:5],
                 cand_pca[[j]]$eigenvect[,1:5],
                 scores = 'sites', scale = T, permutations = 1000)
    
  }
}
}
### Determine which comparisons are significant
prot_sig <- list()
for (i in seq_along(prot)){
  prot_sig[i] <- print(prot[[i]]$signif)
  } 
prot_sig <- unlist(prot_sig)
  
{### Sum of square differences between Procrustese differences
procrat_ss <- list()
prot_ss <- list()
idx = 0                                         # start index to 0
for (i in seq_along(cand)){                       # for each list of candidates
  #not_i = seq_along(cand)[seq_along(cand)!=i]     # every alternative list of candidates 
  #for (j in not_i){
  for (j in i:length(cand_name)){
    idx <- idx + 1                            # increase index
    procrat_ss[idx] <- procrat[[idx]]$ss
    prot_ss[idx] <- prot[[idx]]$ss
    }
}
procrat_ss <- unlist(procrat_ss)
prot_ss <- unlist(prot_ss)
}
{### Quartile for Most and Least Different according to scaled Sum Squares 
procrat_summary <- data.frame(procrat_name,procrat_ss,prot_ss, prot_sig)
x <- quantile(prot_ss, prob = seq(0, 1, length = 11), type = 5)
q10 <- procrat_summary[prot_ss<x[2],1]
q20 <- procrat_summary[prot_ss<x[3],1]
q90 <- as.numeric(rownames(procrat_summary[prot_ss>x[10],]))
q80 <- procrat_summary[prot_ss>x[9],1]
}

## Procrustese Risduals for most differentiated comparisons suitable for QGIS
for (i in seq_along(cand_name)){
    plot(procrat[[i]], kind = 2)
    x <- cbind(pc_plot[[1]][,2:3],residuals(procrat[[i]]),
               residuals(procrat[[i]])/max(residuals(procrat[[i]])))
    colnames(x) <- c("lat","lon","residual","res.scaled")
    write.csv(x,paste0(curr_dir,'/procat_',procrat_name[i],'.csv'))
  }
  
  ## Draw
  for (i in seq_along(cand)){
    par(mfrow=c(1,1))
    levels(snp_admix_data[[1]]$pop) <- c("East","West")
    lineage <- snp_admix_data[[1]]$pop
    bg = c("darkslateblue","gold")
    ## axes 1 & 2
    jpeg(paste0(curr_dir,"/pca_",cand_name[i],".jpeg"), width = 300, height = 300)
    plot(tab[[i]][,c("EV1","EV2")], type="n", scaling=3, xlab = "PC1", ylab = "PC2")
    #points(snp.rda, display="species", pch=20, cex=0.7, col="gray32", scaling=3)           # the SNPs
    points(tab[[i]][,c("EV1","EV2")], 
           display="pop", pch=21, cex=1.3, col="gray32", scaling=3,
           bg=bg[lineage]) # the wolves
    #text(snp.rda, scaling=3, display="bp", col="#0868ac", cex=1)                           # the predictors
    legend("bottomleft", legend=levels(lineage), bty="n", col="gray32", pch=21, cex=1, pt.bg=bg)
    title(main=cand_name[i],cex.lab=0.75)
    mtext(paste0("SS: ",signif(prot_ss[i], digits = 3)),
          side=3, line=0.5, cex.lab=3, adj = 1)
    mtext(paste0("SNPs: ",length(cand[[i]])),
          side=3, line=0.5, cex.lab=3, adj = 0)
    dev.off()
  }
  
  
## write results to file
write.csv(procrat_summary, paste0(curr_dir,'/procrustese.csv'))


## Cumulative Procrustese distance from Neutrality per location
procrat_name[1:length(cand)]
loc_res <- list()
loc_procrat <- list() # summary object
x <- list()
# all the residual and scaled residuals per sample
for (i in 1:length(cand)){
  x[[i]] <- cbind.data.frame(
    pc_plot[[i]][,1:3],
    env_raw[env_raw$vaucher %in% pc_plot[[i]][,1],]$site,
    residuals(procrat[[i]]),
    residuals(procrat[[i]])/max(residuals(procrat[[i]]))
    )
  colnames(x[[i]]) <- c("vaucher","lat","lon","site",
                   "residual","res.scaled")
 # Unique locations per sample
  loc_name = unique(x[[i]]$site)
  z <- list() # location residual
  for (l in seq_along(loc_name)){
    z[l] = mean(x[[i]][x[[i]]$site==loc_name[l],]$residual)
    }
  loc_res[[i]] <- unlist(z)
}
## Combine into single file
loc_procrat <- cbind.data.frame(
  unique(x[[1]]$site),
  unique(x[[1]]$lat),
  unique(x[[1]]$lon)
)
## Single object with residual for each site
for (i in 2:length(cand)){
  loc_procrat <- cbind.data.frame(
    loc_procrat,
    loc_res[[i]]
    )
}
## Mean residual accorss methods
loc_res_mean <- list()
for (i in 1:nrow(loc_procrat)){
  loc_res_mean[i] <- mean(sapply(loc_procrat[-1:-3], "[[", i))
}
loc_procrat[ncol(loc_procrat)+1] <- unlist(loc_res_mean)
## Join variables
colnames(loc_procrat) <- c("site","lat","lon",
                      cand_name[2:length(cand)],
                      "mean")
# write to file
write.csv(loc_procrat, paste0('data/final/residual_table.csv'))


# clean-up
rm(snp_admix_data, snp_df, snp_geno, snp_mx, structure_trace_data,
   #snp_cand, cand, cand_pca, procrat, procrat_list, prot,
   bayescan_results, bayescan_trace_data, 
   bayenv_results, bayescan_outputs, 
   mod, msod_outlier_loci, pcadapt_outlier_loci, selestim_outlier_loci, snp_cand,
   pcadapt_results, pcadapt_inital_run, selestim_results, selestim_trace_data
   )


# save session
session::save.session(session_path("16"), compress = "xz")
