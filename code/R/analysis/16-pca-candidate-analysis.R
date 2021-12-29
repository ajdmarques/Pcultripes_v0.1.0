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
  snp_admix_data[[1]] <- snp_admix_data[[1]][snp_admix_data[[1]]$ind.names %in% rownames(Sp_mat[[1]])]
  load('data/final/results.rda')
  #bayenv_results <- readRDS('data/final/bayenv_cand.rds')
  
## create directory to save Venn diagrams results
unlink("data/intermediate/pca-adapt", force = TRUE, recursive = TRUE)
dir.create("data/intermediate/pca-adapt", showWarnings = FALSE,
           recursive = TRUE)
## create directory to save LFMM results
curr_dir <- paste0("data/intermediate/pca-adapt/",
                   spp_parameters$species_atlas_names)
dir.create(curr_dir, showWarnings = FALSE, recursive = TRUE)

{## define envi.cluster colours
snp_clust <- snp_metadata[[1]][snp_metadata[[1]]$order_id %in% snp_admix_data[[1]]$ind.names,]
  snp_clust$order_id <- as.numeric(snp_clust$order_id)
  snp_clust <- snp_clust[order(snp_clust$order_id),]
    snp_clust$cluster_id <- as.factor(snp_clust$cluster_id)
  levels(snp_clust$cluster_id) <- rda_parameters$env
  bg <- rda_parameters$col
  }
{## read matrix
#snp_mx <- readRDS(snp_mx_path)
  # Update to incldue SNP matrix from GF
  Sp_mat <- readRDS("~/Documents/NGC/dart-adaptation-revised/data/intermediate/gradientforest-results/snp_matrix.rds")
  # candidtaes to Keep 4 = Bayenv2, 6 = RDA, 7 = Neutral
  cand_name = c(
    "Bayescan", "RDA", "Neutral"
  )
  x = match(cand_name,test_names)
  snp_mx <- list()
  for (i in seq_along(x)){
    snp_mx[[i]] <- as.matrix(Sp_mat[[x[i]]])
    }
}
## limit input data to candidates alleles 
snp_cand <- lapply(seq_along(snp_mx), function(i) {
  ## RDA requires complete data frames (i.e., no missing data). 
  ## Using the most common genotype at each SNP across all individuals.
  apply(snp_mx[[i]], 2, function(x)
          replace(x, is.na(x), as.numeric(names(which.max(table(x))))))
})
## generate a single column data frame
snp_df <- lapply(seq_along(snp_cand), function(i) {
  t(snp_cand[[i]])
})
## save and reload as geno object
lapply(seq_along(snp_cand), function(i) {
  write.geno(snp_df[[i]], output.file = paste0(curr_dir,'/snp_',cand_name[i],'.geno'))
  })  
snp_geno <- lapply(seq_along(snp_cand), function(i) {
    read.geno(paste0(curr_dir,'/snp_',cand_name[i],'.geno'))
  })
## generate .gds file
cand_gds <- lapply(seq_along(snp_cand), function(i) {
  SNPRelate::snpgdsCreateGeno(gds.fn = paste0(curr_dir,'/snp_cand_',cand_name[i],'.gds'), genmat = snp_geno[[i]],
                              sample.id = rownames(snp_cand[[i]]), snp.id = colnames(snp_cand[[i]]),
                              snpfirstdim = TRUE)
  snpgdsOpen(paste0(curr_dir,'/snp_cand_',cand_name[i],'.gds'))
  })
  ## perform PCA
  cand_pca <- lapply(seq_along(snp_cand), function(i) {
    snpgdsPCA(cand_gds[[i]], snp.id = colnames(snp_cand[[i]]),
              sample.id = rownames(snp_cand[[i]]), num.thread = 2)
})
  ## variance proportion (%)
  pc.percent <- lapply(seq_along(snp_cand), function(i) {
    cand_pca[[i]]$varprop*100
  })
{## Plot OA results
    ## Get sample id
    sample.id <- read.gdsn(index.gdsn(cand_gds[[1]], "sample.id"))
    ## Get population information
    pop_code <- snp_admix_data[[1]]$pop
    #pop_code <- snp_oa$pop
      levels(pop_code) <- c("East","West")
      lineage <- c(21,22)
  }
## Make a data.frame
tab <- lapply(seq_along(snp_cand), function(i) {
  data.frame(sample.id = cand_pca[[i]]$sample.id,
             pop = factor(pop_code)[match(cand_pca[[i]]$sample.id, sample.id)],
             clust = snp_clust$cluster_id,
             EV1 = cand_pca[[i]]$eigenvect[,1],    # the first eigenvector
             EV2 = cand_pca[[i]]$eigenvect[,2],    # the second eigenvector
             EV3 = cand_pca[[i]]$eigenvect[,3],    # the second eigenvector
             stringsAsFactors = FALSE)
})
  lapply(seq_along(snp_cand), function(i) {
  write.csv(tab[[i]], paste0(curr_dir,'/cand_eigenvecter_',cand_name[i],'.csv'))
  })
## Print 1st PC value to spatially visualize described variation.
pc_plot <- lapply(seq_along(snp_cand), function(i) {
    cbind.data.frame(tab[[i]]$sample.id,sambada_env$lat,sambada_env$lon,tab[[i]]$EV1, tab[[i]]$EV2)
  })
  for (i in seq_along(pc_plot)){
    colnames(pc_plot[[i]]) <- c("vaucher","lat","lon","pc","pc2")
    write.csv(pc_plot[[i]],paste0(curr_dir,'/PC_',cand_name[i],'.csv'))
  }

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

procrat <- list()                               # list of Proc distances
idx = 0                                         # start index to 0
for (i in seq_along(snp_cand)){                 # for each list of candidates
  #not_i = seq_along(cand)[seq_along(cand)!=i]  # every alternative list of candidates 
  #for (j in not_i){
  for (j in i:length(cand_name)){
    idx <- idx + 1                              # increase index
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
for (i in seq_along(snp_cand)){                       # for each list of candidates
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
for (i in seq_along(cand_name)){                       # for each list of candidates
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
# comparisons of interest
procrat_name # Select items to include
#k = c(3,5)

for (i in seq_along(procrat_name)){
    plot(procrat[[i]], kind = 1)
    x <- cbind(pc_plot[[1]][,2:3],residuals(procrat[[i]]),
               residuals(procrat[[i]])/max(residuals(procrat[[i]])))
    colnames(x) <- c("lat","lon","residual","res.scaled")
#    write.csv(x,paste0(curr_dir,'/procat_',procrat_name[i],'.csv'))
  }
  
  ## Draw
  for (i in seq_along(cand_name)){
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
    #mtext(paste0("SS: ",signif(prot_ss[i], digits = 3)),
    #      side=3, line=0.5, cex.lab=3, adj = 1)
    mtext(paste0("SNPs: ",length(rownames(snp_df[[i]]))),
          side=3, line=0.5, cex.lab=3, adj = 0)
    dev.off()
  }
  
  
## write results to file
write.csv(procrat_summary, paste0(curr_dir,'/procrustese.csv'))


## Cumulative Procrustese distance from Neutrality per location
loc_res <- list()
loc_procrat <- list() # summary object
x <- list()
## all the residual and scaled residuals per sample
for (i in seq_along(k)){
  loc_res[[i]] <- cbind.data.frame(
    pc_plot[[i]][,1:3],
    env_raw[env_raw$vaucher %in% pc_plot[[i]][,1],]$site,
    residuals(procrat[[k[i]]]),
    residuals(procrat[[k[i]]])/max(residuals(procrat[[k[i]]]))
    )
  colnames(loc_res[[i]]) <- c("vaucher","lat","lon","site","residual","scaled.residual")
}

## Mean residual accorss methods
mean_procrat <- list()
for (i in seq_along(loc_res)){
mean_procrat[[i]] <- aggregate.data.frame(loc_res[[i]], 
                                     by = list(loc_res[[i]]$site), FUN = mean)
write.csv(mean_procrat[[i]], paste0(curr_dir,"/",
                                    procrat_name[k[i]], "_residual.csv"))
}

# write to file
write.csv(mean_procrat, paste0('data/final/residual_table.csv'))


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
