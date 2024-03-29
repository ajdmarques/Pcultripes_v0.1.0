source('code/R/functions/session_path.R')
## restore session
session::restore.session(session_path("16"))

## load parameters
mclust_parameters <- "code/parameters/mclust.toml" %>%
  RcppTOML::parseTOML() %>%
  `[[`(MODE)

## packages
library(mclust)

# load functions
#source("code/R/functions/*.R")

## read data
eig <- list()
for (i in seq_along(cand_name)){
  eig[[i]] <- read.csv(paste0("data/intermediate/pca-adapt/pelcul/cand_eigenvecter_",
                              cand_name[i],".csv"))
  eig[[i]] <- eig[[i]][,-1]
  #z <- c(107150,107151) # Remove samples from Cuidade Real
  #eig[[i]] <- eig[[i]][!eig[[i]]$sample.id %in% z,]
}
## MSOD causes script to crash there it will be removed
#eig <- eig[-8]
#cand_name <- cand_name[-8]
#mclust_parameters$G <- mclust_parameters$G[-8]

## environmental variables
env_data <- read.csv('data/intermediate/vif-results/env.csv', row.names = 1)
#env_data <- env_data[!duplicated(env_data$.),]
env_data <- env_data[env_data$. %in% eig[[i]]$sample.id,]


## create directory to save Venn diagrams results
unlink("data/intermediate/mclust-adapt/", force = TRUE, recursive = TRUE)
dir.create("data/intermediate/mclust-adapt/", showWarnings = FALSE,
           recursive = TRUE)
## create directory to save LFMM results
curr_dir <- paste0("data/intermediate/mclust-adapt/",
                   spp_parameters$species_atlas_names)
dir.create(curr_dir, showWarnings = FALSE, recursive = TRUE)

## define lineage 
par(mfrow=c(1,1))
lineage <- lapply(seq_along(eig), function(i) {
  x <- eig[[i]]$pop
})
## restrict data to variables
X <- lapply(seq_along(eig), function(i) {
  x <- eig[[i]][,4:6]
})
## correlation among classifications
for (i in seq_along(X)) {
  clPairs(X[[i]], lineage[[i]])
}
## Best BIC values for Gaussian clustering
BIC <- lapply(seq_along(X), function(i) {
  mclustBIC(X[[i]], G = seq(from=1,to=5))
   #BIC <- lapply(seq_along(1:7), function(i) {
  ## Force BIC to calculate set K
  #mclustBIC(X[[i]], G = mclust_parameters$k)
})
for (i in seq_along(BIC)){
  print(summary(BIC[[i]]))
  plot(BIC[[i]])
}
## Gaussian finite mixture model fitted by EM algorithm 
mod1 <- lapply(seq_along(BIC), function(i) {
  Mclust(X[[i]], x = BIC[[i]], G = seq(from=1,to=5))
})
par(mfrow=c(1,1))
for (i in seq_along(mod1)){
  summary(mod1[[i]], parameters = TRUE)  
  plot(mod1[[i]], what = "classification")
  title(main = cand_name[[i]], sub = "classification") 
  table(lineage[[i]], mod1[[i]]$classification)  
  plot(mod1[[i]], what = "uncertainty")
  title(main = cand_name[[i]], sub = "uncertainty")
}
## Bootstrap sequential LRT for the number of mixture components
ICL <- lapply(seq_along(X), function(i) {
  mclustICL(X[[i]], G = seq(from=1,to=5))
})
for (i in seq_along(ICL)){
  summary(ICL[[i]])
  plot(ICL[[i]])
} 

## Bootstrap sequential LRT for the number of mixture components
LRT <- lapply(seq_along(X), function(i) {
  mclustBootstrapLRT(X[[i]], 
                     modelName = substr(names(summary(ICL[[i]])[1]), 1, 3),
                     #maxG = mclust_parameters$G[i],
                     maxG = best_k+1,
                     nboot = mclust_parameters$nboot)
})
LRT

## EM algorithm is used by mclust for maximum likelihood estimation. 
## Initialisation of EM is performed using the partitions obtained from agglomerative hierarchical clustering. 
## LRT used to define the number of clusters. 
## Hscaled SVD transformation;
BIC <- lapply(seq_along(X), function(i) {
  hc1 <- hc(X[[i]], modelName = 'VVV', use = "SVD")
  BIC1 <- mclustBIC(X[[i]], G = 1:max(LRT[[i]]$G), initialization = list(hcPairs = hc1)) # default 
  (hc2 <- hc(X[[i]], modelName = 'VVV', use = "VARS"))
  BIC2 <- mclustBIC(X[[i]], G = 1:max(LRT[[i]]$G), initialization = list(hcPairs = hc2))
  (hc3 <- hc(X[[i]], modelName = 'EEE', use = "SVD"))
  BIC3 <- mclustBIC(X[[i]], G = 1:max(LRT[[i]]$G), initialization = list(hcPairs = hc3))
  mclustBICupdate(BIC1, BIC2, BIC3)
})

for (i in seq_along(BIC)){
  print(summary(BIC[[i]]))
  plot(BIC[[i]])
} 

## Clustering
## (Optional) Redefine Best BIC values based on Hierarchical Clustering 
mod1 <- lapply(seq_along(BIC), function(i) {
  Mclust(X[[i]], x = BIC[[i]])
})
## Assign clusters
mod1dr <- lapply(seq_along(mod1), function(i) {
  MclustDR(mod1[[i]], lambda = 0.05)
})

## Plot ordinations 
## with clusters for PC1 vs PC2/PC3 
for (i in seq_along(mod1dr)){
  # Save file as JPEG
  jpeg(paste0(curr_dir,'/mclust_',cand_name[i],'.jpeg'), width = 500, height = 500)
  # Display PC1~PC2 and PC1~PC3
  par(mfrow=c(1,2))
  ## PC1 ~ PC2
  #plot(mod1dr[[i]]$x, type="n", scaling=3, xlab = "PC1", ylab = "PC2")
  plot(mod1dr[[i]]$x[,1] ~ mod1dr[[i]]$x[,2], type="n", 
       #scaling=3, 
       xlab = "PC1", ylab = "PC2", 
       axes = T, ann = F)
  axis(1)
  axis(2)
  #points(snp.rda, display="species", pch=20, cex=0.7, col="gray32", scaling=3)           # the SNPs
  points(mod1dr[[i]]$x[,1] ~ mod1dr[[i]]$x[,2],
         #mod1dr[[i]]$x, 
         #display="pop", 
         pch=21, cex=1.5, col="gray32", 
         #scaling=3,
         bg=mclust_parameters$bg[mod1dr[[i]]$classification])
  ## Array details
  title(main=cand_name[i],
        cex.lab=1, cex.main = 1 
        #outer = T
        )
  mtext(paste0("SS: ",signif(prot_ss[i], digits = 3)),
        side=3, line=0.5, cex = 1, cex.lab= 1, adj = 1 
        #outer = T
        )
  mtext(paste0("SNPs: ",length(cand[[i]])),
        side=3, line=0.5, cex = 1, cex.lab=3, adj = 0 
        #outer = T
        )
  ## PC1 ~ PC3
  #plot(mod1dr[[i]]$x, type="n", scaling=3, xlab = "PC1", ylab = "PC2")
  plot(mod1dr[[i]]$x[,1] ~ mod1dr[[i]]$x[,3], type="n", 
       scaling=3, 
       xlab = "PC1", ylab = "PC3", 
       axes = T, ann = F)
  axis(1)
  axis(4)
  #points(snp.rda, display="species", pch=20, cex=0.7, col="gray32", scaling=3)           # the SNPs
  points(mod1dr[[i]]$x[,1] ~ mod1dr[[i]]$x[,3],
         #mod1dr[[i]]$x, 
         display="pop", pch=21, cex=1.5, col="gray32", scaling=3,
         bg=mclust_parameters$bg[mod1dr[[i]]$classification])
  #text(snp.rda, scaling=3, display="bp", col="#0868ac", cex=1)                           # the predictors
  #legend("topright", legend=levels(mod1dr[[i]]$class), 
  #       bty="n", col="gray32", pch=21, cex=1, pt.bg= mclust_parameters$bg)
  dev.off()
} 
## join cluster assignments to environmental data
cand_clust <- lapply(seq_along(mod1), function(i) {
  cbind(env_data[,1:3], 
        eig[[i]]$pop, eig[[i]]$clust,
        mod1[[i]]$classification, mod1[[i]]$uncertainty)
})
for (i in seq_along(cand_clust)){
  colnames(cand_clust[[i]]) <- c('vaucher','lat','lon','lineage','env_clust','cand_clust','uncertainty')
  write.csv(cand_clust[[i]], paste0(curr_dir,'/cand_clust_',cand_name[[i]],'.csv'))
}

{# Generate summary table
lat <- unique(cand_clust[[1]]$lat)
pie_loc <- list()
k1 <- list()
k2 <- list()
k3 <- list()
k4 <- list()
k5 <- list()
for (i in 1:length(cand_clust)){
  for (l in 1:length(lat)){
    k1[l] = sum(cand_clust[[i]][cand_clust[[i]]$lat==lat[l],"cand_clust"] == 1)
    k2[l] = sum(cand_clust[[i]][cand_clust[[i]]$lat==lat[l],"cand_clust"] == 2)
    k3[l] = sum(cand_clust[[i]][cand_clust[[i]]$lat==lat[l],"cand_clust"] == 3)
    k4[l] = sum(cand_clust[[i]][cand_clust[[i]]$lat==lat[l],"cand_clust"] == 4)
    k5[l] = sum(cand_clust[[i]][cand_clust[[i]]$lat==lat[l],"cand_clust"] == 5)
  }
  pie_loc[[i]]<- cbind.data.frame(
    unique(cand_clust[[i]]$lat), unique(cand_clust[[1]]$lon),
    unlist(k1), unlist(k2), unlist(k3), unlist(k4), unlist(k5)
    )
  colnames(pie_loc[[i]]) <- c("lat","lon","k1","k2","k3","k4","k5")
  # write to file
  write.csv(pie_loc[[i]], paste0(curr_dir,'/',cand_name[i],'_mclust_pie.csv'))
}
#
}

# clean-up
rm()

# save session
session::save.session(session_path("17"), compress = "xz")
