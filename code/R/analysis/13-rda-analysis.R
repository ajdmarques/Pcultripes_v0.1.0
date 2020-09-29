## source: https://popgen.nescent.org/2018-03-27_RDA_GEA

source('code/R/functions/session_path.R')
## restore session
session::restore.session(session_path("12"))

## Load packages
## -------------
library(psych)    ## Used to investigate correlations among predictors
library(vegan)    ## Used to run RDA

## load parameters
rda_parameters <- "code/parameters/rda.toml" %>%
  RcppTOML::parseTOML() %>%
  `[[`(MODE)

## create directory to save LFMM results
unlink("data/intermediate/rda", force = TRUE, recursive = TRUE)
dir.create("data/intermediate/rda", showWarnings = FALSE,
           recursive = TRUE)
## create directory to save LFMM results
curr_dir <- paste0("data/intermediate/rda/",
                   spp_parameters$species_atlas_names[i])
dir.create(curr_dir, showWarnings = FALSE, recursive = TRUE)

{## read data
  snp_admix_data <- readRDS(snp_ol_admix_path)
  snp_mx <- as.matrix(snp_admix_data[[1]])
  #snp_gi_data <- dartR::gl2gi(snp_admix_data[[1]])
  #dartR::gl2plink(snp_admix_data[[1]], '/snp_admix.raw', curr_dir)
  #snp_plink_data <- read.PLINK(paste0(curr_dir,'/snp_admix.raw'))
  ## read environmental data
  env_data <- read.csv('data/intermediate/vif-results/env.csv', row.names = 1)
  env_data <- env_data[!duplicated(env_data$.),]
  env_data <- env_data[env_data$. %in% snp_admix_data[[1]]$ind.names,]
  colnames(env_data[4:length(colnames(env_data))]) <- rda_parameters
  ## Confirm that genotypes and environmental data are in the same order
  identical(as.factor(rownames(snp_mx)), as.factor(env_data[,1]))
  #identical(rownames(snp_gi_data$tab), env_data[,1]) 
  ## order according to genotype id
  snp_mx <- snp_mx[order(as.numeric(rownames(snp_mx))),]
  env_data <- env_data[match(rownames(snp_mx),env_data$.),]
  #snp_gi_data <- snp_gi_data$tab[order(rownames(snp_gi_data$tab)),]
  #env_data <- env_data[match(rownames(snp_gi_data),env_data$.),]
}

{## RDA requires complete data frames (i.e., no missing data). 
  sum(is.na(snp_mx))
## Using the most common genotype at each SNP across all individuals.
snp_inp <- apply(snp_mx, 2, function(x)
  replace(x, is.na(x), as.numeric(names(which.max(table(x))))))
  sum(is.na(snp_inp)) # No NAs
## Perform RDA
snp.rda <- rda(snp_inp ~ ., data = env_data[,4:length(colnames(env_data))], scale=T)
## adj.r.squared is % variaiton described by RDA 
  RsquareAdj(snp.rda)
## eigenvalues of constrained axes reflect the variance explained by each canonical axis:
  summary(eigenvals(snp.rda, model = "constrained"))
  screeplot(snp.rda)
## Test RDA for significance using formal tests using F-statistics
signif.full <- anova.cca(snp.rda, parallel=getOption("mc.cores")) # default is permutation=999
  signif.full
## check Variance Inflation Factors for the predictor variables
  vif.cca(snp.rda) # retain values < 10 or better < 5
}

{## Plot RDA
## Apply lineage data
  levels(snp_admix_data[[1]]$pop) <- c("East","West")
lineage <- snp_admix_data[[1]]$pop
bg <- c("#ff7f00","#1f78b4")
  par(mfrow=c(1,1))
  ## axes 1 & 2
  jpeg(paste0(curr_dir,"/lineage","_1-2",".jpg"), width = 500, height = 500)
  plot(snp.rda, type="n", scaling=3)
  points(snp.rda, display="species", pch=20, cex=0.7, col="gray32", scaling=3)           # the SNPs
  points(snp.rda, display="sites", pch=21, cex=1.3, col="gray32", scaling=3, bg=bg[lineage]) # the wolves
  text(snp.rda, scaling=3, display="bp", col="#0868ac", cex=1)                           # the predictors
  legend("topleft", legend=levels(lineage), bty="n", col="gray32", pch=21, cex=1, pt.bg=bg)
  dev.off()
  ## axes 1 & 3
  jpeg(paste0(curr_dir,"/lineage","_1-3",".jpg"), width = 500, height = 500)
  plot(snp.rda, type="n", scaling=3, choices=c(1,3))
  points(snp.rda, display="species", pch=20, cex=0.7, col="gray32", scaling=3, choices=c(1,3))           # the SNPs
  points(snp.rda, display="sites", pch=21, cex=1.3, col="gray32", scaling=3, bg=bg[lineage], choices=c(1,3)) # the wolves
  text(snp.rda, scaling=3, display="bp", col="#0868ac", cex=1, choices=c(1,3))                           # the predictors
  legend("topleft", legend=levels(lineage), bty="n", col="gray32", pch=21, cex=1, pt.bg=bg)
  dev.off()
}

{## Apply cluster data
  snp_metadata[[1]] <- snp_metadata[[1]][!duplicated(snp_metadata[[1]]$order_id),]
  snp_metadata[[1]] <- snp_metadata[[1]][snp_metadata[[1]]$order_id %in% snp_admix_data[[1]]$ind.names,]
  snp_metadata[[1]] <- snp_metadata[[1]][match(snp_admix_data[[1]][order(snp_admix_data[[1]]$ind.names),]$ind.names,
                                               snp_metadata[[1]]$order_id),]
  snp_metadata[[1]] <- snp_metadata[[1]][order(as.numeric(snp_metadata[[1]]$order_id)),]
  snp_clust <- as.factor(snp_metadata[[1]]$cluster_id)
  levels(snp_clust) <- rda_parameters$env
  #levels(snp_clust) <- c('PT.N.Coast','Alicante','Barcelona','Orense','Madrid',
   #                      'Salamanca','Salamanca.PT','Zaragosa','PT.N.Interior','Murcia',
    #                     'Cuidade Real','Catalan.Coast','Extramendura','Cadiz Coast','PT.S.Interior','PT.S.Coast')
  bg <- rda_parameters$col
  #bg <- c('#a60ce3','#a6cee3','#1f78b4','#b2df8a','#33a02c',
   #       '#fb0299','#fb9a99','#e31a1c','#fdbf6f','#ff7f00',
    #      '#44f665','#cab2d6','#6a3d9a','#ffff99','#b1f628','#b15928')
  par(mfrow=c(1,1))
  ## axes 1 & 2
  jpeg(paste0(curr_dir,"/cluster","_1-2",".jpg"), width = 800, height = 800)
  plot(snp.rda, type="n", scaling=3)
  points(snp.rda, display="species", pch=20, cex=0.7, col="gray32", scaling=3)           # the SNPs
  points(snp.rda, display="sites", pch=21, cex=1.3, col="gray32", scaling=3, bg=bg[snp_clust]) # the wolves
  text(snp.rda, scaling=3, display="bp", col="#0868ac", cex=1)                           # the predictors
  legend("topleft", legend=levels(snp_clust), bty="n", col="gray32", pch=21, cex=1, pt.bg=bg)
  dev.off()
  ## axes 1 & 3
  jpeg(paste0(curr_dir,"/cluster","_1-3",".jpg"), width = 800, height = 800)
  plot(snp.rda, type="n", scaling=3, choices=c(1,3))
  points(snp.rda, display="species", pch=20, cex=0.7, col="gray32", scaling=3, choices=c(1,3))           # the SNPs
  points(snp.rda, display="sites", pch=21, cex=1.3, col="gray32", scaling=3, bg=bg[snp_clust], choices=c(1,3)) # the wolves
  text(snp.rda, scaling=3, display="bp", col="#0868ac", cex=1, choices=c(1,3))                           # the predictors
  legend("topleft", legend=levels(snp_clust), bty="n", col="gray32", pch=21, cex=1, pt.bg=bg)
  dev.off()
}

{## Species scores for the first three constrained axes
load.rda <- vegan::scores(snp.rda, choices=c(1:3), display="species")  
  jpeg(paste0(curr_dir,"/loadings.jpg"), width = 800, height = 300)
  par(mfrow=c(1,3))
  hist(load.rda[,1], main="Loadings on RDA1")
  hist(load.rda[,2], main="Loadings on RDA2")
  hist(load.rda[,3], main="Loadings on RDA3") 
  dev.off()}

{## Define function for outliers, where x = vector of loadings and z = number of standard deviations
outliers <- function(x,z){
  lims <- mean(x) + c(-1, 1) * z * sd(x)     # find loadings +/-z sd from mean loading     
  x[x < lims[1] | x > lims[2]]               # locus names in these tails
}
  cand1 <- outliers(load.rda[,1],3)
  cand2 <- outliers(load.rda[,2],3)
  cand3 <- outliers(load.rda[,3],3)
  ncand <- length(cand1) + length(cand2) + length(cand3)
  ncand
  ## Organize results in data frame with the axis, SNP name, loading, & correlation with each predictor:
  cand1 <- cbind.data.frame(rep(1,times=length(cand1)), names(cand1), unname(cand1))
  cand2 <- cbind.data.frame(rep(2,times=length(cand2)), names(cand2), unname(cand2))
  cand3 <- cbind.data.frame(rep(3,times=length(cand3)), names(cand3), unname(cand3))
  colnames(cand1) <- colnames(cand2) <- colnames(cand3) <- c("axis","snp","loading")
cand <- rbind(cand1, cand2, cand3)
  #cand$snp <- as.character(cand$snp)
  ## Add correlations of each candidate SNP with the eight environmental predictors:
  foo <- matrix(nrow=(ncand), ncol=(length(colnames(env_data))-3))  
  colnames(foo) <- colnames(env_data[,4:length(colnames(env_data))])
for (i in 1:length(cand$snp)) {
  nam <- cand[i,2]
  snp.gen <- snp_inp[,nam]
  foo[i,] <- apply(env_data[,4:length(colnames(env_data))],2,function(x) cor(x,snp.gen))
}
cand <- cbind.data.frame(cand,foo)  
  head(cand)
  ## remove duplicates
  length(cand$snp[duplicated(cand$snp)])
  foo <- cbind(cand$axis, duplicated(cand$snp)) 
  ## number of duplicates per RDA
  table(foo[foo[,1]==1,2])
  table(foo[foo[,1]==2,2])
  table(foo[foo[,1]==3,2])
## remove duplicates
cand <- cand[!duplicated(cand$snp),]
rda_cand <- cand
  write.table(noquote(rda_cand), paste0(curr_dir,'/rda_candidates'), sep = "\t")
}

{# Which of the predictors each candidate SNP is most strongly correlated with:
  #cand = cand[,-length(colnames(cand))]
  c.pred <- length(colnames(cand))+1
  c.corr <- length(colnames(cand))+2
for (i in 1:length(cand$snp)) {
  bar <- cand[i,]
  cand[i,c.pred] <- names(which.max(abs(bar[4:(c.pred-1)]))) # gives the variable
  cand[i,c.corr] <- max(abs(bar[4:(c.corr-2)]))              # gives the correlation
}
  # rename columns
  colnames(cand)[c.pred] <- "predictor"
  colnames(cand)[c.corr] <- "correlation"
  # frequency of predictors among outliers
  table(cand$predictor)
}

{## Plot the SNPs
  sel <- cand$snp
  env <- cand$predictor
  bg <- RColorBrewer::brewer.pal(length(colnames(rda_cand)),"Set3")
  for (i in 4:length(colnames(rda_cand)) ){
    env[env == colnames(rda_cand)[i]] <- bg[i]
  }
  #env[env=="bio14"] <- '#e41a1c'
  #env[env=="bio2"] <- '#377eb8'
  #env[env=="bio3"] <- '#4daf4a'
  #env[env=="bio8"] <- '#984ea3'
  #env[env=="bio9"] <- '#ff7f00'
  #env[env=="elevate"] <- '#ffff33'
  ## color by predictor:
  col.pred <- rownames(snp.rda$CCA$v) # pull the SNP names
for (i in 1:length(sel)) {           # color code candidate SNPs
  foo <- match(sel[i],col.pred)
  col.pred[foo] <- env[i]
}
  col.pred[grep("/",col.pred)] <- '#f1eef6' # non-candidate SNPs
  empty <- col.pred
  empty[grep("#f1eef6",empty)] <- rgb(0,1,0, alpha=0) # transparent
  empty.outline <- ifelse(empty=="#00FF0000","#00FF0000","gray32")
}
{par(mfrow=c(1,1))
  # axes 1 & 2
  jpeg(paste0(curr_dir,"/cand_snp","_1-2",".jpg"), width = 800, height = 800)
  plot(snp.rda, type="n", scaling=3, xlim=c(-1,1), ylim=c(-1,1))
  points(snp.rda, display="species", pch=21, cex=1, col="gray32", bg=col.pred, scaling=3)
  points(snp.rda, display="species", pch=21, cex=1, col=empty.outline, bg=empty, scaling=3)
  text(snp.rda, scaling=3, display="bp", col="#0868ac", cex=1)
  legend("bottomright", legend=colnames(rda_cand)[4:length(colnames(rda_cand))],
         bty="n", col="gray32", pch=21, cex=1, pt.bg=bg)
  dev.off()
  # axes 1 & 3
  jpeg(paste0(curr_dir,"/cand_snp","_1-3",".jpg"), width = 800, height = 800)
  plot(snp.rda, type="n", scaling=3, xlim=c(-1,1), ylim=c(-1,1), choices=c(1,3))
  points(snp.rda, display="species", pch=21, cex=1, col="gray32", bg=col.pred, scaling=3, choices=c(1,3))
  points(snp.rda, display="species", pch=21, cex=1, col=empty.outline, bg=empty, scaling=3, choices=c(1,3))
  text(snp.rda, scaling=3, display="bp", col="#0868ac", cex=1, choices=c(1,3))
  legend("bottomright", legend=colnames(rda_cand)[4:length(colnames(rda_cand))],
         bty="n", col="gray32", pch=21, cex=1, pt.bg=bg)
  dev.off()
}

## save results
snp_mx_path <- 'data/intermediate/snp_pop_admix.matrix'
saveRDS(snp_mx, snp_mx_path)

# cleanup
rm(snp_admix_data, snp_gi_data, snp_gi_imp, snp_scale, snp_mx, load.rda, snp_inp, snp.rda)

# save session
session::save.session(session_path("13"), compress = "xz")
