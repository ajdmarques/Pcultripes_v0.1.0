#browseVignettes()
library(gradientForest)
require(gradientForest)
library(ggplot2)
library(reshape2)

source('code/R/functions/session_path.R')
## restore session
session::restore.session(session_path("13")) 

##
# Path to store results
gf_path <- "~/Documents/NGC/dart-adaptation-revised/data/intermediate/gradientforest-results/"
# Define current directory
curr_dir <- "~/Documents/NGC/dart-adaptation-revised/data/intermediate/gradientforest-results/"
# Load environmental data post VIF
env_path <- "~/Documents/NGC/dart-adaptation-revised/data/intermediate/vif-results/env.csv"
env_data <- read.csv(env_path)
#snp_adaptive_path <- "~/Documents/NGC/dart-adaptation-revised/data/intermediate/snp_neutral.rds"
#snp_adaptive_data <- readRDS(snp_adaptive_path)
snp_ol_admix_path <- "~/Documents/NGC/dart-adaptation-revised/data/intermediate/snp_ol_admix.rds"
snp_ol_admix_data <- readRDS(snp_ol_admix_path)
z <- c(107150,107151) # Remove samples from Cuidade Real
snp_ol_admix_data[[1]] <- snp_ol_admix_data[[1]][!snp_ol_admix_data[[1]]$ind.names %in% z]
# list of outliers identified by EA methods
outlier_path <- "~/Documents/NGC/dart-adaptation-revised/data/intermediate/ea_outlier.rds"
outlier_data <- readRDS(outlier_path)
for (i in seq_along(outlier_data)){
  outlier_data[[i]] <- outlier_data[[i]][!outlier_data[[i]]$ind.names %in% z]
}
# Names of tests
test_names = c("Bayescan", "SelEstim", "PCAdapt", "Bayenv2", "LFMM", "RDA", "Neutral")
for (i in seq_along(outlier_data)){
  print(paste0(test_names[i], " contains ", outlier_data[[i]]$n.loc, " SNPs"))
}
# samples included in environmental data
if (length(snp_ol_admix_data[[1]]$ind.names) > length(env_data$.)) {
  snp_ol_admix_data[[1]] = snp_ol_admix_data[[1]][snp_ol_admix_data[[1]]$ind.names %in% 
                                                    as.character(env_data$.)]       
} else {
  env_data = env_data[as.character(env_data$.) %in% 
                        snp_ol_admix_data[[1]]$ind.names,] 
}

# Convert genlight to matrix of alleles
outlier_matrix <- list()
for (i in seq_along(outlier_data)) {
  outlier_matrix[[i]] <- as.matrix(outlier_data[[i]])
  ## Using the most common genotype at each SNP across all individuals.
  outlier_matrix[[i]] <- apply(outlier_matrix[[i]], 2, function(x)
    replace(x, is.na(x), as.numeric(names(which.max(table(x))))))
  print(paste0(i, "Number of NA data ",sum(is.na(outlier_matrix[[i]])))) # No NAs
  outlier_matrix[[i]][order(as.numeric(rownames(outlier_matrix[[i]]))),]
}
# Sort env by SNP
env_data <- env_data[match(rownames(outlier_matrix[[i]]),env_data$.),]
# Introduce data to script
Sp_mat = list()
for (i in seq_along(outlier_data)) {
  Sp_mat[[i]] = outlier_matrix[[i]]
}
saveRDS(Sp_mat, paste0(curr_dir,"snp_matrix.rds"), compress = "xz")
#MODE = "debug"
if (MODE == "debug"){
  for (i in seq_along(Sp_mat)){
    if(dim(Sp_mat[[i]])[2] > 100)
      Sp_mat[[i]] = Sp_mat[[i]][,sort(sample(dim(Sp_mat[[i]])[2], 100))]
  }
  #env_data <- env_data[,1:5]
} else {
  for (i in seq_along(Sp_mat)){
    if(dim(Sp_mat[[i]])[2] > 1000)
      Sp_mat[[i]] = Sp_mat[[i]][,sort(sample(dim(Sp_mat[[i]])[2], 1000))]
  }
}
# To make SNP data compatible with gradientForest, SNPs must be renamed
# Example "44453974-64-G/A" -> "44453974.64.G.A"
for (i in seq_along(Sp_mat)){
  colnames(Sp_mat[[i]]) <- gsub('[-\\/]','.', colnames(Sp_mat[[i]]))
}


Phys_site = env_data[,-1:-2]
Phys_site[,c("bio14","bio8","bio9","bio2")] = Phys_site[,c("bio14","bio8","bio9","bio2")]/10
####

## create directory to save Venn diagrams results
unlink("data/intermediate/gradientforest-results/", force = TRUE, recursive = TRUE)
dir.create("data/intermediate/gradientforest-results/", showWarnings = FALSE,
           recursive = TRUE)
## create directory to save gradientForest results
curr_dir <- paste0("data/intermediate/gradientforest-results/",
                   spp_parameters$species_atlas_names,"/")
dir.create(curr_dir, showWarnings = FALSE, recursive = TRUE)

# load the site by species matrix, object=Sp_mat
# includes transformed biomass of 110 species at 197 sites
#.libPaths()
#gradientForest_doc <- "/home/adam/R/x86_64-pc-linux-gnu-library/3.6/gradientForest/doc/" 
#load(paste0(gradientForest_doc,"GZ.sps.mat.Rdata"))
for (i in seq_along(Sp_mat)){
  print(dim(Sp_mat[[i]]))
}

# load the site physical data, object=Phys_site
# includes 28 predictors at 197 sites
#load(paste0(gradientForest_doc,"GZ.phys.site.Rdata"))
dim(Phys_site)

# gradientForestis  a  wrapper  function  that  callsextendedForest, 
# a modified version of randomForest.
# extendedForest assesses the importance of each variable for prediction accuracy; 
# information that is further collated and processed by gradientForest.
nSites = list()
nSpecs = list()
for (i in seq_along(Sp_mat)){
  nSites[[i]] <- dim(Sp_mat[[i]])[1]
  nSpecs[[i]] <- dim(Sp_mat[[i]])[2]
}
# set depth of conditional permutation
lev = list()
for (i in seq_along(Sp_mat)){
  lev[[i]] <- floor(log2(nSites[[i]]*0.368/2)) 
}
lev

# The summary shows the number of species with positive R2
# ie. species  that  could  be  predicted  to  any  extent  by  the  available  predictor.
# The returned  object  is  a  list  containing  
# the  data,  predictor  importances,  species R2’s  
# and  other information described in the html help pages under Value.
gf = list()
for (i in seq_along(Sp_mat)){
  gf[[i]] <- gradientForest(cbind(Phys_site,Sp_mat[[i]]), 
                            predictor.vars=colnames(Phys_site), 
                            response.vars=colnames(Sp_mat[[i]]), 
                            ntree=500, transform = NULL, compact=T, nbin = 201, 
                            maxLevel=lev, corr.threshold=0.5)
}
gf
names(gf[[i]])

# The first is the predictor overall importance plot.  
# This show the mean accuracy importance and the mean importance weighted by species R2.
for (i in seq_along(Sp_mat)){
  plot(gf[[i]],plot.type="O")
}

# Write to dataframe.
importance_df <- matrix(nrow = length(Sp_mat), ncol = length(Phys_site))
rownames(importance_df) <- test_names
colnames(importance_df) <- names(Phys_site)
for (i in seq_along(Phys_site)){
  for (j in seq_along(Sp_mat)){
    importance_df[j,i] <- importance(gf[[j]])[names(Phys_site)[i]]
  }
}
# Identify outliers in each predictor
jpeg(paste0(curr_dir, "boxplot_importance.jpeg"))
OutVals <- boxplot(importance_df, 
                   xlab = "predictor", ylab = "R2 weighted importance")$out
dev.off()
# Identify if outlier are significant
# a Dixon test was used as it is suitable for n < 25 
dixon_p <- list()
for (i in seq_along(Phys_site)){
  dixon_p[[i]] <- paste("Predictor:", names(Phys_site)[i], 
                        "pvalue:", outliers::dixon.test(importance_df[,i])$p.value,
                        "outlier test:", names(outliers::dixon.test(importance_df[,i])$statistic))
}
dixon_p

# Identify the predictors to visualise
bioclim_name <- list("Latitude","Longitude","Elevation",
                     "Isothermality (Isotherm)", "Precipitation of Driest Month (PDM)",
                     "Mean Temperature of Wettest Quarter (MTWQ)",
                     "Mean Temperature of Driest Quarter (MTWQ)",
                     "Mean Diurnal Range (MDR)")
indx = 0
p = 0.05
important_plots <- list()
xlab_plot <- list()
for (i in seq_along(dixon_p)){
  if (outliers::dixon.test(importance_df[,i])$p.value < p){
    indx = indx + 1
    important_plots[[indx]] <- names(Phys_site)[i]
    xlab_plot[[indx]] <- bioclim_name[[i]]
  }
}

# Identify geospatial predictors
gps_predictors <- list()
for (i in seq_along(names(Phys_site)[1:2])){
  gps_predictors[[i]] <- names(Phys_site)[1:2][i]
}

# Identify climatic predictors
clim_predictors <- list()
for (i in seq_along(names(Phys_site)[-1:-2])){
  clim_predictors[[i]] <- names(Phys_site)[-1:-2][i]
}

# Load customized function
source("~/Documents/NGC/dart-adaptation-revised/code/R/functions/my.species.cumulative.plot.R")

######
### Results for Geospatial predictors
######
# Create a new object to save results to
cplot <- list()
for (i in seq_along(Sp_mat)){
  cplot[[i]] <- my.species.cumulative.plot(gf[[i]], imp.vars = gps_predictors)
}

# Create a list for each table 
cumulative_plots <- list()
for (j in seq_along(gps_predictors)){
  cumulative_plots[[j]] <- as.data.frame(cplot[[1]][gps_predictors[[j]]])
  for (i in seq_along(Sp_mat)[-1]){
    x <- as.data.frame(cplot[[i]][gps_predictors[[j]]])
    cumulative_plots[[j]] <- merge(cumulative_plots[[j]], x,
                                   by = names(cumulative_plots[[j]])[1], 
                                   all = T)
  }
  names(cumulative_plots[[j]]) <- c(gps_predictors[[j]], test_names)
}

# Basic line plot
d <- list()
p <- list()
for (j in seq_along(gps_predictors)){
  # melt variable into a ggplot compatible data frame
  d[[j]] <- melt(cumulative_plots[[j]], id.vars=gps_predictors[j])
  # reassigna  generic title for data frame
  names(d[[j]]) <- c("x", "test", "value")
  # Plot everything on the same plot
  p[[j]] <- ggplot(d[[j]], aes(x, value, col = test)) +
    geom_point() +
    stat_smooth() +
    labs(x = bioclim_name[1:2][j],
         y = "Cumulative R2")
  # save plot
  ggsave(paste0(curr_dir,"r2-",gps_predictors[[j]],".png"), 
         plot = last_plot())
}

# Save output
saveRDS(cplot, paste0(curr_dir,"/cumulative_geospatial.rds"), compress = "xz")

######
### Results for Climatic predictors
######
# Create a new object to save results to
cplot <- list()
for (i in seq_along(Sp_mat)){
  cplot[[i]] <- my.species.cumulative.plot(gf[[i]], imp.vars = clim_predictors)
}

# Create a list for each table 
cumulative_plots <- list()
for (j in seq_along(clim_predictors)){
  cumulative_plots[[j]] <- as.data.frame(cplot[[1]][clim_predictors[[j]]])
  for (i in seq_along(Sp_mat)[-1]){
    x <- as.data.frame(cplot[[i]][clim_predictors[[j]]])
    cumulative_plots[[j]] <- merge(cumulative_plots[[j]], x,
                                   by = names(cumulative_plots[[j]])[1], 
                                   all = T)
  }
  names(cumulative_plots[[j]]) <- c(clim_predictors[[j]], test_names)
}

# Basic line plot
d <- list()
p <- list()
for (j in seq_along(clim_predictors)){
  # melt variable into a ggplot compatible data frame
  d[[j]] <- melt(cumulative_plots[[j]], id.vars=clim_predictors[j])
  # reassigna  generic title for data frame
  names(d[[j]]) <- c("x", "test", "value")
  # Plot everything on the same plot
  p[[j]] <- ggplot(d[[j]], aes(x, value, col = test)) +
    geom_point() +
    stat_smooth() +
    labs(x = bioclim_name[-1:-2][j],
         y = "Cumulative R2")
  # save plot
  ggsave(paste0(curr_dir,"r2-",clim_predictors[[j]],".png"), 
         plot = last_plot())
}

# Save output
saveRDS(cplot, paste0(curr_dir,"/cumulative_climatic.rds"), compress = "xz")

######
### Results for significant predictors according to Dixon's P
######
# Create a new object to save results to
cplot <- list()
for (i in seq_along(Sp_mat)){
  cplot[[i]] <- my.species.cumulative.plot(gf[[i]], imp.vars = important_plots)
}

# Create a list for each table 
cumulative_plots <- list()
for (j in seq_along(important_plots)){
  cumulative_plots[[j]] <- as.data.frame(cplot[[1]][important_plots[[j]]])
  for (i in seq_along(Sp_mat)[-1]){
    x <- as.data.frame(cplot[[i]][important_plots[[j]]])
    cumulative_plots[[j]] <- merge(cumulative_plots[[j]], x,
                                   by = names(cumulative_plots[[j]])[1], 
                                   all = T)
  }
  names(cumulative_plots[[j]]) <- c(important_plots[[j]], test_names)
}

# Basic line plot
d <- list()
p <- list()
for (j in seq_along(important_plots)){
  # melt variable into a ggplot compatible data frame
  d[[j]] <- melt(cumulative_plots[[j]], id.vars=important_plots[j])
  # reassigna  generic title for data frame
  names(d[[j]]) <- c("x", "test", "value")
  # Plot everything on the same plot
  p[[j]] <- ggplot(d[[j]], aes(x, value, col = test)) +
    geom_point() +
    stat_smooth() +
    labs(x = xlab_plot[j], 
         cex.lab = .1,
         y = " ") 
  }

# Cumulative figure
#box_p <- stack(as.data.frame(importance_df))
box_p <- melt(t(as.data.frame(importance_df)))
bp <- ggplot(box_p, aes(x = Var1, y = value)) +
  geom_boxplot() +
  geom_dotplot(aes(x = Var1, y = value, fill = Var2), dotsize=1.5,
               binaxis='y', stackdir='center', stackgroups = TRUE, binpositions = "all") +
  labs(x="Predictor", y = "R2 weighted importance") +
  theme_classic()
bp
# <http://www.sthda.com/english/articles/24-ggpubr-publication-ready-plots/81-ggplot2-easy-way-to-mix-multiple-graphs-on-the-same-page/>
#jpeg("~/Documents/NGC/dart-adaptation-revised/data/final/cumulative_R2.jpeg")
ggpubr::ggarrange( 
  bp, labels = "a", label.x = .08,
  legend = F, nrow = 2, heights = c(1, 2),
  ggpubr::ggarrange(
    p[[2]],
    p[[1]],
    p[[3]],
    p[[4]],
  #rremove("x.text"),
  labels = c("b", "c", "d", "e"), # Assign the necessary number of titles
  common.legend = T, 
  ncol = ceiling(length(important_plots)/2),
  nrow = ceiling(length(important_plots)/2)))
#dev.off()
ggsave(paste0(curr_dir,"important-variables.png"), 
       plot = last_plot(), height = 7, width = 7.5)

# Save output
saveRDS(cplot, paste0(curr_dir,"/cumulative_importance.rds"), compress = "xz")

# Clean environment
rm(bayescan_outlier_loci, pcadapt_outlier_loci, selestim_outlier_loci,
   bayenv_results, selestim_trace_data, snp_ol_admix_data, 
   cand1, cand2, cand3, gf, outlier_data, outlier_matrix, 
   structure_trace_data)

# save session
session::save.session(session_path("14"), compress = "xz")
