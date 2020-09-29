source('code/R/functions/session_path.R')
# restore session
session::restore.session(session_path("10"))

# load parameters
sambada_parameters <- "code/parameters/sambada.toml" %>%
  RcppTOML::parseTOML() %>%
  `[[`(MODE)

# packages
library(LEA)
library(dplyr)

# load functions
source("code/R/functions/*.R")

# read data
snp_pop_admix_path <- 'data/intermediate/snp_admix_pop.rds'
snp_admix_data <- readRDS(snp_pop_admix_path)
# read environmental data
env_data <- read.csv('data/intermediate/vif-results/env.csv', row.names = 1)
env_data <- env_data[!duplicated(env_data$.),]

# create directory to save SamBada results
unlink("data/intermediate/sambada-results", force = TRUE, recursive = TRUE)
dir.create("data/intermediate/sambada-results", showWarnings = FALSE,
           recursive = TRUE)
# create directory to save Bayenv results
curr_dir <- paste0("data/intermediate/sambada-results/",
                   spp_parameters$species_atlas_names[i])
dir.create(curr_dir, showWarnings = FALSE, recursive = TRUE)

# Generate input files
{
# Write genotypes in the lfmm format
  dartR::gl2plink(x = snp_admix_data[[1]], 
         outfile = "snp_plink.ped", 
         outpath = curr_dir)
  x <- read.csv(paste0(curr_dir,"/snp_plink.ped"), header = T, row.names = 1)
  x[(x == -9)] = 9
  x <- x[order(as.numeric(rownames(x))),]
  write.lfmm(x, paste0(curr_dir,"/genotypes.lfmm"))
# Create a pcaProject object: pc
  pc = pca(paste0(curr_dir,"/genotypes.lfmm"), scale = F)
# Import eigenvectors from first two PCs
  eigenvect <- read.delim(file = paste0(curr_dir,"/genotypes.pca/genotypes.eigenvectors"), 
                           sep = " ", header = F )
# match environmental data to genotypes
  z <- env_data[env_data$. %in% as.numeric(snp_admix_data[[1]]$ind.names),]
  z <- z[order(z[,1]),]
# append eigenvectors to environmental data
sambada_env <- z[,3] %>%
  cbind(z[,2]) %>%
  cbind(z[,{4:length(colnames(env_data))}]) %>%
  cbind(eigenvect[,1])
  colnames(sambada_env)[length(colnames(env_data))] <- "pop"
  colnames(sambada_env)[1:2] <- c("longitude","latitude")
  sambada_env_path <- paste0(curr_dir,"/sambada_env")
  write.table(sambada_env[], sambada_env_path, row.names = z[,1], sep = "\t")
  # need to add "id" column to 
  system(noquote(paste0("sed -i '1 s/^/\"id\" /' ",sambada_env_path)), intern = TRUE)
# convert to matrix file for samabada
sambada_snp <- as.matrix(snp_admix_data[[1]])
  sambada_snp <- sambada_snp[order(as.numeric(rownames(sambada_snp))),]
  sambada_snp_path <- paste0(curr_dir,"/snp_data")
  write.table(sambada_snp, sambada_snp_path, sep = "\t")
  system(noquote(paste0("sed -i '1 s/^/\"id\" /' ", sambada_snp_path)), intern = TRUE)
}

# Initiate sambada in bash
sambada_params <- {
  # set working directory 
  species_path <- getwd()
    setwd(curr_dir)
  # create parameter file
  system(noquote(paste0("echo ","'",
                        sambada_parameters$shape,
                        "'", " >> param.txt")), intern = TRUE)
  system(noquote(paste0("echo ","'",
                        sambada_parameters$header,
                        "'", " >> param.txt")), intern = TRUE)
  system(noquote(paste0("echo ","'",
                        sambada_parameters$delim, " \" \"",
                        "'", " >> param.txt")), intern = TRUE)
  system(noquote(paste0("echo ","'",
                        sambada_parameters$n_env, length(colnames(sambada_env))+1,
                        "'", " >> param.txt")), intern = TRUE)
  system(noquote(paste0("echo ","'", 
                        sambada_parameters$n_mark, length(colnames(sambada_snp))+1,
                        "'", " >> param.txt")), intern = TRUE)
  system(noquote(paste0("echo ","'", 
                        sambada_parameters$n_ind, length(row.names(sambada_snp)),
                        "'", " >> param.txt")), intern = TRUE)
  system(noquote(paste0("echo ","'", 
                        sambada_parameters$ind_id, " \"id\"",
                        "'", " >> param.txt")), intern = TRUE)
  system(noquote(paste0("echo ","'", 
                        sambada_parameters$spatial,
                        "'", " >> param.txt")), intern = TRUE)
  system(noquote(paste0("echo ","'", 
                        sambada_parameters$autocor,
                        "'", " >> param.txt")), intern = TRUE)
  system(noquote(paste0("echo ","'", 
                        sambada_parameters$pop,
                        "'", " >> param.txt")), intern = TRUE)
  system(noquote(paste0("echo ","'", 
                        sambada_parameters$dim,
                        "'", " >> param.txt")), intern = TRUE)
  system(noquote(paste0("echo ","'", 
                        sambada_parameters$save,
                        "'", " >> param.txt")), intern = TRUE)
  setwd(species_path)
  # read parameter file
  sambada_params_path <- paste0(curr_dir,"/param.txt")
  read.delim(sambada_params_path)
}

# Run SamBada 
{  
  # path to sambada
  sambada_path <- 'code/sambada/binaries/sambada'
  system(noquote(paste0("cp ", "./", sambada_path, " ./",curr_dir)), intern = TRUE)  
  # set directory for output
  setwd(curr_dir)
  # initiate Samabda
  system(paste0("./sambada param.txt sambada_env snp_data"), intern = TRUE)
  # retain models that control for population
  system(paste0("head -1 snp_data-Out-2 > results.out"), intern = TRUE)
  system(paste0("grep 'pop' snp_data-Out-2 >> results.out"), intern = TRUE)
# save results
sambada_results <- read.table("results.out", header = T, sep = " ")
  setwd(species_path)
  saveRDS(sambada_results, 'data/final/sambada_cand.rds', compress = 'xz')
}

{
# Assign alpha with Bonferroni correction
pVal.threshold <- sambada_parameters$alpha/
  (length(colnames(sambada_env))*
     length(colnames(sambada_snp)))
# Computing the p-values
p.values <- cbind(sambada_results, pvalG=pchisq(sambada_results$GscorePop, 1, lower.tail=F), 
                 pvalWald=pchisq(sambada_results$WaldScorePop, 1, lower.tail=F))
  # Selecting significant models for the G score
  p.values[p.values$pvalG<pVal.threshold,]
  # Selecting significant models for the Wald score
  p.values[p.values$pvalWald<pVal.threshold,]
  # add outlier data to results
  sambada_results$G.outlier <- ifelse(sambada_results$Marker %in% 
                                    p.values[p.values$pvalG<pVal.threshold,]$Marker,
                                      "TRUE", "NA")
  sambada_results$Wald.outlier <- ifelse(sambada_results$Marker %in% 
                                      p.values[p.values$pvalWald<pVal.threshold,]$Marker,
                                    "TRUE", "NA")
}

# cleanup
rm(snp_admix_data)

# save session
session::save.session(session_path("11"), compress = "xz")
