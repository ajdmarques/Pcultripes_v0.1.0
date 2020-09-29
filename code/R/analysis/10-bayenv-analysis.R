source('code/R/functions/session_path.R')
# restore session
session::restore.session(session_path("07"))

# load parameters
bayenv_parameters <- "code/parameters/bayenv.toml" %>%
  RcppTOML::parseTOML() %>%
  `[[`(MODE)

# load functions
source("code/R/functions/gl2bayenv.R")

# read data
snp_ol_data <- readRDS(snp_ol_path)
# read environmental data
env_data <- read.csv('data/intermediate/vif-results/env.csv', row.names = 1)

# create directory to save bayenv results
unlink("data/intermediate/bayenv-results", force = TRUE, recursive = TRUE)
dir.create("data/intermediate/bayenv-results", showWarnings = FALSE,
           recursive = TRUE)
# create directory to save Bayenv results
curr_dir <- paste0("data/intermediate/bayenv-results/",
                   spp_parameters$species_atlas_names[i])
dir.create(curr_dir, showWarnings = FALSE, recursive = TRUE)

# write Gelight into Bayenv output
bayenv_input <- lapply(seq_along(snp_ol_data), function(i) {
  ## save data in Bayenv format
  gl2bayenv(snp_ol_data[[i]], paste0(curr_dir, "/data.txt"))
  read.delim(paste0(curr_dir, "/data.txt"), header = F, sep = "\t")
})

# generate covariance matrix
bayenv_matrix <- lapply(seq_along(snp_ol_data), function(i) {
  # call 'bayenv2' in bash
  bayenv2_path <- normalizePath(paste0('code/bayenv/','bayenv2')) 
  # call input file
  bayenv_input_path <- normalizePath(paste0(curr_dir, "/data.txt"))
  # call parameters p, k, r
  p = length(colnames(bayenv_input[[1]]))
  k = bayenv_parameters$matrix_k
  r = .Random.seed[i]
  # define output
  matrix_out_path <- normalizePath(paste0(getwd(),"/",curr_dir, "/matrix.out"))
  matrix_in_path <- normalizePath(paste0(getwd(),"/",curr_dir, "/matrix.in"))
  # initiate bayenv
  system(paste0(bayenv2_path, 
    " -i ", bayenv_input_path,
    " -p ", p,
    " -k ", k,
    " -r ", r,
    " > ", matrix_out_path, ""), intern = TRUE)
  # write matrix tail to input
  system(paste0("tail -n ", p+1, " ", matrix_out_path, 
                " > ", matrix_in_path))
  # read output
  x <- read.delim(matrix_in_path, header = F, sep = "\t")
  x <- x[1:p]
})

# generate environmental input
bayenv_env <- lapply(seq_along(snp_ol_data), function(i) {
  # match variables to genotypes
  x <- env_data[env_data$. %in% snp_ol_data[[1]]$ind.names,]
    x <- x[order(x[,1]),] 
    colnames(x)[1] <- "Sample-ID"
  # order snp and environment data
  y <- dartR::gl2demerelate(snp_ol_data[[1]])
    y <- y[,c("Sample-ID", "Population")]
    y <- y[order(y[,1]),]
    # join pop factor to environment table
    y <- merge(x,y, by = "Sample-ID")
  # scale variables then aggregate variables by population
  t(aggregate(scale(y[,4:(length(colnames(y))-1)]), by = list(y$Population), FUN = mean))
})
# write environment input file
write.table(bayenv_env[[1]][-1,], paste0(curr_dir,"/bayenv_env.txt"), 
            sep = '\t', row.names = F, col.names = F, quote = F)

# move 'calc_bf.sh' and call
calc_bf_path <- normalizePath(paste0('code/bayenv/','calc_bf.sh'))
system(paste0("cp ", calc_bf_path, " ", getwd(), "/", curr_dir),
       intern = TRUE)
# move 'bayenv2' and call
bayenv2_path <- normalizePath(paste0('code/bayenv/','bayenv2'))
system(paste0("cp ", bayenv2_path, " ", getwd(), "/", curr_dir),
       intern = TRUE)

# run bayenv analysis
bayenv_results <- lapply(seq_along(snp_ol_data), function(i) {
  # call programs
  calc_bf_path <- normalizePath(paste0(curr_dir,'/calc_bf.sh'))
  bayenv2_path <- normalizePath(paste0(curr_dir,'/bayenv2'))
  # call bayen input
  bayenv_input_path <- normalizePath(paste0(curr_dir, "/data.txt"))
  # call environment input
  bayenv_env_path <- normalizePath(paste0(curr_dir, "/bayenv_env.txt"))
  # call matrix
  bayenv_matrix_path <- normalizePath(paste0(curr_dir, "/matrix.in"))
  # call paraeters p, k, n
  p = length(colnames(bayenv_input[[1]]))
  k = bayenv_parameters$matrix_k
  n = length(rownames(bayenv_env[[1]]))-1
  # initiate bayenv
  species_path <- getwd()
  setwd(paste0(getwd(), "/", curr_dir))
  system(paste0("./calc_bf.sh", " ",
                "data.txt", " ",
                "bayenv_env.txt", " ",
                "matrix.in", " ",
                p, " ",
                k, " ",
                n), intern = TRUE)
  setwd(species_path)
  # read input
  read.table(paste0(curr_dir, "/bf_environ.bayenv_env.txt"), 
             header = F, 
             sep = "\t", 
             row.names = snp_ol_data[[1]]$loc.names, 
             colClasses = c(rep("NULL", 1),rep("numeric", n)))
  saveRDS(bayenv_results, 'data/final/bayenv_cand.rds', compress = 'xz')
  }) 

# interpret results
{
  # perform PCA on variables
  PC <- prcomp(bayenv_results[[1]])
  jpeg(paste0(curr_dir, "/bayenv_pc.jpeg")) 
  plot(PC, type = 'l')
  dev.off()
  # detect outlier based on first 2 PCs
  mod <- lm(PC$x[,1] ~ PC$x[,2], data = as.data.frame(x = PC$x, df1 = 6, df2 = 6))
  # creat object 'bayenv_results' with outliers defined
  bf_outlier <- car::outlierTest(mod)
  bf_outlier <- t(as.data.frame(do.call(rbind, bf_outlier)))
  # project outliers
  jpeg(paste0(curr_dir, "/bayenv_manhattan.jpeg"))
  plot(PC$x[,1],
       pch = 19, 
       cex = .5, 
       xlab = "SNP", ylab = "PC1",
       col=ifelse(rownames(PC$x) %in% rownames(bf_outlier), "blue", "grey"),
       main = "Bayenv")
  dev.off()
  # add outlier data to results
  bayenv_results[[1]]$outlier <- ifelse(rownames(bayenv_results[[1]]) %in% rownames(bf_outlier),
                                        "TRUE", "NA")
  colnames(bayenv_results[[1]]) <- c(paste0("PC",1:(length(rownames(bayenv_env[[1]]))-1)),
                                     "outlier")
}

# cleanup
rm(snp_ol_data)

# save session
session::save.session(session_path("10"), compress = "xz")
