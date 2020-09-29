source('code/R/functions/session_path.R')
## restore session
session::restore.session(session_path("11"))

## load parameters
lfmm_parameters <- "code/parameters/lfmm.toml" %>%
  RcppTOML::parseTOML() %>%
  `[[`(MODE)

## packages
library(lfmm)
library(adegenet)
library(LEA)

## load functions
#source("code/R/functions/*.R")

{## read data
snp_admix_data <- readRDS(snp_ol_admix_path)
## read environmental data
env_data <- read.csv('data/intermediate/vif-results/env.csv', row.names = 1)
env_data <- env_data[!duplicated(env_data$.),]
env_data <- env_data[env_data$. %in% snp_admix_data[[1]]$ind.names,]
}

## create directory to save LFMM results
unlink("data/intermediate/lfmm-results", force = TRUE, recursive = TRUE)
dir.create("data/intermediate/lfmm-results", showWarnings = FALSE,
           recursive = TRUE)
## create directory to save LFMM results
curr_dir <- paste0("data/intermediate/lfmm-results/",
                   spp_parameters$species_atlas_names[i])
dir.create(curr_dir, showWarnings = FALSE, recursive = TRUE)

{## Convert Genlight into matrix object
snp_mx <- as.matrix(snp_admix_data[[1]])
snp_scale <- apply(snp_mx, 2, function(x)
  ## Some scaling values are null, corresponding alleles are removed.
  replace(x, is.na(x), as.numeric(names(which.max(table(x)))))) 
}  
{## set the plot so that all figures per variable are joined
par(mfrow=c(1, 1))
## generate pca and select the number of axese, e.g. 3
PC <- prcomp(snp_scale)
plot(PC$sdev[1:20]^2, xlab = 'PC', ylab = "Variance explained")
points(3,PC$sdev[3]^2, type = "h", lwd = 3, col = "blue")
}

{## Run LFMM
  ## ridge estimates minimize regularized least-squares problem with an L2 penalty
  y = snp_scale
  ## Change environmnetal variable 
  x <- env_data[, 4:length(colnames(env_data))] 
  ## scale environmental variable
  x <- apply(x, 2, scale)
  ## set the number of principle components based on broken stick
  k = sum(PC$sdev[1:20]^2 >= mean(PC$sdev[1:20]^2))
  ## run cross-validation to determine suitable lambda
  curr_cv <-
    lfmm::lfmm_ridge_CV(Y = y, X = x, K = k,
                        n.fold.row = lfmm_parameters$row, 
                        n.fold.col = lfmm_parameters$col, 
                        lambdas =  lfmm_parameters$lambda) %>% 
    dplyr::group_by(lambda) %>%
    dplyr::summarise(error_mean = mean(err), error_sd = sd(err)) %>%
    dplyr::ungroup() %>%
    tibble::as_tibble()
  ## Fit an LFMM, i.e, compute B, U, V estimates
  snp_lfmm <- lfmm_ridge(Y = y, 
                         X = x,
                         lambda = curr_cv$lambda[which.min(curr_cv$error_mean)],
                         it.max = lfmm_parameters$itr,
                         K = k)
  ## performs association testing using the fitted model
  p.values <- lfmm_test(Y = y,
                        X = x,
                        lfmm = snp_lfmm,
                        calibrate = "gif") 
  pvalues <- p.values$calibrated.pvalue
  pVal.threshold <- lfmm_parameters$alpha
  al = apply(!pvalues < pVal.threshold, 1, any)
  av = apply(pvalues < pVal.threshold, 2, sum)
  ## Display quartiles
  par(mfrow=c(2, length(colnames(pvalues))/2))
  for (i in colnames(pvalues)){
    qqplot(rexp(length(pvalues[,i]), rate = log(10)),
         -log10(pvalues), xlab = "Expected quantile",
         pch = 19, cex = .4)
    abline(0,1)}
  ## identify candidate loci according to pvalue
  f1 <- function(num) {
    m <- if_else(num <= pVal.threshold,true = TRUE, false = FALSE)
    return(m)
  }
  sig.set <-
    data.frame(
    snp_admix_data[[1]]$loc.names, 
    apply(pvalues, 2, f1))
  colnames(sig.set) <- c('snp', colnames(pvalues))
  write.csv(sig.set, paste0(curr_dir,"/lfmm_results.csv"))
  ## Manhattan plot
  #par(mfrow=c(2, length(colnames(pvalues))/2))
  #for (i in colnames(pvalues)){
  #  plot(-log10(pvalues[,i]),
  #       pch = 19,
  #       cex = .2,
  #       xlab = "SNP", ylab = "-Log P",
  #       ylim = c(0, max(-log10(pvalues))),
  #       col = "grey",
  #       main = i)
  #  points(match(sig.set[sig.set[,i]==T,'snp'], rownames(pvalues)), 
  #         -log10(pvalues[which(pvalues[,i] <= pVal.threshold)]), 
  #         pch = 19, cex = 0.5,
  #         type = "p", 
  #         col = "blue")}
}

# cleanup
rm(snp_admix_data, snp_admix_pop_data, snp_lfmm, snp_ol_data, snp_pop_prob_data, snp_raw_pop_data,
   y, bayescan_results, p.values, PC, pcadapt_inital_run, pcadapt_results, pvalues, sambada_results,
   sambada_snp, selestim_results, snp_mx, snp_scale, test)

# save session
session::save.session(session_path("12"), compress = "xz")
