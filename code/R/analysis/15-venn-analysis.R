source('code/R/functions/session_path.R')
# restore session
#session::restore.session(session_path("14"))
session::restore.session(session_path("14"))

# load parameters
venn_parameters <- "code/parameters/venn.toml" %>%
  RcppTOML::parseTOML() %>%
  `[[`(MODE)

# packages
library(VennDiagram)

# load functions
#source("code/R/functions/*.R")

{# read data
#  load('data/final/results.rda')
#  bayenv_results <- readRDS('data/final/bayenv_cand.rds')
#  bayescan_outlier_loci <-readRDS('data/intermediate/bayescan-results/bayescan_outlier_loci.rds')
}

# create directory to save Venn diagrams results
unlink("data/intermediate/venn-results", force = TRUE, recursive = TRUE)
dir.create("data/intermediate/venn-results", showWarnings = FALSE,
           recursive = TRUE)
# create directory to save LFMM results
curr_dir <- paste0("data/intermediate/venn-results/",
                   spp_parameters$species_atlas_names[i])
dir.create(curr_dir, showWarnings = FALSE, recursive = TRUE)

# Load outlier data
Sp_mat <- readRDS("~/Documents/NGC/dart-adaptation-revised/data/intermediate/gradientforest-results/snp_matrix.rds")

# Names of tests
test_names = c("Bayescan", "SelEstim", "PCAdapt", "Bayenv2", "LFMM", "RDA", "Neutral")
test_cand <- list()
for (i in seq_along(test_names)){
  test_cand[[i]] <- colnames(Sp_mat[[i]])
}
str(test_cand)

# Call list of outliers
#  {BS <- bayescan_outlier_loci
#    BS <- BS[[1]][BS[[1]]$outlier==TRUE,]$loc_name}
    #BS <- bayescan_outlier_loci[[1]][bayescan_outlier_loci[[1]]$outlier==TRUE,]$loc_name
#  {SE <- selestim_outlier_loci
#    SE <- SE[[1]][SE[[1]]$outlier==TRUE,]$loc_name}
    #SE <- selestim_outlier_loci[[1]][selestim_outlier_loci[[1]]$outlier==TRUE,]$loc_name 
#  PA <- pcadapt_outlier_loci[[1]][pcadapt_outlier_loci[[1]]$outlier==TRUE,]$loc_name
#  {bayenv_results <- readRDS("data/final/bayenv_cand.rds")
#    BE <- rownames(bayenv_results[bayenv_results$outlier==TRUE,]) 
#      BE <- as.character(BE)}
  #{sambada_results <- readRDS("data/final/sambada_cand.rds")
   # SB <- sambada_results[sambada_results$G.outlier==TRUE & sambada_results$Wald.outlier==TRUE,]$Marker 
    # SB <- as.character(SB)}
#  {LF <- lapply(seq_along(sig.set[,-1])+1, function(i){
#      sig.set[sig.set[,i]==T,]$snp }) 
#    LF <- as.character(unlist(LF))
#    LF <- unique(LF)}
#  {RD <- cand$snp
#    RD <-as.character(RD)}
  #MM <- msod_outlier_loci[[1]][msod_outlier_loci[[1]]$outlier==T,]$loc_name
  {DC <- read.csv("data/intermediate/minotaur/pelcul/minotaur_data.csv")
    DC <- as.character(DC[DC$dcms_outlier == T,]$snp)}
  
{# Generate Venn diagrams
  VennDiagram::venn.diagram(x = test_cand[1:3],
                            category.names = test_names[1:3],
                            filename = paste0(curr_dir,'/venn-OA.png'),
                            output = TRUE,
                            imagetype="png" ,
                            #height = 350 , 
                            #width = 350 , 
                            #resolution = 500,
                            #compression = "lzw",
                            #lwd = 2,
                            #lty = 'blank',
                            #lwd = c(length(SE)/length(c(SE, BS, PA))+.5,
                            #        length(BS)/length(c(SE, BS, PA))+.5,
                            #        length(PA)/length(c(SE, BS, PA))+.5),
                            #fill = c('blue', 'purple', 'green'),
                            #cex = 0.5,
                            #fontface = "bold",
                            #fontfamily = "sans",
                            #cat.cex = 0.4,
                            #cat.fontface = "bold",
                            #cat.default.pos = "outer",
                            #cat.pos = c(-20, 20, 180),
                            #cat.dist = c(0.09, 0.2, 0.05),
                            #cat.fontfamily = "sans",
                            #rotation = 1
  )
  VennDiagram::venn.diagram(x = test_cand[4:6],
                            category.names = test_names[4:6],
                            filename = paste0(curr_dir,'/venn-EAA.png'),
                            output = TRUE,
                            imagetype="png" ,
                            #height = 540 , 
                            #width = 560 , 
                            #resolution = 500,
                            #compression = "lzw",
                            #lwd = 2,
                            #lty = 'blank',
                            #fill = c('red', 'orange', 'yellow', 'pink'),
                            #cex = 0.5,
                            #fontface = "bold",
                            #fontfamily = "sans",
                            #cat.cex = 0.4,
                            #cat.fontface = "bold",
                            #cat.default.pos = "outer",
                            cat.pos = c(180, 180, 180),
                            cat.dist = c(0.02, 0.02, 0.04),
                            #cat.fontfamily = "sans",
                            #rotation = 1
  )
  # Method comparisons
  VennDiagram::venn.diagram(x = list(unique(unlist(test_cand[1:3])), 
                                     unique(unlist(test_cand[4:6])),
                                     DC
                                     ),
                            category.names = c("OA" , 
                                               "EAA",  
                                               "DCMS"
                                               ),
                            filename = paste0(curr_dir,'/venn-methods.png'),
                            output = TRUE,
                            imagetype="png" ,
                            #height = 500 , 
                            #width = 500 , 
                            #resolution = 500,
                            #compression = "lzw",
                            #lwd = 2,
                            #lty = 'blank',
                            #fill = c('purple', 'red', 'grey'),
                            #cex = .5,
                            #fontface = "bold",
                            #fontfamily = "sans",
                            #cat.cex = 0.4,
                            #cat.fontface = "bold",
                            #cat.default.pos = "outer",
                            #cat.pos = c(0, 0),
                            #cat.dist = c(0.08, 0.08),
                            #cat.fontfamily = "sans",
                            #rotation = 1
  )
}
'''  
# Overlaps with DCMS
  VennDiagram::venn.diagram(x = list(SE, BS, PA, DC),
                            category.names = c("SelEstim" , "Bayescan" , "PCAdapt", "DCMS"),
                            filename = paste0(curr_dir,'/venn-OA-DCMS.png'),
                            output = TRUE,
                            imagetype="png"
                            )
  VennDiagram::venn.diagram(x = list(BE, LF, RD, DC),
                            category.names = c("Bayenv2" , "LFMM", "RDA", "DCMS"),
                            filename = paste0(curr_dir,'/venn-EA-DCMS.png'),
                            output = TRUE,
                            imagetype="png"
  )
  VennDiagram::venn.diagram(x = list(unique(c(BS,SE,PA)), 
                                     unique(c(BE,LF,RD)),
                                     DC),
                            category.names = c("OA" , "EAA", "DCMS"),
                            filename = paste0(curr_dir,'/venn-DCMS.png'),
                            output = TRUE,
                            imagetype="png"
  )
  
{# proportion of shared SNPs among OA
oa_prop = c(sum(SE %in% PA)/((length(SE)+length(PA))-sum(SE %in% PA)),
             sum(SE %in% BS)/((length(SE)+length(BS))-sum(SE %in% BS)),
             sum(BS %in% PA)/((length(BS)+length(PA))-sum(BS %in% PA)))
# proportion of shared SNPs among EAA
eaa_prop = c(sum(BE %in% SB)/((length(BE)+length(SB))-sum(BE %in% SB)),
              sum(BE %in% LF)/((length(BE)+length(LF))-sum(BE %in% LF)),
              sum(BE %in% RD)/((length(BE)+length(RD))-sum(BE %in% RD)),
              sum(SB %in% LF)/((length(SB)+length(LF))-sum(SB %in% LF)),
              sum(SB %in% RD)/((length(SB)+length(RD))-sum(SB %in% RD)),
              sum(LF %in% RD)/((length(LF)+length(RD))-sum(LF %in% RD)))
# proportion of shared SNPs between methods
x_prop = c(sum(SE %in% BE)/((length(SE)+length(BE))-sum(SE %in% BE)),
            sum(SE %in% SB)/((length(SE)+length(SB))-sum(SE %in% SB)),
            sum(SE %in% LF)/((length(SE)+length(LF))-sum(SE %in% LF)),
            sum(SE %in% RD)/((length(SE)+length(RD))-sum(SE %in% RD)),
            sum(BS %in% BE)/((length(BS)+length(BE))-sum(BS %in% BE)),
            sum(BS %in% SB)/((length(BS)+length(SB))-sum(BS %in% SB)),
            sum(BS %in% LF)/((length(BS)+length(LF))-sum(BS %in% LF)),
            sum(BS %in% RD)/((length(BS)+length(RD))-sum(BS %in% RD)),
            sum(PA %in% BE)/((length(PA)+length(BE))-sum(PA %in% BE)),
            sum(PA %in% SB)/((length(PA)+length(SB))-sum(PA %in% SB)),
            sum(PA %in% LF)/((length(PA)+length(LF))-sum(PA %in% LF)),
            sum(PA %in% RD)/((length(PA)+length(RD))-sum(PA %in% RD)))
# Two-sample Kolmogorov-Smirnov test for difference in proportions
ks.test(x = x_prop, y = c(oa_prop,eaa_prop), alternative = 'less')
ks.test(x = oa_prop, y = eaa_prop, alternative = 'less')
}
'''
#clean-up
rm(bayescan_results, bayescan_trace_data, pcadapt_results, sambada_results, 
   selestim_results, selestim_trace_data, structure_trace_data)
    
# save session
session::save.session(session_path("15"), compress = "xz")
  
