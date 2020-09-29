source('code/R/functions/session_path.R')
# restore session
session::restore.session(session_path("07"))

# packages
library(psych) #Calls: pairs.panels
library(car) #Calls: vif
library(plyr) #Calls: rbind.fill
library(dplyr) #Calls: %>%

# load parameters
vif_parameters <- "code/parameters/vif.toml" %>%
  RcppTOML::parseTOML() %>%
  `[[`(MODE)

# create directory to save vif results
unlink("data/intermediate/vif-results", force = TRUE, recursive = TRUE)
dir.create("data/intermediate/vif-results", showWarnings = FALSE,
           recursive = TRUE)
curr_dir <- paste0("data/intermediate/vif-results/")

# read environmental variable
env_raw <- read.csv(paste0('data/raw/', vif_parameters$env_data), header = T)

# Fit a linear model to the data
mod <- lm(elevate ~ ., data = env_raw[,3:length(colnames(env_raw))])
plot(mod, which=1:2)

# check IF variables are aliased
alias(mod) 
# if aliased, remove dependant variable
env_raw[,"bio6"] <- NULL
env_raw[,vif_parameters$lat] <- NULL
env_raw[,vif_parameters$lon] <- NULL

# remove aliased variable
mod <- lm(elevate ~ ., data = env_raw[,3:length(colnames(env_raw))])

# Generate scatterplots with Pearson correlations
pairs.panels(env_raw[,3:length(colnames(env_raw))])

# calculate variance inflation factors
vif(mod)
# Choose VIF cutoff
cutoff=2
# Create function to sequentially drop the variable with the largest VIF until VIF > cutoff
flag=TRUE
viftable=data.frame()
while(flag==TRUE) {
  vfit=vif(mod)
  viftable=rbind.fill(viftable,as.data.frame(t(vfit)))
  if(max(vfit)>cutoff) { mod=
    update(mod,as.formula(paste(".","~",".","-",names(which.max(vfit))))) }
  else { flag=FALSE } }

# Look at the final model
print(mod)

# And show the order in which variables were dropped
print(viftable)

# And associated VIFs
print(vfit)
vif_names <- names(vfit)

# restore raw data
env_raw <- read.csv(paste0('data/raw/', vif_parameters$env_data), header = T)
# join columns for retention
vif_env <- env_raw[,"vaucher"] %>%
  cbind(env_raw[c("lat","lon","elevate")]) %>%
  cbind(env_raw[, colnames(env_raw) %in% vif_names[vif_names != c("lat", "lon")]])
  
# write environment data to file
vif_env_path <- "data/intermediate/vif-results/env.csv"
write.csv(vif_env, vif_env_path)

