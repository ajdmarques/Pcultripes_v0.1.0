source('code/R/functions/session_path.R')
# restore session
session::restore.session(session_path("07"))

# save session
session::save.session("data/final/results.rda", compress = "xz")
