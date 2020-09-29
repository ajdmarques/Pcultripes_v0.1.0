#' Extract references from an Rmarkdown file
#'
#' This function extracts reference codes (e.g. @r1) from an RMarkdown file and
#' returns them.
#'
#' @param x \code{character} file path.
#'
#' @return \code{character} vector with references codes.
extract_refs <- function(x) {
  # load file
  txt <- readLines(x)
  # remove YAML header
  txt <- txt[-1 * seq(1, max(which(txt == "---")))]
  # collapse text
  txt <- paste(txt, collapse = " ")
  # extract text in square brackets
  r <- gsub("[\\[\\]]", "", regmatches(txt, gregexpr("\\[.*?\\]", txt))[[1]])
  # remove square brackets
  r <- gsub("[", "", r, fixed = TRUE)
  r <- gsub("]", "", r, fixed = TRUE)
  # split r by semi-colon and space delimiters
  r <- unlist(strsplit(r, ";"))
  r <- unlist(strsplit(r, " "))
  # subset words that start with @ symbol
  r <- r[startsWith(r, "@")]
  # return result
  r
}
