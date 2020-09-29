#' Export DArT SNP data to a BayeScan file
#'
#' This function is used to export a \code{\link[adegenet]{genlight}} object
#' as a BayeScan input data file.
#'
#' @param x \code{\link[adegenet]{genlight}} object.
#'
#' @param path \code{character} output file path.
#'
#' @param spid_path \code{character} PGDSpider converter template file.
#'
#' @param pgdspider_path \code{character} PGDSpider program file path.
#'
#' @details Note that this function requires the PGDSpider (version 2.1.1.5)
#'   program to be installed.
#'
#' @return invisible \code{logical} \code{TRUE} indicating success.
gl2bayescan <- function(x, path,
                        spid_path = "code/templates/convert_str_to_bs.spid",
                        pgdspider_path = paste0("code/pgdspider/",
                                                 "PGDSpider_2.1.1.5",
                                                 "/PGDSpider2-cli.jar")) {
  # assert that arguments are valid
  assertthat::assert_that(inherits(x, "genlight"),
                          assertthat::is.string(path),
                          file.exists(spid_path),
                          file.exists(pgdspider_path))
  assertthat::assert_that(!is.null(x@pop),
                          msg = "argument to x has no population data")
  # save snp data as a structure file
  curr_data_path <- tempfile(fileext = ".str")
  withr::with_dir(dirname(curr_data_path),
    dartR::gl2structure(x,
      addcolumns = data.frame(PopData = as.character(as.integer(x@pop))),
      outfile = basename(curr_data_path)))
  # save spid conversion format file
  ## load default settings from template
  curr_spid <- readLines(spid_path)
  pos <- grep("STRUCTURE_PARSER_NUMBER_LOCI_QUESTION=", curr_spid)
  curr_spid[pos] <- paste0("STRUCTURE_PARSER_NUMBER_LOCI_QUESTION=",
                           sprintf("%i", length(x@loc.names)))
  ## save modified spid conversion file
  curr_spid_path <- tempfile(fileext = ".spid")
  writeLines(curr_spid, curr_spid_path)
  # run pgdspider
  pgdspider_path <- normalizePath(pgdspider_path)
  if (!startsWith(path, "/"))
    path <- paste0(getwd(), "/", path)
  withr::with_dir(dirname(pgdspider_path),
  system(paste0("java -Xmx1024m -Xms512m -jar ", basename(pgdspider_path),
         " -inputfile \"", curr_data_path, "\" ",
         "-inputformat STRUCTURE ",
         "-outputfile \"", path, "\" ",
         "-outputformat GESTE_BAYE_SCAN ",
         "-spid \"", curr_spid_path, "\""), intern = TRUE))
  # return success
  invisible(TRUE)
}
