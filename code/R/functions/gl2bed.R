#' Export DArT SNP data to a .beed (plink) file
#'
#' This function is used to export a \code{\link[adegenet]{genlight}} object
#' as a .bed (plink) input data file.
#'
#' @param x \code{\link[adegenet]{genlight}} object.
#'
#' @param path \code{character} output file path.
#'
#' @param x \code{\link[adegenet]{genlight}} object.
#'
#' @param path \code{character} output file path.
#'
#' @param spid_path \code{character} PGDSpider converter template file.
#'
#' @param pgdspider_path \code{character} PGDSpider program file path.
#'
#' @param plink_path \code{character} plink executable file path.
#'
#' @details Note that this function requires the PGDSpider (version 2.1.1.5)
#'   program to be installed. It also requires the plink software program
#'   (version 1.07) to be installed.
#'
#' @return invisible \code{logical} \code{TRUE} indicating success.
gl2bed <- function(x, path,
                   spid_path = "code/templates/convert_str_to_ped.spid",
                   pgdspider_path = paste0("code/pgdspider/",
                                           "PGDSpider_2.1.1.5",
                                           "/PGDSpider2-cli.jar"),
                   plink_path = "code/plink/plink-1.07-x86_64/plink") {
  # assert arguments are valid
  assertthat::assert_that(inherits(x, "genlight"),
                          assertthat::is.string(path),
                          assertthat::is.string(spid_path),
                          assertthat::is.string(pgdspider_path),
                          assertthat::is.string(plink_path),
                          file.exists(spid_path),
                          file.exists(pgdspider_path),
                          file.exists(plink_path))
  assertthat::assert_that(!is.null(x@pop),
                          msg = "argument to x has no population data")
  # define tempfile names
  curr_data_path <- tempfile(fileext = ".str")
  curr_ped_path <- tempfile(fileext = ".ped")
  # save snp data as a structure file
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
   PED_WRITER_MAP_FILE_QUESTION=
   pos <- grep("PED_WRITER_MAP_FILE_QUESTION=", curr_spid)
  curr_spid[pos] <- paste0("PED_WRITER_MAP_FILE_QUESTION=",
                           raster::extension(curr_ped_path, ".map"))
  ## save modified spid conversion file
  curr_spid_path <- tempfile(fileext = ".spid")
  writeLines(curr_spid, curr_spid_path)
  # run pgdspider to convert structure file to ped file
  pgdspider_path <- normalizePath(pgdspider_path)
  withr::with_dir(dirname(pgdspider_path),
  system(paste0("java -Xmx1024m -Xms512m -jar ", basename(pgdspider_path),
         " -inputfile \"", curr_data_path, "\" ",
         "-inputformat STRUCTURE ",
         "-outputfile \"", curr_ped_path, "\" ",
         "-outputformat PED ",
         "-spid \"", curr_spid_path, "\""), intern = TRUE))
  # run plink to convert ped file to bed file
  if (!startsWith(path, "/"))
    path <- paste0(getwd(), "/", path)
  system(paste0(normalizePath(plink_path), " ",
                "--file ", tools::file_path_sans_ext(curr_ped_path), " ",
                "--noweb ",
                "--maf 0.0 ",
                "--make-bed ",
                "--out ", tools::file_path_sans_ext(path)))
  # delete plink log file if needed
  if (file.exists("plink.log"))
    unlink("plink.log")
  # return invisible success
  invisible(TRUE)
}
