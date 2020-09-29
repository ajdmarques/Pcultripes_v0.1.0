#' Make pie chart grobs
#'
#' This function makes grobs for rendering pie charts on a map.
#'
#' @param x \code{\link[sf]{sf}} object. Each row should be for a different
#'   sample and each column should correspond to the proportion of admixture.
#'
#' @param width \code{numeric} width of the pie charts on the plot.
#'
#' @param xlim \code{numeric} x-axis limits.
#'
#' @param ylim \code{numeric} y-axis limits.
#'
#' @param force \code{numeric} buffer used for repelling the viewports.
#'
#' @return \code{list} objects. Each element in the \code{list} object
#'   corresponds to a different row in the argument to \code{x}. Within
#'   the list for each \code{x}, the
make_pie_chart_grobs <- function(x, width, xlim, ylim, force = 1) {
  # assert arguments are valid
  assertthat::assert_that(inherits(x, "sf"), assertthat::is.number(width),
                          isTRUE(width > 0),
                          is.numeric(xlim), length(xlim) == 2,
                          is.numeric(ylim), length(ylim) == 2,
                          all(is.finite(xlim)), all(is.finite(ylim)),
                          assertthat::is.number(force))
  # create vector with population colors
  pop_colors <- setNames(
    suppressWarnings(RColorBrewer::brewer.pal(ncol(x) - 1, "Set1")),
    setdiff(names(x), "geometry"))
  pop_colors <- pop_colors[seq_len(ncol(x) - 1)]
  # make pie chart ggplot2 canvas
  pie_chart_canvas <-
    ggplot2::ggplot() +
    ggplot2::theme_classic() +
    ggplot2::coord_polar("y", start = 0) +
    ggplot2::theme(axis.ticks = ggplot2::element_blank(),
                   axis.text = ggplot2::element_blank(),
                   axis.title = ggplot2::element_blank(),
                   axis.line = ggplot2::element_blank(),
                   axis.ticks.length = ggplot2::unit(0, "null"),
                   axis.ticks.margin = ggplot2::unit(0, "null"),
                   panel.grid = ggplot2::element_blank(),
                   legend.position = "none",
                   legend.margin = ggplot2::unit(0, "null"),
                   plot.margin = ggplot2::unit(c(0, 0, 0, 0), "null"),
                   panel.border = ggplot2::element_blank(),
                   plot.background = ggplot2::element_blank(),
                   panel.background = ggplot2::element_blank()) +
    ggplot2::scale_fill_manual(values = pop_colors)
  # make pie chart grobs
  g <- plyr::llply(seq_len(nrow(x)), function(i) {
    pie_chart_canvas +
    ggplot2::geom_col(ggplot2::aes(x = group, y = proportion,
                                   fill = population),
                      x %>%
                      as.data.frame() %>%
                      select(-geometry) %>%
                      dplyr::filter(seq_len(nrow(.)) == i) %>%
                      tidyr::gather(population, proportion) %>%
                      dplyr::mutate(group = "a"))
  })
  # preliminary calculations
  ## create centroids/boxes in sf coordinate system
  grob_centroids <-
    x %>%
    as("Spatial") %>%
    as.data.frame() %>%
    dplyr::rename(x = coords.x1, y = coords.x2) %>%
    dplyr::select(x, y)
  grob_boxes <-
    x %>%
    as("Spatial") %>%
    as.data.frame() %>%
    dplyr::rename(x = coords.x1, y = coords.x2) %>%
    dplyr::mutate(id = seq_len(nrow(.))) %>%
    plyr::ddply("id", function(x) {
      data.frame(x1 = x$x - width, y1 = x$y - width,
                 x2 = x$x + width, y2 = x$y + width)
    }) %>%
    plyr::rbind.fill() %>%
    dplyr::select(x1, y1, x2, y2)
  ## redefine xlim and ylim if the extent including the boxes is larger
  if (min(grob_boxes$x1) < xlim[1]) xlim[1] <- min(grob_boxes$x1) * 0.9
  if (min(grob_boxes$y1) < ylim[1]) ylim[1] <- min(grob_boxes$y1) * 0.9
  if (max(grob_boxes$x2) > xlim[2]) xlim[2] <- max(grob_boxes$x2) * 1.1
  if (max(grob_boxes$y2) > ylim[2]) ylim[2] <- max(grob_boxes$y2) * 1.1
  ## check that limits are correct
  assertthat::assert_that(all(unlist(grob_boxes[, c("x1", "x2")]) > xlim[1]),
                          all(unlist(grob_boxes[, c("x1", "x2")]) < xlim[2]),
                          all(unlist(grob_boxes[, c("y1", "y2")]) > ylim[1]),
                          all(unlist(grob_boxes[, c("y1", "y2")]) < ylim[2]))
  ## create centroids/boxes in 0-1 coordinate system
  sclim <- xlim
  if (abs(diff(ylim)) < abs(diff(xlim))) sclim <- ylim
  width_sc <- scales::rescale(sclim[1] + width, from = sclim) -
              scales::rescale(sclim[1], from = sclim)
  ## create rescaled centroids and box extents
  grob_centroids_sc <- grob_centroids %>%
                       dplyr::mutate(x = scales::rescale(x, from = xlim),
                                     y = scales::rescale(y, from = ylim))
  grob_boxes_sc <-
    grob_centroids_sc %>%
    dplyr::mutate(x1 = x - width_sc, y1 = y - width_sc,
                  x2 = x + width_sc, y2 = y + width_sc) %>%
    dplyr::select(x1, y1, x2, y2)
  ## double check that ggrepel is given sane input
  assertthat::assert_that(all(unlist(grob_centroids_sc) <= 1),
                          all(unlist(grob_centroids_sc) >= 0),
                          all(unlist(grob_boxes_sc) <= 1),
                          all(unlist(grob_boxes_sc) >= 0))
  ## calculate hypotenuse
  h <- as.numeric(dist(matrix(c(xlim, ylim), ncol = 2)))
  # compute data for repelling pie charts
  ## calculate centroids
  r <- ggrepel:::repel_boxes(
    data_points = as.matrix(grob_centroids_sc),
    point_padding_x = 0, point_padding_y = 0,
    boxes = as.matrix(grob_boxes_sc), xlim = c(0, 1), ylim = c(0, 1),
    hjust = 0.5, vjust = 0.5, force = force * 1e-6,
    maxiter = 2000, direction = "both") %>%
    setNames(c("x", "y")) %>%
    as.matrix()
  ## double check that ggrepel gives sane output
  assertthat::assert_that(nrow(r) == nrow(grob_centroids), ncol(r) == 2,
                          assertthat::noNA(c(r)),
                          all(c(r) <= 1), all(c(r) >= 0),
                          all(is.finite(c(r))))
  ## rescale new centroids back to original coordinate system
  new_grob_centroids <-
    r %>%
    as.data.frame() %>%
    setNames(c("x", "y")) %>%
    dplyr::mutate(x = scales::rescale(x, from = c(0, 1), to = xlim),
                  y = scales::rescale(y, from = c(0, 1), to = ylim))
  # create output
  out <- lapply(seq_len(nrow(grob_centroids)), function(i) {
    # initialize output list
    o <- list()
    # grob
    o$grob <- ggplot2::ggplotGrob(g[[i]])
    # grob dimensions
    o$dims <- data.frame(x = new_grob_centroids[i, 1],
                         y = new_grob_centroids[i, 2],
                         xmin = new_grob_centroids[i, 1] - width,
                         xmax = new_grob_centroids[i, 1] + width,
                         ymin = new_grob_centroids[i, 2] - width,
                         ymax = new_grob_centroids[i, 2] + width)
    # add original dimensions
    o$dims0 <- data.frame(x = grob_centroids[i, 1],
                          y = grob_centroids[i, 2],
                          xmin = grob_centroids[i, 1] - width,
                          xmax = grob_centroids[i, 1] + width,
                          ymin = grob_centroids[i, 2] - width,
                          ymax = grob_centroids[i, 2] + width)
    # create line
    o$line <- data.frame(x = c(new_grob_centroids[i, 1],
                               grob_centroids[i, 1]),
                         y = c(new_grob_centroids[i, 2],
                               grob_centroids[i, 2]))
    # return output
    o
  })
  # return output
  out
}
