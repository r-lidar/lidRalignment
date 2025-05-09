#' Brute Force Point Cloud Registration
#'
#' This function performs a brute-force registration of a moving point cloud to a reference point cloud
#' by iterating over a grid of rotation and translation parameters to find the best alignment.
#'
#' @param ref A `LAS` object containing the reference point cloud.
#' @param mov A `LAS` object containing the moving point cloud to be aligned.
#' @param res Numeric. The resolution for voxel-based downsampling (default: 0.5).
#' @param max_offset Numeric. The maximum translation offset in X and Y directions (default: 4). Reducing
#' the value decrease computation time but reduces the search range.
#' @param ... Additional parameters, including `debug = TRUE` to enable visualization.
#' @param verbose logical
#'
#' @return A 4x4 transformation matrix that best aligns `mov` with `ref`, with an attribute `rms`
#' indicating the root mean square (RMS) registration error.
#'
#' @details The function performs coarse registration by evaluating transformations over a large parameter grid,
#' followed by a fine registration step with a refined grid. RMS error is used to assess the alignment quality.
#' Optionally, debug visualizations can be displayed using `rgl`.
#' @noMd
brute_force_registration <- function(ref, mov, res = 0.5, max_offset = 8, verbose = TRUE, ...)
{
  p = list(...)

  strategy_ref = attr(ref, "strategy")
  strategy_mov = attr(mov, "strategy")
  strategy = NULL

  if (is.null(strategy_ref) || is.null(strategy_mov))
    warning("Input point clouds do not have an attribute 'strategy' that record the registration strategy")
  else if (strategy_mov != strategy_ref)
    stop("Input point clouds do not have the same attribute 'strategy' that record the registration strategy")
  else
    strategy = strategy_ref

  X <- Y <- Z <- . <- NULL

  vref = as.matrix(lidR::decimate_points(ref, lidR::random_per_voxel(res)@data[, .(X,Y,Z)]))
  umov = as.matrix(lidR::decimate_points(mov, lidR::random_per_voxel(res)@data[, .(X,Y,Z)]))

  #x = plot(decimate_points(ref, random_per_voxel(res)), pal = "yellow", size = 2)
  #plot(decimate_points(mov, random_per_voxel(res)), add = x, pal = "red", size = 2)

  # Coarse registration
  angles <- c(0, seq(-180, 180, by = 2) * pi / 180)
  dx <- c(0, seq(-max_offset, max_offset, by = 1))
  dy <- c(0, seq(-max_offset, max_offset, by = 1))
  param_grid <- expand.grid(angle = angles, dx = dx, dy = dy, dz = 0)

  if (!is.null(strategy) && strategy == "chm-dtm")
  {
    dtm_ref = lidR::pixel_metrics(ref, ~min(Z), 1)
    dtm_mov = lidR::pixel_metrics(mov, ~min(Z), 1)

    Z0 = terra::extract(dtm_mov, cbind(0,0))$V1
    Z = terra::extract(dtm_ref, param_grid[,2:3])$V1
    param_grid$dz = Z-Z0
  }

  param_grid <- as.matrix(param_grid)

  results = rms_scan_grid(vref, umov, param_grid)
  best_params <- results[which.min(results$rms),]

  rmsi = results$rms[1]

  if (isTRUE(p$debug))
  {
    col_gradient <- grDevices::colorRampPalette(c("blue", "green", "yellow", "red"))(100)
    color_vector = results$dy
    col_index <- as.integer(100 * (color_vector - min(color_vector)) / (max(color_vector) - min(color_vector)))

    rgl::open3d()
    rgl::plot3d(results$angle*180/pi, results$dx, results$rms,
                col = col_gradient[col_index],
                xlab = 'Angle',
                ylab = 'Translation X',
                zlab = 'RMS')
    rgl::spheres3d(best_params$angle*180/pi, best_params$dx, best_params$rms, radius = 3, col = "black")
  }

  a = best_params$angle*180/pi

  # Fine registration
  angles <- seq(a-3, a+2, by = 1) * pi / 180
  dx <- seq(best_params$dx-1, best_params$dx+1, by = .25)
  dy <- seq(best_params$dy-1, best_params$dy+1, by = .25)
  dz <- seq(best_params$dz-0.5, best_params$dz+0.5, by = 0.1)
  param_grid <- expand.grid(angle = angles, dx = dx, dy = dy, dz = dz)
  param_grid <- as.matrix(param_grid)

  results = rms_scan_grid(vref, umov, param_grid)
  best_params <- results[which.min(results$rms),]


  if (isTRUE(p$debug))
  {
    col_gradient <- grDevices::colorRampPalette(c("blue", "green", "yellow", "red"))(100)
    color_vector = results$dy
    col_index <- as.integer(100 * (color_vector - min(color_vector)) / (max(color_vector) - min(color_vector)))

    rgl::points3d(results$angle*180/pi, results$dx, results$rms,
                col = col_gradient[col_index],
                xlab = 'Angle',
                ylab = 'Translation X',
                zlab = 'RMS')
    rgl::spheres3d(best_params$angle*180/pi, best_params$dx, best_params$rms, radius = 3, col = "blue")
  }

  M = rotation_matrix(best_params$angle*180/pi)
  M[1,4] = best_params$dx
  M[2,4] = best_params$dy
  M[3,4] = best_params$dz

  attr(M, "rms_init") = round(rmsi, 4)
  attr(M, "rms") = round(best_params$rms, 4)

  if (verbose) rtm_info(M)

  return(M)
}

icp = function(vref, umov, overlap = 90, tz_only = FALSE, rz_only = TRUE, tolerance = 1e-6)
{
  vref = as.matrix(vref@data)
  umov = as.matrix(umov@data)
  cpp_icp(umov, vref, tz_only = tz_only, rz_only = rz_only, overlap = overlap, tolerance = tolerance)
}

