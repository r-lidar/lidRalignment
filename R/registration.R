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
#'
#' @examples
#' @export
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

  rotate_z <- function(pc, angle)
  {
    rotation_matrix <- matrix(c(cos(angle), -sin(angle), 0,
                                sin(angle), cos(angle), 0,
                                0, 0, 1), nrow = 3, byrow = TRUE)
    t(rotation_matrix %*% t(pc))
  }

  translate_xyz <- function(pc, dx, dy, dz)
  {
    pc + matrix(c(dx, dy, dz), nrow = nrow(pc), ncol = 3, byrow = TRUE)
  }

  compute_rms <- function(params, vref, umov)
  {
    angle <- params["angle"]
    dx <- params["dx"]
    dy <- params["dy"]
    dz <- params["dz"]

    rotated <- rotate_z(umov, angle)
    translated <- translate_xyz(rotated, dx, dy, dz)

    #rgl::plot3d(vref, col = "blue")
    #rgl::points3d(translated, col = "red")

    rms_val <- compare_alignment(vref, translated)

    list(angle = angle, dx = dx, dy = dy, dz = dz, rms = rms_val)
  }

  vref = as.matrix(sf::st_coordinates(lidR::decimate_points(ref, lidR::random_per_voxel(res))))
  umov = as.matrix(sf::st_coordinates(lidR::decimate_points(mov, lidR::random_per_voxel(res))))


  #x = plot(decimate_points(ref, random_per_voxel(res)), pal = "yellow", size = 2)
  #plot(decimate_points(mov, random_per_voxel(res)), add = x, pal = "red", size = 2)

  # Coarse registration
  angles <- c(0, seq(-180, 180, by = 2) * pi / 180)
  dx <- c(0, seq(-max_offset, max_offset, by = 1))
  dy <- c(0, seq(-max_offset, max_offset, by = 1))
  param_grid <- expand.grid(angle = angles, dx = dx, dy = dy, dz = 0)

  if (!is.null(strategy) && strategy == "chm-dtm")
  {
    dtm_ref = pixel_metrics(ref, ~min(Z), 1)
    dtm_mov = pixel_metrics(mov, ~min(Z), 1)

    Z0 = terra::extract(dtm_mov, cbind(0,0))$V1
    Z = terra::extract(dtm_ref, param_grid[,2:3])$V1
    param_grid$dz = Z-Z0
  }


  n = lidR::get_lidr_threads()
  groups <- cut(seq(1, nrow(param_grid)), breaks = n*8, labels = FALSE)
  param_grid <- split(param_grid, groups)

  results <- pbapply::pblapply(param_grid, function(x, vref, umov)
  {
    ans = apply(x, 1, compute_rms, vref = vref, umov = umov)
    ans = data.table::rbindlist(ans)
  }, vref = vref, umov = umov, cl = n)
  results <- data.table::rbindlist(results)
  best_params <- results[which.min(results$rms),]

  rmsi = results$rms[1]

  if (!is.null(p$debug))
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

  n = lidR::get_lidr_threads()
  groups <- cut(seq(1, nrow(param_grid)), breaks = n*4, labels = FALSE)
  param_grid <- split(param_grid, groups)

  results <- pbapply::pblapply(param_grid, function(x, vref, umov)
  {
    ans = apply(x, 1, compute_rms, vref = vref, umov = umov)
    ans = data.table::rbindlist(ans)
  }, vref = vref, umov = umov, cl = n)
  results <- data.table::rbindlist(results)
  best_params <- results[which.min(results$rms),]


  if (!is.null(p$debug))
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

#' CloudCompare ICP Registration
#'
#' This function performs point cloud registration using the CloudCompare ICP (Iterative Closest Point) algorithm.
#' CloudCompare must be installed to use this function. You can download it from: https://www.danielgm.net/cc/.
#'
#' @param vref A `LAS` containing the reference point cloud
#' @param umov A `LAS` containing the moving point cloud
#' @param min_error_diff Numeric. The minimum error difference for the ICP algorithm to terminate (default: 1e-5).
#' @param overlap Numeric. to specify the percentage of (final) overlap (integer number between 10 and 100 - default = 100)
#' @param rot Character. 'XYZ' or 'X' or 'Y' or 'Z' or 'NONE' to constrain the rotation around a given axis
#' @param skip_txy Logical. Block XY translation.
#' (or no rotation at all)
#' @param cc Character. The path to the CloudCompare executable
#' @param verbose logical
#'
#' @return A 4x4 transformation matrix resulting from the ICP registration.
#'
#' @details The function exports both point clouds to temporary `.xyz` files, runs CloudCompare ICP
#' via a system call, retrieves the transformation matrix from the generated file, and then deletes
#' temporary files.
#' @export
#' @md
icp = function(vref, umov, min_error_diff = 1e-5, overlap = 90, rot = "Z", skip_txy = FALSE, skip_tz = FALSE, cc = find_cloudcompare(), verbose = TRUE)
{
  . <- X <- Y <- Z <- NULL

  on.exit(file.remove(fmatrix))

  stopifnot(overlap <= 100, overlap >= 10)
  rot = match.arg(rot, c('XYZ','X','Y','Z','NONE'))

  ovref = paste0(tempdir(),"/vref.xyz")
  oumov = paste0(tempdir(),"/umov.xyz")

  data.table::fwrite(vref@data[, .(X,Y,Z)], ovref, sep = " ", col.names = FALSE)
  data.table::fwrite(umov@data[, .(X,Y,Z)], oumov, sep = " ", col.names = FALSE)

  skip_txy = if(skip_txy) " -SKIP_TX -SKIP_TY " else ""
  skip_tz  = if(skip_tz) " -SKIP_TZ " else ""

  tool = paste0("-ICP -MIN_ERROR_DIFF ", min_error_diff, " -OVERLAP ", overlap, " -ROT ", rot, skip_txy, skip_tz)
  args = paste0("-SILENT -O ", oumov, " -O ", ovref)
  cmd  = paste(cc, args, tool)

  system(cmd, ignore.stdout = TRUE, ignore.stderr = TRUE, timeout = 60)

  fmatrix <- tools::file_path_sans_ext(oumov)
  fmatrix <- paste0(fmatrix, "_REGISTRATION_MATRIX")
  dir     <- dirname(fmatrix)
  pattern <- paste0("^", basename(fmatrix))
  fmatrix <- list.files(dir, pattern = pattern, full.names = TRUE)

  M = data.table::fread(fmatrix)
  M <- as.matrix(M)

  if (verbose) rtm_info(M)

  return(M)
}

#' Find CloudCompare executable path
#'
#' This function attempts to locate the CloudCompare executable based on the operating system.
#' On Windows, it queries the registry to find the installation path. On macOS and Linux (Flathub version),
#' it returns an error.
#'
#' @return A character string with the path to the CloudCompare executable.
#' @export
find_cloudcompare = function()
{
  if (Sys.info()["nodename"] == "FFGG-1009803")
  {
    return("/home/jr/Logiciels/CloudCompare/build/qCC/CloudCompare")
  }

  os = get_os()

  if (os == "linux")
  {
    stop("The Linux version of CloudCompare is a Flathub package and cannot be called from the command line. Please provide the path to a compiled version of CloudCompare.")
  }

  if (os == "windows")
  {
    path = "C:/Program Files/CloudCompare/CloudCompare.exe"
    if (!file.exists(path)) stop(paste("Cannot find CloudCompare in", path))
    path = "\"C:/Program Files/CloudCompare/CloudCompare.exe\""
    return(path)
  }
  if (os == "osx")
  {
    stop("Finding CloudCompare on macOS is not yet supported. Please provide the path to CloudCompare.")
  }

  stop("OS not detected")
}

#' Compare the Alignment of Two Point Clouds
#'
#' Computes the root mean square error (RMSE) of the nearest neighbor distances between two point clouds.
#' This function measures how well `pc2` aligns with `pc1` by finding the closest points in `pc1` for
#' each point in `pc2` and computing the RMSE for half of the point cloud. It uses the 50% of the
#' most aligned points to account for non-overlapping features. Thus, the returned value is not
#' a true RMSE but rather an arbitrary alignment score.
#'
#' @param pc1 A matrix of XYZ coordinates or a `LAS` object from the `lidR` package.
#' @param pc2 A matrix of XYZ coordinates or a `LAS` object from the `lidR` package.
#' @noMd
compare_alignment <- function(pc1, pc2)
{
  if (methods::is(pc1, "LAS")) pc1 = as.matrix(sf::st_coordinates(pc1))
  if (methods::is(pc2, "LAS")) pc2 = as.matrix(sf::st_coordinates(pc2))
  K = RANN::nn2(pc1, pc2, k = 1)
  D = K$nn.dists
  D = sort(as.numeric(D))
  n = length(D)/2
  sqrt(mean(D[1:n]^2))
}


