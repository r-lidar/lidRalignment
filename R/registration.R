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
#' \dontrun{
#' ref <- readLAS("reference.las")
#' mov <- readLAS("moving.las")
#' M <- brute_force_registration(ref, mov, res = 0.5, max_offset = 4, debug = TRUE)
#' }
#'
#' @export
brute_force_registration <- function(ref, mov, res = 0.5, max_offset = 4, verbose = TRUE, ...)
{
  p = list(...)

  rotate_z <- function(pc, angle)
  {
    rotation_matrix <- matrix(c(cos(angle), -sin(angle), 0,
                                sin(angle), cos(angle), 0,
                                0, 0, 1), nrow = 3, byrow = TRUE)
    t(rotation_matrix %*% t(pc))
  }

  translate_xy <- function(pc, dx, dy)
  {
    pc + matrix(c(dx, dy, 0), nrow = nrow(pc), ncol = 3, byrow = TRUE)
  }

  rms <- function(pc1, pc2)
  {
    K = RANN::nn2(pc1, pc2, k = 1)
    D = K$nn.dists
    D = sort(as.numeric(D))
    n = length(D)/2
    sqrt(mean(D[1:n]^2))
  }

  compute_rms <- function(params)
  {
    angle <- params["angle"]
    dx <- params["dx"]
    dy <- params["dy"]

    rotated <- rotate_z(umov, angle)
    translated <- translate_xy(rotated, dx, dy)

    #rgl::plot3d(vref, col = "blue")
    #rgl::points3d(translated, col = "red")

    rms_val <- rms(vref, translated)

    list(angle = angle, dx = dx, dy = dy, rms = rms_val)
  }

  vref = as.matrix(sf::st_coordinates(lidR::decimate_points(ref, lidR::random_per_voxel(res))))
  umov = as.matrix(sf::st_coordinates(lidR::decimate_points(mov, lidR::random_per_voxel(res))))

  #x = plot(decimate_points(ref, random_per_voxel(res)), pal = "yellow", size = 2)
  #plot(decimate_points(mov, random_per_voxel(res)), add = x, pal = "red", size = 2)

  # Coarse registration
  angles <- seq(-180, 180, by = 2) * pi / 180
  dx <- seq(-max_offset, max_offset, by = 1)
  dy <- seq(-max_offset, max_offset, by = 1)
  param_grid <- expand.grid(angle = angles, dx = dx, dy = dy)

  results <- pbapply::pbapply(param_grid, 1, compute_rms)
  results <- data.table::rbindlist(results)
  best_params <- results[which.min(results$rms),]

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

    color_vector = results$dy
    col_index <- as.integer(100 * (color_vector - min(color_vector)) / (max(color_vector) - min(color_vector)))
    rgl::open3d()
    rgl::plot3d(results$angle*180/pi, results$dy, results$rms,
                xlab = 'Angle',
                ylab = 'Translation Y',
                zlab = 'RMS')
    rgl::spheres3d(best_params$angle*180/pi, best_params$dy, best_params$rms, radius = 3, col = "black")
  }

  a = best_params$angle*180/pi

  # Fine registration
  angles <- seq(a-2, a+2, by = 1) * pi / 180
  dx <- seq(best_params$dx-1, best_params$dx+1, by = .25)
  dy <- seq(best_params$dy-1, best_params$dy+1, by = .25)
  param_grid <- expand.grid(angle = angles, dx = dx, dy = dy)

  results <- pbapply::pbapply(param_grid, 1, compute_rms)
  results <- data.table::rbindlist(results)
  best_params <- results[which.min(results$rms),]


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
    rgl::spheres3d(best_params$angle*180/pi, best_params$dx, best_params$rms, radius = 0.1, col = "black")

    color_vector = results$dy
    col_index <- as.integer(100 * (color_vector - min(color_vector)) / (max(color_vector) - min(color_vector)))
    rgl::open3d()
    rgl::plot3d(results$angle*180/pi, results$dy, results$rms,
                xlab = 'Angle',
                ylab = 'Translation Y',
                zlab = 'RMS')
    rgl::spheres3d(best_params$angle*180/pi, best_params$dy, best_params$rms, radius = 0.1, col = "black")
  }

  M = rotation_matrix(best_params$angle*180/pi)
  M[1,4] = best_params$dx
  M[2,4] = best_params$dy

  attr(M, "rms") = best_params$rms

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
icp = function(vref, umov, min_error_diff = 1e-5, overlap = 90, rot = "Z", skip_txy = FALSE, cc = find_cloudcompare(), verbose = TRUE)
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

  tool = paste0("-ICP -MIN_ERROR_DIFF ", min_error_diff, " -OVERLAP ", overlap, " -ROT ", rot, skip_txy)
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
    path = "\"C:/Program Files/CloudCompare/CloudCompare.exe\""
    if (!file.exists(path)) stop(paste("Cannot find CloudCompare in", path))
    return(path)
  }
  if (os == "osx")
  {
    stop("Finding CloudCompare on macOS is not yet supported. Please provide the path to CloudCompare.")
  }

  stop("OS not detected")
}

