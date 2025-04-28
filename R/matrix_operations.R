#' Create a Translation Matrix
#'
#' Generates a 4x4 homogeneous transformation matrix for translation.
#'
#' @param tx Translation along the X-axis.
#' @param ty Translation along the Y-axis.
#' @param tz Translation along the Z-axis.
#' @return A 4x4 transformation matrix.
#' @export
translation_matrix <- function(tx, ty, tz)
{
  matrix(c(
    1, 0, 0, tx,
    0, 1, 0, ty,
    0, 0, 1, tz,
    0, 0, 0, 1
  ), nrow = 4, byrow = TRUE)
}

#' Create a Rotation Matrix
#'
#' Generates a 4x4 homogeneous transformation matrix for rotation around X, Y, or Z.
#'
#' @param angle Rotation angle in degrees.
#' @param axis Axis of rotation ('x', 'y', or 'z').
#' @return A 4x4 transformation matrix.
#' @export
rotation_matrix <- function(angle, axis = "z")
{
  angle <- angle * pi / 180  # Convert to radians
  c <- cos(angle)
  s <- sin(angle)

  M <- switch(axis,
              x = matrix(c(
                1,  0,  0, 0,
                0,  c, -s, 0,
                0,  s,  c, 0,
                0,  0,  0, 1
              ), nrow = 4, byrow = TRUE),
              y = matrix(c(
                c,  0,  s, 0,
                0,  1,  0, 0,
                -s,  0,  c, 0,
                0,  0,  0, 1
              ), nrow = 4, byrow = TRUE),
              z = matrix(c(
                c, -s,  0, 0,
                s,  c,  0, 0,
                0,  0,  1, 0,
                0,  0,  0, 1
              ), nrow = 4, byrow = TRUE),
              stop("Invalid axis. Use 'x', 'y', or 'z'.")
  )
  return(M)
}

#' Apply Translation to a LAS Object
#'
#' Translates a LAS point cloud by a given offset.
#'
#' @param las A `LAS` object.
#' @param tx Translation along X-axis (default 0).
#' @param ty Translation along Y-axis (default 0).
#' @param tz Translation along Z-axis (default 0).
#' @return A modified LAS object.
#' @noRd
translate_las <- function(las, tx = 0, ty = 0, tz = 0)
{
  las@data$X <- las@data$X + tx
  las@data$Y <- las@data$Y + ty
  las@data$Z <- las@data$Z + tz
  h = lidR::header(las)
  h[["X offset"]] = ceiling(min(las$X))
  h[["Y offset"]] = ceiling(min(las$Y))
  h[["Z offset"]] = ceiling(min(las$Z))
  lidR::quantize(las@data[["X"]], h[["X scale factor"]], h[["X offset"]])
  lidR::quantize(las@data[["Y"]], h[["Y scale factor"]], h[["Y offset"]])
  lidR::quantize(las@data[["Z"]], h[["Z scale factor"]], h[["Z offset"]])
  return(lidR::las_update(las))
}

#' Apply a Transformation Matrix to a LAS Object
#'
#' Applies a 4x4 transformation matrix to a LAS point cloud. The function can handle both `LAS` objects
#' and LAS file paths, transforming the point cloud data accordingly.
#'
#' @param las A `LAS` object or a character string representing the file path to a LAS file.
#' @param M A 4x4 transformation matrix.
#' @param crs Optional CRS to assign to the transformed LAS. Defaults to `sf::NA_crs_`.
#' @return A transformed `LAS` object. If a file path is provided for `las`, a new LAS file is created and returned.
#' @export
#' @md
transform_las <- function(las, M, crs = sf::NA_crs_)
{
  stopifnot(dim(M) == c(4, 4))

  if (methods::is(las, "LAS"))
  {
    coords <- as.matrix(sf::st_coordinates(las))
    coords <- cbind(coords, 1)
    transformed <- coords %*% t(M)

    xyz = as.data.frame(transformed[, 1:3])
    names(xyz) <- c("X", "Y", "Z")
    h = lidR::header(las)
    h[["X offset"]] = ceiling(min(xyz[["X"]]))
    h[["Y offset"]] = ceiling(min(xyz[["Y"]]))
    h[["Z offset"]] = ceiling(min(xyz[["Z"]]))
    lidR::quantize(xyz[["X"]], h[["X scale factor"]], h[["X offset"]])
    lidR::quantize(xyz[["Y"]], h[["Y scale factor"]], h[["Y offset"]])
    lidR::quantize(xyz[["Z"]], h[["Z scale factor"]], h[["Z offset"]])
    las_new <- lidR::LAS(xyz, h)

    lidR::st_crs(las_new) <- crs

    return(las_new)
  }

  if (is.character(las))
  {
    stopifnot(file.exists(las))

    ofile = paste0(tools::file_path_sans_ext(las), "_registered.laz")

    mat = lasR::load_matrix(M)
    reg = lasR::transform_with(mat)
    crs = lasR::set_crs(crs$wkt)
    io  = lasR::write_las(ofile)
    pipeline = mat + reg + crs + io
    ans = lasR::exec(pipeline, on = las, progress = TRUE)
    return(ans$write_las)
  }

  stop("Invalid input 'las'")
}

#' Extract Translation and Rotation from a Transformation Matrix
#'
#' Prints translation and rotation angles from a 4x4 transformation matrix.
#'
#' @param M A 4x4 transformation matrix.
#' @return None (prints results).
#' @export
rtm_info <- function(M)
{
  translation <- M[1:3, 4]
  names(translation) <- c("X", "Y", "Z")

  R <- M[1:3, 1:3]
  theta_x <- atan2(R[3,2], R[3,3])
  theta_y <- atan2(-R[3,1], sqrt(R[3,2]^2 + R[3,3]^2))
  theta_z <- atan2(R[2,1], R[1,1])

  angles <- round(c(theta_x, theta_y, theta_z) * 180 / pi, 3)
  names(angles) <- c("X", "Y", "Z")

  cat("--- Transformation Matrix Info ---\n")
  cat("Translation:\n")
  for (i in names(translation)) {
    cat(sprintf(" %-2s: %1.3f m\n", i, translation[i]))
  }

  cat("Rotation:\n")
  for (i in names(angles)) {
    cat(sprintf(" %-2s: %1.3f\u00B0\n", i, angles[i]))
  }
  cat("----------------------------------\n")
}

#' Combine Multiple Transformation Matrices
#'
#' This function takes any number of matrices and computes their combined
#' transformation through matrix multiplication. The matrices should be provided
#' in the order they are applied.
#'
#' @param ... A variable number of matrices to be combined. These matrices should
#'   be compatible for matrix multiplication (i.e., the number of columns of one
#'   matrix should match the number of rows of the next).
#'
#' @return A matrix representing the combined transformation. The output is the
#'   result of multiplying all input matrices in reverse order.
#'
#' @examples
#' M1 <- matrix(c(1, 0, 0, 1), nrow = 2)  # Identity matrix (no transformation)
#' M2 <- matrix(c(0, -1, 1, 0), nrow = 2) # 90-degree rotation matrix
#' M3 <- matrix(c(1, 0, 3, 1), nrow = 2)  # Translation matrix (move by 3 on x-axis)
#'
#' M <- combine_transformations(M1, M2, M3)
#' @export
combine_transformations <- function(...)
{
  matrices <- list(...)
  result <- matrices[[1]]
  for (i in 2:length(matrices)) result <- matrices[[i]] %*% result
  return(result)
}