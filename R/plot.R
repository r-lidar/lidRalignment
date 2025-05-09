#' Visualize Alignment of Point Clouds
#'
#' This function visualizes the alignment of a moving (`mov`) point cloud with a reference (`ref`) point cloud.
#' If a transformation matrix (`M`) is provided, the function applies the transformation and displays
#' the translation and rotation.
#'
#' @param ref,mov LAS objects
#' @param M A 4x4 transformation matrix (optional). If provided, the moving dataset is transformed accordingly.
#' @param size Numeric, point size for visualization. Default is 2.
#'
#' @return A 3D plot of the reference and transformed moving point clouds. If a transformation is applied,
#' translation and rotation are also visualized with lines.
#'
#' @noRd
show_alignment = function(ref, mov, M = NULL, size = 2)
{
  . <- X <- Y <- Z <- NULL

  trans = !is.null(M)
  if (trans)
  {
    icentroid <- colMeans(mov@data[, .(X,Y,Z)])

    mov = transform_las(mov, M)
    tx <- M[1, 4]
    ty <- M[2, 4]
    tz <- M[3, 4]
    theta <- atan2(M[2, 1], M[1, 1])

    #cat(sprintf("Translation: (%.2f, %.2f, %.2f)\n", tx, ty, tz))
    #cat(sprintf("Rotation (Z-axis): %.2f degrees\n", theta*(180 / pi)))

    fcentroid <- colMeans(mov@data[, .(X,Y,Z)])
    p = c(25, 0)
    a = seq(0, theta, length.out = 25)
    arcx = p[1]*cos(a) - p[2]*sin(a)
    arcy = p[1]*sin(a) + p[2]*cos(a)

    arcx = arcx + fcentroid[1]
    arcy = arcy + fcentroid[2]

    trx = c(icentroid[1], icentroid[1] +tx)
    try = c(icentroid[2], icentroid[2] +ty)
  }

  x = lidR::plot(ref, pal = "yellow", size = size)
  lidR::plot(mov, add = x, pal = "red", size = size)

  if (trans)
  {
    rgl::lines3d(arcx-x[1],arcy-x[2],-1, col = "blue", lwd = 5)
    rgl::lines3d(trx-x[1],try-x[2],-1, col = "blue", lwd = 5)
  }
}