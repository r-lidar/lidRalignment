#' Extract Features from a LiDAR Point Cloud
#'
#' This function extracts features from a LiDAR point cloud based on the specified strategy in order to
#' retain points that are likely matching in order to perform a registration of the point clouds.
#'
#' @param las A `LAS` object containing the LiDAR point cloud data.
#' @param strategy Character. The extraction strategy to use. `"chm-dtm"` or `"mls-tls"`. See details.
#' @param verbose logical.
#' @return A `LAS` object containing the extracted features.
#'
#' @details
#' - The `"chm-dtm"` strategy extracts a digital terrain model (DTM) and a canopy height model (CHM),
#'   and merges them into a new `LAS` object.
#' - The `"trunks"` strategy applies a series of filtering and processing steps including noise removal,
#'   smoothing, decimation, anisotropy filtering, and connected component filtering in order to extract
#'   trees in the first 4 meters of the point clouds
#' @noRd
#' @md
extract_features = function(las, strategy = "chm-dtm", verbose = TRUE)
{
  .N <- N <- hag <- anisotropy <- clusterID <- NULL

  strategy = match.arg(strategy, c("chm-dtm", "trunks"))

  if (strategy == "chm-dtm")
  {
    if(verbose) cat(" Computing Digital Terrain Model...\n")
    dtm = suppressMessages(lidR::rasterize_terrain(las, res = 0.5))

    if(verbose) cat(" Computing Canopy Heigh Model...\n")
    chm = lidR::rasterize_canopy(las, 0.25)

    dtm_points = terra::as.data.frame(dtm, xy = TRUE)
    dtm_points$Classification = 2L
    chm_points = terra::as.data.frame(chm, xy = TRUE)
    chm_points$Classification = 1L
    points = rbind(dtm_points, chm_points)
    names(points) = c('X', 'Y', 'Z', "Classification")

    h = rlas::header_create(points)
    h[["X offset"]] = ceiling(min(las$X))
    h[["Y offset"]] = ceiling(min(las$Y))
    h[["Z offset"]] = ceiling(min(las$Z))
    lidR::quantize(points[["X"]], h[["X scale factor"]], h[["X offset"]])
    lidR::quantize(points[["Y"]], h[["Y scale factor"]], h[["Y offset"]])
    lidR::quantize(points[["Z"]], h[["Z scale factor"]], h[["Z offset"]])

    res = lidR::LAS(points, h)
    res = lidR::height_above_ground(res, dtm)
    res@header@VLR = list()

    lidR::st_crs(res) = lidR::st_crs(las)
    attr(res, "strategy") = "chm-dtm"
    return(res)
  }

  if (strategy == "trunks")
  {
    #las = readLAS("~/Téléchargements/las_ref.las")
    #las = readLAS("~/Téléchargements/las_mov.las")

    if(verbose) cat("  Computing Digital Terrain Model...\n")
    dtm = suppressMessages(lidR::rasterize_terrain(las, res = 0.1))
    las = lidR::height_above_ground(las, dtm = dtm, algorithm = lidR::tin())

    if(verbose) cat("  Slicing between 0.15 and 3 m...\n")
    las = lidR::filter_poi(las, hag > 0.15, hag < 3)

    if(verbose) cat("  Decimating with 3 cm resolution...\n")
    vox = lidR:::C_voxel_id(las, 0.03)
    vox = duplicated(vox)
    las = las[!vox]

    if(verbose) cat("  Smoothing 3D...\n")
    las = smooth3d(las, 0.05, ncpu = lidR::get_lidr_threads(), progress = FALSE)

    if(verbose) cat("  Filtering anisotropic features\n")
    las = compute_anisotropy(las, k = 20)
    las = lidR::filter_poi(las, anisotropy > 0.9)

    if(verbose) cat("  Comuting connected component...\n")
    las = lidR::connected_components(las, 0.05, 100)
    las = lidR::filter_poi(las, clusterID > 0)

    if(verbose) cat("  Keeping only 50 most important objects...\n")
    cluster_size = las@data[, .N, by = clusterID]
    data.table::setorder(cluster_size, -N)

    if (nrow(cluster_size) > 50)
    {
      cluster_size = cluster_size[1:50]
      las = lidR::filter_poi(las, clusterID %in% cluster_size$clusterID)
    }

    las = lidR::remove_lasattribute(las, "hag")
    las = lidR::remove_lasattribute(las, "anisotropy")
    las = lidR::remove_lasattribute(las, "clusterID")

    #dtm_points = terra::as.data.frame(dtm, xy = TRUE)
    #names(dtm_points) = c('X', 'Y', 'Z')
    #points = las@data[, .(X,Y,Z)]
    #las = LAS(rbind(points, dtm_points))
    #plot(las, color = "anisotropy", breaks = "quantile", legend = T)

    return(las)
  }
}

compute_anisotropy = function(las, k = 0)
{
  eigen = lidR::point_eigenvalues(las, k = k, metrics = T)
  las = lidR::add_lasattribute_manual(las, eigen$anisotropy, "anisotropy", "anisotropy", "float")
  return(las)
}

smooth3d = function(las, radius = 0.05, weight = NULL, ncpu = 8, progress = TRUE)
{
  ti = Sys.time()

  if (is.function(radius)) radius = radius(las$Z)
  if (length(radius) == 1L) radius = rep(radius, lidR::npoints(las))
  stopifnot(length(radius) == lidR::npoints(las))
  stopifnot(all(radius > 0))
  if (is.null(weight)) weight = rep(1, lidR::npoints(las))
  stopifnot(length(weight) == lidR::npoints(las), radius > 0)

  z = cpp_smooth3d(las, radius = radius, weight = weight, ncpu = ncpu, progress, FALSE)

  tf = Sys.time()
  #print(difftime(tf,ti))

  las$X = z$X
  las$Y = z$Y
  las$Z = z$Z
  lidR::las_quantize(las)
  las
}

