get_os <- function()
{
  sysinf <- Sys.info()
  if (!is.null(sysinf)){
    os <- sysinf['sysname']
    if (os == 'Darwin')
      os <- "osx"
  } else { ## mystery machine
    os <- .Platform$OS.type
    if (grepl("^darwin", R.version$os))
      os <- "osx"
    if (grepl("linux-gnu", R.version$os))
      os <- "linux"
  }
  tolower(os)
}

#' Find CloudCompare executable path
#'
#' This function attempts to locate the CloudCompare executable based on the operating system.
#' On Windows, it queries the registry to find the installation path. On macOS and Linux (Flathub version),
#' it returns an error. If the path is not found automatically use `register_cloudcompare()`
#'
#' @return A character string with the path to the CloudCompare executable.
#' @noRd
find_cloudcompare <- function()
{
  if (exists("cc", envir = .ccpath_env, inherits = FALSE))
  {
    path <- get("cc", envir = .ccpath_env)
    if (!is.null(path) && nzchar(path))
    {
      if (!file.exists(path))
        stop(paste("Cannot find CloudCompare in", path))

      return(path)
    }
  }

  if (Sys.info()["nodename"] == "FFGG-1009803") {
    return("/home/jr/Logiciels/CloudCompare/build/qCC/CloudCompare")
  }

  os <- get_os()

  if (os == "linux") {
    stop("Linux version of CloudCompare is a Flathub package. Please provide a compiled path using register_cloudcompare().")
  }

  if (os == "windows") {
    path <- "C:/Program Files/CloudCompare/CloudCompare.exe"
    if (!file.exists(path)) {
      stop(paste("Cannot find CloudCompare in", path))
    }
    return(shQuote(path))
  }

  if (os == "osx") {
    stop("Finding CloudCompare on macOS is not yet supported. Please provide the path using register_cloudcompare().")
  }

  stop("OS not detected")
}

#' @param path A string with full path to the CloudCompare binary
#' @rdname find_cloudcompare
register_cloudcompare <- function(path)
{
  assign("cc", path, envir = .ccpath_env)
}

.ccpath_env <- new.env(parent = emptyenv())

#' CloudCompare ICP Registration
#'
#' This function performs point cloud registration using the CloudCompare ICP (Iterative Closest Point) algorithm.
#' CloudCompare must be installed to use this function. You can download it from: https://www.danielgm.net/cc/.
#'
#' @param vref A `LAS` containing the reference point cloud or nx3 matrix
#' @param umov A `LAS` containing the moving point cloud or nx3 matrix
#' @param min_error_diff Numeric. The minimum error difference for the ICP algorithm to terminate (default: 1e-5).
#' @param overlap Numeric. to specify the percentage of (final) overlap (integer number between 10 and 100 - default = 100)
#' @param rot Character. 'XYZ' or 'X' or 'Y' or 'Z' or 'NONE' to constrain the rotation around a given axis
#' @param skip_txy,skip_tz Logical. Block XY or Z translations.
#' (or no rotation at all)
#' @param cc Character. The path to the CloudCompare executable
#' @param verbose logical
#'
#' @return A 4x4 transformation matrix resulting from the ICP registration.
#'
#' @details The function exports both point clouds to temporary `.xyz` files, runs CloudCompare ICP
#' via a system call, retrieves the transformation matrix from the generated file, and then deletes
#' temporary files.
#' @md
#' @noRd
cc_icp = function(vref, umov, min_error_diff = 1e-5, overlap = 90, rot = "Z", skip_txy = FALSE, skip_tz = FALSE, cc = find_cloudcompare(), verbose = TRUE)
{
  . <- X <- Y <- Z <- NULL

  on.exit(file.remove(fmatrix))

  stopifnot(overlap <= 100, overlap >= 10)
  rot = match.arg(rot, c('XYZ','X','Y','Z','NONE'))

  ovref = paste0(tempdir(),"/vref.xyz")
  oumov = paste0(tempdir(),"/umov.xyz")

  if (is.matrix(vref)) vref = data.table::as.data.table(vref)
  if (is.matrix(umov)) umov = data.table::as.data.table(umov)
  if (methods::is(vref, "LAS")) vref = vref@data[, .(X,Y,Z)]
  if (methods::is(umov, "LAS")) umov = umov@data[, .(X,Y,Z)]

  data.table::fwrite(vref, ovref, sep = " ", col.names = FALSE)
  data.table::fwrite(umov, oumov, sep = " ", col.names = FALSE)

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
