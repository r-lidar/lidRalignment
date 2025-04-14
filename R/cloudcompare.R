#' Find CloudCompare executable path
#'
#' This function attempts to locate the CloudCompare executable based on the operating system.
#' On Windows, it queries the registry to find the installation path. On macOS and Linux (Flathub version),
#' it returns an error. If the path is not found automatically use `register_cloudcompare()`
#'
#' @return A character string with the path to the CloudCompare executable.
#' @export
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
