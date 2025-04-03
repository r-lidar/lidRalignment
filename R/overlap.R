#' Automatically adjust overlap parameter for ICP
#'
#' @param overlap default value
#' @param r plot radius
#' @param M brute force registration matrix
#' @export
adjust_overlap = function(overlap, r, M)
{
  d = sqrt(sum(M[1:2,4]^2))
  ans = circle_overlap(r, d)
  ans = floor(ans / 10) * 10
  ans = min(overlap, ans)
  ans = max(ans, 10)
  ans
}


circle_overlap <- function(r, d)
{
  if (d >= 2 * r) { return(0) }  # No overlap
  if (d == 0) { return(100) }    # Full overlap
  A <- 2 * r^2 * acos(d / (2 * r)) - (d / 2) * sqrt(4 * r^2 - d^2)
  return(A/(pi*r^2)*100)
}