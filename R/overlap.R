#' Automatically adjust overlap parameter for ICP
#'
#' @param overlap default value
#' @param r plot radius
#' @param M brute force registration matrix
#' @noRd
adjust_overlap = function(overlap, ref, mov, M)
{
  href = sf::st_convex_hull(ref)
  hmov = sf::st_convex_hull(transform_las(mov, M))
  sf::st_crs(hmov) = sf::st_crs(href)
  poverlap = sf::st_intersection(href, hmov)
  a = sf::st_area(poverlap)
  A = sf::st_area(hmov)
  ans = a/A*100

  #d = sqrt(sum(M[1:2,4]^2))
  #ans = circle_overlap(r, d)
  ans = floor(ans / 10) * 10
  ans = min(overlap, ans)
  ans = ans-10
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