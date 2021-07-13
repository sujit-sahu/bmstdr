#' Values of three covariates for 100 grid locations in New York
#' averaged over the 62 days during the months of July and August, 2006.
#' @source See the NYgrid data  set in spTimer package.
#'  Each data row is the mean of the available covariate
#'  values at the 100 grid locations during the months of July and August in 2006.
#'
#'  @format A data frame with 100 rows and 8 columns:
#' \describe{
#'   \item{s.index}{site index (1 to 28)}
#'   \item{Longitude}{Longitude of the site}
#'   \item{Latitude}{Latitude of the site}
#'   \item{utmx}{UTM X-coordinate of the site}
#'   \item{utmy}{UTM Y-coordinate of the site}
#'   \item{xmaxtemp}{Average maximum temperature (degree celsius) at the site over 62 days in  July and August, 2006}
#'   \item{xwdsp}{Average windspeed (nautical mile per hour) over 62 days in  July and August, 2006}
#'   \item{xwdsp}{Average relative humidity over 62 days in  July and August, 2006}
#' }
#usethis::use_data(gridnyspatial, overwrite = TRUE)
#' @examples
#' \dontrun{
#'  gridnyspatial
#'  summary(gridnyspatial)
#' }
"gridnyspatial"
