#' Average ozone concentration values and three covariates from 28 sites in New York.
#'
#' @source See the NYdata set in spTimer package.
#'  Each data row is the mean of the available daily 8-hour maximum average ozone concentrations
#'  in parts per billion (ppb) at each of the
#'  28 sites. The daily values are for the months of July and August in 2006.
#'
#'  @format A data frame with 28 rows and 9 columns:
#' \describe{
#'   \item{s.index}{site index (1 to 28)}
#'   \item{Longitude}{Longitude of the site}
#'   \item{Latitude}{Latitude of the site}
#'   \item{utmx}{UTM X-coordinate of the site}
#'   \item{utmy}{UTM Y-coordinate of the site}
#'   \item{yo3}{Average ozone concentration value (ppb) at the site over 62 days in July and August, 2006}
#'   \item{xmaxtemp}{Average maximum temperature (degree celsius) at the site over 62 days in  July and August, 2006}
#'   \item{xwdsp}{Average windspeed (nautical mile per hour) over 62 days in  July and August, 2006}
#'   \item{xrh}{Average relative humidity over 62 days in  July and August, 2006}
#' }
# usethis::use_data(nyspatial, overwrite = TRUE)
#' @examples
#' \dontrun{
#'  nyspatial
#'  summary(nyspatial)
#'  pairs(nyspatial[, 6:9])
#' }
"nyspatial"
