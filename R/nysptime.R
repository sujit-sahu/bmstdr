#' Daily 8-hour maximum ozone concentration values and three covariates from 28 sites in New York for
#' the 62 days during the months of July and August, 2006.
#'
#' @source It is the same as the NYdata set in the spTimer package
#' with two added columns providing the UTM X- and Y- coordinates.
#'  Each data row is  for a particular site and day as detailed below.
#'
#'  @format A data frame with 1736 rows and 12 columns:
#' \describe{
#'   \item{s.index}{site index (1 to 28)}
#'   \item{Longitude}{Longitude of the site}
#'   \item{Latitude}{Latitude of the site}
#'   \item{utmx}{UTM X-coordinate of the site}
#'   \item{utmy}{UTM Y-coordinate of the site}
#'   \item{Year}{This is 2006 for all the rows}
#'   \item{Month}{Month taking values 7 for July and 8 for August}
#'   \item{Day}{Day taking values 1 to 31}
#'   \item{y8hrmax}{Daily 8-hour maximum ozone concentration value}
#'   \item{xmaxtemp}{Maximum temperature (degree Celsius)}
#'   \item{xwdsp}{Windspeed (nautical mile per hour)}
#'   \item{xrh}{Relative humidity}
#' }
# usethis::use_data(nysptime, overwrite = TRUE)
#' @examples
#' \dontrun{
#'  nysptime
#'  summary(nysptime[, 9:12])
#' }
"nysptime"

#' Values of three covariates for 100 grid locations in New York for
#' the 62 days during the months of July and August, 2006.
#'
#' @source It is the same data set as NYgrid data set in the spTimer package
#' with two added columns providing the UTM X- and Y- coordinates.
#'  Each data row is  for a particular grid site and day as detailed below.
#'
#'  @format A data frame with 6200 rows and 11 columns:
#' \describe{
#'   \item{s.index}{site index (1 to 100)}
#'   \item{Longitude}{Longitude of the site}
#'   \item{Latitude}{Latitude of the site}
#'   \item{utmx}{UTM X-coordinate of the site}
#'   \item{utmy}{UTM Y-coordinate of the site}
#'   \item{Year}{This is 2006 for all the rows}
#'   \item{Month}{Month taking values 7 for July and 8 for August}
#'   \item{Day}{Day taking values 1 to 31}
#'   \item{xmaxtemp}{Maximum temperature (degree Celsius)}
#'   \item{xwdsp}{Windspeed (nautical mile per hour)}
#'   \item{xrh}{Relative humidity}
#' }
# usethis::use_data(gridnysptime, overwrite = TRUE)
#' @examples
#' \dontrun{
#'  nysptime
#'  summary(gridnysptime[, 9:11])
#' }
"gridnysptime"
