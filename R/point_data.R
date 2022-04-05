#' Average air pollution values from 28 sites in New York.
#'
#' @source This is obtained by calculating site-wise averages 
#' of the NYdata set in the 'spTimer' package
#' \insertCite{spTimer;textual}{bmstdr}.
#'  Each data point is the mean of the available daily 8-hour maximum average ozone concentrations
#'  in parts per billion (ppb) at each of the
#'  28 sites. The daily values are for the month of July and August in 2006.
#' @format A vector with 28 real values.
## library(spTimer)
## attach(NYdata)
## ydata <- as.vector(tapply(o8hrmax, INDEX=s.index,  FUN=mean, na.rm=TRUE))
ydata <- c(45.08065,  48.34758,  54.65677,  44.95774,  43.92903,  55.37726,  51.08589,
44.13661,  44.82355,  41.24903,  45.86645,  45.88557,  48.93226,  46.31767,  49.96500,
48.18387,  45.64532,  39.87500,  54.25150,  57.59550,  48.31339,  46.14145,  52.99596,
48.16083,  44.35355,  39.18065,  52.80557,  52.39814)
# usethis::use_data(ydata, overwrite = TRUE)
#' @references
#' \insertAllCited{}
#' @examples
#'  summary(ydata)
"ydata"
#' Average ozone concentration values and three covariates from 28 sites in New York.
#'
#' @source See the NYdata set in spTimer package, \insertCite{spTimer;textual}{bmstdr}.
#'  Each data row is the mean of the available daily 8-hour maximum average ozone 
#'  concentrations in parts per billion (ppb) at each of the
#'  28 sites. The daily values are for the months of July and August in 2006.
#'  @format A data frame with 28 rows and 9 columns:
#' \describe{
#'   \item{s.index}{site index (1 to 28)}
#'   \item{Longitude}{Longitude of the site}
#'   \item{Latitude}{Latitude of the site}
#'   \item{utmx}{UTM X-coordinate of the site}
#'   \item{utmy}{UTM Y-coordinate of the site}
#'   \item{yo3}{Average ozone concentration value (ppb) at the site over 62 days in July and August, 2006}
#'   \item{xmaxtemp}{Average maximum temperature (degree Celsius) at the site over 62 days in  July and August, 2006}
#'   \item{xwdsp}{Average wind speed (nautical mile per hour) over 62 days in  July and August, 2006}
#'   \item{xrh}{Average relative humidity over 62 days in  July and August, 2006}
#' }
# usethis::use_data(nyspatial, overwrite = TRUE)
#' @references
#' \insertAllCited{}
#' @examples
#'  head(nyspatial)
#'  summary(nyspatial)
#'  pairs(nyspatial[, 6:9])
"nyspatial"
#' Daily 8-hour maximum ozone concentration values and three covariates from 28 sites in New York for
#' the 62 days during the months of July and August, 2006.
#'
#' @source It is the same as the NYdata set in the spTimer package, 
#' \insertCite{spTimer;textual}{bmstdr}, 
#' with two added columns providing the UTM X- and Y- coordinates.
#'  Each data row is  for a particular site and a day as detailed below.
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
#'   \item{xwdsp}{wind speed (nautical mile per hour)}
#'   \item{xrh}{Relative humidity}
#' }
# usethis::use_data(nysptime, overwrite = TRUE)
#' @references
#' \insertAllCited{}
#' @examples
#' head(nysptime)
#' summary(nysptime[, 9:12])
"nysptime"

#' Values of three covariates for 100 grid locations in New York for
#' the 62 days during the months of July and August, 2006.
#' @source It is the same data set as NYgrid data set in the spTimer package, 
#' \insertCite{spTimer;textual}{bmstdr}.
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
#'   \item{xwdsp}{wind speed (nautical mile per hour)}
#'   \item{xrh}{Relative humidity}
#' }
# usethis::use_data(gridnysptime, overwrite = TRUE)
#' @references
#' \insertAllCited{}
#' @examples
#'  head(gridnysptime)
#'  summary(gridnysptime[, 9:11])
"gridnysptime"
#' Values of three covariates for 100 grid locations in New York
#' averaged over the 62 days during the months of July and August, 2006.
#' @source See the NYgrid data  set in spTimer package, 
#' \insertCite{spTimer;textual}{bmstdr}.
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
#'   \item{xmaxtemp}{Average maximum temperature (degree Celsius) at the site over 62 days in  July and August, 2006}
#'   \item{xwdsp}{Average wind speed (nautical mile per hour) over 62 days in  July and August, 2006}
#'   \item{xwdsp}{Average relative humidity over 62 days in  July and August, 2006}
#' }
#usethis::use_data(gridnyspatial, overwrite = TRUE)
#' @references
#' \insertAllCited{}
#' @examples
#'  summary(gridnyspatial)
"gridnyspatial"
#' Values of three covariates for 100 grid locations in New York for
#' the 62 days during the months of July and August, 2006.
#'
#' @source It is the same data set as NYgrid data set in the spTimer package, 
#' 
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
#'   \item{xwdsp}{wind speed (nautical mile per hour)}
#'   \item{xrh}{Relative humidity}
#' }
# usethis::use_data(gridnysptime, overwrite = TRUE)
#' @references
#' \insertAllCited{}
#' @examples
#'  summary(gridnysptime[, 9:11])
"gridnysptime"

#' The color palette used to draw maps to illustrate the package 
#' bmstdr, see  \insertCite{Sahubook;textual}{bmstdr}
#' It has the values in order: dodgerblue4, dodgerblue2, firebrick2,  
#' firebrick4 and purple.     
# usethis::use_data(colpalette, overwrite = TRUE)
#' @examples
#'  colpalette
#' @references
#' \insertAllCited{}
"colpalette"
