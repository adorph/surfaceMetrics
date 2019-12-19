#' @description Calculate surface root mean square for an input landscape using a moving window analysis
#' using a matrix to specify the area around the focal cell (spatial scale) that you want to be used to make
#' the calculations. Returns a new raster layer.
#'
#' @param landscape Raster* Layer
#' @param moving_window A matrix (the moving window) e.g. 3 x 3 matrix with values of 1. The matrix foes not need
#'  to be square, but the sides do need to be odd numbers. Equivalent to scale of measurement.
#' @param filename (optional) Character. A filename for the new raster output.
#'
#' @details
#' This function calculates surface root mean square values for each cell by applying a moving window of a specified
#' size to the area surrounding each cell. The function uses the equation outlined by Stout et al. (1993). The scale at
#' which the metrics area calculated for each focal cell will be determined by the resolution of the raster.
#'
#' @examples
#' \dontrun{
#' # Load a raster of your landscape into R
#' file <- system.file("external/test.grd", package="raster")
#' mylandscape <- raster(file)
#'
#' # Create a matrix of the area in the landscape you want to use to calculate your variable
#' mywindow <- matrix(1, nrow = 5, ncol = 5)
#'
#' #Apply the function
#' mySq <- moving_Sq(landscape, mywindow)
#' plot(mySq)
#'
#' #Outputing it to a file
#' mySq2 <- moving_Sq(landscape, mywindow, filename="./mySq.tiff")
#' }
#'
#' @references
#' Stout, K.J., Sullivan, W.P., Mainsah, E., Luo, N., Mathia, T., Zahouani, H., 1993. The Development of Methods
#' for the Characterisation of Roughness in Three Dimensions. Publication no. EUR 15178 EN of the Commission of
#' the European Communities Dissemination of Scientific and Technical Knowledge Unit, Luxembourg.

moving_Sq <- function(landscape, moving_window, filename=''){
  require(raster)

  stopifnot(is.matrix(moving_window))
  d <- dim(moving_window)
  if (prod(d) == 0) { stop('ncol and nrow of moving_window must be > 0') }
  if (min(d %% 2) == 0) { stop('moving_window must have uneven sides') }

  newrast <- focal(landscape, moving_window, fun=function(x){
    MN <- nrow(moving_window) * ncol(moving_window) #Get the moving window dimensions
    vals <- as.matrix(x) #Convert it to a matrix
    mu <- mean(vals) #Find the mean
    #Calculate the root mean square error
    Sq <- sqrt(sum((vals - mu)^2, na.rm=T)/MN)
    return(Sq)
  })
  if (filename  != '') {
    writeRaster(newrast, filename)
  }
  return(newrast)
}
