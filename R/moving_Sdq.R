#' @description Calculate surface root mean square gradient for an input landscape through a moving window analysis
#' using a matrix to specify the area around the focal cell (spatial scale) that you want to be used to make
#' the calculations. This function calls a C++ function through the Rcpp package. Returns a new raster layer.
#'
#' @param landscape Raster* Layer
#' @param moving_window A matrix (the moving window) e.g. 3 x 3 matrix with values of 1. The matrix foes not need
#' to be square, but the sides do need to be odd numbers. Equivalent to scale of measurement.
#' @param filename (optional) Character. A filename for the new raster output.
#'
#' @details
#' This function calculates surface root mean square gradient values for each cell by applying a moving window of
#' a specified size to the area surrounding each cell. The function uses the equation outlined by Stout et al. (1993).
#' The scale at which the metrics area calculated for each focal cell will be determined by the resolution of the raster.
#' It is slow to run for large rasters.
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
#' mySa <- moving_Sdq(landscape, mywindow)
#' plot(mySdq)
#'
#' #Outputing it to a file
#' mySdq2 <- moving_Sdq(landscape, mywindow, filename="./mySdq.tiff")
#' }
#'
#' @references
#' Stout, K.J., Sullivan, W.P., Mainsah, E., Luo, N., Mathia, T., Zahouani, H., 1993. The Development of Methods
#' for the Characterisation of Roughness in Three Dimensions. Publication no. EUR 15178 EN of the Commission of
#' the European Communities Dissemination of Scientific and Technical Knowledge Unit, Luxembourg.


moving_Sdq <- function(landscape, moving_window, filename=''){
  require(raster)
  require(Rcpp)

  stopifnot(is.matrix(moving_window))
  d <- dim(moving_window)
  if (prod(d) == 0) { stop('ncol and nrow of moving_window must be > 0') }
  if (min(d %% 2) == 0) { stop('moving_window must have uneven sides') }

  #Load a C++ function through Rcpp to quickly compute the cell adjacencies
  cppFunction('NumericMatrix shiftSdq(const NumericMatrix& m1){
              NumericMatrix out(m1.nrow(),m1.ncol());

              for(size_t i=0; i <= (m1.nrow() - 1); i++) {
                for(size_t j=0; j <= (m1.ncol() - 1); j++) {
                  if (std::isnan(m1(i,j))) {
                    out(i,j) = std::nan("");
                  } else {
                    out(i,j) = pow((m1(i,j) - m1(i-1,j)), 2) + pow((m1(i,j) - m1(i,j-1)), 2);
                  }
                }
              }
              return out;
            }')

  #Convert the input raster to a matrix for the C++ function
  mrast <- as.matrix(landscape)
  #Run the C++ function
  out <- shiftSdq(mrast)
  #Convert the matrix output to a raster with the same extent and crs as the input raster
  temprast <- raster(out)
  extent(temprast) <- extent(landscape)
  crs(temprast) <- crs(landscape)

  #Run the moving window analysis on this raster output
  newrast <- raster:: focal(temprast, moving_window, fun= function(x) {

    vals <- as.matrix(x) #Convert it to a matrix

    M <- nrow(moving_window) #get the raster x dimension
    N <- ncol(moving_window) #get the raster y dimension

    #Calculate the root mean square gradient
    Sdq <- sqrt(sum(vals, na.rm=T)/((M-1)*(N-1)))
    return(Sdq)
  })

  frast <- crop(newrast, landscape)

  if (filename  != '') {
    writeRaster(frast, filename, format="GTiff", overwrite=T)
  }
  return(frast)
}
