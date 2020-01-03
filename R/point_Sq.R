#' @description This function takes an input raster that has been clipped to the area
#' around a site that you want to be measured and calculates surface root mean square for
#' that area. It returns a single value for that layer.
#'
#' @param landscape Raster* Layer clipped to the scale you want to measure
#'
#' @details
#' This function calculates a surface root mean square value for supplied raster. The function
#' uses the equation outlined by Stout et al. (1993).
#'
#' @examples
#' \dontrun{
#' # Load a raster of your landscape into R
#' file <- system.file("external/test.grd", package="raster")
#' mylandscape <- raster(file)
#'
#' # Crop the raster to a section of the landscape you want to calculate surface
#' #kurtosis for
#' croppedland <- crop(mylandscape, extent(mylandscape, 67, 77, 30,40))
#'
#' #Apply the function
#' mySq <- point_Sq(croppedland)
#' mySq
#'
#' @references
#' Stout, K.J., Sullivan, W.P., Mainsah, E., Luo, N., Mathia, T., Zahouani, H., 1993. The Development of Methods
#' for the Characterisation of Roughness in Three Dimensions. Publication no. EUR 15178 EN of the Commission of
#' the European Communities Dissemination of Scientific and Technical Knowledge Unit, Luxembourg.

point_Sq <- function(landscape){
  require(raster)

  vals <- as.matrix(landscape) #Convert it to a matrix
  MN <- dim(landscape)[1]*dim(landscape)[2] #Multiply the raster x dimension by the raster y dimension

  #Calculate the root mean square error
  Sq <- sqrt(sum((vals - mean(vals, na.rm=T))^2, na.rm=T)/MN)
  return(Sq)
}
