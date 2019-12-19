#' @description Calculate surface root mean square gradient for an input landscape using a moving window analysis
#' using a matrix to specify the area around the focal cell (spatial scale) that you want to be used to make
#' the calculations. Returns a new raster layer.
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

  stopifnot(is.matrix(moving_window))
  d <- dim(moving_window)
  if (prod(d) == 0) { stop('ncol and nrow of moving_window must be > 0') }
  if (min(d %% 2) == 0) { stop('moving_window must have uneven sides') }

  dx <- res(landscape)[1]
  dy <- res(landscape)[2]

  out <- NULL #assign an object for your results data

  #Calculate new cell value based on equation in SPIP user guide and based on whether the call has raster boundaries that allow for the calculation of this value
  for(i in 1:(nrow(landscape))){
    for(j in 1:(ncol(landscape))){

      #Assign the position of the cell in the matrix
      position <- c(i,j)
      rowPosition <- position[1]
      colPosition <- position[2]

      xcoord <- xFromCol(landscape, colPosition)
      ycoord <- yFromRow(landscape, rowPosition)

      #Extract value of the cell in question:
      cell <- landscape[rowPosition, colPosition]

      #Check if a cell has an NA value, if it does keep it as an NA
      if(is.na(cell)) {val1 <- NA }
      else {
        #check that the cell has adjacent cells that can be used to do the calculation:
        if(((rowPosition - 1) > 0) & ((colPosition - 1 ) > 0)) {
          #extract value of the cell to the left:
          cellLeft <- landscape[rowPosition, (colPosition-1)]
          #extract the value of the cell above:
          cellUpper <- landscape[(rowPosition-1), colPosition]
          #calculate the new value of the cell:
          val1 <- ((cell-cellLeft)/dx)^2 + ((cell-cellUpper)/dy)^2
        } else { val1 <- NA } #if there is not adjacent cells assign that cell an NA value

      }
      out <- rbind(out, as.data.frame(cbind(x=xcoord, y=ycoord, val1)))
      remove(position, rowPosition, colPosition, xcoord, ycoord, val1)
    }}

  temprast <- rasterFromXYZ(out)

  newrast <- raster:: focal(temprast, moving_window, fun= function(x) {

    vals <- as.matrix(x) #Convert it to a matrix

    M <- nrow(moving_window) #get the raster x dimension
    N <- ncol(moving_window) #get the raster y dimension

    #Calculate the root mean square gradient!
    Sdq <- sqrt(sum(vals, na.rm=T)/((M-1)*(N-1)))
    return(Sdq)
  })

  if (filename  != '') {
    writeRaster(newrast, filename)
  }
  return(newrast)
}
