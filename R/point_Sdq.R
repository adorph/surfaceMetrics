#' @description This function takes an input raster that has been clipped to the area
#' around a site that you want to be measured and calculates surface root mean square gradient for
#' that area. It returns a single value for that layer.
#'
#' @param landscape Raster* Layer clipped to the scale you want to measure
#'
#' @details
#' This function calculates a surface root mean square gradient for the supplied raster.
#' The function uses the equation outlined by Stout et al. (1993).
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
#' mySdq <- point_Sdq(croppedland)
#' mySdq
#'
#' @references
#' Stout, K.J., Sullivan, W.P., Mainsah, E., Luo, N., Mathia, T., Zahouani, H., 1993. The Development of Methods
#' for the Characterisation of Roughness in Three Dimensions. Publication no. EUR 15178 EN of the Commission of
#' the European Communities Dissemination of Scientific and Technical Knowledge Unit, Luxembourg.

point_Sdq <- function(landscape){
  require(raster)

  mrast <- as.matrix(landscape) #Change the raster to matrix format

  M <- dim(landscape)[1] #get the raster x dimension
  N <- dim(landscape)[2] #get the raster y dimension

  #Check the resolution of the cells for calculations below
  dx <- res(landscape)[1]
  dy <- res(landscape)[2]

  out <- matrix(nrow=nrow(mrast), ncol=ncol(mrast)) #assign an object for your results data

  #Calculate new cell value based on equation in SPIP user guide and based on whether the call has raster boundaries that allow for the calculation of this value
  for(i in 1:(nrow(mrast))){
    for(j in 1:(ncol(mrast))){

      #Assign the position of the cell in the matrix
      position = c(i,j)
      rowPosition <- position[1]
      colPosition <- position[2]

      #Check if adjacent cell positions we need for the calculation are out of raster boundaries
      upperBound1 = (rowPosition - 1) > 0
      leftBound1 = (colPosition - 1 ) > 0

      #Extract value of the cell in question:
      cell <- mrast[rowPosition, colPosition]

      #Check if a cell has an NA value, if it does keep it as an NA
      if(is.na(cell)){
        out[i,j] <- NA

      } else {
        #check that the cell has adjacent cells that can be used to do the calculation:
        if(upperBound1 & leftBound1) {
          #extract value of the cell to the left:
          cellLeft <- mrast[rowPosition, (colPosition-1)]
          #extract the value of the cell above:
          cellUpper <- mrast[(rowPosition-1), colPosition]
          #calculate the new value of the cell:
          out[i,j] <- ((cell-cellLeft)/dx)^2 + ((cell-cellUpper)/dy)^2

        } else {
          #if there is not adjacent cells assign that cell an NA value - this means that we lose two of the edges in the matrix. If I were to do this as a moving window then you wouldn't lose as many values... Laters problem
          out[i,j] <- NA
        }
      }
    }
  }

  #Calculate te root mean square gradient!
  Sdq <- sqrt(sum(out, na.rm=T)/((M-1)*(N-1)))

  return(Sdq)
}
