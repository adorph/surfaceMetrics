sm.Sdq <- function(raster){
  require(raster)
  require(tidyr)
  mrast <- as.matrix(raster) #Change the raster to matrix format

  M <- dim(raster)[1] #get the raster x dimension
  N <- dim(raster)[2] #get the raster y dimension

  #Check the resolution of the cells for calculations below
  dx <- res(raster)[1]
  dy <- res(raster)[2]

  out <- NULL #assign an object for your results data

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
        val1 <- NA

      } else {
        #check that the cell has adjacent cells that can be used to do the calculation:
        if(upperBound1 & leftBound1) {
          #extract value of the cell to the left:
          cellLeft <- mrast[rowPosition, (colPosition-1)]
          #extract the value of the cell above:
          cellUpper <- mrast[(rowPosition-1), colPosition]
          #calculate the new value of the cell:
          val1 <- ((cell-cellLeft)/dx)^2 + ((cell-cellUpper)/dy)^2

        } else {
          #if there is not adjacent cells assign that cell an NA value - this means that we lose two of the edges in the matrix. If I were to do this as a moving window then you wouldn't lose as many values... Laters problem
          val1 <- NA
        }
      }
      res <- cbind(rowPosition, colPosition, val1)
      out <- rbind(out, as.data.frame(res))
    }
  }
  #Convert the data frame back into a matrix format:
  data_wide <- spread(out, colPosition, val1)
  data_wide <- data.frame(data_wide[,-1], row.names=data_wide[,1])
  sampmatrix <- as.matrix(data_wide)

  #Add all the values in the sample matrix together:
  sum.mat <- sum(sampmatrix, na.rm=T)

  #Calculate te root mean square gradient!
  Sdq <- sqrt(sum.mat/((M-1)*(N-1)))

  return(Sdq)
}
