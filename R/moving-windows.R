remove(list=ls())

require(raster)
xrast <-raster("C:/Users/adorph/Documents/GIS/GIS_Otways/SDRIVE_StudyArea/DEM/dem/hdr.adf")

#Define your raster
xr <- crop(xrast, extent(xrast, 1, 5, 1, 5))
#Define your moving window size
moving_window <- matrix(1, nrow = 5, ncol = 5)

moving_Sku <- function(landscape, moving_window, filename=''){
  #"landscape" is your whole raster.
  #"moving window" is the area around each pixel that you want to take into account in the calculation of
  #your metric (spatial scale), given as a matrix. It has to have an odd number of sides or will error.
  #You also have the option to write the output raster to a file using "filename"

  require(raster)

  stopifnot(is.matrix(moving_window))
  d <- dim(moving_window)
  if (prod(d) == 0) { stop('ncol and nrow of moving_window must be > 0') }
  if (min(d %% 2) == 0) { stop('moving_window must have uneven sides') }

  t <- focal(landscape, moving_window, fun=function(x){
    MN <- nrow(moving_window) * ncol(moving_window) #Get the moving window dimensions
    vals <- as.matrix(x) #Convert it to a matrix
    mu <- mean(vals) #Find the mean
    #Calculate the root mean square error
    Sq <- sqrt(sum((vals - mu)^2, na.rm=T)/MN)
    #calculate surface kurtosis
    Sku <- (sum((vals - mu)^4, na.rm=T)/(MN*(Sq^4)))
    return(Sku)
  })
  if (filename  != '') {
    writeRaster(t, filename)
  }
  return(t)
}

plot(moving_Sku(xr, moving_window))

moving_Sdq <- function(landscape, moving_window, filename=''){
  require(raster)

  stopifnot(is.matrix(moving_window))
  d <- dim(moving_window)
  if (prod(d) == 0) { stop('ncol and nrow of moving_window must be > 0') }
  if (min(d %% 2) == 0) { stop('moving_window must have uneven sides') }

  dx <- res(x)[1]
  dy <- res(x)[2]

  out <- NULL #assign an object for your results data

  #Calculate new cell value based on equation in SPIP user guide and based on whether the call has raster boundaries that allow for the calculation of this value
  for(i in 1:(nrow(landscape))){
    for(j in 1:(ncol(landscape))){

      #Assign the position of the cell in the matrix
      position <- c(i,j)
      rowPosition <- position[1]
      colPosition <- position[2]

      xcoord <- rasterToPoints(landscape)[i,1]
      ycoord <- rasterToPoints(landscape)[j,1]

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
      out <- as.data.frame(cbind(xcoord, ycoord, val1))
    }}


  newrast <- raster:: focal(out, moving_window, fun= function(x) {

    vals <- as.matrix(x) #Convert it to a matrix

    M <- nrow(moving_window) #get the raster x dimension
    N <- ncol(moving_window) #get the raster y dimension

    #Calculate the root mean square gradient!
    Sdq <- sqrt(sum(vals, na.rm=T)/((M-1)*(N-1)))
    return(Sdq)
    })

  # if (filename  != '') {
  #   writeRaster(newrast, filename)
  # }
  return(newrast)
}

as.matrix(moving_Sdq(xr, moving_window))

