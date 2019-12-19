sm.Sa <- function(x, useMean = TRUE) {
  require(raster)
  d1 <- getValues(x) #extract values as vector from raster as spatial location does not matter
  mu <- mean(d1, na.rm=T) #calculate mean of raster values
  MN <- dim(x)[1]*dim(x)[2] #Multiply the raster x dimension by the raster y dimension

  out <- NULL #set up an empty object to recieve values
  #If useMean=TRUE calculate the cell value minus the mean value for all the cells in the matrix.
  #If useMean=FALSE the mean will not be subtracted.
  for(i in 1:length(d1)){
    out[i] <- ifelse(useMean==TRUE, d1[i]- mu, d1[i])
  }

  sumval <- sum(out) #Add together all the cell values from the output matrix
  Sa <- (sumval/MN) #Calculate the root mean square value
  return(Sa)
}
