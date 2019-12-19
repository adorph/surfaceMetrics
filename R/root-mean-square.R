sm.Sq <- function(x, useMean=TRUE){
  require(raster)
  d1 <- (getValues(x)) #extract values from raster
  mu <- mean(d1, na.rm=T) #calculate mean of raster values
  MN <- dim(x)[1]*dim(x)[2] #Multiply the raster x dimension by the raster y dimension

  out <- NULL #set up an empty matrix to recieve values
  #If useMean=TRUE calculate the cell value minus the mean value for all the cells in the matrix and square this number.
  #If useMean=FALSE the mean will not be subtracted.
  for(i in 1:length(d1)){
    out[i] <- ifelse(useMean==TRUE, (d1[i]- mu)^2, (d1[i])^2)
  }
  Sq <- sqrt(sum(out, na.rm=T)/MN) #Add together all the cell values from the output matrix and calculate the root mean square value
  return(Sq)
}
