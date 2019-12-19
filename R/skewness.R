sm.Ssk <- function(x, useMean=TRUE) {
  require(raster)
  source("C:/Users/adorph/Documents/R/surfaceMetrics/R/root-mean-square.R")

  d1 <- getValues(x) #extract values from raster
  mu <- mean(d1, na.rm=T) #calculate mean of raster values
  MN <- dim(x)[1]*dim(x)[2] #Multiply the raster x dimension by the raster y dimension

  Sq <- ifelse(useMean==TRUE, sm.Sq(x, useMean=TRUE), sm.Sq(x, useMean=FALSE))

  out <- NULL #set up an empty matrix to recieve values
  for(i in 1:length(d1)) { #calculate the cell value minus the mean value for all the cells in the matrix and square this number
    out[i] <- ifelse(useMean==TRUE, (d1[i]- mu)^3, (d1[i])^3)
  }

  Ssk <- (sum(out, na.rm=T)/(MN*(Sq^3))) #Add together all the cell values from the second output matrix and calculate surface skewness
  return(Ssk)
}
