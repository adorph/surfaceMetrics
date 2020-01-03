require(raster)
inrast <- raster("C:/Users/adorph/Documents/2_Otways_EnvironmentalHeterogeneity/0-gradient-surface-analysis/raster_reproj/aridity.tif")
# inrast <- raster("C:/Users/adorph/Documents/2_Otways_EnvironmentalHeterogeneity/0-gradient-surface-analysis/raster_clipped/aridityscale1020_site1.tif")

intSdq6rasteraddition <- function(x){
  require(parallel)
  require(raster)

  #Calculate new cell value based on equation in SPIP user guide and based on whether the call has raster boundaries that allow for the calculation of this value
  parallelforeachcell <- function(i){
    require(raster)

    #Change the raster to matrix format
    mrast <- as.matrix(x)

    # M <- dim(x)[1] #get the raster x dimension
    # N <- dim(x)[2] #get the raster y dimension

    #Assign the position of the cell within the raster:
    cellpos <- which(mrast==mrast[i], arr.ind = T)
    #Assign the position of the cell in the matrix
    rowPosition <- cellpos[1]
    colPosition <- cellpos[2]

    #Check the resolution of the cells for calculations below
    dx <- res(x)[1]
    dy <- res(x)[2]

    #Check if adjacent cell positions we need for the calculation are out of raster boundaries
    upperBound3 = (rowPosition - 3) > 0
    leftBound3 = (colPosition - 3) > 0
    lowerBound3 = (rowPosition + 3) <= nrow(mrast)
    rightBound3 = (colPosition + 3) <= ncol(mrast)

    #Extract value of the cell in question:
    cell <- mrast[rowPosition, colPosition]

    #Check if a cell has an NA value, if it does keep it as an NA
    if(is.na(cell)){
      val1 <- NA
    } else {
      #check that the cell has adjacent cells that can be used to do the calculation:
      if(upperBound3 & leftBound3 & rightBound3 & lowerBound3) {

        #extract value of the cells to the left:
        cellLeft1 <- mrast[rowPosition, (colPosition-1)]
        cellLeft2 <- mrast[rowPosition, (colPosition-2)]
        cellLeft3 <- mrast[rowPosition, (colPosition-3)]
        #extract value of the cells to the right:
        cellRight1 <- mrast[rowPosition, (colPosition+1)]
        cellRight2 <- mrast[rowPosition, (colPosition+2)]
        cellRight3 <- mrast[rowPosition, (colPosition+3)]
        #extract the value of the cells above:
        cellUpper1 <- mrast[(rowPosition-1), colPosition]
        cellUpper2 <- mrast[(rowPosition-2), colPosition]
        cellUpper3 <- mrast[(rowPosition-3), colPosition]
        #extract the value of the cells above:
        cellLower1 <- mrast[(rowPosition+1), colPosition]
        cellLower2 <- mrast[(rowPosition+2), colPosition]
        cellLower3 <- mrast[(rowPosition+3), colPosition]
        #calculate the new value of the cell:
        val1 <- sqrt(((-cellLeft3 + 9*cellLeft2 - 45*cellLeft1 +
                  45*cellRight1-9*cellRight2 + cellRight3) /
                  60*(dx*3))^2 +
                ((-cellUpper3 + 9*cellUpper2 - 45*cellUpper1 +
                  45*cellLower1 - 9*cellLower2 + cellLower3) /
                  60*(dx*3))^2)
      } else {
          val1 <- NA
          }
    }
    return(val1)
  }

  cl <- makeCluster(detectCores())
  clusterEvalQ(cl, {library(raster)}) #reads in the libraries you may need within the clusters
  clusterExport(cl, varlist = c("x"), envir=environment())
  r_val <- clusterApply(cl, 1:ncell(x), parallelforeachcell)
  stopCluster(cl)

  out <- matrix(nrow=x@nrows, ncol=x@ncols, unlist(r_val))

  #Add all the values in the sample matrix together:
  # sum.mat <- sum(sampmatrix, na.rm=T)

  #Calculate te root mean square gradient!
  # Sdq <- sqrt(sum.mat/((M-1)*(N-1)))

  return(out)
}

aridity <- intSdq6rasteraddition(inrast)
write.csv(as.data.frame(aridity), "./aridity_intermediate.csv")
