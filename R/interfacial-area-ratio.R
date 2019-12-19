#INCOMPLETE
#In development

require(raster)
# inrast <- raster("C:/Users/adorph/Documents/2_Otways_EnvironmentalHeterogeneity/0-gradient-surface-analysis/raster_reproj/aridity.tif")
x <- raster("C:/Users/adorph/Documents/2_Otways_EnvironmentalHeterogeneity/0-gradient-surface-analysis/raster_clipped/demscale1020_site154.tif")


mrast <- as.matrix(x)

dx <- res(x)[1]
dy <- res(x)[2]

calcAij  <- function(i){
  position <- which(mrast == mrast[i], arr.ind = T)

  rowPosition <- position[1]
  colPosition <- position[2]

  lowerBound <- rowPosition + 1 <= nrow(mrast)
  rightBound <- colPosition + 1 <= ncol(mrast)

  if(lowerBound & rightBound) {
    cell <- mrast[rowPosition, colPosition]
    cellLower <- mrast[(rowPosition + 1), colPosition]
    cellRight <- mrast[rowPosition, (colPosition + 1)]
    cellDiag <- mrast[(rowPosition + 1), (colPosition + 1)]

    Aij <- ((sqrt(dy^2 + (cell - cellRight)^2) +
    sqrt(dy^2 + (cellDiag - cellLower)^2)) *
    (sqrt(dx^2 + (cell - cellLower)^2) +
     sqrt(dx^2 + (cellRight - cellDiag)^2)))/4
  } else {
    Aij <- NA
  }
  # res <- cbind(rowPosition, colPosition, Aij)
  return(Aij)
}

matAij <- sapply(1:ncell(mrast), calcAij, simplify = TRUE)
M <- dim(x)[1]
N <- dim(x)[2]

Sdr <- ((sum(matAij, na.rm=T) - (M-1)*(N-1)*dx*dy)/((M-1)*(N-1)*dx*dy))*100
