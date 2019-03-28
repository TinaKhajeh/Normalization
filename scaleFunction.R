setwd("/Users/Tina/Desktop/Thesis Code")
# first scale data in 0 to 30:
#function scale between 0 yo value
 
scaleFunction <- function(matrix, value){
  indexNa = is.na(matrix)
  matrix[indexNa] = 8#put the temparory value for genes that have NA values(This will not affect data and will be replace at the end
  minValue = min(matrix)
  maxValue = max(matrix)
  scaledMatrix = ((matrix-minValue)/(maxValue-minValue))*value
  scaledMatrix[indexNa] = NA
  return (scaledMatrix)
}
