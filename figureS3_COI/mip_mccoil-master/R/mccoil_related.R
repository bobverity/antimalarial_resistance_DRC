
assignGTforREALMcCOIL <- function(wsaf, err = 0.1){

  GT <- matrix(NA, dim(wsaf)[1], dim(wsaf)[2])

  GT <- ifelse(wsaf > 1-err, 1,
               ifelse(wsaf < 0+err, 0,
                      0.5))
  GT[is.na(GT)] <- -1
  return(GT)

}
