#' Random walk on Z^2
#' 
#' Generates a random walk in 2-dimensional integer space
#' 
#' @param n Integer specifying the number of steps of the generated walk.
#' 
#' @return Vector of coordinates for the walk
#' 
#' @examples
#' random_walk(n = 10)
#' 
#' @importFrom stats runif
#' 
#' @export
random_walk <- function(n){
  walk <- vector("list", n+1)
  walk[[1]] <- c(0, 0)
  
  for (i in 2:(n+1)){
    # Step: back or forth
    step <- ifelse(runif(1) > 0.5, 1, -1)
    # Direction: vertical of sideways. It's 50/50.
    direction <- ifelse(runif(1) > 0.5, 1, 2)
    still <- setdiff(c(1, 2), direction)
    
    walk[[i]][direction] <- walk[[i-1]][direction] + step
    walk[[i]][still] <- walk[[i-1]][still]
  }
  return(walk)
}

