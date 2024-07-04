burstiness_KJ <- function(x){
  out <- burstiness_GB(x)
  B <- out[1]
  itv <- which(x > 0,arr.ind=TRUE)[,1]
  itv_dist <- diff(itv)
  r <- sd(itv_dist) / mean(itv_dist)
  if (r > sqrt(sum(x)-1)){
    print("#########################################")
    print("##OUT OF BOUNDARY FROM KIM & JO'S PAPER##")
    print("#########################################")
  } else {
    n <- sum(x)
    A <- (sqrt(n+1)-sqrt(n-1)+ B*(sqrt(n+1)+sqrt(n-1)))/
      (sqrt(n+1)+sqrt(n-1)-2+B*(sqrt(n+1)-sqrt(n-1)-2))
  }
  return(A)
}