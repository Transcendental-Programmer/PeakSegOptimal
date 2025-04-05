### Compute the weighted Normal (squared error) loss function
NormalLoss <- structure(function(observation, seg.mean, weight=1){
  stopifnot(is.numeric(observation))
  stopifnot(is.numeric(seg.mean))
  stopifnot(is.numeric(weight))
  n.data <- length(observation)
  
  if(length(seg.mean) == 1){
    seg.mean <- rep(seg.mean, n.data)
  }
  if(length(weight) == 1){
    weight <- rep(weight, n.data)
  }
  
  stopifnot(n.data == length(seg.mean))
  stopifnot(n.data == length(weight))
  if(any(weight < 0)){
    stop("NormalLoss undefined for negative weight")
  }
  
  # Calculate weighted squared error loss
  loss <- weight * (observation - seg.mean)^2
  sum(loss)
}, ex=function(){
  NormalLoss(1, 1)
  NormalLoss(0, 1)
  NormalLoss(c(1, 2, 3), c(1, 2, 3))
  NormalLoss(c(1, 2, 3), 2)
})
