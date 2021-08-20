#.........................................................
# Carl Schmertmann
#   created 6 Aug 2021
#   edited  6 Aug 2021
#
# standalone function to produce quadratic spline fits
# for an age range, given the 4 parameter values (R,alpha,P,H)
#
# This version of the function returns a list with $x and $y
# components, suitable for plotting directly
#.........................................................

QS = function(min_age = 12,
              max_age = 55,
              delta   = .25,
              R       = .20,
              alpha   = 15,
              P       = 28,
              H       = 34) {
  
  age     = seq(from = min_age+delta/2, 
                to   = max_age-delta/2, 
                by   = delta)
  
  w     = min(.75, .25+.025*(P-alpha))
  beta  = max (  min(50, 4*H-3*P) , (4*H-P)/3 )
  
  D     = P-20          # delay index
  C     = (P+50)/2 - H  # control index
  
  knot    = vector("numeric",5)
  knot[1] = alpha
  knot[2] = alpha + w*(P-alpha)
  knot[3] = P
  knot[4] = (P+H)/2
  knot[5] = (H+beta)/2
  
  A      = matrix(NA,5,5)
  target = vector("numeric",5)
  
  # target value at P=1
  A[1,]     = (pmax(P-knot, 0))^2
  target[1] = 1
  
  # target value at H=1/2
  A[2,]     = (pmax(H-knot, 0))^2
  target[2] = 1/2
  
  # target value at beta=0
  A[3,]     = (pmax(beta-knot, 0))^2
  target[3] = 0
  
  # target slope at P=0
  A[4,]     = 2*pmax(P-knot, 0) 
  target[4] = 0
  
  # target slope at beta=0
  A[5,]     = 2*pmax(beta-knot, 0) 
  target[5] = 0
  
  # calculate thetas
  theta = solve(A,target)
  
  tmp = (pmax( outer(age,knot,"-") , 0))^2
  
  y = (age >= alpha) * (age <= beta) * R * (tmp %*% theta)  # ASFRs
  
  return(list(x = age,
              y = as.vector(y) ))
  
} #QS

