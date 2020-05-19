zeromean_isotropic_pr <- function( s, r ) {
  # zero-mean, isotropic : angular marginalization
  
  pr = r/s^2 * exp(-r^2/(2*s^2))
  
  result = data.frame(r,pr)
  return(result)
}


####################################################################################
zeromean_anisotropic_diagonal_pr<- function( sx, sy, r ) {
  # zero-mean, anisotropic, diagonal : angular marginalization
  
  a  = (sy^2 + sx^2)/(2*sx*sy)^2
  b  = (sx^2 - sy^2)/(2*sx*sy)^2
  pr = r/(sx*sy) * exp(-a*r^2) * besselI(abs(-b*r^2),0)
  
  result = data.frame(r,pr)
  return(result)
}


####################################################################################
zeromean_anisotropic_nondiagonal_pr <- function( sx, sy, p, r ) {
  # zero-mean, anisotropic, nondiagonal : angular marginalization
  
  sxt = sqrt( (sx^2+sy^2)/2 + sqrt((sx^2+sy^2)^2/4 - (sx^2*sy^2 - p^2*sx^2*sy^2)) )
  syt = sqrt( (sx^2+sy^2)/2 - sqrt((sx^2+sy^2)^2/4 - (sx^2*sy^2 - p^2*sx^2*sy^2)) )
  a   = (sxt^2 + syt^2)/(2*sxt*syt)^2
  b   = (sxt^2 - syt^2)/(2*sxt*syt)^2
  pr  = r/(sxt*syt) * exp(-a*r^2) * besselI(abs(-b*r^2),0)
  
  result = data.frame(r,pr)
  return(result)
}


####################################################################################
nonzeromean_isotropic_pr <- function( mx, my, s, r ) {
  # non-zero-mean, isotropic : angular marginalization
  
  a  = sqrt(mx^2 + my^2)
  pr = r/s^2 * exp(-(r^2 + a^2)/(2*s^2)) * besselI(abs((r*a)/s^2),0)
  
  result = data.frame(r,pr)
  return(result)
}


####################################################################################
nonzeromean_anisotropic_diagonal_pr <- function( mx, my, sx, sy, r ) {
  # non-zero-mean, anisotropic, diagonal : angular marginalization
  
  a   = 1/(sx*sy) * exp(-(mx^2*sy^2 + my^2*sx^2)/(2*sx^2*sy^2))
  b   = (sx^2-sy^2)/(4*sx^2*sy^2)
  c   = sqrt((mx/sx^2)^2 + (my/sy^2)^2)
  psi = atan2( (my*sx^2), (mx*sy^2) )
  
  d   = rep(0, length(r))
  for (k in 1 : 100) { # truncated series
    d = d + (besselI(abs(b*r^2),k) * (sign(b*r^2)^k)) * (besselI(abs(c*r),2*k) * (sign(c*r)^(2*k))) * cos(2*k*psi)
  }
  
  pr  = a*r * exp(-(r^2*(sx^2+sy^2))/(4*sx^2*sy^2)) * (besselI(abs(b*r^2),0) * besselI(abs(c*r),0) + 2*d)
  
  if( max(diff(pr)) > 0.1 ) {
    print( 'WARNING: possible numeric instability in nonzeromean_anisotropic_diagonal_pr' )
  }
  
  result = data.frame(r,pr)
  return(result)
  
}


####################################################################################
nonzeromean_anisotropic_nondiagonal_pr <- function(mx, my, sx, sy, p, r) {
  # non-zero-mean, anisotropic, nondiagonal : angular marginalization
  
  w   = 1/2*atan2((2*p*sx*sy), (sx^2-sy^2))
  mxt = mx*cos(w) + my*sin(w)
  myt = -mx*sin(w) + my*cos(w)
  sxt = sqrt( (sx^2+sy^2)/2 + sqrt((sx^2+sy^2)^2/4 - (sx^2*sy^2 - p^2*sx^2*sy^2)) )
  syt = sqrt( (sx^2+sy^2)/2 - sqrt((sx^2+sy^2)^2/4 - (sx^2*sy^2 - p^2*sx^2*sy^2)) )
  
  a   = 1/(sxt*syt) * exp(-(mxt^2*syt^2 + myt^2*sxt^2)/(2*sxt^2*syt^2))
  b   = (sxt^2-syt^2)/(4*sxt^2*syt^2)
  c   = sqrt((mxt/sxt^2)^2 + (myt/syt^2)^2)
  psi = atan2( (myt*sxt^2), (mxt*syt^2) )
  
  d   = rep(0, length(r))
  for (k in 1 : 100) { # truncated series
    d = d + (besselI(abs(b*r^2),k) * (sign(b*r^2)^k)) * (besselI(abs(c*r),2*k) * (sign(c*r)^(2*k))) * cos(2*k*psi)
  }
  
  pr  = a*r * exp(-(r^2*(sxt^2+syt^2))/(4*sxt^2*syt^2)) * (besselI(abs(b*r^2),0) * besselI(abs(c*r),0) + 2*d)
  
  if( max(diff(pr)) > 0.1 ) {
    print( 'WARNING: possible numeric instability in nonzeromean_anisotropic_diagonal_pr' )
  }
  
  result = data.frame(r,pr)
  return(result)
  
}
