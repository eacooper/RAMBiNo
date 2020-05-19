zeromean_isotropic_pt <- function( s, t ) {
  # zero-mean, isotropic : radial marginalization
  
  pt = 1/(2*pi)*rep(1,length(t))
  
  result = data.frame(t,pt)
  return(result)
}


####################################################################################
zeromean_anisotropic_diagonal_pt<- function( sx, sy, t ) {
  # zero-mean, anisotropic, diagonal : radial marginalization
  
  a  = 2*pi*sx*sy
  b  = (cos(t)^2)/(2*sx^2) + (sin(t)^2)/(2*sy^2)
  pt = 1./(2*a*b)
  
  result = data.frame(t,pt)
  return(result)
}


####################################################################################
zeromean_anisotropic_nondiagonal_pt <- function( sx, sy, p, t ) {
  # zero-mean, anisotropic, nondiagonal : radial marginalization
  
  a  = 2*pi*sx*sy*sqrt(1-p^2)
  b  = 1/(2*(1-p^2)) * (cos(t)^2/(sx^2) + sin(t)^2/(sy^2) - 2*p*sin(t) * cos(t)/(sx*sy))
  pt = 1./(2*a*b)
  
  result = data.frame(t,pt)
  return(result)
}


####################################################################################
nonzeromean_isotropic_pt <- function( mx, my, s, t ) {
  # non-zero-mean, isotropic : radial marginalization
  
  a  = 1/s * sqrt(mx^2 + my^2)
  b  = 1/s * (mx*cos(t) + my*sin(t))
  pt = 1/(sqrt(2*pi)) * dnorm(a) * (1 + b*pnorm(b)/dnorm(b))
  
  result = data.frame(t,pt)
  return(result)
}


####################################################################################
nonzeromean_anisotropic_diagonal_pt <- function( mx, my, sx, sy, t ) {
  # non-zero-mean, anisotropic, diagonal : radial marginalization
  
  a  = cos(t)^2/sx^2 + sin(t)^2/sy^2
  b  = sqrt(mx^2/sx^2 + my^2/sy^2)
  c  = (cos(t)*mx/sx^2 + sin(t)*my/sy^2) / sqrt(a)
  
  pt = 1/(a*sqrt(2*pi*sx^2*sy^2)) * dnorm(b) * (1 + c*pnorm(c)/dnorm(c))
  
  result = data.frame(t,pt)
  return(result)
  
}


####################################################################################
nonzeromean_anisotropic_nondiagonal_pt <- function(mx, my, sx, sy, p, t) {
  # non-zero-mean, anisotropic, nondiagonal : radial marginalization
  
  c  = 1 / (sx*sy*sqrt(1-p^2))
  a  = c^2 * (sy^2*cos(t)^2 - p*sx*sy*sin(2*t) + sx^2*sin(t)^2)
  #b  = generateN(mx,my,0,0,sx,sy,p)
  b = 1/(2*pi*sx*sy*sqrt(1-p^2)) * exp( -(1/(2*(1-p^2))*((0-mx)^2/(sx^2) + (0-my)^2/(sy^2) - 2*p*(0-mx)*(0-my)/(sx*sy))) )
  d  = (c^2/sqrt(a)) * (mx*sy*(sy*cos(t) - p*sx*sin(t)) + my*sx*(sx*sin(t) - p*sy*cos(t)))
  
  pt = 1/a * (b + c*d*pnorm(d)*dnorm((c*(mx*sin(t) - my*cos(t)))/sqrt(a)) )
  
  result = data.frame(t,pt)
  return(result)
  
}
