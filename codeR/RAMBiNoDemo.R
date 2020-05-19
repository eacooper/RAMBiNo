#   TITLE: 
#     RAMBiNo
#   
#   DESCRIPTION: 
#   This toolbox provides 12 functions for marginalizing a bivariate normal
#   distribution of increasingly more general form. Functions ending with
#   "_pt" marginalize over radius and return the 1-D probability (pt) as a
#   function of angle (t) in radians. Functions ending with "_pr"
#   marginalize over angle and return the 1-D probability (pr) as a
#   function of radius (r). 
#   
#   The main input parameters to these function are the bivariate normal:
#   mean (mx, my) 
#   standard deviation (s for isotropic, sx, sy, for anisotropic)
#   covariance (p) 
#
#   Input parameters for angle (t, for *_pt)
#   and radius (r, for *_pr) must also be passed for which the 1-D probability
#   will be evaluated.
#
#   Each function returns these parameters back in a dataframe called "result", 
#   along with the probability evaluated at each point (pt or pr)
#   
#   zero-mean, isototropic: 
#     zeromean_isotropic_pt(s)
#     zeromean_isotropic_pr(s)
#   zero-mean, anisototropic, diagonal
#     zeromean_anisotropic_diagonal_pt(sx, sy) 
#     zeromean_anisotropic_diagonal_pr(sx, sy)
#   zero-mean, anisototropic, non-diagonal
#     zeromean_anisotropic_nondiagonal_pt(sx, sy, p) 
#     zeromean_anisotropic_nondiagonal_pr(sx, sy, p)
#   non-zero-mean, isototropic: 
#      nonzeromean_isotropic_pt(mx, my, s)
#      nonzeromean_isotropic_pr(mx, my, s)
#   non-zero-mean, anisototropic, diagonal
#     nonzeromean_anisotropic_diagonal_pt(mx, my, sx, sy) 
#     nonzeromean_anisotropic_diagonal_pr(smx, my, x, sy)
#   non-zero-mean, anisototropic, non-diagonal
#     nonzeromean_anisotropic_nondiagonal_pt(mx, my, sx, sy, p) 
#     nonzeromean_anisotropic_nondiagonal_pr(mx, my, sx, sy, p)
#   
#   NOTE:
#     Because the functions nonzeromean_anisotropic_diagonal_pr(), and
#   nonzeromean_anisotropic_nondiagonal_pr(), rely on Bessel functions,
#   they can become unstable for large values. A warning is printed if
#   the marginalizations are thought to be numericaly unstable.
#   
#   AUTHORS:
#     Emily A. Cooper and Hany Farid
#   
#   CITATION:
#     E.A. Cooper and H. Farid. A Toolbox for the Radial and Angular Marginalization of Bivariate
#   Normal Distribution, <<<TBD>>>
#     
#     DATE:
#     April 20, 2020


# load functions
source("analytic_functions_pr.R")
source("analytic_functions_pt.R")


####################################################################################
# zero-mean, isotropic
s = 2

# angular marginalization
r      = seq(from = 0, to = 5*s, by = 5*s/1000)
result = zeromean_isotropic_pr( s, r )

plot(result$r,result$pr,ylim=c(0, 1),type='l')
title('zero-mean, isotropic p(r)')

# radial marginalization
t = seq(from = -180, to = 179, by = 359/1000) * pi/180
result = zeromean_isotropic_pt( s, t )

plot(result$t,result$pt,ylim=c(0, 1),type='l')
title('zero-mean, isotropic p(t)')

####################################################################################
# zero-mean, anisotropic, diagonal
sx = 3
sy = 2

# angular marginalization
r      = seq(from = 0, to = 5*(max(sx,sy)), by = 5*(max(sx,sy))/1000)
result = zeromean_anisotropic_diagonal_pr(sx, sy, r)

plot(result$r,result$pr,ylim=c(0, 1),type='l')
title('zero-mean, anisotropic, diagonal p(r)')

# radial marginalization
t = seq(from = -180, to = 179, by = 359/1000) * pi/180
result = zeromean_anisotropic_diagonal_pt( sx, sy, t )

plot(result$t,result$pt,ylim=c(0, 1),type='l')
title('zero-mean, anisotropic, diagonal p(t)')

####################################################################################
# zero-mean, anisotropic, nondiagonal
sx = 3
sy = 2
p =  0.75

# angular marginalization
r      = seq(from = 0, to = 5*(max(sx,sy)), by = 5*(max(sx,sy))/1000)
result = zeromean_anisotropic_nondiagonal_pr(sx, sy, p, r)

plot(result$r,result$pr,ylim=c(0, 1),type='l')
title('zero-mean, anisotropic, nondiagonal p(r)')

# radial marginalization
t = seq(from = -180, to = 179, by = 359/1000) * pi/180
result = zeromean_anisotropic_nondiagonal_pt(sx, sy, p, t)

plot(result$t,result$pt,ylim=c(0, 1),type='l')
title('zero-mean, anisotropic, nondiagonal p(t)')

####################################################################################
# non-zero-mean, isotropic
mx = 1.5
my = -1.5
s = 2

# angular marginalization
r      = seq(from = 0, to = 5*s, by = 5*s/1000)
result = nonzeromean_isotropic_pr(mx, my, s, r)

plot(result$r,result$pr,ylim=c(0, 1),type='l')
title('non-zero-mean, isotropic p(r)')

# radial marginalization
t = seq(from = -180, to = 179, by = 359/1000) * pi/180
result = nonzeromean_isotropic_pt(mx, my, s, t)

plot(result$t,result$pt,ylim=c(0, 1),type='l')
title('non-zero-mean, isotropic p(t)')

####################################################################################
# non-zero-mean, anisotropic, diagonal
mx = 1.5
my = -1.5
sx = 3
sy = 2


# angular marginalization
r      = seq(from = 0, to = 5*(max(sx,sy)), by = 5*(max(sx,sy))/1000)
result = nonzeromean_anisotropic_diagonal_pr( mx, my, sx, sy, r )

plot(result$r,result$pr,ylim=c(0, 1),type='l')
title('non-zero-mean, anisotropic, diagonal p(r)')

# radial marginalization
t = seq(from = -180, to = 179, by = 359/1000) * pi/180
result = nonzeromean_anisotropic_diagonal_pt( mx, my, sx, sy, t )

plot(result$t,result$pt,ylim=c(0, 1),type='l')
title('non-zero-mean, anisotropic, diagonal p(t)')


####################################################################################
# non-zero-mean, anisotropic, nondiagonal
mx = 1.5
my = -1.5
sx = 3
sy = 2
p = 0.75


# angular marginalization
r      = seq(from = 0, to = 5*(max(sx,sy)), by = 5*(max(sx,sy))/1000)
result = nonzeromean_anisotropic_nondiagonal_pr( mx, my, sx, sy, p, r )

plot(result$r,result$pr,ylim=c(0, 1),type='l')
title('non-zero-mean, anisotropic, nondiagonal p(r)')

# radial marginalization
t = seq(from = -180, to = 179, by = 359/1000) * pi/180
result = nonzeromean_anisotropic_nondiagonal_pt( mx, my, sx, sy, p, t )

plot(result$t,result$pt,ylim=c(0, 1),type='l')
title('non-zero-mean, anisotropic, nondiagonal p(t)')