% generate a bivariate normal distribution
function[N] = generateN( x, y, mx, my, sx, sy, p )
N = 1/(2*pi*sx*sy*sqrt(1-p^2)) * exp( -(1/(2*(1-p^2))*((x-mx).^2/(sx^2) + (y-my).^2/(sy^2) - 2*p*(x-mx).*(y-my)/(sx*sy))) );
