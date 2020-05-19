% zero-mean, anisotropic, nondiagonal : radial marginalization
function[t,pt] = zeromean_anisotropic_nondiagonal_pt( sx, sy, p, t )
if( ~exist( 't', 'var' ) )
    t = [-180 : 359/1000 : 179] * pi/180; % angle (radians)
end
a  = 2*pi*sx*sy*sqrt(1-p^2);
b  = 1/(2*(1-p^2)) * (cos(t).^2/(sx^2) + sin(t).^2/(sy^2) - 2*p*sin(t).*cos(t)/(sx*sy));
pt = 1./(2*a*b);
