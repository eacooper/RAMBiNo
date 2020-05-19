% zero-mean, anisotropic, diagonal : radial marginalization
function[t,pt] = zeromean_anisotropic_diagonal_pt( sx, sy, t )
if( ~exist( 't', 'var' ) )
    t = [-180 : 359/1000 : 179] * pi/180; % angle (radians)
end
a  = 2*pi*sx*sy;
b  = (cos(t).^2)/(2*sx^2) + (sin(t).^2)/(2*sy^2);
pt = 1./(2*a*b);
