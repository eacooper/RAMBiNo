% nonzero-mean, anisotropic, diagonal : radial marginalization
function[t,pt] = nonzeromean_anisotropic_diagonal_pt( mx, my, sx, sy, t )
if( ~exist( 't', 'var' ) )
    t = [-180 : 359/1000 : 179] * pi/180; % angle (radians)
end
a  = cos(t).^2/sx^2 + sin(t).^2/sy^2;
b  = sqrt(mx^2/sx^2 + my^2/sy^2);
c  = (cos(t)*mx/sx^2 + sin(t)*my/sy^2) ./ sqrt(a);
pt = 1./(a*sqrt(2*pi*sx^2*sy^2)) .* normpdf(b) .* (1 + c.*normcdf(c)./normpdf(c));