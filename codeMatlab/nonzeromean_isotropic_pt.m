% nonzero-mean, isotropic : radial marginalization
function[t,pt] = nonzeromean_isotropic_pt( mx, my, s, t )
if( ~exist( 't', 'var' ) )
    t = [-180 : 359/1000 : 179] * pi/180; % angle (radians)
end
a  = 1/s * sqrt(mx^2 + my^2);
b  = 1/s * (mx*cos(t) + my*sin(t));
pt = 1./(sqrt(2*pi)) .* normpdf(a) .* (1 + b.*normcdf(b)./normpdf(b));