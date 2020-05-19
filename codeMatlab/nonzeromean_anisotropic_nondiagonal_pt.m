% non-zero-mean, anisotropic, nondiagonal : radial marginalization
function[t,pt] = nonzeromean_anisotropic_nondiagonal_pt( mx, my, sx, sy, p, t )
if( ~exist( 't', 'var' ) )
    t = [-180 : 359/1000 : 179] * pi/180; % angle (radians)
end
c  = 1 / (sx*sy*sqrt(1-p^2));
a  = c^2 * (sy^2*cos(t).^2 - p*sx*sy*sin(2*t) + sx^2*sin(t).^2);
b  = generateN(mx,my,0,0,sx,sy,p);
d  = (c^2./sqrt(a)) .* (mx*sy*(sy*cos(t) - p*sx*sin(t)) + my*sx*(sx*sin(t) - p*sy*cos(t)));
pt = 1./a .* (b + c*d.*normcdf(d).*normpdf((c*(mx*sin(t) - my*cos(t)))./sqrt(a)) );
