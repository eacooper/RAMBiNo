% non-zero-mean, anisotropic, nondiagonal : angular marginalization
function[r,pr] = nonzeromean_anisotropic_nondiagonal_pr( mx, my, sx, sy, p, r )
if( ~exist( 'r', 'var' ) )
    s = max(sx,sy);
    r = [0 : 5*s/1000 : 5*s]; % radius: default is 0 to 5 x [std.dev.(s)]
end
w   = 1/2*atan2((2*p*sx*sy), (sx^2-sy^2));
mxt = mx*cos(w) + my*sin(w);
myt = -mx*sin(w) + my*cos(w);
sxt = sqrt( (sx^2+sy^2)/2 + sqrt((sx^2+sy^2)^2/4 - (sx^2*sy^2 - p^2*sx^2*sy^2)) );
syt = sqrt( (sx^2+sy^2)/2 - sqrt((sx^2+sy^2)^2/4 - (sx^2*sy^2 - p^2*sx^2*sy^2)) );
a   = 1/(sxt*syt) * exp(-(mxt^2*syt^2 + myt^2*sxt^2)/(2*sxt^2*syt^2));
b   = (sxt^2-syt^2)/(4*sxt^2*syt^2);
c   = sqrt((mxt/sxt^2)^2 + (myt/syt^2)^2);
psi = atan2( (myt*sxt^2), (mxt*syt^2) );
d   = zeros( size(r) );
for k = 1 : 100 % truncated series
    d = d + (besseli(k,b*r.^2) .* besseli(2*k,c*r) * cos(2*k*psi));
end
pr = a*r .* exp(-(r.^2*(sxt^2+syt^2))/(4*sxt^2*syt^2)) .* (besseli(0,b*r.^2) .* besseli(0,c*r) + 2*d);
if( max(diff(pr)) > 0.1 )
    fprintf( 'WARNING: possible numeric instability in nonzeromean_anisotropic_nondiagonal_pr\n' );
end

