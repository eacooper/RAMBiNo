% nonzero-mean, anisotropic, diagonal : angular marginalization
function[r,pr] = nonzeromean_anisotropic_diagonal_pr( mx, my, sx, sy, r )
if( ~exist( 'r', 'var' ) )
    s = max(sx,sy);
    r = [0 : 5*s/1000 : 5*s]; % radius: default is 0 to 5 x [std.dev.(s)]
end
a   = 1/(sx*sy) * exp(-(mx^2*sy^2 + my^2*sx^2)/(2*sx^2*sy^2));
b   = (sx^2-sy^2)/(4*sx^2*sy^2);
c   = sqrt((mx/sx^2)^2 + (my/sy^2)^2);
psi = atan2( (my*sx^2), (mx*sy^2) );
d   = zeros( size(r) );
for k = 1 : 100 % truncated series
    d = d + (besseli(k,b*r.^2) .* besseli(2*k,c*r) * cos(2*k*psi));
end
pr = a*r .* exp(-(r.^2*(sx^2+sy^2))/(4*sx^2*sy^2)) .* (besseli(0,b*r.^2) .* besseli(0,c*r) + 2*d);
if( max(diff(pr)) > 0.1 )
    fprintf( 'WARNING: possible numeric instability in nonzeromean_anisotropic_diagonal_pr\n' );
end
