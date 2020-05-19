% zero-mean, anisotropic, nondiagonal : angular marginalization
function[r,pr] = zeromean_anisotropic_nondiagonal_pr( sx, sy, p, r )
if( ~exist( 'r', 'var' ) )
    s = max(sx,sy);
    r = [0 : 5*s/1000 : 5*s]; % radius: default is 0 to 5 x [std.dev.(s)]
end
sxt = sqrt( (sx^2+sy^2)/2 + sqrt((sx^2+sy^2)^2/4 - (sx^2*sy^2 - p^2*sx^2*sy^2)) );
syt = sqrt( (sx^2+sy^2)/2 - sqrt((sx^2+sy^2)^2/4 - (sx^2*sy^2 - p^2*sx^2*sy^2)) );
a   = (sxt^2 + syt^2)/(2*sxt*syt)^2;
b   = (sxt^2 - syt^2)/(2*sxt*syt)^2;
pr  = r/(sxt*syt) .* exp(-a*r.^2) .* besseli(0,-b*r.^2);

