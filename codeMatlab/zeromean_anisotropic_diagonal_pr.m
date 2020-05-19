% zero-mean, anisotropic, diagonal : angular marginalization
function[r,pr] = zeromean_anisotropic_diagonal_pr( sx, sy, r )
if( ~exist( 'r', 'var' ) )
    s = max(sx,sy);
    r = [0 : 5*s/1000 : 5*s]; % radius: default is 0 to 5 x [std.dev.(s)]
end
a  = (sy^2 + sx^2)/(2*sx*sy)^2;
b  = (sx^2 - sy^2)/(2*sx*sy)^2;
pr = r/(sx*sy) .* exp(-a*r.^2) .* besseli(0,-b*r.^2);
