% nonzero-mean, isotropic : angular marginalization
function[r,pr] = nonzeromean_isotropic_pr( mx, my, s, r )
if( ~exist( 'r', 'var' ) )
    r = [0 : 5*s/1000 : 5*s]; % radius: default is 0 to 5 x [std.dev.(s)]
end
a  = sqrt(mx^2 + my^2);
pr = r/s^2 .* exp(-(r.^2 + a^2)/(2*s^2)) .* besseli(0,(r*a)/s^2);
