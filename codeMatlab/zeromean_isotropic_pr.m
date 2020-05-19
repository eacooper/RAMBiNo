% zero-mean, isotropic : angular marginalization
function[r,pr] = zeromean_isotropic_pr( s, r )
if( ~exist( 'r', 'var' ) )
    r = [0 : 5*s/1000 : 5*s]; % radius: default is 0 to 5 x [std.dev.(s)]
end
pr = r/s^2 .* exp(-r.^2/(2*s^2));