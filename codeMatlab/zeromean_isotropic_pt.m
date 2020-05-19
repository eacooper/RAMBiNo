% zero-mean, isotropic : radial marginalization
function[t,pt] = zeromean_isotropic_pt( s, t )
if( ~exist( 't', 'var' ) )
    t = [-180 : 359/1000 : 179] * pi/180; % angle (radians)
end
pt = 1/(2*pi)*ones(size(t));