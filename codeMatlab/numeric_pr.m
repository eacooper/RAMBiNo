% marginalization over angle (returns radius R and p(r))
function[R,pr] = numeric_pr( mx, my, sx, sy, p, T, R )
pr = zeros(1,length(R));
C  = 2*pi*R; % circumference
k = 1;
for r = R
    x = r*cos(T);
    y = r*sin(T);
    pr(k) = mean( C(k) .* generateN( x, y, mx, my, sx, sy, p ) );
    k = k + 1;
end
