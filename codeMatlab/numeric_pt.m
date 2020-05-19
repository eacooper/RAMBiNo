%  marginalization over radius (returns angle T and p(t))
function[T,pt] = numeric_pt( mx, my, sx, sy, p, T, R )
pt = zeros(1,length(T));
C  = 2*pi*R; % circumference
k = 1;
for t = T
    x = R*cos(t);
    y = R*sin(t);
    pt(k) = sum( C .* generateN( x, y, mx, my, sx, sy, p ) );
    k = k + 1;
end
pt = length(pt) * 1/(2*pi) * pt/sum(pt); % mean of distribution should be 1/(2*pi)

