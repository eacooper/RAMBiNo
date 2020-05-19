%{

TITLE: 
    RAMBiNo

DESCRIPTION: 
    This toolbox provides 12 functions for marginalizing a bivariate normal
    distribution of increasingly more general form. Functions ending with
    "_pt" marginalize over radius and return the 1-D probability (pt) as a
    function of angle (t) in radians. Functions ending with "_pr"
    marginalize over angle and return the 1-D probability (pr) as a
    function of radius (r). 

    The input parameters to this function are the bivariate normal mean
    (mx, my), standard deviation (s for isotropic, sx, sy, for
    anisotropic), and covariance (p). An optional parameter angle (t, for *_pt)
    and radius (r, for *_pr) may be passed for which the 1-D probability
    will be evaluated.

        zero-mean, isototropic: 
            [t,pt] = zeromean_isotropic_pt(s)
            [r,pr] = zeromean_isotropic_pr(s)
        zero-mean, anisototropic, diagonal
            [t,pt] = zeromean_anisotropic_diagonal_pt(sx, sy) 
            [r,pr] = zeromean_anisotropic_diagonal_pr(sx, sy)
        zero-mean, anisototropic, non-diagonal
            [t,pt] = zeromean_anisotropic_nondiagonal_pt(sx, sy, p) 
            [r,pr] = zeromean_anisotropic_nondiagonal_pr(sx, sy, p)
        non-zero-mean, isototropic: 
            [t,pt] = nonzeromean_isotropic_pt(mx, my, s)
            [r,pr] = nonzeromean_isotropic_pr(mx, my, s)
        non-zero-mean, anisototropic, diagonal
            [t,pt] = nonzeromean_anisotropic_diagonal_pt(mx, my, sx, sy) 
            [r,pr] = nonzeromean_anisotropic_diagonal_pr(smx, my, x, sy)
        non-zero-mean, anisototropic, non-diagonal
            [t,pt] = nonzeromean_anisotropic_nondiagonal_pt(mx, my, sx, sy, p) 
            [r,pr] = nonzeromean_anisotropic_nondiagonal_pr(mx, my, sx, sy, p)

NOTE:
    Because the functions nonzeromean_anisotropic_diagonal_pr(), and
    nonzeromean_anisotropic_nondiagonal_pr(), rely on Bessel functions,
    they can become unstable for large values. A warning is printed if
    the marginalizations are thought to be numericaly unstable.

AUTHORS:
    Emily A. Cooper and Hany Farid

CITATION:
    E.A. Cooper and H. Farid. A Toolbox for the Radial and Angular Marginalization of Bivariate
    Normal Distribution, <<<TBD>>>

DATE:
    April 20, 2020
%}

% --------------------------------------------
function[] = RAMBiNoDemo()

RAMBiNo( 0, 0, 2, 2, 0, 1 ); % zero-mean, isotropic
RAMBiNo( 0, 0, 3, 2, 0, 2 ); % zero-mean, anisotropic, diagonal
RAMBiNo( 0, 0, 3, 2, 0.75, 3 ); % zero-mean, anisotropic, non-diagonal

RAMBiNo( 1.5, -1.5, 2, 2, 0, 4 ); % non-zero-mean, isotropic
RAMBiNo( 1.5, -1.5, 3, 2, 0, 5 ); %  non-zero-mean, anisotropic, diagonal
RAMBiNo( 1.5, -1.5, 3, 2, 0.75, 6 ); % non-zero-mean, anisotropic, non-diagonal

% --------------------------------------------
function[] = RAMBiNo( mx, my, sx, sy, p, num )

% generate 2-D bivariate normal
[x,y] = meshgrid( -15:0.02:15, -15:0.02:15 );
N     = generateN( x, y, mx, my, sx, sy, p );

% analytic marginalization
if(num == 1 )
    s = sx; % sx = sy
    [t,pt] = zeromean_isotropic_pt( s );
    [r,pr] = zeromean_isotropic_pr( s );
elseif( num == 2 )
    [t,pt] = zeromean_anisotropic_diagonal_pt( sx, sy );
    [r,pr] = zeromean_anisotropic_diagonal_pr( sx, sy );
elseif( num == 3 )
    [t,pt] = zeromean_anisotropic_nondiagonal_pt( sx, sy, p );
    [r,pr] = zeromean_anisotropic_nondiagonal_pr( sx, sy, p );
elseif( num == 4 )
    s = sx; % sx = sy
    [t,pt] = nonzeromean_isotropic_pt( mx, my, s );
    [r,pr] = nonzeromean_isotropic_pr( mx, my, s );
elseif( num == 5 )
    [t,pt] = nonzeromean_anisotropic_diagonal_pt( mx, my, sx, sy );
    [r,pr] = nonzeromean_anisotropic_diagonal_pr( mx, my, sx, sy );
elseif( num == 6 )
    [t,pt] = nonzeromean_anisotropic_nondiagonal_pt( mx, my, sx, sy, p );
    [r,pr] = nonzeromean_anisotropic_nondiagonal_pr( mx, my, sx, sy, p );
else
    return;
end

% numeric marginalization
[t_numeric,pt_numeric] = numeric_pt( mx, my, sx, sy, p, t, r ); % marginalize over radius
[r_numeric,pr_numeric] = numeric_pr( mx, my, sx, sy, p, t, r ); % marginalize over angle

% display 2-D and 1-D analytic and numeric distributions
sp = (num-1)*3 + 1;
figure(1); pos = get(gcf,'Position'); set( gcf, 'Position', [pos(1) pos(2) 650 800] );
subplot(6,3,sp); imagesc(N); axis image off; colorbar; colormap gray; set( gca, 'Fontsize', 12 );
title( 'g(x,y)' );

subplot(6,3,sp+1); h = plot( t, pt, '-', t_numeric, pt_numeric, '--' ); axis square; set( gca, 'Fontsize', 12 );
set( h(1), 'LineWidth', 3 ); set( h(2), 'LineWidth', 3 );
title( 'p(\theta)' );
axis( [min(t) max(t) 0 1.2*max(pt_numeric)] );

subplot(6,3,sp+2); h = plot( r, pr, '-', r_numeric, pr_numeric, '--' ); axis square; set( gca, 'Fontsize', 12 );
set( h(1), 'LineWidth', 3 ); set( h(2), 'LineWidth', 3 );
title( 'p(r)' );
axis( [min(r) max(r) 0 1.2*max(pr_numeric)] );
drawnow;