% Approximate dy/dxi(x,nu).
%
%   dxiy = ApproximateDXiY(y_cells, b, yGrid, wcr)
%
% Inputs:
%       y_cells: y(x,nu) values, stored as a cell array
%                  where y_cells{b} stores the values when
%                  nu is at the center of box b.
%             b: box number to approximate dy/dxi at
%         yGrid: y grid
%           wcr: WedgeCreator used for the boxes.
%
% Output:
%          dxiy: n x n x i1 x ... x in array.
%                 dxiy(i,j,:,...,:) stores the approximate values of
%                 dy_i/dxi_j at (x,nu) for each x in the domain grid
%                 for the nu of the given box b.
%
function dxiy = ApproximateDXiY(y_cells, b, yGrid, wcr)

% Get the adjoining nu values (nu_p, nu_m) for b + 1, b - 1
% (the weird mods handle wraparound when b = 1 or b = #boxes.)
boxes = wcr.Count();
b_p = mod(b,boxes) + 1;
b_m = mod(b-2,boxes) + 1;
nu_p = wcr.Center(b_p);
nu_m = wcr.Center(b_m);

n = length(nu_p);
assert(n == 2, 'Approximating d\xi/dy is only implemented for n = 2.');

% Calculate the difference in angles between nu_p and nu_m,
%  including the sign.
dth = acos(dot(nu_p,nu_m));
boxDir = (-nu_m(2)*nu_p(1)+nu_m(1)*nu_p(2)) > 0;
dth = dth * boxDir;

% Calculate dy/dtheta numerically on the grid edges.
% Use PeriodicDiff in case the y grid is periodic.
dthy = PeriodicDiff(y_cells{b_p}, y_cells{b_m}, yGrid)*(1/dth);
% Interpolate dy/dtheta to gridpoints.
dthy = InterpolateCornersToCenters(dthy, n, 1);
    
% Calculate dy/dxi as dy/dtheta * dtheta/dxi,
%  In general, dtheta/dxi(x,xi) = xi^\perp / |xi|^2,
%  but here we assume xi = nu, so dtheta/dxi = nu^\perp.
nu = wcr.Center(b);
dxiy(:,1,:,:) = -nu(2)*dthy;
dxiy(:,2,:,:) = +nu(1)*dthy;

end