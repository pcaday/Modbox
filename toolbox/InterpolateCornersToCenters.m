% Linearly interpolate a function defined on the corners of a grid
%  to a function on the center of the grid.
%
% For instance, function values at the x's in this 2D example
%  would be interpolated to values at the o's:
%
%       x     x     x     x
%          o     o     o
%       x     x     x     x
%          o     o     o
%       x     x     x     x
%
%
%   fc = InterpolateCornersToCenters(fe, n)
%   fc = InterpolateCornersToCenters(fe, n, odims)
%
%  odims: number of dimensions in the codomain of the function
%           (# of dimensions of the values).
%         Default 0 (scalar function).
%
%     fe: values for the corners
%          (J1 x ... x Jk x I1 x ... x In array, where k = odims)
%     fc: values for the centers
%          (J1 x ... x Jk x (I1-1) x ... x (In-1) array.
%      n: number of dimensions in the domain
%
function fc = InterpolateCornersToCenters(fe, n, odims)
    % odims defaults to 0 (scalar function)
    if nargin < 3, odims = 0; end
    
    % Resize the input from
    %   J1 x ... x Jk x I1 x ... x In
    % to
    %   M x I1 x ... x In
    % where M is the product J1*J2*...*Jk
    sz = size(fe);
    assert(length(sz) <= odims+n, 'Too many dimensions in input array.');
    
    sz = [sz ones(1,n)];
    nsz = [prod(sz(1:odims)) sz(odims+1:odims+n) 1];
    fe = reshape(fe, nsz);
    
    switch n
        case 0
            fc = fe;
        case 1      % untested
            fc = 0.5 * (fe(:,1:end-1) + fe(:,2:end));
        case 2
            fc = 0.25 * (...
                fe(:,1:end-1,1:end-1) + ...
                fe(:,2:end,  1:end-1) + ...
                fe(:,1:end-1,2:end  ) + ...
                fe(:,2:end,  2:end));
        case 3      % untested.
            fc = 0.125 * (...
                fe(:,1:end-1,1:end-1,1:end-1) + ...
                fe(:,2:end,  1:end-1,1:end-1) + ...
                fe(:,1:end-1,2:end,  1:end-1) + ...
                fe(:,2:end,  2:end,  1:end-1) + ...
                fe(:,1:end-1,1:end-1,2:end)   + ...
                fe(:,2:end,  1:end-1,2:end)   + ...
                fe(:,1:end-1,2:end,  2:end)   + ...
                fe(:,2:end,  2:end,  2:end));
    end
    
    % Reshape fc back to the correct size
    sz = [sz(1:odims) max(0,sz(odims+1:odims+n) - 1)];
    fc = reshape(fc, sz);
end