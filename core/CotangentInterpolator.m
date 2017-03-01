% Interpolate in (x,nu) space
%
% First, create the interpolator:
%     I = CotangentInterpolator(v, xGrid, wcr)
% Then you can call it:
%    Iv = I(x, xi)
%
%     v:  I1 x ... x In x M array of (x,nu) values   [M = #boxes]
% xGrid:  Grid object for the x grid.
%   wcr:  WedgeCreator for this set of boxes.
%     I:  interpolator function

% x, xi:  n x J1 x ... x Jr array of x and xi values
%    Iv:  J1 x ... x Jr array of interpolated values
function I = CotangentInterpolator(v, xGrid, wcr)

switch wcr.grid.dims
    case 2
        % Retrieve angles for the boxes.
        %
        % ASSUMPTION: the boxes are ordered by angle from 0 to 2*pi.
        nus = wcr.Centers();
        nuAngles = mod(atan2(nus(2,:), nus(1,:)), 2*pi);
        
        % Handle wraparound by duplicating data at the beginning
        % and end.
        nuAngles = [nuAngles(end)-2*pi nuAngles nuAngles(1:2)+2*pi];
        %v = cat(3, v(:,:,end), v, v(:,:,1:2));
        v = v(:,:,[end 1:end 1:2]);
        
        Ig = griddedInterpolant([xGrid.gridVectors nuAngles], v, 'linear');
        I = @(x,xi) CotangentInterpolator2D(Ig,x,xi);
    otherwise
        error('Only 2D supported currently.');
end

end

function Iv = CotangentInterpolator2D(Ig,x,xi)
th = mod(atan2(xi(2,:), xi(1,:)), 2*pi);
Iv = Ig(x(1,:), x(2,:), th);
sz = size(x);
Iv = reshape(Iv, [sz(2:end) 1]);
end