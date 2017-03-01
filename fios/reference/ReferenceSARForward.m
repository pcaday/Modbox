function Rf = ReferenceSARForward(f, stGrid, gamma, c)
% Compute the circular Radon transform Rf of a function f
%  along a given curve.
% The Radon transform is computed directly by approximating line
%  integrals around circles centered at points on the curve.
%
%  Rf = ReferenceSARForward(f, stGrid, gamma, c)
%
% Inputs
%         f: Input function (type Function)
%    stGrid: Grid to compute R_\gamma f on (type Grid)
%     gamma: Curve (function handle: @(s) xVals),
%              or a curveInfo structure.
%            xVals is 2 x i1 x ... x ik, where i1 x ... x ik is the
%              size of the input s.
%         c: Sound speed [scalar] (default 1)
%
% Outputs
%        Rf: Circular Radon transform (type Function)
%

% Substitute default sound speed if none given.
if nargin < 4, c = 1; end

if isstruct(gamma), gamma = gamma.paramF; end

CheckTypes({f,stGrid,gamma,c},...
    {'Function','Grid','function_handle','scalar'});

% Approximate spacing for samples along the circles.
ms = min(f.grid.ds) / 2;

stGV = stGrid.gridVectors();
sGV = stGV{1};
tGV = stGV{2};
p = gamma(sGV);
sz = size(stGrid);
Rfvals = zeros(sz);

for si = 1:sz(1)
    for ti = 1:sz(2)
        % Calculate radius for this circle.
        t = tGV(ti);
        r = t / c;
        
        % Calculate number of angles (nth) to sample
        %  at along this circle.
        nth = max(ceil(2*pi*r / ms), 4);
        dth = (2*pi) / nth;                 % dth: angle spacing
        th = (1:nth) * dth;                 % th: vector of angles
        xy = [p(1,si) + r*cos(th); ...
              p(2,si) + r*sin(th)];         % xy: coords of sample points
        v = f.SampleAt(xy);                 % v: sample values
        Rfvals(si,ti) = sum(v)*r*dth;       % Numerically integrate the
                                            %  sample values.
    end
end

Rf = Function.WithValues(stGrid, Rfvals);

end