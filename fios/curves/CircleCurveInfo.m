function curveInfo = CircleCurveInfo(rad, thetas)
    if nargin < 1, rad = 0.5; end
    if nargin < 2, thetas = [0 2*pi]; end
    
    curveInfo.lrTestF = @(x) CircleLRTest(x,rad);
    curveInfo.tangentF = @(x) CircleTangent(x);
    curveInfo.changeXToS = @(x) CircleS(x);
    curveInfo.paramF = @(s) CircleParam(s, rad);
    curveInfo.definingF = @(x) CircleDefiningFunction(x,rad);
    curveInfo.curvatureF = @(x) CircleCurvature(x,rad);
    
    curveInfo.sRange = thetas;        % Range of s parameter values
    curveInfo.hasInside = true;
end

% Check whether points are on the left or right of the curve.
% Also doubles as an inside test for closed curves like the circle.
function lr = CircleLRTest(x, rad)
    lr = CircleDefiningFunction(x,rad) < 0;
end

% Get the tangent(s) to the curve at given point(s)
function t = CircleTangent(x)
    t = zeros(size(x), class(x));
    t(1,:) = -x(2,:);
    t(2,:) =  x(1,:);
end

% Convert from x coordinates to a coordinate s along the curve.
% The choice of the parameterization s is up to the curve.
%
% Here we're taking the theta value (as in polar coordinates) for s.
function s = CircleS(x)
    s = mod(atan2(coord(x,2), coord(x,1)), 2*pi);
end

% Parameterizing function for the curve
%  Given s coordinate values, returns the corresponding x coordinates.
function x = CircleParam(s, rad)
    x = rad * fcat(cos(s), sin(s));
end

% Defining function for the curve.
%
% The curve is the zero set of the defining function.
% Note only full circles are currently implemented...
function v = CircleDefiningFunction(x,rad)
    v = coord(x,1).^2 + coord(x,2).^2 - rad*rad;
end

% Curvature function for the curve.
%
% For the circle, constant 1/r.
function c = CircleCurvature(x,rad)
    sz = size(x);
    c = repmat(1/rad, [sz(2:end) 1]);
end