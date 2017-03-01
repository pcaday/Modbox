% A CurveInfo structure describing the parabola y = ax^2.
%
function curveInfo = ParabolaCurveInfo(a)
    if nargin < 1, a = 0.5; end
    
    curveInfo.lrTestF = @(x) ParabolaLRTest(x,a);
    curveInfo.tangentF = @(x) ParabolaTangent(x, a);
    curveInfo.changeXToS = @(x) ParabolaS(x);
    curveInfo.paramF = @(s) ParabolaParam(s,a);
    curveInfo.definingF = @(x) ParabolaDefiningFunction(x,a);
    curveInfo.curvatureF = @(x) ParabolaCurvature(x,a);
    curveInfo.hasInside = true;
    
    curveInfo.sRange = [-inf inf];        % Range of s parameter values
end

% Check whether points are on the left or right of the curve.
% Also doubles as an inside test for closed curves like the circle.
function lr = ParabolaLRTest(x,a)
    lr = (ParabolaDefiningFunction(x,a) > 0);
end

% Get the tangent(s) to the curve at given point(s)
function t = ParabolaTangent(x, a)
    t = zeros(size(x), class(x));
    t(1,:) = 1;
    t(2,:) = 2*a*x(1,:);
end

% Convert from x coordinates to a coordinate s along the curve.
% The choice of the parameterization s is up to the curve.
%
% Here we're taking the x1 coordinate.
function s = ParabolaS(x)
    s = coord(x,1);
end

function x = ParabolaParam(s, a)
    x = fcat(s, a*s.*s);
end

% Defining function for the curve.
%
% The curve is the zero set of the defining function.
%
function v = ParabolaDefiningFunction(x,a)
    v = coord(x,2) - a*coord(x,1).^2;
end



% Curvature function for the curve.
%
% For the parabola, the curvature is
%                2a
%   c(y) =  ------------
%            1 + 4a^2y^2
function c = ParabolaCurvature(x,a)
    y = coord(x,2);
    c = 2*a ./ (1 + 4*a*a*y.*y);
end