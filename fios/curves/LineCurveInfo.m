function curveInfo = LineCurveInfo(origin, direction)
    if nargin < 1, origin = [0; 0]; end
    if nargin < 2, direction = [0; 1]; end
    
    % Make the direction and origin column vectors; also make the direction
    %  a unit vector.
    direction = direction(:) / norm(direction);
    origin = origin(:);
    
    curveInfo.lrTestF = @(x) LineLRTest(x,origin,direction);
    curveInfo.tangentF = @(x) LineTangent(x,direction);
    curveInfo.changeXToS = @(x) LineS(x,origin,direction);
    curveInfo.paramF = @(s) LineParam(s,origin,direction);
    curveInfo.definingF = @(x) LineDefiningFunction(x,origin,direction);
    curveInfo.curvatureF = @(x) LineCurvature(x);
    curveInfo.hasInside = false;
    curveInfo.sRange = [-inf inf];        % Range of s parameter values
end


function lr = LineLRTest(x,origin,direction)
    lr = LineDefiningFunction(x,origin,direction) < 0;
end

function t = LineTangent(x,direction)
    sz = size(x);
    t = repmat(direction, [1 sz(2:end)]);
end

function x = LineParam(s,origin,direction)
    x = zeros([2 size(s)]);
    x(1,:) = origin(1) + direction(1)*s;
    x(2,:) = origin(2) + direction(2)*s;
end

function s = LineS(x,origin,direction)
    s = (coord(x,1)-origin(1))*direction(1) ...
       +(coord(x,2)-origin(2))*direction(2);
end

function a = LineDefiningFunction(x,origin,direction)
    a = (coord(x,1)-origin(1))*+direction(2) ...
       +(coord(x,2)-origin(2))*-direction(1);
end

function c = LineCurvature(x)
    sz = size(x);
    c = zeros([sz(2:end) 1]);
end