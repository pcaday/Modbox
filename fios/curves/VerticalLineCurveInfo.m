function curveInfo = VerticalLineCurveInfo(line_x)
    if nargin < 1, line_x = 0; end
    
    curveInfo.lrTestF = @(x) VLineLRTest(x,line_x);
    curveInfo.tangentF = @(x) VLineTangent(x);
    curveInfo.changeXToS = @(x) VLineS(x);
    curveInfo.paramF = @(s) VLineParam(s, line_x);
    curveInfo.definingF = @(x) VLineDefiningFunction(x);
    curveInfo.curvatureF = @(x) VLineCurvature(x);
    curveInfo.hasInside = false;
    curveInfo.sRange = [-inf inf];        % Range of s parameter values
end


function lr = VLineLRTest(x, line_x)
    lr = fsqueeze(x(1,:,:) < line_x);
end

function t = VLineTangent(x)
    t = zeros(size(x), class(x));
    t(2,:) = 1;
end

function x = VLineParam(s, line_x)
    x = fcat(repmat(line_x, size(s)),s);
end

function s = VLineS(x)
    s = fsqueeze(x(2,:,:));
end

function a = VLineDefiningFunction(x)
    a = coord(x,1);
end

function c = VLineCurvature(x)
    sz = size(x);
    c = zeros([sz(2:end) 1]);
end