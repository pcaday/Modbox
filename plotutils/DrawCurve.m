% Draw a curve described by a curveInfo structure.
%
%   DrawCurve(curveInfo)
%   DrawCurve(curveInfo, [smin smax])
%   DrawCurve(curveInfo, [smin smax],...)
%   h = DrawCurve(...)
%
% smin and smax, if given, are the range of s-values for the curve,
%  otherwise defaults are taken from the curveInfo structure, if specified
%  there.
% h is the handle to the drawn curve.
%
% Extra arguments after [smin smax] are passed to plot.
%  e.g., you can call DrawCurve(curveInfo, [smin smax], '--r') to draw the
%  curve in a red dashed line.
%
function h = DrawCurve(curveInfo, sRange, varargin)
    if nargin < 2 || isempty(sRange),
        smin = curveInfo.sRange(1);
        smax = curveInfo.sRange(2);
        if ~isfinite(smin), smin = -1; end  % Substitute defaults for infinite curves
        if ~isfinite(smax), smax = 1; end
        if smin >= smax, smax = smin + 2; end   % Avoid empty s range after substituting defaults.
    else
        smin = sRange(1);
        smax = sRange(2);
    end
    
    x = curveInfo.paramF(linspace(smin,smax,101).');    
    hh = plot(x(1,:).', x(2,:).', varargin{:});
    
    if nargout > 0, h = hh; end
end