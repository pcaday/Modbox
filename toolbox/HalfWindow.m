% Window function that ramps smoothly from 0 to 1
%
%   v = HalfWindow(u, off, on, kind)
%
% The window function ramps smoothly from u = 0 at off
%  to u = 1 at on.
%
% off can be larger or smaller than on:
%
%  1            --|--       --|--
%              /                 \
%             /         OR        \
%  0    ---|--                     --|---
%
%         off     on          on    off
%
% u can be a matrix or a scalar; v will be the same size as u.
%
% kind can be one of:
%    'smooth'             Smooth (to all orders) cutoff. Symmetric
%                           about u = 0.5.
%    'l2normalized'       Modified version of 'smooth' such that
%                           w(u)^2 + w(1-u)^2 = 1 for all u.
%    'linear'             Linear ramp function.
%                           w(u) = u.
%    'sine'               Sine-like ramp function, C^2 continuous
%                           w(u) = u - sin(2*pi*u)/(2*pi) for 0 <= u <= 1.
%    'cubic'              Cubic ramp function, C^1 continuous
%                           w(u) = -2u^3 + 3u^2           for 0 <= u <= 1.
%
function v = HalfWindow(u, off, on, kind)
v = HalfWindow01((u - off) / (on - off), kind);
end


% Window function v(u) that ramps smoothly from 0
%  at u = 0 to 1 at u = 1.
% It is 0 for u <= 0 and 1 for u >= 1,
%  and has the property
%             v(u)^2 + v(1-u)^2 = 1.
%
% u can be a matrix or a scalar.
%
function v = HalfWindow01(u, kind)
switch kind
    case 'l2normalized'
        % C^infinity
        v0 = 1 ./ (1 + exp(-1 ./ (1-u)) ./ exp(-1 ./ u));       % unnormalized window
        vflip = 1 ./ (1 + exp(-1 ./ u) ./ exp(-1 ./ (1-u)));
        
        v = v0 ./ sqrt(v0.^2 + vflip.^2);
    case 'smooth'
        % C^infinity
        v = 1 ./ (1 + exp(-1 ./ (1-u)) ./ exp(-1 ./ u));
    case 'linear'
        % C^0
        v = u;
    case 'sine'
        % C^2
        v = u - sin(2*pi*u)/(2*pi);
    case 'cubic'
        % C^1
        v = -2*u.^3 + 3*u.^2;
end

v(u >= 1) = 1;
v(u <= 0) = 0;
v(isnan(u)) = 0;        % NaN's come if off = on in HalfWindow,
                        %  and u = off = on. Arbitrarily assign 0.

end