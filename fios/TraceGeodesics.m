% Approximate geodesics using raytracing.
%
%  [ys,etas] = TraceGeodesics(metric,x,xi,sgrid)
%  [y,eta] = TraceGeodesics(metric,x,xi,sgrid,'final')
%
% Inputs
%       metric - A Metric object representing the Riemannian metric to use.
%            x - n x anything x ... x anything array of initial points
%           xi - n x anything x ... x anything array of initial directions
%                    (same size as x, or a single vector to use the same
%                     direction for all initial points.)
%        sgrid - A Grid object representing the coordinate along the 
%                 geodesics. Its grid spacing is used to determine step
%                 size when computing geodesics.
%
% Outputs
%  If the 'finals' flag is included, only the final covectors are returned:
%            y - n x anything x ... x anything array of final points
%          eta - n x anything x ... x anything array of final directions
%
%  Otherwise, the whole geodesic is returned:
%           ys - Function object whose value at each s is the set of
%                 points of all geodesics at that s.
%         etas - Function object whose value at each s is the set of
%                 directions of all geodesics at that s.
%
% Standard RK4 is used for time-stepping the geodesic equation.

function [y,eta] = TraceGeodesics(metric,x,xi,sGrid,flag)

if nargin < 5,
    flagFinalOnly = false;
else
    flagFinalOnly = strcmp(flag, 'final');
end

saveYs   = ~flagFinalOnly && nargout > 0;
saveEtas = ~flagFinalOnly && nargout > 1;


% Subroutine for computing dv/dt for the geodesic equation.
    function dvdt = GetDvdt(x,v)
        dvdt = zeros(size(v));
        if isConformallyFlat2D,
            % Special case for 2D conformally flat metrics
            %  (for speed).
            dl1 = dlambda{1}(x);
            dl2 = dlambda{2}(x);
            a = v(1,:).^2 - v(2,:).^2;
            b = 2*v(1,:).*v(2,:);
            dl1 = dl1(:);
            dl2 = dl2(:);
            a = a(:);
            b = b(:);
            dvdt(1,:) = -a.*dl1 - b.*dl2;
            dvdt(2,:) =  a.*dl2 - b.*dl1;
        else
            % General case, using Christoffel symbols.
            for k = 1:n
                for i = 1:n
                    for j = 1:n
                        dvdt(k,:) = dvdt(k,:) - v(i,:) .* v(j,:) .* ...
                            metric.Christoffel(i,j,k,x(:,:)).';
                    end
                end
            end
        end
    end

ds = sGrid.ds;
nt = sGrid.ns;
n = metric.dims;

isConformallyFlat2D = (n == 2) && isa(metric, 'EuclideanConformalMetric');
if isConformallyFlat2D,
    dlambda = metric.dlambda;
end

% Expand xi if necessary (expand to the size of x)
xi = repmattosize(xi, size(x));

% xc, vc: current x and v (velocity) values.
vc = metric.UnitVector(x, xi);
xc = x;

if saveYs,   y_   = zeros(numel(x), nt); end
if saveEtas, eta_ = zeros(numel(x), nt); end

% Loop over time steps
for ts = 1:nt
    % Solve the geodesic equation with RK4.
    
    % Stage 1
    q1x = ds * vc;
    q1v = ds * GetDvdt(xc,vc);
    
    % Stage 2
    vm = vc + 0.5*q1v;
    q2x = ds * vm;
    q2v = ds * GetDvdt(xc + 0.5*q1x, vm);
    
    % Stage 3
    vm = vc + 0.5*q2v;
    q3x = ds * vm;
    q3v = ds * GetDvdt(xc + 0.5*q2x, vm);
    
    % Stage 4
    vm = vc + q3v;
    q4x = ds * vm;
    q4v = ds * GetDvdt(xc + q3x, vm);
    
    % Step
    xc = xc + 1/6 * (q1x + 2*q2x + 2*q3x + q4x);
    vc = vc + 1/6 * (q1v + 2*q2v + 2*q3v + q4v);
    
    % Save
    if saveYs,   y_(:,ts)   = xc(:); end
    if saveEtas, eta_(:,ts) = vc(:); end
end

if flagFinalOnly,
    y   = xc;
    eta = vc;
else
    shape = [size(x) nt];
    if saveYs,   y   = Function.WithValues(sGrid, reshape(y_,shape));   end
    if saveEtas, eta = Function.WithValues(sGrid, reshape(eta_,shape)); end
end

end