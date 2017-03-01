function xyg = GridToXYTGrid(g, Tmax, dt)
    % Create an XYTGrid for the reflection.
    if nargin < 2, Tmax = 1.5 * norm(g.maxes - g.mins); end
    if nargin < 3, dt = 0.05; end
    
    nt = ceil(Tmax/dt)+1;
    xyg = XYTGrid(g.mins(1), g.maxes(1), g.mins(2), g.maxes(2), g.ns(1), g.ns(2), Tmax, nt);
end