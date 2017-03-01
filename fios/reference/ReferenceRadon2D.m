% Direct, brute-force Radon transform
function Rf = ReferenceRadon2D(f, yGrid)
    Rfvals = zeros(size(yGrid));
    gvs = yGrid.gridVectors();
    
    dx = min(f.grid.ds);        % Distance between sample points along the
                                % hyperplanes (lines)
    % Find furthest point from the origin.
    furthest = max(abs(yGrid.maxes), abs(yGrid.mins));
    % Get its distance from the origin.
    diam = hypot(furthest(1),furthest(2));
    % Prepare vector of points to sample along each line.
    p = linspace(-diam, diam, ceil((2*diam)/dx))';
    
    % Loop over s and theta (not very efficient!)
    for thi = 1:size(yGrid,2)
        th = gvs{2}(thi);               % Current theta value
        om = [cos(th) sin(th)];         % Unit vector perpendicular to line
        omp = [-sin(th) cos(th)];       % Unit vector parallel to line.
        for si = 1:size(yGrid,1)
            s = gvs{1}(si);             % Current s value.
            % Sample the function along this line.
            samples = f.SampleAt(bsxfun(@plus, s*om, p*omp).');
            % Integrate the samples to approximate Rf(s,theta)
            Rfvals(si,thi) = sum(samples(~isnan(samples))) * dx;
        end
    end
    
    Rf = Function.WithValues(yGrid, Rfvals);
end