% Take the difference of two functions taking values on
%  a (potentially) periodic grid.
%
%   d = PeriodicDiff(f1,f2,grid)
%
%   f1,f2: vector-valued functions (n x anything x ... x anything arrays)
%    grid: a grid object
%
%       d: the difference between f1 and f2
function d = PeriodicDiff(f1,f2,grid)

d = f1 - f2;
isPeriodic = grid.periodic;
p = (grid.maxes-grid.mins) .* isPeriodic(:).';
for i = 1:grid.dims
    if isPeriodic(i),
        d(i,:) = mod(d(i,:), p(i));
        d(i,:) = d(i,:) - p(i) * (d(i,:) > p(i)/2);
    end
end

end