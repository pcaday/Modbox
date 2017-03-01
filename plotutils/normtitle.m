%  Add a title to the plot consisting of a given string and
%   the L^2 norm of the given function (plain array or Function object)
%
%  normtitle(string, f)
%
function normtitle(string, f)
    if isa(f, 'Function')
        dx = prod(f.grid.ds);
        f = f.f;
    else
        dx = 1;
    end
    normF = norm(f(:)) * sqrt(dx);
    if ~isempty(string)
        title(sprintf('%s (L^2 norm = %g)', string, normF));
    else
        title(sprintf('L^2 norm = %g', normF));        
    end
end