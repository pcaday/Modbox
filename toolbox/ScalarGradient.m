% Take the gradient of a Function object f.
%
%  dfs = ScalarGradient(f)
%
% f is a scalar Function object.
% dfs is a length-n cell array (n = #dims of f), where dfs{i} is
%  a Function object representing df/dx_i.
function dfs = ScalarGradient(f)
    assert(isequal(f.m, 1), 'f must be a scalar function.');
    
    % Create a dummy codomain grid 
    codomainGrid = Grid(0, 1, 2);
    
    g = f.grid;
    n = g.dims;
    f = reshape(f.f, [1 size(f.f)]);
    
    df = Gradient(f, g, codomainGrid);
    dfs = cell(1,n);
    for i = 1:n
        dfs{i} = Function.WithValues(g, fsqueeze(coord(df,i)));
    end
end