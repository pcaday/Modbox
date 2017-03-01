% Compute the gradient of a vector-valued function.
%
% df = Gradient(f, domainGrid, codomainGrid)
%
% f is an n x I1 x I2 x ... x Ik array, representing an n-vector-
%  valued function on an I1 x I2 x ... x Ik grid.
% domainGrid, codomainGrid are Grid objects for the function's
%  domain and codomain.
%
% df is a n x k x I1 x I2 x ... x Ik array, containing numerical
%  derivatives of f with respect to each coordinate
%  (so df(i,j,:,...,:) contains df_i/dx_j, where x_j is the j'th
%   coordinate and f_i is the i'th coordinate function of f)
%
function df = Gradient(f, domainGrid, codomainGrid)
ds = domainGrid.ds;
ns = domainGrid.ns;
dims = domainGrid.dims;
df = zeros([codomainGrid.dims dims size(domainGrid)], class(f));

subAll(1:dims) = {':'};

for i = 1:dims
    n = ns(i);
    
    if domainGrid.periodic(i),        
        % Prepare subscripts for periodic centered differences
        sub1 = subAll;        sub1{i} = [2:n 1]; 
        sub2 = subAll;        sub2{i} = [n   1:n-1];
        df(:,i,subAll{:}) = PeriodicDiff(f(:,sub1{:}), f(:,sub2{:}), codomainGrid) ...
            ./ (2*ds(i));
    else
        % Prepare subscripts for centered differences for the interior
        %  points, and one-sided differences for the boundary points.
        sub1 = subAll;        sub1{i} = [2:n n]; 
        sub2 = subAll;        sub2{i} = [1   1:n-1];
        df(:,i,subAll{:}) = PeriodicDiff(f(:,sub1{:}), f(:,sub2{:}), codomainGrid) ...
            ./ (2*ds(i));
        % For the boundary points, need to multiply by 2.
        sub1{i} = [1 n]; 
        df(:,i,sub1{:}) = df(:,i,sub1{:}) * 2;
    end
end

end