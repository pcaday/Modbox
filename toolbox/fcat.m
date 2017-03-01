% fcat: "First cat"
%
%  B = fcat(A1,...,An)
%
% Add a new singleton dimension at the beginning of A1 ... An
% and concatenate the arrays along the new first dimension.

function B = fcat(varargin)
    n = ndims(varargin{1});
    B = permute(cat(n+1, varargin{:}), [n+1 1:n]);
end