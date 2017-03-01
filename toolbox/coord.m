% coord(p,i)
%
% Same as p(i,:,...,:) but removes the singleton first dimension.
%
function v = coord(p,i)
    sz = size(p);
    v = reshape(p(i,:), [sz(2:end) 1]);
end