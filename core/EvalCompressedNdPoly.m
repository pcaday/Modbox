% Evalute the polynomials given by (H,I), with dimension d, at a
% point or a vector of points:
%
%  vals = EvalCompressedNdPoly(H,I,X)
%
%  or
%
%  vals = EvalCompressedNdPoly(H,I,X1,X2,...,Xd);
%
% In the first form, X is an n x d matrix of points; each row of X
% represents one point, and vals is a n x h matrix, where h is the number
% of spherical harmonics.
%
% In the second form, X1,...,Xd are arrays of coordinates (of the same
% dimension), and vals is a matrix of the same dimensions as X1...Xd.
%
% 


function vals = EvalCompressedNdPoly(h,I,varargin)
    d = size(I,2);
    if length(varargin) == d,
        % We have d arguments. Reshape them into columns and concatenate
        % them into a d x m matrix.
        %
        % Also remember their old shape, so we can restore it later.
        %
        % Do this case first, so that if d == 1, we follow this route.
        finalShape = size(varargin{1});
        icols = cellfun(@(Xi) Xi(:), varargin, 'UniformOutput', false);
        inputs = [icols{:}];
    elseif length(varargin) == 1,
        % We have 1 argument. Check that it's a matrix with d columns.
        inputs = varargin{1};
        finalShape = [size(inputs,1) 1];
        assert(size(inputs,2)==d, 'EvalHarmonicsBasisAt:badInput',...
            'Expected a matrix with %d columns.', d);
    else
        error('EvalHarmonicsBasisAt:wrongNumInputs',...
            'Expected %d inputs or a matrix with %d columns.', d, d);
    end
    
    % Evaluate!
    % Take the log of the inputs, and replace any -Inf's by the largest
    % possible negative number (so if they are multiplied by zero entries
    % in I, we get 0 instead of NaN's).
    loginputs = log(inputs)';
    loginputs(isinf(loginputs)) = -realmax;
    % Then evaluate each monomial in I on the inputs by multiplying I by
    % the log of the inputs and exponentiating, and then get the polynomial
    % values by multiplying by h.
    vals = h'*exp(I*loginputs);
    % Reshape :)
    vals = reshape(vals, finalShape);
    
    % Remove imaginary parts if everything should be real.
    if isreal(inputs) && isreal(h)
        vals = real(vals);
    end
end
