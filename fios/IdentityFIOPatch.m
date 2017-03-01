classdef IdentityFIOPatch < FIOPatch
    methods
        function patch = IdentityFIOPatch(parent)
            patch.parent = parent;
            patch.dims = parent.dims;
            patch.degree = 0;
        end
        
        % y = fio.CanonicalTransformationY(x,xi)
        %
        %   Apply the canonical transformation to the covectors (x,xi)
        %  and return the y values.
        %
        %  x is an n x anything x ... x anything array.
        %  xi can be an n-vector (same xi for all x's), or an array with
        %    the same size as x.
        %
        %  The output y is an array of the same dimensions as x.
        function y = CanonicalTransformationY(patch,x,xi) %#ok<INUSL,INUSD>
            y = x;
        end
        
        % [dxy, dxiy] = fio.DCanonicalTransformationY(x,xi)
        %
        %   Get the first derivatives of the canonical transformation, with
        %  respect to x (dxy) and xi (dxiy).
        %
        %  x and xi are as in CanonicalTransformation,
        %  while dxy, dxiy are n x n x anything x ... x anything arrays.
        %
        %   dxy(i,j,...) is the derivative of y_i with respect to x_j.
        %
        %  Does not need to be implemented; if it is not, the derivatives
        %   will be approximated using finite differences.
        function [dxy, dxiy] = DCanonicalTransformationY(patch,x,xi) %#ok<INUSD>
            sz = size(x);
            sz(1) = [];
            sz(end+1) = 1;
            
            n = patch.dims;
            
            dxiy = zeros([n n sz], class(x));
            dxy = repmat(eye(n, class(x)), [1 1 sz]);
        end

        % ps = fio.PrincipalSymbol(x,xi)
        %
        %   Evaluate the principal symbol at given (x,xi) values.
        %
        %  x is an n x anything x ... x anything array.
        %  xi can be an n-vector (same xi for all x's), or an array with
        %    the same size as x.
        %
        %  The output ps is an anything x ... x anything array (same as x
        %   but without the leading dimension).
        %
        %  By default, the principal symbol is identically 1.
        %
        % function ps = PrincipalSymbol(patch,x,xi)
     end
end