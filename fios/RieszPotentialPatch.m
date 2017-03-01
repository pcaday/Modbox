classdef RieszPotentialPatch < FIOPatch
    properties
        alpha
    end
    
    methods
        function patch = RieszPotentialPatch(alpha,dims)
            if nargin < 1, alpha = 1; end
            if nargin < 2, dims = 2; end
            patch.dims = dims;
            patch.degree = -alpha;
            patch.alpha = alpha;
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
        function ps = PrincipalSymbol(patch,x,xi) %#ok<INUSD>
            sz = size(x);
            sz(1) = [];
            sz(end+1) = 1;
            
            % Principal symbol is constant: (2pi)^(-alpha).
            psval = (2*pi)^(-patch.alpha);
            psval = cast(psval, class(x));
            ps = repmat(psval, sz);
        end
     end
end