classdef IdentityDiffeo < Diffeo
    properties
        dims
    end
    
    
    methods
        function diffeo = IdentityDiffeo(dims)
            diffeo.dims = dims;
        end
        
        %  n = diffeo.Dimension()
        %
        %  Return the number of dimensions of the domain and codomain
        %  of the diffeomorphism.
        function n = Dimension(diffeo)
            n = diffeo.dims;
        end
        
        %  xt = diffeo.ForwardTransform(x)
        %
        %   Apply the diffeomorphism to the points x, where
        %  x is an n x anything x ... x anything array.
        %
        %   Returns an array xt of the transformed points, of the same
        %  dimensions as x.
        %
        %
        %  [xt,xit] = diffeo.ForwardTransform(x,xi)
        %
        %   Push forward the covectors (x,xi) by the diffeomorphism,
        %  where x and xi are both n x anything x ... x anything arrays.
        %
        %   Returns arrays (xt,xit) of the pushforwards, of the same
        %  dimensions as (x,xi).
        function [xt,xit] = ForwardTransform(diffeo, x, xi) %#ok<INUSL>
            xt = x;
            if nargin == 3,
                xit = xi;
            end
        end
        
        %  x = diffeo.InverseTransform(xt)
        %
        %   Apply the inverse of the diffeomorphism to the points xt, where
        %  xt is an n x anything x ... x anything array.
        %
        %   Returns an array xt of the transformed points, of the same
        %  dimensions as x.
        %
        %
        %  [x,xi] = diffeo.InverseTransform(xt,xit)
        %
        %   Pull back the covectors (xt,xit) by the diffeomorphism,
        %  where xt and xit are both n x anything x ... x anything arrays.
        %
        %   Returns arrays (x,xi) of the pullbacks, of the same dimensions
        %  as (xt,xit).
        function [x,xi] = InverseTransform(diffeo, xt, xit) %#ok<INUSL>
            x = xt;
            if nargin == 3,
                xi = xit;
            end
        end
        
        %  D = diffeo.DInverseCotangent(x,xi)
        %
        %   Find the first derivative of the pullback map on covectors
        %
        function [dxdxt, dxidxt, dxidxit] = DInverseCotangent(diffeo,x,xi) %#ok<INUSD>
            sz = size(x);
            n = diffeo.dims;
            zero = zeros([n sz], class(x));
            id = zero;
            for i = 1:n
                id(i,i,:) = 1;
            end
            dxdxt = id;
            dxidxit = id;
            dxidxt = zero;
        end
        
        %  D = diffeo.DInverseCotangentT(xt,xit)
        %
        %   Find the first derivative of the pullback map on covectors
        %   Like DInverseCotangent, but in codomain coordinates.
        %
        function [dxdxt, dxidxt, dxidxit] = DInverseCotangentT(diffeo,xt,xit)
            % No difference!
            [dxdxt,dxidxt,dxidxit] = diffeo.DInverseCotangent(xt,xit);
        end
        
        
        
        % newGrid = diffeo.TransformGrid(originalGrid)
        %
        %   Given a grid, return a grid of similar size that contains the
        %  image of the input grid under the diffeomorphism.
        %
        % newGrid and originalGrid are Grid objects.
        function newGrid = TransformGrid(diffeo, originalGrid)        %#ok<INUSL>
            newGrid = originalGrid;
        end
    end
end