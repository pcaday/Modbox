classdef Diffeo < handle
    methods (Abstract)
        %  n = diffeo.Dimension()
        %
        %  Return the number of dimensions of the domain and codomain
        %  of the diffeomorphism.
        n = Dimension(diffeo)
        
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
        [xt,xit] = ForwardTransform(diffeo, x, xi)   
        
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
        [x,xi] = InverseTransform(diffeo, xt, xit)
        
        %  D = diffeo.DInverseCotangent(x,xi)
        %
        %   Find the first derivative of the pullback map on covectors
        %
        [dxdxt, dxdxit, dxidxt, dxidxit] = DInverseCotangent(diffeo, x, xi)
        
        %  D = diffeo.DInverseCotangentT(xt,xit)
        %
        %   Find the first derivative of the pullback map on covectors
        %   Like DInverseCotangent, but in codomain coordinates.
        %
        [dxdxt, dxdxit, dxidxt, dxidxit] = DInverseCotangentT(diffeo, xt, xit)
        
        % newGrid = diffeo.TransformGrid(originalGrid)
        %
        %   Given a grid, return a grid of similar size that contains the
        %  image of the input grid under the diffeomorphism.
        %
        % newGrid and originalGrid are Grid objects.
        tGrid = TransformGrid(diffeo, originalGrid)       
        
    end
end