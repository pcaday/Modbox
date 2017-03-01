classdef Bend2DDiffeo < Diffeo
    properties
        i           % The coordinate that will be bent (1 or 2)
        alpha       % The amount of bending
        center      % The center for bending (length-2 vector)
    end
    
    methods
        function diffeo = Bend2DDiffeo(i,alpha,center)
            if nargin < 2, alpha = 0.5; end
            if nargin < 3, center = [0; 0]; end
            diffeo.i = i;
            diffeo.alpha = alpha;
            diffeo.center = center(:);
        end
        
        %  n = diffeo.Dimension()
        %
        %  Return the number of dimensions of the domain and codomain
        %  of the diffeomorphism.
        function n = Dimension(diffeo) %#ok<MANU>
            n = 2;
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
        function [xt,xit] = ForwardTransform(diffeo, x, xi)
            alpha_ = diffeo.alpha;
            center_ = diffeo.center;
            i_ = diffeo.i;
            j_ = 3 - i_;
            xt = x;
            xt(i_,:,:) = xt(i_,:,:) + alpha_*(xt(j_,:,:) - center_(j_)).^2;
            if nargout >= 2,
                % Do the covectors too.
                xit = repmattosize(xi, size(x));
                xit(j_,:,:) = xit(j_,:,:)...
                    - xit(i_,:,:) .* (2*alpha_*(xt(j_,:,:) - center_(j_)));
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
        function [x,xi] = InverseTransform(diffeo, xt, xit)
            alpha_ = diffeo.alpha;
            center_ = diffeo.center;
            i_ = diffeo.i;
            j_ = 3 - i_;
            x = xt;
            x(i_,:,:) = x(i_,:,:) - alpha_*(x(j_,:,:) - center_(j_)).^2;
            if nargout >= 2,
                % Do the covectors too.
                xi = repmattosize(xit, size(xt));
                xi(j_,:,:) = xi(j_,:,:)...
                    + xi(i_,:,:) .* (2*alpha_*(x(j_,:,:) - center_(j_)));
            end
        end
        
        %  D = diffeo.DInverseCotangent(x,xi)
        %
        %   Find the first derivative of the pullback map on covectors
        %
        function [dxdxt, dxidxt, dxidxit] = DInverseCotangent(diffeo,x,xi)
            alpha_ = diffeo.alpha;
            center_ = diffeo.center;
            i_ = diffeo.i;
            j_ = 3 - i_;
            
            id = zeros([2 size(x)], class(x));
            id(1,1,:) = 1;
            id(2,2,:) = 1;
            
            diff = 2*alpha_*(x(j_,:,:) - center_(j_));
            
            dxdxt = id;
            dxdxt(i_,j_,:,:) = -diff;
            
            dxidxit = id;
            dxidxit(j_,i_,:,:) = diff;
            
            dxidxt = zeros([2 size(x)], class(x));
            dxidxt(j_,j_,:) = 2*alpha_*xi(i_,:);
        end
        
        %  D = diffeo.DInverseCotangentT(xt,xit)
        %
        %   Find the first derivative of the pullback map on covectors
        %   Like DInverseCotangent, but in codomain coordinates.
        %
        function [dxdxt, dxidxt, dxidxit] = DInverseCotangentT(diffeo,xt,xit)
            % Same as DInverseCotangent!
            [dxdxt,dxidxt,dxidxit] = diffeo.DInverseCotangent(xt,xit);
        end
        
        % newGrid = diffeo.TransformGrid(originalGrid)
        %
        %   Given a grid, return a grid of similar size that contains the
        %  image of the input grid under the diffeomorphism.
        %
        % newGrid and originalGrid are Grid objects.
        function newGrid = TransformGrid(diffeo, originalGrid)
            i_ = diffeo.i;
            j_ = 3 - i_;
            alpha_ = diffeo.alpha;
            center_ = diffeo.center;
            
            % If the original grid is periodic in the coordinate we're
            %  bending, no need to change it. Just let it wrap around!
            if originalGrid.periodic(i_)
                newGrid = originalGrid.copy;
                return;
            end
            
            % Find maxDeform, the maximum value of |x_i - x~_i|
            maxDeform = alpha_ * max(...
                ([originalGrid.mins(j_) originalGrid.maxes(j_)] ...
                  - center_(j_)) .^ 2);
            
            % Copy the xGrid and expand it appropriately.
            newGrid = originalGrid.copy();
            if maxDeform > 0,
                newGrid.maxes(i_) = newGrid.maxes(i_) + maxDeform;
            else
                newGrid.mins(i_) = newGrid.mins(i_) + maxDeform;
            end
            
            % Keep the same spacing in the x~ grid as was in the x grid.
            newGrid.setDS(originalGrid.ds);
        end
    end
end