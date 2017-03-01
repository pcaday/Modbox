classdef Radon2DPatch < FIOPatch
    properties
        sign
    end
    
    methods
        function patch = Radon2DPatch(sign, parent)
            patch.sign = sign;
            
            patch.dims = 2;
            patch.degree = -1/2;
            patch.parent = parent;
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
        function y = CanonicalTransformationY(patch,x,xi) 
            s = x(1,:,:).*xi(1,:,:) + x(2,:,:).*xi(2,:,:);
            s = (patch.sign * s) ./ hypot(xi(1,:,:), xi(2,:,:));
            th = atan2(xi(2,:,:), xi(1,:,:));
            if patch.sign == -1,
                th = th + pi;
            end
            th = mod(th, 2*pi);
            if numel(th) == 1,
                th = repmat(th, size(s));
            end
            y = cat(1,s,th);
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
        function [dxy, dxiy] = DCanonicalTransformationY(patch,x,xi) 
            % Debugging: disables using exact derivatives.
%             [dxy,dxiy] = DCanonicalTransformationY@FIOPatch(patch,x,xi);
%             return;
             
            
            valueClass = class(x);
            
            sz = size(x);
            sz = [sz(2:end) 1];
            dxy = zeros([2 2 sz], valueClass);
            dxiy = zeros([2 2 sz], valueClass);
            
            absXi = hypot(xi(1,:), xi(2,:));
            xi(1,:) = xi(1,:) ./ absXi;
            xi(2,:) = xi(2,:) ./ absXi;

            xiA2 = zeros(size(xi), valueClass);
            xiA2(1,:) = xi(1,:) ./ absXi;
            xiA2(2,:) = xi(2,:) ./ absXi;

            xi = xi * patch.sign;
            xdotxip = x(1,:).*-xi(2,:) + x(2,:).*xi(1,:);
            
            dxy(1,1,:) = xi(1,:);
            dxy(1,2,:) = xi(2,:);            
            
            dxiy(1,1,:) = -xiA2(2,:) .* xdotxip;
            dxiy(1,2,:) = +xiA2(1,:) .* xdotxip;
            dxiy(2,1,:) = -xiA2(2,:);
            dxiy(2,2,:) = +xiA2(1,:);
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
        function ps = PrincipalSymbol(patch,x,xi) %#ok<INUSL>
			% The principal symbol for the 2D Radon transform is sqrt(2*pi/|xi|).
			absXi = fsqueeze(hypot(xi(1,:,:), xi(2,:,:)));
            ps = sqrt(2*pi ./ absXi);
            
            % If xi is a single vector, expand ps into a constant function
            %  with the same value for every x.
            sz = size(x);
            ps = repmattosize(ps, [sz(2:end) 1]);
        end
    end
end