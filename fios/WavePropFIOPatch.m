classdef WavePropFIOPatch < FIOPatch
    properties
        t
    end
    
    methods
        function patch = WavePropFIOPatch(t)
            if nargin < 1, t = 0.2; end
            
            patch.dims = 2;
            patch.degree = 0;
            
            patch.t = t;
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
            absXi = sqrt(xi(1,:).^2+xi(2,:).^2);
            dOverAbsXi = patch.t ./ absXi;
            
            y = zeros(size(x));
            y(1,:) = x(1,:) + xi(1,:) .* dOverAbsXi;
            y(2,:) = x(2,:) + xi(2,:) .* dOverAbsXi;
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
            % Disable exact derivatives (for testing):
            % [dxy,dxiy] = DCTYQuickFiniteDifferences(patch,x,xi);
            % return;
            
            sz = size(x);
            sz(1) = [];
            sz(end+1) = 1;
            
            absXi = sqrt(xi(1,:).^2+xi(2,:).^2);
            dAbsXiN3 = patch.t ./ (absXi.*absXi.*absXi);
            
            dxiy = zeros([2 2 sz], class(x));
            dxiy(1,1,:) =  xi(2,:).*xi(2,:).*dAbsXiN3;
            dxiy(1,2,:) = -xi(1,:).*xi(2,:).*dAbsXiN3;
            dxiy(2,1,:) =  dxiy(1,2,:);
            dxiy(2,2,:) =  xi(1,:).*xi(1,:).*dAbsXiN3;
            
            dxy = repmat(eye(2, class(x)), [1 1 sz]);
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
        function ps = PrincipalSymbol(patch,x,xi)
            global GlobalProcessNewMaslov
            global GlobalProcessMaslovSign
            
            % Get constant 1 principal symbol
            ps = PrincipalSymbol@FIOPatch(patch,x,xi);
            
            % Multiply by Maslov factor
            if GlobalProcessNewMaslov,
                if patch.t < 0,
                    ps = ps * (1i * GlobalProcessMaslovSign);
                end
            else
                ps = ps * exp(-pi*1i/4 * sign(patch.t) * GlobalProcessMaslovSign);
            end
        end
     end
end