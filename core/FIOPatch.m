classdef FIOPatch < handle
    properties
        degree      % Order of the symbol (= degree of FIO)
        parent      % Parent FIO object.
        dims        % Number of dimensions in the domain and codomain
    end
    
    methods
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
        function y = CanonicalTransformationY(patch,x,xi) %#ok<INUSD,STOUT>
            error('Either the CanonicalTransformationY or CanTransYAndPrincipalSymbol methods must be implemented.');
        end
        
        % [dxy, dxiy] = fio.DCanonicalTransformationY(x,xi)
        %
        %   Get the first derivatives of the canonical transformation, with
        %  respect to x (dxy) and xi (dxiy).
        %
        %  x and xi are as in CanonicalTransformation,
        %  while dxy, dxiy are n x n x anything x ... x anything arrays.
        %
        %   dxy(i,j,...) is the derivative of y_i with respect to x_j;
        %    similarly for dxiy.
        %
        %  Does not need to be implemented; if it is not, the derivatives
        %   will be approximated using finite differences.
        function [dxy, dxiy] = DCanonicalTransformationY(patch,x,xi)
            [dxy, dxiy] = patch.DCTYQuickFiniteDifferences(x,xi);
        end

        % [dxy, dxiy] = fio.DCTYQuickFiniteDifferences(x,xi)
        %
        %   Get the first derivatives of the canonical transformation, with
        %  respect to x (dxy) and xi (dxiy).
        %
        %   Unlike DCanonicalTransformationY,
        %       x is assumed to be an array of grid points, while
        %      xi is assumed to be a single vector.
        %
        %  dxy, dxiy are n x n x anything x ... x anything arrays.
        %
        %  dxy(i,j,...) is the derivative of y_i with respect to x_j;
        %   similarly for dxiy.
        %
        %  Does not need to be implemented; if it is not, the derivatives
        %   will be approximated using finite differences.
        function [dxy, dxiy] = DCTYQuickFiniteDifferences(patch,x,xi)
            e = 1e-3;
            n = patch.dims;
            if n ~= 2,
                warning('FIOPatch:DCTYQFDnot2D',...
                    'Currently only implemented for 2D');
                dxy = [];
                dxiy = [];
                return;
            end

            % xirot = [1 -e; e 1] * xi;
            xirot = zeros(size(xi));
            xirot(1,:) = xi(1,:) - e*xi(2,:);
            xirot(2,:) = xi(2,:) + e*xi(1,:);
   
            y = patch.CanonicalTransformationY(x,xi);
            yrot = patch.CanonicalTransformationY(x,xirot);

            dxiyPerp = PeriodicDiff(yrot, y, patch.parent.outLogicalGrid)...
                        * (1. / e);
            invNormXi2 = 1 ./ sum(xi.^2, 1);
            dxiy(:,1,:,:) = bsxfun(@times, dxiyPerp, invNormXi2) * -xi(2);
            dxiy(:,2,:,:) = bsxfun(@times, dxiyPerp, invNormXi2) * +xi(1);
            

            dxy = zeros([n size(x)], class(x));
            
%             if isequal(size(xi), [patch.dims 1]),
%                 h1 = x(1,2,1) - x(1,1,1);
%                 h2 = x(2,1,2) - x(2,1,1);
% 
%                 dxy(:,1,:,:) = PeriodicDiff(y(:,[2:end end],:), ...
%                                             y(:,[1 1:end-1],:), ...
%                                             patch.parent.outLogicalGrid) ...
%                                 * (1 / (2*h1));
%                 dxy(:,1,[1 end],:) = dxy(:,1,[1 end],:) * 2;
% 
%                 dxy(:,2,:,:) = PeriodicDiff(y(:,:,[2:end end]), ...
%                                             y(:,:,[1 1:end-1]), ...
%                                             patch.parent.outLogicalGrid) ...
%                                 * (1 / (2*h2));
%                 dxy(:,2,:,[1 end]) = dxy(:,2,:,[1 end]) * 2;
%             else
                xm = x;
                xm(1,:,:) = xm(1,:,:) + e;
                yshift = patch.CanonicalTransformationY(xm, xi);
                dxy(:,1,:,:) = PeriodicDiff(yshift, y, ...
                                            patch.parent.outLogicalGrid) ...
                                * (1 / e);
                xm = x;
                xm(2,:,:) = xm(2,:,:) + e;
                yshift = patch.CanonicalTransformationY(xm, xi);
                dxy(:,2,:,:) = PeriodicDiff(yshift, y, ...
                                            patch.parent.outLogicalGrid) ...
                                * (1 / e);
%             end

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
        function ps = PrincipalSymbol(patch,x,xi) %#ok<INUSL,INUSD>
            sz = size(x);
            ps = ones([sz(2:end) 1], class(x));
        end
        
        % [y,a] = fio.CanTransYAndPrincipalSym(x,xi)
        %
        %   Apply the canonical transformation and compute principal symbol
        %  simultaneously. (Combination of methods CanonicalTransformationY
        %  and PrincipalSymbol)
        %
        %   Does not need to be implemented; if is not, the
        %  CanonicalTransformationY and PrincipalSymbol methods will be
        %  used instead.
        %
        %  x is an n x anything x ... x anything array.
        %  xi can be an n-vector (same xi for all x's), or an array with
        %    the same size as x.
        %
        %  y: y coordinates of the input covectors after applying the
        %      canonical transformation
        %          (n x anything x ... x anything)
        % ps: principal symbol
        %          (anything x ... x anything)
        function [y,ps] = CanTransYAndPrincipalSym(patch,x,xi) %#ok<INUSD>
            y = [];
            ps = [];
        end
     end
end