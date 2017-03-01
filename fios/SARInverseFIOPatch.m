classdef SARInverseFIOPatch < FIOPatch
    properties
        curveInfo
        side
    end
    
    methods
        function patch = SARInverseFIOPatch(curveInfo,side)
            patch.dims = 2;
            patch.degree = 1/2;
            
            patch.curveInfo = curveInfo;
            patch.side = side;
        end
        
        % y = fio.CanonicalTransformationY(x,xi)
        %
        %   Apply the canonical transformation.
        %
        %  x is an n x anything x ... x anything array.
        %  xi can be an n-vector (same xi for all x's), or an array with
        %    the same size as x.
        %
        %  y: y coordinates of the input covectors after applying the
        %      canonical transformation
        %          (n x anything x ... x anything)
        %
        % Here x = (r,s); xi = (rho,sigma), and "y" is what we would
        %  normally call "x" (because the x and y coordinates are
        %  reversed.)
        function y = CanonicalTransformationY(patch,x,xi) 
            % Note this is not at all optimized for speed.
            rho = xi(2,:);
            sig = xi(1,:);
            as = sqrt(rho.*rho - sig.*sig) .* sign(rho) * patch.side;
            
            % The domain of the canonical transformation
            %  is |rho| >= |sigma|; outside of this the canonical
            %  transformation is NaN's
             
            if isequal(size(x),size(xi)),
                as(abs(rho) < abs(sig)) = NaN;
            elseif abs(rho) < abs(sig),
                as(:) = NaN;
            end
            
            r = x(2,:);
            s = x(1,:);
            
            p = patch.curveInfo.paramF(s);          % Points on curves
            t = patch.curveInfo.tangentF(p);        % Tangent vector
            n = zeros(size(t));
            n(1,:) = -t(2,:);
            n(2,:) = +t(1,:);                       % Normal vector
            
            v = zeros(size(t));                     % Compute
            v(1,:) = sig.*t(1,:) + as.*n(1,:);      %  original singularity.
            v(2,:) = sig.*t(2,:) + as.*n(2,:);
            
            y = zeros(size(x));
            y(1,:) = p(1,:) + v(1,:) ./ rho .* r;   % Move in that direction
            y(2,:) = p(2,:) + v(2,:) ./ rho .* r;   %  for a distance of r.
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
        %  Here we are using a principal symbol of 1 on the domain of
        %   the canonical transformation (|rho| >= |sigma|) and NaN
        %   (undefined) elsewhere.
        function ps = PrincipalSymbol(patch,x,xi) %#ok<INUSL>
            sz = size(x);
            
            if isequal(size(x), size(xi)),
                rho = xi(2,:); sigma = xi(1,:);
                ps = ones([sz(2:end) 1], class(x));
                ps(abs(rho) < abs(sigma)) = nan;
            else
                rho = xi(2); sigma = xi(1);
                if abs(rho) >= abs(sigma),
                    ps = ones([sz(2:end) 1], class(x));
                else
                    ps =  nan([sz(2:end) 1], class(x));
                end            
            end
        end
    end
end