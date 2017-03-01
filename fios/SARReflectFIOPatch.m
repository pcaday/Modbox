classdef SARReflectFIOPatch < FIOPatch
    properties
        curveInfo
        sign
        maskInside
        
        % Parameters for the intersection point-finding algorithm
        hopMax = 8
        hopMin = 1e-5
        divisions = 4
    end
    
    methods
        function patch = SARReflectFIOPatch(curveInfo,sign,maskInside,parent)
            patch.dims = 2;
            patch.degree = 0;
            patch.parent = parent;
            patch.maskInside = maskInside;
            
            patch.curveInfo = curveInfo;
            patch.sign = sign;
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
        function y = CanonicalTransformationY(patch,x,xi)
            global GlobalSARReflectVerbosity
            
            ip = patch.FindIntersections(x,xi);
            
            % Get tangents to the curve at the intersection points
            %  and normalize them.
            t = patch.curveInfo.tangentF(ip);
            t = bsxfun(@rdivide, t, hypot(t(1,:), t(2,:)));
            
            % Reflect the original x's across the tangent lines
            %  to get the y's.
            d = 2 * (- t(2,:).*(ip(1,:)-x(1,:))...
                     + t(1,:).*(ip(2,:)-x(2,:)));
            y = zeros(size(x));
            y(1,:) = x(1,:) + d .* -t(2,:);
            y(2,:) = x(2,:) + d .* t(1,:);
            
            if GlobalSARReflectVerbosity == 3,
                nonnan = find(~isnan(ip(1,:)));
                if ~isempty(nonnan),
                    wh = nonnan(ceil(rand(1) * length(nonnan)));
                    scatter(x(1,wh), x(2,wh), 'rx');
                    hold on
                    plot(x(1,wh)+[0 xi(1)*0.3], x(2,wh)+[0 xi(2)*0.3], 'r-');
                    scatter(ip(1,wh), ip(2,wh), 'bo');
                    plot(ip(1,wh)+t(1,wh)*[-0.3 0.3], ip(2,wh)+t(2,wh)*[-0.3 0.3], 'b-');
                    scatter(y(1,wh), y(2,wh), 'gx');
                    hold off
                    axis equal
                    title('Example reflection');
                    pause
                end
            end
            
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
%        function [dxy, dxiy] = DCanonicalTransformationY(patch,x,xi) %#ok<INUSD>


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
            ps = patch.PrincipalSymbol@FIOPatch(x,xi);
            
            % Find the Maslov factor (quick & dirty)
            %  Ideally, we would reuse the information from the
            %  previous canonical transformation calculations.
            
            % Flatten and normalize xi.
            xi = reshape(xi, size(xi,1), []);
            xi = bsxfun(@rdivide, xi, hypot(xi(1,:), xi(2,:)));
            
            % Get the intersection points.
            ip = patch.FindIntersections(x,xi);
            
            % Get the curvatures at each.
            c = patch.curveInfo.curvatureF(ip);
            
            % Get the dot product xi . gamma'(s)
            gp = patch.curveInfo.tangentF(ip);
            dp = xi(1,:).*-gp(2,:) + xi(2,:).*gp(1,:);
            dp = fsqueeze(dp);
            
            % The signature in the Maslov factor is
            %   sgn(xi . (c(s) gamma'(s)))
            sg = sign(dp.*c);               %#ok<CPROP>
            ps(:) = ps(:) .* exp(pi*1i/4*sg);
        end
        
        
        
        function ip = FindIntersections(patch, x, xi)
            global GlobalSARReflectVerbosity
            
            xi = xi * patch.sign;
            defF = patch.curveInfo.definingF;
            divisions_ = patch.divisions;
            
            h = patch.hopMax;
            
            % Flatten and normalize xi.
            xi = reshape(xi, size(xi,1), []);
            xi = bsxfun(@rdivide, xi, hypot(xi(1,:), xi(2,:)));
            
            % Flatten x and initialize xLB and xUB, which will hold
            %  the coordinates at each end of our search intervals.
            xLB = reshape(x, size(x,1), []);
            xUB = bsxfun(@plus, xLB, xi * h);
            
            % Create vector of divisions+1 equally spaced points
            %  between 0 and 1.
            fracs = (0:divisions_) / divisions_;
            fracs = reshape(fracs, [1 1 length(fracs)]);
            
            % Record number of x values
            nx = size(xLB,2);
            
            
            while h > patch.hopMin,
                % Debugging output: show search intervals.
                if GlobalSARReflectVerbosity == 2,
                    % Select a subset of the whole grid to display
                    i = size(x,2);
                    j = size(x,3);
                    ii = round(linspace(1,i,15));
                    jj = round(linspace(1,j,15));
                    [ii,jj] = ndgrid(ii,jj);
                    ll = sub2ind([size(x,2) size(x,3)], ii(:), jj(:));
                    plot([xLB(1,ll); xUB(1,ll)], [xLB(2,ll); xUB(2,ll)]);
                    title('Search intervals');
                    pause(0.5)
                end
                
                % Compute equally spaced points between xLB and xUB.
                xc = bsxfun(@plus, bsxfun(@times, 1-fracs, xLB), bsxfun(@times, fracs, xUB));
                % Get values of defining function at all these points
                v = defF(xc);
                % Look for sign changes.
                signV = sign(v); %#ok<CPROP>
                signCh = signV(:,2:end) ~= signV(:,1:end-1);
                hasSignCh = any(signCh,2);
                [~,sci] = max(signCh,[],2);
                % If no sign changes, look for the nearest values to zero.
                [~,nearest] = min(abs(v),[],2);
                newLB = sci;
                newUB = sci + 1;
                newLB(~hasSignCh) = max(1,            nearest(~hasSignCh) - 1);
                newUB(~hasSignCh) = min(divisions_+1, nearest(~hasSignCh) + 1);
                
                % Pick out the new upper and lower bounds.
                nxList = (1:nx).';
                lbInds = sub2ind([nx divisions_+1], nxList, newLB);
                ubInds = sub2ind([nx divisions_+1], nxList, newUB);
                xLB = xc(:,lbInds);
                xUB = xc(:,ubInds);
                
                h = h / divisions_;
            end
            
            % Average the final lower and upper bounds to get
            %  the approximate intersection points.
            xLB = 0.5*(xLB+xUB);

            % Debugging output: display found intersection points
            if GlobalSARReflectVerbosity == 2,
                scatter(xLB(1,hasSignCh), xLB(2,hasSignCh), 'bo');
                hold on
                scatter(x(1,~hasSignCh), x(2,~hasSignCh), 'ro');
                arrow(x(:,1), x(:,1)+xi(:,1), 'EdgeColor', 'g');         % from File Exchange
                hold off
                title('Intersection points (blue = found, red = not found)');
                pause
            end
            
            % If we did not find a sign change by the last step
            %  (and therefore didn't find a intersection point), set the
            %  intersection point to NaN.
            xLB(:,~hasSignCh) = NaN;
            
            % For curves with an inside, mask out the outside/inside
            %  according to the maskInside setting.
            if patch.maskInside && patch.curveInfo.hasInside,
                inside = patch.curveInfo.lrTestF(x);
                if patch.maskInside == 1,
                    xLB(:,~inside) = NaN;
                elseif patch.maskInside == -1,
                    xLB(:,inside) = NaN;
                end
            end

            % Rename xLB --> ip to reflect the fact that it now
            %  holds intersection points.
            ip = xLB;     
            
            if GlobalSARReflectVerbosity > 0,
                fprintf('SARReflect: %d/%d intersection points found\n',...
                    nnz(~isnan(ip))/2, numel(ip)/2);
            end
        end
    end
end