classdef GeodesicXRay2DPatch < FIOPatch
    properties
        metric          % Metric object
        dM              % Curve info structure for dM
        direction       % Either +1 or -1: determines direction the
                        %  geodesic is traced.
        
        dt = 0.05       % Time step for numerically solving the geodesic ODE.
        tMax = 5        % Maximum time to trace geodesics before giving up
                        %  if they do not hit the boundary.
    end
    
    methods
        % Create a GeodesicXRay2DPatch, an FIOPatch for the
        %  2D geodesic X-ray transform.
        function patch = GeodesicXRay2DPatch(metric, dMCurveInfo, direction)
            patch.dims = 2;
            patch.degree = 0;
            
            patch.metric = metric;
            patch.dM = dMCurveInfo;
            patch.direction = direction;
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
            global GlobalGeoXRayVerbosity
            
            
            % Local versions of properties
            dt_ = patch.dt;
            tMax_ = patch.tMax;
            m = patch.metric;
            insideTest_ = patch.dM.lrTestF;
            
            % Expand xi if necessary (expand to the size of x)
            xi = repmattosize(xi, size(x));
            
            % xc, vc: current x and v (velocity) values.
            vc = m.UnitVector(x, xi);
            vc = m.Perpendicular(x, vc);
            vc = vc * patch.direction;
            xc = x;
            
            % Initialize out, a logical array which keeps track of which
            %  geodesics are still contained completely in M.
            out = ~insideTest_(xc);

            % Initialize xf, vf, the final x-coordinates and velocities of
            %  the geodesic (lying somewhere on dM).
            xf = nan(size(x), class(x));
            vf = nan(size(x), class(x));
            
            % Loop over time steps
            for ts = 1:ceil(tMax_/dt_)
                % Solve the geodesic equation.
                
                % Propagate x half a step forward
                xm = xc + vc * dt_ * 0.5;
                
                % Calculate dv/dt, using the geodesic equation ODE.
                dvdtm = zeros(size(vc));
                for k = 1:2
                    for i = 1:2
                        for j = 1:2
                            dvdtm(k,:) = dvdtm(k,:) - vc(i,:) .* vc(j,:) .* ...
                                m.Christoffel(i,j,k,xm(:,:)).';
                        end
                    end
                end
                
                % Propagate v a half step, and a full step forward
                vm = vc + dvdtm * dt_ * 0.5;
                vc = vc + dvdtm * dt_;
                
                % Propagate x a full step forward with the half-step v
                xn = xc + vm * dt_;
                
                % Update which geodesics have escaped M; record their
                %  escape locations in xf.
                newOut = ~insideTest_(xn) & ~out;
                xf(:,newOut) = 0.5*(xc(:,newOut) + xn(:,newOut));
                vf(:,newOut) = vm(:,newOut);
                out = out | newOut;
                
                xc = xn;
                                
                % Check for geodesics leaving M.
                % If every geodesic has escaped, we're done!
                if all(out), break; end
            end
            
            % At this point, xf and vf contain the final positions and
            %  velocities of all the geodesics starting inside the domain
            %  which escaped. For other starting points, xf and vf are NaN.
            
            % Get s values for the exit locations
            y(1,:,:) = patch.dM.changeXToS(xf);
            
            % Get angles (alpha values) at each escape location.
            % Here we're assuming the curve is oriented CCW.
            gp = patch.dM.tangentF(xf);
            y(2,:,:) = m.AngleBetween(xf, vf, gp);
            
            if GlobalGeoXRayVerbosity > 1,
                % Informational output
                fprintf('Geodesic tracing: breaking after time %f\n', ts*dt_);

                l = m.Length(xc,vc);
                l(isnan(l)) = [];
                fprintf('Geodesic tracing: min/max lengths %f/%f\n', min(l(:)), max(l(:)));
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
        % [dxy, dxiy] = DCanonicalTransformationY(patch,x,xi) %#ok<INUSD>

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