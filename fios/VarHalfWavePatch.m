classdef VarHalfWavePatch < FIOPatch
    properties
        metric          % Metric object
        t               % Time to propagate geodesics into the future
        
        dt = 0.02       % Time step for numerically solving the geodesic ODE.
    end
    
    methods
        function patch = VarHalfWavePatch(dims, metric, t)
            n = dims;
            
            patch.dims = n;
            patch.degree = 0;
            
            patch.metric = metric;
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
            global GlobalHalfWaveVerbosity
            
            % Subroutine for computing dv/dt for the geodesic equation.
            function dvdt = GetDvdt(x,v)
                dvdt = zeros(size(v));
                if isConformallyFlat2D,
                    % Special case for 2D conformally flat metrics
                    %  (for speed).
                    dl1 = dlambda{1}(x);
                    dl2 = dlambda{2}(x);
                    a = v(1,:).^2 - v(2,:).^2;
                    b = 2*v(1,:).*v(2,:);
                    dl1 = dl1(:);
                    dl2 = dl2(:);
                    a = a(:);
                    b = b(:);
                    dvdt(1,:) = -a.*dl1 - b.*dl2;
                    dvdt(2,:) =  a.*dl2 - b.*dl1;
                else
                    % General case, using Christoffel symbols.
                    for k = 1:n
                        for i = 1:n
                            for j = 1:n
                                dvdt(k,:) = dvdt(k,:) - v(i,:) .* v(j,:) .* ...
                                    m.Christoffel(i,j,k,x(:,:)).';
                            end
                        end
                    end
                end
            end
            
            % Local versions of properties
            dt_ = patch.dt;
            t_ = patch.t;
            m = patch.metric;
            n = patch.dims;
            
            isConformallyFlat2D = (n == 2) && isa(m, 'EuclideanConformalMetric');
            if isConformallyFlat2D,
                dlambda = m.dlambda;
            end
            
            % Adjust dt so it goes evenly into t
            nt = ceil(t_/dt_);
            dt_ = t_/nt;
            
            % Expand xi if necessary (expand to the size of x)
            xi = repmattosize(xi, size(x));
            
            % xc, vc: current x and v (velocity) values.
            vc = m.UnitVector(x, xi);
            xc = x;
            
            if GlobalHalfWaveVerbosity > 0,
                % Informational output
                fprintf('Integrating to time %.2f in %d steps.', t_, nt);
            end

            if GlobalHalfWaveVerbosity > 1 && n == 2,
                % Get sound speed picture and set up axis
                ax = [min(xc(1,:)) max(xc(1,:)) min(xc(2,:)) max(xc(2,:))];
                ax = 2 * ax;
                axis(ax);
                g = Grid([ax(1) ax(3)], [ax(2) ax(4)], [64 64]);
                ss = Function.WithHandle(g, @(x) exp(-m.lambda(x)));
            end

            % Loop over time steps
            for ts = 1:nt
                % Solve the geodesic equation with RK4.
                
                % Stage 1
                q1x = dt_ * vc;
                q1v = dt_ * GetDvdt(xc,vc);
                
                % Stage 2
                vm = vc + 0.5*q1v;
                q2x = dt_ * vm;
                q2v = dt_ * GetDvdt(xc + 0.5*q1x, vm);
                
                % Stage 3
                vm = vc + 0.5*q2v;
                q3x = dt_ * vm;
                q3v = dt_ * GetDvdt(xc + 0.5*q2x, vm);
                
                % Stage 4
                vm = vc + q3v;
                q4x = dt_ * vm;
                q4v = dt_ * GetDvdt(xc + q3x, vm);
                
                % Step
                xc = xc + 1/6 * (q1x + 2*q2x + 2*q3x + q4x);
                vc = vc + 1/6 * (q1v + 2*q2v + 2*q3v + q4v);
                
                if GlobalHalfWaveVerbosity > 1 && n == 2,
                    % Select a subset of the whole grid to display
                    is = size(x,2);
                    js = size(x,3);
                    ii = round(linspace(1,is,10));
                    jj = round(linspace(1,js,10));
                    [ii,jj] = ndgrid(ii,jj);
                    
                    % Plot geodesics.
                    ll = sub2ind([size(x,2) size(x,3)], ii(:), jj(:));
                    
                    % Freeze the axis limits after the first step.
                    cla
                    ss.plot;
                    colormap copper
                    l = 0.15;
                    hold on
                    plot([xc(1,ll); xc(1,ll)+l*vc(1,ll)], [xc(2,ll); xc(2,ll)+l*vc(2,ll)], 'b-');
                    
                    drawnow
                end
                
            end
            
            
            if GlobalHalfWaveVerbosity > 0,
                % Informational output
                l = m.Length(xc,vc);
                l(isnan(l)) = [];
                fprintf(' (min/max |v| %f/%f)\n', min(l(:)), max(l(:)));
            end
            
            y = xc;
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
        function ps = PrincipalSymbol(patch,x,xi)
            ps = patch.PrincipalSymbol@FIOPatch(x,xi);
            ps = ps * exp(-pi*1i/4);
        end
     end
end