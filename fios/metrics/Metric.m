classdef Metric < handle
    properties
        dims                    % # of dimensions
        g                       % The metric, an n x n cell array of
                                %  function handles.
        gInverse                % Inverse of the metric, n x n cell array
                                %  of function handles.
        dg                      % First derivatives of the metric,
                                %  n x n x n cell array of function handles
                                %
                                %  dg{k,i,j} = d_k(g_ij)
                                %
                                % (derivative of g_ij w.r.t. x_k)
        Rm                      % Function handle to compute full curvature
                                %  tensor.
    end
    
    methods
        % Compute a Christoffel symbol.
        %
        %   values = metric.Christoffel(i,j,k,x)
        %
        % computes \Gamma_{ij}^k
        %
        %  i,j,k: indices
        %      x: n x anything matrix of points (first dimension =
        %                coordinate #)
        function v = Christoffel(metric,i,j,k,x)
            gInverse_ = metric.gInverse; % Use underscores for local copies
            dg_ = metric.dg;
            v = 0;
            for l = 1:metric.dims
                v = v + gInverse_{k,l}(x) .* ...
                    (dg_{i,j,l}(x)+dg_{j,i,l}(x)...
                    -dg_{l,i,j}(x));
            end
            v = v / 2;
        end
        
        function v = ScalarCurvature(metric,x)
            gInverse_ = metric.gInverse;
            n = metric.dims;
            Rm_ = metric.Rm;
            v = 0;
            % Trace over the inside and outside pair of indices
            %  of Rm.
            for i = 1:n
                for j = 1:n
                    for k = 1:n
                        for l = 1:n
                            v = v + gInverse_{i,l}(x)...
                                .* gInverse_{j,k}(x)...
                                .* Rm_(i,j,k,l,x);
                        end
                    end
                end
            end
        end
        
        function v = Ricci(metric,i,j,x)
            gInverse_ = metric.gInverse;
            n = metric.dims;
            Rm_ = metric.Rm;
            v = 0;
            % Trace over the outside pair of indices
            for k = 1:n
                for l = 1:n
                    v = v + gInverse_{k,l}(x)...
                        .* Rm_(k,i,j,l,x);
                end
            end
        end
        
        % Compute a pullback metric, given:
        %
        %        phi: map to pull back by (function handle)
        %    dphi{i}: cell array of function handles of 1st derivatives
        % d2phi{i,j}: cell array of function handles of 2nd derivatives
        pm = Pullback(m, phi, dphi, d2phi)
        
        
        % Calculate lengths of the vector(s) v, based at point(s) x.
        function l = Length(m, x, v)
            assert(isequal(size(x), size(v)));
            
            n = m.dims;
            xr = reshape(x, n, []);
            l = zeros(size(x(1,:)));
            
            for i = 1:n
                for j = 1:i-1
                    l = l + 2 * v(i,:) .* v(j,:) .* m.g{i,j}(xr).';
                end
                l = l + v(i,:).^2 .* m.g{i,i}(xr).';
            end
            
            l = sqrt(l);
            l = preshape(l,x);
        end
        
        % Calculate inner product of the vector(s) v, w based at point(s) x
        function p = InnerProduct(m, x, v, w)
            assert(all(size(x) == size(v)));
            
            n = m.dims;
            xr = reshape(x, n, []);
            p = zeros(size(x(1,:)));
            
            for i = 1:m.dims
                for j = 1:(i-1)
                    p = p + (v(i,:).*w(j,:) + v(j,:).*w(i,:))...
                            .* m.g{i,j}(xr).';
                end
                p = p + v(i,:).*w(i,:) .* m.g{i,i}(xr).';
            end
            
            p = preshape(p,x);
        end
        
        
        % Make the vector(s) v unit vectors
        function w = UnitVector(m, x, v)
            l = m.Length(x,v);
            % Divide each component by the associated length.
            w = v ./ frep(l, m.dims);
        end
        
        % Compute perpendicular vector(s) w to the given vector(s) v
        %  w.r.t. to the metric, at given point(s) (2D only).
        function w = Perpendicular(m, x, v)
            n = m.dims;
            assert(n == 2, 'Metric:dimNot2',...
                'The manifold must be 2-dimensional');
            assert(all(size(x) == size(v)));

            x = reshape(x, n, []);
            g11 = m.g{1,1}(x);
            g12 = m.g{1,2}(x);
            g22 = m.g{2,2}(x);
            
            d = (g11.*g22 - g12.*g12).^(-1/2);
            w = zeros(size(v));
            w(1,:) = (-g12.*v(1,:).' - g22.*v(2,:).') .* d;
            w(2,:) = (g11.*v(1,:).' + g12.*v(2,:).') .* d;
        end
        
        
        % Compute unit vector(s) with given angle(s) from the vector d/dx1
        %  w.r.t. to the metric, at given point(s) (2D only).
        %
        function v = UnitVectorOfAngle(m, x, alpha)
            assert(m.dims == 2, 'Metric:dimNot2',...
                'The manifold must be 2-dimensional');
            assert(all(size(x) == size(alpha)));
            
            g11 = m.g{1,1}(x);
            g12 = m.g{1,2}(x);
            detG = g11 .* m.g{2,2}(x) - g12 .* g12;
            
            v = frep(cos(alpha),2) .* fcat(g11^(-1/2), zeros(size(alpha)))...
                + frep(sin(alpha) .* (g11.*detG)^{-1/2},2) .* fcat(-g12, g11);
        end
        
        % Compute the angle(s) between the vector(s) v and w,
        %  based at point(s) x (2D only).
        function v = AngleBetween(m, x, v, w)
            ip = m.InnerProduct(x,v,w) ./ (m.Length(x,v) .* m.Length(x,w));
            % Take cross product to get which side of v that
            %  w is on. (Euclidean is good enough for our purposes.)
            side = (-v(2,:).*w(1,:) + v(1,:).*w(2,:)) < 0;
            side = reshape(side, size(ip));
            v = acos(ip) .* (1-2*side) + (2*pi.*side);
        end
        
        % Rotate the vector(s) v with base point(s) x by angle(s) alpha
        %  (2D only).
        function w = Rotate(m, x, v, alpha)
            n = m.dims;
            assert(n == 2, 'Metric:dimNot2',...
                'The manifold must be 2-dimensional');
            w = frep(cos(alpha),n) .* v + ...
                frep(sin(alpha),n) .* m.Perpendicular(x,v);
        end
        
        
        % Determinant of the metric tensor entries in local coordinates
        function det = Det(m, x)
            assert(m.dims == 2, 'Metric:dimNot2',...
                'Det only implemented for 2D manifolds.');
            
            det = m.g{1,1}(x) .* m.g{2,2}(x) - m.g{1,2}(x) .^ 2;
        end
    end
end