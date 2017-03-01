classdef EuclideanConformalMetric < Metric
    properties
        lambda          % Scaling function (function handle).
                        %  The metric will be exp(2*lambda) times the flat
                        %  metric.
                        %
                        %  lambda = @(x) vals
                        %
                        %   x is an n x anything x ... x anything array.
                        %   vals should be an array having the same size
                        %    as x, minus its first dimension

        dlambda         % First derivatives of scaling function as
                        %  a length-n cell array of function handles
                        %
                        %  dlambda{i} = @(x) vals
                        %    Here,      i = variable to differentiate in
                        %               x = points to evaluate at.
                        %
                        %   x is an n x anything x ... x anything array.
                        %    vals should be an array having the same size
                        %    as x, minus its first dimension

        d2lambda        % Second derivatives of scaling function as
                        %  an n x n cell array of function handles.
                        %
                        %  dlambda{i,j} = @(x) vals
                        %    Here,    i,j = variables to differentiate in
                        %               x = points to evaluate at.
                        %
                        %   x is an n x anything x ... x anything array.
                        %    vals should be an array having the same size
                        %    as x, minus its first dimension

    end
    
    
    methods
        function m = EuclideanConformalMetric(dims, lambda, dlambda, d2lambda)
            % Constructor:
            %
            %    m = EuclideanConformalMetric(dims, lambda, dlambda, d2lambda)
            %
            % Creates a metric conformal to the Euclidean metric:
            %    m = exp(2*lambda)*e,
            % where e is the Euclidean metric.
            %
            %     dims: # of dimensions
            %   lambda: function handle representing lambda (L).
            %  dlambda: First derivatives of lambda:  dlambda{i} = dlambda/dx_i
            % d2lambda: Second derivatives of lambda: d2lambda{i,j} = d^2 lambda/dx_i*dx_j
            %
            %
            %    m = EuclideanConformalMetric(Lstruct)
            % 
            % The arguments can also be given in a struct
            %  (Lstruct) with fields:
            %      dims
            %    lambda
            %   dlambda
            %  d2lambda
            %
            
            % Allow for no arguments (sort of a default constructor)
            
            % If the arguments are passed in a structure...
            haveStruct = false;
            if nargin > 0,
                if isstruct(dims),
                    Lstruct = dims;
                    dims =      Lstruct.dims;
                    lambda =    Lstruct.lambda;
                    dlambda =   Lstruct.dlambda;
                    d2lambda =  Lstruct.d2lambda;
                    haveStruct = true;
                end
                
                m.dims = dims;
            end
            
            if nargin > 1 || haveStruct,
                m.lambda = lambda;
                m.dlambda = dlambda;
                m.d2lambda = d2lambda;
                
                m.g = cell(dims);
                m.gInverse = cell(dims);
                
                zero = @(x) 0;
                % Generate the metric
                for i = 1:dims
                    for j = 1:dims
                        if i ~= j,
                            m.g{i,j} = zero;
                            m.gInverse{i,j} = zero;
                        else
                            m.g{i,j}        = @(x) exp(2*lambda(x));
                            m.gInverse{i,j} = @(x) exp(-2*lambda(x));
                        end
                    end
                end
                
                % Generate first derivatives of the metric
                m.dg = cell([dims dims dims]);
                for i = 1:dims
                    for j = 1:dims
                        for k = 1:dims
                            if i ~= j,
                                m.dg{k,i,j} = zero;
                            else
                                m.dg{k,i,j} = @(x) 2*dlambda{k}(x).*exp(2*lambda(x));
                            end
                        end
                    end
                end
                
                % Assign the curvature tensor finder
                m.Rm = @(varargin) m.CurvatureTensor(varargin{:});
            end
        end

    end
    
    methods (Static)
        function m = DiscreteMultiplier(c)
        % Construct a conformally flat metric g = c^2*e, where
        %  c is given by its function values.
        %
        % m = DiscreteMultiplier(c)
        %
        % Input 
        %    c: Function object
        % Output  
        %    m: EuclideanConformalMetric object representing the
        %        metric c^2*e.
            
            % Dimension count
            dims = c.grid.dims;
            
            % Create (empty) metric object
            m = EuclideanConformalMetric(dims);
            
            % Create c^2
            c2 = c.copy();
            c2.f = c.f .* c.f;
            
            % Create c^(-2)
            c2Inverse = c2.copy();
            c2Inverse.f = 1 ./ c2.f;
            
            % Create lambda = ln(c)
            lamb = c.copy();
            lamb.f = log(c.f);
            m.lambda = @(x) lamb.SampleAt(x);
            
            % Find the gradient of c.
            dcs = ScalarGradient(c);
            
            % Create dlambda.  d(lambda)/di = (dc/di) / c.
            m.dlambda = cell(1,dims);
            for i = 1:dims
                dl = dcs{i}.copy();
                dl.f = dcs{i}.f ./ c.f;
                m.dlambda{i} = @(x) dl.SampleAt(x);
            end
            
            % Avoid waste; create function handles!
            zero_fh =      @(x) 0;
            c2_fh =        @(x) c.SampleAt(x);
            c2Inverse_fh = @(x) c2Inverse.SampleAt(x);
            
            % Find the gradient of c^2
            %  d(c^2)/di = 2*c*(dc/di).
            dc2s_fh = cell(1,dims);
            for i = 1:dims
                dc2s = dcs{i}.copy();
                dc2s.f = 2 * dcs{i}.f .* c.f;
                dc2s_fh{i} = @(x) dc2s.SampleAt(x);
            end
            
            % Generate the metric
            for i = 1:dims
                for j = 1:dims
                    if i ~= j,
                        m.g{i,j} = zero_fh;
                        m.gInverse{i,j} = zero_fh;
                    else
                        m.g{i,j}        = c2_fh;
                        m.gInverse{i,j} = c2Inverse_fh;
                    end
                end
            end
            
            % Generate first derivatives of the metric
            m.dg = cell([dims dims dims]);
            for i = 1:dims
                for j = 1:dims
                    for k = 1:dims
                        if i ~= j,
                            m.dg{k,i,j} = zero_fh;
                        else
                            m.dg{k,i,j} = dc2s_fh{k};
                        end
                    end
                end
            end
            
            % Curvature tensor not implemented...
        end
    end
        
    methods
        function vals = CurvatureTensor(m, i, j, k, l, x)
            d2l = m.d2lambda;
            dl = m.dlambda;
            vals = 0;
            if j == l, vals = vals + d2l{i,k}(x) - dl{i}(x).*dl{k}(x); end
            if j == k, vals = vals - d2l{i,l}(x) + dl{i}(x).*dl{l}(x); end
            if i == l, vals = vals - d2l{j,k}(x) + dl{j}(x).*dl{k}(x); end
            if i == k, vals = vals + d2l{j,l}(x) - dl{j}(x).*dl{l}(x); end
            
            if (i == l && j == k) || (i == k && j == l),
                normGradL = 0;
                for p = 1:m.dims
                    normGradL = normGradL + dl{p}(x).^2;
                end
                
                if (i == l && j == k), vals = vals - normGradL; end
                if (i == k && j == l), vals = vals + normGradL; end
            end
            
            vals = vals .* exp(2*m.lambda(x));
        end
        
        function v = Christoffel(metric,i,j,k,x)
            dl = metric.dlambda;
            v = 0;
            
            if j == k, v = v + dl{i}(x); end
            if i == k, v = v + dl{j}(x); end
            if i == j, v = v - dl{k}(x); end
        end
    end
end