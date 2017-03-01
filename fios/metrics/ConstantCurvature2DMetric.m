classdef ConstantCurvature2DMetric < EuclideanConformalMetric
    properties
        kappa
    end
    
    
    methods
        function m = ConstantCurvature2DMetric(kappa)
            if kappa == 0,
                zero = @(varargin) 0;
                lambda = zero;
                dlambda = {zero zero};
                d2lambda = {zero zero; zero zero};
            else
                lambda = @(x) real(preshape(log(2/kappa) - log(sum(x.^2,1) + 1/kappa), x));
                dlambda = {...
                    @(x) preshape(-2*x(1,:) ./ (sum(x(:,:).^2,1) + 1/kappa), x)...
                    @(x) preshape(-2*x(2,:) ./ (sum(x(:,:).^2,1) + 1/kappa), x)};
                d2lambda = {...
                    @(x) preshape(2*(x(1,:).^2-x(2,:).^2-1/kappa) ./ (sum(x(:,:).^2,1)+1/kappa).^2, x) ...
                    @(x) preshape(4*x(1,:).*x(2,:) ./ (sum(x(:,:).^2,1)+1/kappa).^2, x); ...
                    @(x) preshape(4*x(1,:).*x(2,:) ./ (sum(x(:,:).^2,1)+1/kappa).^2, x) ...
                    @(x) preshape(2*(-x(1,:).^2+x(2,:).^2-1/kappa) ./ (sum(x(:,:).^2,1)+1/kappa).^2, x)};
            end
            
            m = m@EuclideanConformalMetric(2, lambda, dlambda, d2lambda);
            m.kappa = kappa;
        end
    end
    
end