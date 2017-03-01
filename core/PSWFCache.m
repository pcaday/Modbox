classdef PSWFCache < handle & matlab.mixin.Copyable
    properties
        cData = {}          % 1 x N cell array holding the cache data.
                            %  cData{n} holds data for dimension n.
                            %  cData{n} is an array of structs:
                            %      cData{n}(i).c        -- c value
                            %      cData{n}(i).lambdas  -- array of
                            %                              eigenvalues
                            %      cData{n}(i).psis     -- cell array of
                            %                              function handles
                            %                              for generalized
                            %                              PSWFs.
                            
        ratio = 2           % Ratio between successive c values in cache
                            %   in a given dimension.
        
        pswfTol = 1e-8      % Tolerance value for PSWFs. Only generalized
                            %  PSWFs with eigenvalues at least this big
                            %  are stored (see GeneralizedPSWF.m)
                            
        pswfK = 20          % K value for PSWFs. This is the number of
                            %  basis functions for the radial part of the
                            %  PSWFs (see PSWFRadial.m)
        
        useBessel = true    % Use Bessel functions instead of Jacobi
                            %  polynomials to numerically evaluate the
                            %  PSWFs. Slower, but more accurate and
                            %  the PSWFs work for x > 1.
        
        samples = 500       % # of samples for radial PSWF functions
                            %  (see PSWFRadial.m)
    end
    
    methods
        function ClearCache(cache)
            % Clear the cache
            cache.cData = {};
        end
        
        function [roundedC, psis, lambdas] = GetPSWFs(cache, n, c)
            % Retrieve a PSWF from the cache.
            %
            %  [roundedC, psis, lambdas] = cache.GetPSWFs(n, c);
            %
            % Inputs
            %       n = # of dimensions
            %       c = (minimum) bandwidth
            %
            % Outputs
            %  roundedC = actual bandwidth to use (from cache)
            %      psis = cell array of PSWFs
            %               as function handles.
            %   lambdas = array of eigenvalues corresponding to the PSWFs
            %               w.r.t. to the operator
            %                    Af(y) = \int e^(ix.y) f(x) dx
            
            % First, round c up to the nearest power (positive or negative)
            %  of the cache ratio.
            roundedC = cache.RoundC(c);
            
            % Check if we have cached PSWFs and eigenvalues for this
            %  n and (rounded) c.
            if length(cache.cData) >= n && ~isempty(cache.cData{n})
                existingC = [cache.cData{n}.c];
                cDex = find(existingC == roundedC);
                if ~isempty(cDex)
                    % Yes we do, return them.
                    lambdas = cache.cData{n}(cDex).lambdas;
                    psis = cache.cData{n}(cDex).psis;
                    return
                end
            end
            
            % If not, generate them and add them to the cache.
            [psis, lambdas] = GeneralizedPSWFs(roundedC, n, ...
                cache.pswfTol, cache.pswfK, cache.useBessel, cache.samples);
            
            % Create the struct array for this dimension if necessary.
            %  (commented out -- let MATLAB do this automatically and then
            %    we don't have to worry about problems with "dissimilar
            %    structs.")
            %if length(cache) < n || ~isstruct(cache.cData{n})
            %    cache.cData{n} = struct([]);
            %end
            
            if length(cache.cData) < n || ~isstruct(cache.cData{n})
                newEntry = 1;
            else
                newEntry = length(cache.cData{n}) + 1;
            end
            
            % Add the entry
            cacheStruct.c = roundedC;
            cacheStruct.lambdas = lambdas;
            cacheStruct.psis = psis;
            
            cache.cData{n}(newEntry) = cacheStruct;
        end
        
        
        
        function roundedC = RoundC(cache, c)
            % First, round c up to the nearest power (positive or negative)
            %  of the cache ratio.
            assert(all(isreal(c) & c > 0), ...
                'PSWFCache:ToCachedCValue:badC',...
                'c must be positive');
            r = cache.ratio;
            roundedC = r .^ ceil(log(c) ./ log(r));
        end
    end
    
end