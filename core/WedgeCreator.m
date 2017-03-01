
classdef WedgeCreator < handle
    properties
        grid
    end
        
    methods
        function N = Count(wcr) %#ok<*STOUT,*MANU>
        % N = wedgeCreator.Count()
        %
        %   Return the total number of wedges, N.
        %
        end
        
        function f = Wedge(wcr, i) %#ok<*INUSD>
        % f = wedgeCreator.Wedge(i)
        %
        %   Get the values of the frequency cutoff function f for wedge i.
        %   f will be an array of function values, representing the values
        %   of the function on some grid.
        %
        %   The grid will be specified by the creator of the
        %   WedgeCreator and passed into the WedgeCreator's constructor.
        end
        
        function nu = Center(wcr, i)
        % nu = wedgeCreator.Center(i)
        %
        %   Get a unit vector, nu, pointing toward the center of wedge i.
        %
        end
        
        function nus = Centers(wcr)
        % nus = wedgeCreator.Centers
        %
        %   Get an n x k matrix of unit vectors
        %   pointing toward the center of the wedges.
        %
        %   The i'th column is the vector pointing to wedge i.
        %
        end
        
        function mask = Sum(wcr)
        % mask = wedgeCreator.Sum
        %
        %   Returns the sum of all the wedges (a wedge sum? :p)
        %
        end
        
        function boxes = NearestBoxes(wcr, xi)
        % boxes = wedgeCreator.NearestBoxes(xi)
        %
        %   Find the nearest box numbers to the given xis
        %
        %    xi is an n x i1 x ... x ik array
        % boxes is a  m x i1 x ... x ik array of box numbers
        %
        %   m is determined by the WedgeCreator.
        end
        
        function lpCeil = LowPassCeiling(wcr)
        % xiLP = wedgeCreator.LowPassCeiling()
        %   
        %   Return the largest frequencies not completely handled
        %  by this WedgeCreator. In other words, if the sum of the wedges
        %  does not add up to 1 at xi, then
        %   |xi_i| <= |xiLP_i| for all i.
        end
    end
    
    methods
        function nuBasis = NuBasis(wcr, i)
        % nuBasis = wedgeCreator.NuBasis(i)
        %
        %   Get an orthogonal basis for R^n with nu as its first element.
        %
        % The basis is returned as an n x n matrix with the basis vectors
        %  as its columns.
            nu = wcr.Center(i);
            nu = nu(:);
            nuPerp = null(nu');
            nuBasis = [nu nuPerp];
        end
    end
end