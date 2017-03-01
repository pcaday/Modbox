classdef FIO < handle
    properties
        dims            % Dimension of the domain and codomain (n)
        patches         % Array of FIOPatch objects
        wrapFlag        % Wrapping flag:
                        %   0: no wraparound (nonperiodic grid)
                        %   1: wraparound (periodic grid)
                        %   2: wraparound + automatically unwrap
                        %         coordinates
        outLogicalGrid  % Grid object representing the codomain.
                        %  Only the periodicity is needed, and the width
                        %  of the domain in the coordinates that are
                        %  periodic.
    end
end