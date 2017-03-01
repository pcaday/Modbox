% Check that a list of variables matches up with their corresponding types
%  and raise an assertion if there is a mismatch.
%
% CheckTypes(variable, type, argNum)
%   variable: arbitrary value
%       type: string
%               Types can be
%                   - class names ('double', 'uint64', 'cell', ...)
%                   - 'integer' 
%                   - 'scalar ' + any type
%                             e.g. 'scalar double' checks for a scalar
%                                 double value.
%                   - 'scalar' - synonym for 'scalar numeric'
%     argNum: (optional) argument #, used for error messages (default: 1)
%
%  CheckTypes(variables, types, argNums)
%    variables: cell array of values
%        types: cell array of strings (corresponding types)
%      argNums: cell array of argument numbers, used for error messages
%                 default: {1,2,...}
%
function CheckTypes(variables, types, argumentNums)
    if ~iscell(variables)
        variables = {variables};
    end
    if ~iscell(types)
        types = {types};
    end
    assert(all(size(variables) == size(types)),...
        'CheckTypes:inputSizeMismatch',...
        'Inputs have differing lengths.');
    n = numel(variables);
    if nargin < 3
        argumentNums = num2cell(1:n);       % Make cell array {1 2 ... }
    elseif ~iscell(argumentNums)
        argumentNums = {argumentNums};
    end
    
    for i = 1:n
        variable = variables{i};
        type = types{i};
        if strcmp(type, 'scalar')
            type = 'scalar numeric';
        end
        if strncmp(type, 'scalar ', 7)      % 'scalar ' has 7 characters
            assert(isscalar(variable), 'CheckTypes:notScalar', ...
                'Argument %d is not a scalar', argumentNums{i});
            type(1:7) = [];                 % delete 'scalar ' from start
        end
        switch type
            case 'integer'
                assert(isnumeric(variable), 'CheckTypes:notInteger',...
                    'Argument %d is not an integer', argumentNums{i});
                assert(all(fix(variable) == variable), 'CheckTypes:notInteger',...
                    'Argument %d is not an integer', argumentNums{i});
            otherwise
                assert(isa(variable, type), 'CheckTypes:notOfClass',...
                    'Argument %d is not a %s', argumentNums{i}, type);
        end
    end
end