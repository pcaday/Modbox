% Apply an FIO to a test function, and display the results.
%
%   Au = TestApplyFIO(args)
%   [Au,vars] = TestApplyFIO(args)
%
% The arguments are specified as fields in a structure; see ApplyFIO for
%  a description of the arguments it processes. 
%
% TestApplyFIO also has its own (optional) argument:
%        doPieces: logical value; true displays the contributions to the
%                   result from each coordinate choice separately.

function [Au,vars] = TestApplyFIO(args)

if nargin == 0,
    args.dummy = [];
end

%%% Call ApplyFIO to do the work...
[Au,vars] = ApplyFIO(args);
%%% and ReferenceComparisonPlot to display it.
ReferenceComparisonPlot(vars);

end
