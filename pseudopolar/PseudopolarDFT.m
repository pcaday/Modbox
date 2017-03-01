function fhat = PseudopolarDFT(f, oversampling)
    if nargin < 2 || isempty(oversampling), oversampling = 1; end
    
    % Note that PPFFT expects 2D data in meshgrid-ordering
    %  (y in the first index; x in the second index).
    % Can either transpose here or flip below.
    vals = PPFFT(f.f, oversampling, oversampling);
    % Move parts of fhat around because I can't remember the way the data
    %  from PPFFT is arranged for my life.
    
    N = size(vals, 1) / 2;
    M = size(vals, 2) / 2;
    vals = [vals([N+1:2*N 1], 1:M)...       % top    (y+)
            vals(N+1:-1:1, M+1:2*M)...      % left   (x-)
            vals(N+1:-1:1, 1:M)...          % bottom (y-)
            vals([N+1:2*N 1], M+1:2*M)];    % right  (x+)
    % Swap x and y?
    vals = fliplr(vals);
    % Rotate 45 degrees so the horizontal direction comes first.
    vals = circshift(vals, [0 ceil((M-1)/2)]);
    
    igrid = PseudopolarGrid(f.grid, oversampling);
    fhat = Function.WithValues(igrid, vals);
end


% Pseudopolar output (with meshgrid ordering)
% Swap x and y for ndgrid ordering.
%
%         angles
%
%           |
%           |
%     Y-    |    X-
%           |
%           |
% ----------+---------- 0      radius    <-- row N+1
%           |
%           |
%     Y+    |    X+
%           |
%           |
