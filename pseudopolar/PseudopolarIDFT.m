function f = PseudopolarIDFT(fhat, grid)
    % First, rearrange fhat so it is in the order expected by IPPFFT.
    
    N = size(fhat.f, 1) - 1;
    M = size(fhat.f, 2) / 4;
    
    % Un-rotate by 45 degrees
    vals = circshift(fhat.f, [0 -ceil((M-1)/2)]);
    % Flip to swap x and y.
    vals = fliplr(vals);
    
    % We need to average the zero frequency components (r=0) from
    %  each pair of opposing directions; same for the maximum frequency
    %  components for each pair. Prepare by multiplying all these components
    %  by 0.5.
    
    vals([1 end], :) = vals([1 end], :) * 0.5;
    
    vals = [vals(end, 1:M) + vals(end, 2*M+1:3*M)...        % max frequency - vertical
            vals(end, M+1:2*M) + vals(end, 3*M+1:4*M);...   % max frequency - horizontal
            vals(end-1:-1:2, 2*M+1:3*M)...                  % bottom
            vals(end-1:-1:2, M+1:2*M);...                   % left
            vals(1, 1:M) + vals(1, 2*M+1:3*M)...            % zero frequency - vertical
            vals(1, M+1:2*M) + vals(1, 3*M+1:4*M);...       % zero frequency - horizontal
            vals(2:end-1, 1:M) ...                          % top
            vals(2:end-1, 3*M+1:4*M)];                      % right
    
    % Now call IPPFFT to do the hard work.
    f_rc = IPPFFT(vals);
    f = Function.WithValues(grid, f_rc);
end