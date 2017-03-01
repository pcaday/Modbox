% Enable or disable pseudopolar DFT
%
%   UsePPDFT(true)
%     Enable pseudopolar DFT.
%
%   UsePPDFT(false)
%     Disable pseudopolar DFT and use standard Cartesian DFT.
%
function UsePPDFT(enable)
global GlobalUsePseudopolarFFT
global GlobalDFTRoutine
global GlobalDFTGridRoutine
global GlobalIDFTRoutine

GlobalUsePseudopolarFFT = enable;
if enable,
    GlobalDFTRoutine = @PseudopolarDFT;
    GlobalDFTGridRoutine = @(grid, varargin) PseudopolarGrid(grid, varargin{:});
    GlobalIDFTRoutine = @PseudopolarIDFT;
else
    GlobalDFTRoutine = @(f, varargin) f.DFT(varargin{:});
    GlobalDFTGridRoutine = @(grid, varargin) grid.DFTGrid(varargin{:});
    GlobalIDFTRoutine = @(f, grid) f.IDFT(grid);
end

end