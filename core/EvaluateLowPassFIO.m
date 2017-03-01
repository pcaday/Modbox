function [AfLP, AfLPpatches] = EvaluateLowPassFIO(ApLP, f, outGrid, valueClass, zeroPadding)

global GlobalEvalLPVerbosity
global GlobalDFTRoutine
global GlobalIDFTRoutine

% If ApLP is an array of patches, process each of them in turn.
if numel(ApLP) > 1
    AfLPpatches = cell(size(ApLP));
    AfLP = Function.Zeros(outGrid, valueClass);
    for i = 1:length(ApLP)
        thisAfLP = EvaluateLowPassFIO(ApLP(i), f, outGrid, valueClass, zeroPadding);
        AfLPpatches{i} = thisAfLP;
        AfLP.f = AfLP.f + thisAfLP.f;
    end
    return;
end

% Get some objects we need
localPatch = ApLP.pPatch.localPatch;
globalPatch = localPatch.globalPatch;

% First step: transform f by the diffeomorphism, yielding ft.
inGridT = ApLP.pPatch.inGrid;

xt = inGridT(:);
x = localPatch.diffeo.InverseTransform(xt);
ft = f.Pullback(inGridT, x);

% Second step: Compute the low-pass frequency coefficients
%  by minimizing the 2-norm of the error between the low-passed f
%  and a weighted sum of plane waves with the given frequencies.

% Lowpass f
highPass = ApLP.wedgeCreator.Sum;
ftLPHat = GlobalDFTRoutine(ft, zeroPadding);
ftLPHat.f = ftLPHat.f .* (1 - highPass);
ftLP = GlobalIDFTRoutine(ftLPHat, inGridT);

% Calculate each plane wave and store in matrix PW.
nFreqs = ApLP.nFreqs;
nAngles = ApLP.nAngles;
nRadii = ApLP.nRadii;

x = inGridT(:);
PW = zeros(inGridT.NPoints, nFreqs, valueClass);
for i = 1:nFreqs
    xi = ApLP.freqs(:,i);
    PW(:,i) = exp(1i*x.'*xi);
end

% Find the closest least-squares approximation
coeff = PW\ftLP.f(:);

% If f is real, make sure the Fourier coefficients are symmetric.
if isreal(f.f) || ~any(imag(f.f(:)))
    coeff_flip = conj([coeff(end/2+1:end); coeff(1:end/2)]);
    coeff = 0.5 * (coeff+coeff_flip);
end

if GlobalEvalLPVerbosity > 1,
    recon = PW*coeff;
    recon = Function.WithValues(inGridT, reshape(recon, size(inGridT)));
    err = recon.copy;
    err.f = ftLP.f - recon.f;
    
    subplot 221
    ft.plot
    normtitle('Original after SRD applied', ft)
    colorbar
    
    subplot 222
    ftLP.plot
    normtitle('With lowpass applied', ftLP)
    colorbar
    
    subplot 224
    recon.plot
    normtitle('Reconstructed from low-frequency plane waves', recon);
    colorbar
    
    subplot 223
    err.plot
    normtitle('Error in reconstruction', err);
    colorbar
    
    pause
end

coeff = reshape(coeff, nRadii, nAngles);
PW = reshape(PW, [], nRadii, nAngles);

% Third step: add the weighted linear combination of plane waves
%  for each group of frequencies, and push them forward.
degree = globalPatch.degree;
%wrapFlag = globalPatch.parent.wrapFlag;

AfLP = Function.Zeros(outGrid, valueClass);
AfLP.f = complex(AfLP.f);           % Make sure we get the complex part.

for j = 1:nAngles
    LC = Function.Zeros(inGridT, valueClass);
    for k = 1:nRadii
        xi = ApLP.freqs(:,k,j);
        thisPW = coeff(k,j) * PW(:,k,j) .* (norm(xi) ^ degree);
        LC.f = LC.f + reshape(thisPW, size(LC.f));
    end
    LC.f = LC.f .* ApLP.amps{j};
    
    Pushforward(AfLP, LC, ApLP.warps{j}); % , wrapFlag);
    
    if GlobalEvalLPVerbosity > 2,
        subplot 121
        LC.plot
        normtitle(sprintf('Weighted plane waves for angle %d', j), LC);
        
        subplot 122
        AfLP.plot
        normtitle('Accumulated output', AfLP);
        
        drawnow
    end
end

AfLPpatches = {AfLP};

end