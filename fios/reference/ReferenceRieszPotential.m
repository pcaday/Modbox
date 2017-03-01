function fR = ReferenceRieszPotential(f, alpha)
    fhat = f.DFT();
    xi = fhat.grid.AllPoints;
    xiAbs = sqrt(coord(xi,1).^2 + coord(xi,2).^2);
    fhat.f = fhat.f .* (2*pi*xiAbs).^(-alpha);
    
    % Zero out the zero-frequency component if it would be infinite.
    if alpha > 0,
        subs = DFTCenterSubs(fhat.grid);
        fhat.f(subs{:}) = 0;
    end
    
    fR = fhat.IDFT(f.grid);
end