function fProp = ReferenceWaveProp(f, t)
    fhat = f.DFT();
    xi = fhat.grid.AllPoints;
    xiAbs = sqrt(coord(xi,1).^2 + coord(xi,2).^2);
    fhat.f = fhat.f .* exp(-1i*t*xiAbs);
    fProp = fhat.IDFT(f.grid);
end
