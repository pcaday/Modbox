function [baseArgs, vdArgs] = PVBaseArgs(fio)
    switch (fio)
        case 'radon'
            baseArgs.fio = 'radon';
            baseArgs.n = 128;
            baseArgs.wcrSmoothness = 0.5;
            baseArgs.bendFactor = 3;
            baseArgs.zeroPad = 2;
            baseArgs.GlobalBadDiffeoThreshold = 0.85;
        case 'waveprop'
            baseArgs.fio = 'waveprop';
            baseArgs.n = 64;
            baseArgs.diffeos = 'identity';
            baseArgs.bendFactor = 1;
            baseArgs.GlobalBadDiffeoThreshold = 0.85;
            baseArgs.image = 'dash';
            baseArgs.T = 0.4;
        case 'halfwave1'
            baseArgs.fio = 'halfwave';
            baseArgs.n = 64;
            baseArgs.diffeos = 'identity';
            baseArgs.bendFactor = 1;
            baseArgs.GlobalBadDiffeoThreshold = 0.85;
            baseArgs.image = 'dash';
            baseArgs.T = 0.3;
            baseArgs.periodic = true;
            baseArgs.zeroPad = 1;
            baseArgs.metric = 'soundspeed';
            baseArgs.soundSpeed = 'deep-well';
        case 'halfwave2'
            baseArgs.fio = 'halfwave';
            baseArgs.n = 64;
            baseArgs.diffeos = 'four-plus-identity';
            baseArgs.bendFactor = 1;
            baseArgs.GlobalBadDiffeoThreshold = 0.85;
            baseArgs.image = 'tilted-dash';
            baseArgs.T = 0.6;
            baseArgs.periodic = true;
            baseArgs.zeroPad = 1;
            baseArgs.ratio = 8;
            baseArgs.metric = 'soundspeed';
            baseArgs.soundSpeed = 'deep-well';            
        case 'halfwave-flat'
            baseArgs.fio = 'halfwave';
            baseArgs.n = 64;
            baseArgs.diffeos = 'identity';
            baseArgs.bendFactor = 1;
            baseArgs.GlobalBadDiffeoThreshold = 0.85;
            baseArgs.image = 'dash';
            baseArgs.metric = 'constant';
            baseArgs.periodic = true;
        case 'identity'
            baseArgs.fio = 'identity';
            baseArgs.n = 128;
            baseArgs.diffeos = 'identity';
        case 'riesz'
            baseArgs.fio = 'riesz';
            baseArgs.n = 128;
            baseArgs.diffeos = 'identity';
            baseArgs.alpha = 1;
            baseArgs.image = 'square-and-circle';
        case 'riesz2'
            baseArgs.fio = 'riesz';
            baseArgs.n = 128;
            baseArgs.alpha = 2;
            baseArgs.diffeos = 'identity';
            baseArgs.image = 'circle';
        case 'sarreflect1'
            baseArgs.fio = 'sarreflect';
            baseArgs.n = 64;
            baseArgs.diffeos = 'identity';
            baseArgs.curve = 'parabola';
            baseArgs.image = 'offset-circle';
            baseArgs.maskInside = 1;
        case 'sarreflect2'
            baseArgs.fio = 'sarreflect';
            baseArgs.n = 64;
            baseArgs.diffeos = 'two-plus-identity';
            baseArgs.curve = 'parabola';
            baseArgs.image = 'dash';
            baseArgs.bendFactor = 2;
            baseArgs.maskInside = -1;
        case 'sarreflect3'
            baseArgs.fio = 'sarreflect';
            baseArgs.n = 64;
            baseArgs.diffeos = 'two-plus-identity';
            baseArgs.curve = 'lower-half-circle';
            baseArgs.image = 'low-dash';
            baseArgs.bendFactor = 1;
            baseArgs.maskInside = -1;
        case 'sarreflect4'
            baseArgs.fio = 'sarreflect';
            baseArgs.n = 64;
            baseArgs.diffeos = 'two-plus-identity';
            baseArgs.curve = 'shallow-parabola';
            baseArgs.image = 'slight-low-dash';
            baseArgs.bendFactor = 1;
            baseArgs.maskInside = 0;
        otherwise
            error('PVBaseArgs:unknownFIOSetting', ...
                'Unknown FIO setting: %s', fio);
    end
    
    vdArgs.doErrors = true;
    vdArgs.showArgs = true;
end