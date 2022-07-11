
% So there's linspace and logspace, but what if you want to have something
% that goes up linearly in a power? Like square space or cube space?

% power space makes a vector that runs from LO to HI in NUM steps that are
% spaced by the power specified in POW.

% note this is slightly different from the behaviour of logspace, which I
% always find a bit confusing anyway, in that it goes from the REAL number
% LO to the REAL number HI, rather than from log10(LO) to log10(HI) like
% logspace does.

% when would you use this? Instead of plotting contour maps of logged
% properties, you can now simply power the colormap instead, provided it's
% not diverging, at any power up to 10. If base 10 logging was too much,
% try a power of 2 or 3 etc.

% e.g. powerspace(2,32,4,2) = [2 8 18 32];

function OUT = powerspace(lo,hi,num,pow)

OUT = linspace(lo.^(1/pow),hi.^(1/pow),num).^pow;

end






