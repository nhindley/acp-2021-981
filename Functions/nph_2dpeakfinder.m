
% Very basic 2D peak finder. Peaks are identifed as any point that is TOL
% times larger than the average of all its neighbours.

% EDIT: new strategy - pad the outside edges with NaNs so that you can just do
% all the middle ones.

% fhandle is handle to a function, e.g. @nanmean or @nanmedian. nanmean is
% default. Currently only supporting mean/median.

function varargout = nph_2dpeakfinder(IN,tol,fhandle)

nanpad = nan(size(IN)+2);
nanpad(2:end-1,2:end-1) = IN;
IN = nanpad;
OUT = IN;

sz = size(IN);
[I,J] = meshgrid(1:sz(1),1:sz(2));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Smush middle section around to get mean around each data point

IN_mid = IN(2:end-1,2:end-1);
sz_mid = size(IN_mid);

area_mean = nan([sz-2 8]);

mid = I > 1 & I < sz(1) & J > 1 & J < sz(2);

nudges = [...
    1 1 ; ...
    -1 -1 ; ...
    1 -1 ; ...
    -1 1 ; ...
    0 1 ; ...
    1 0 ; ...
    0 -1 ; ...
    -1 0];

for i = 1:8
    area_mean(:,:,i) = reshape(IN(circshift(circshift(mid,nudges(i,1),1),nudges(i,2),2)),sz_mid);
end

% get mean around each point in central section:
switch nargin
    case 3
        area_mean = fhandle(area_mean,3);
    otherwise
        area_mean = nanmean(area_mean,3);
end

% set threshold:
spikes_mid = IN_mid > (tol * area_mean);
IN_mid(spikes_mid) = area_mean(spikes_mid);

OUT = IN_mid;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OUTPUT

switch nargout
    case {0,1}
        varargout{1} = OUT;
    case 2
        varargout{1} = OUT;
        varargout{2} = spikes_mid;
end

