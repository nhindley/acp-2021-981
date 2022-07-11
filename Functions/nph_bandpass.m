

% nph_bandpass.m
%
% As anyone that knows me will have realised by now, I don't work in
% frequencies.
%
% INPUT bandpass limits in ASCENDING NORMALISED WAVELENGTHS in DATA POINTS.
%
% bpinfo = [stop1 pass1 pass2 stop2]; (ascending for wavelength).
%
% Apply the filter using FILTERED_DATA = filter(OUT,DATA);
% or ZERO_PHASE_FILTERED_DATA = filtfilt(OUT,DATA);
%

function OUT = nph_bandpass(bpinfo,varargin)
% Returns a discrete-time filter object.

% Elliptic Bandpass filter designed using FDESIGN.BANDPASS.

% All frequency values are normalized to 1.

finfo = 2./bpinfo;
% (don't forget the factor of 2, it's normalised relative to the Nyquist freq)

Astop1 = 60;              % First Stopband Attenuation (dB)
Apass  = 0.01;            % Passband Ripple (dB)
Astop2 = 60;              % Second Stopband Attenuation (dB)
match  = 'both';          % Band to match exactly

% Construct an FDESIGN object and call its ELLIP method.
h  = fdesign.bandpass(finfo(4),finfo(3),finfo(2),finfo(1),Astop1, Apass, ...
    Astop2);

OUT = design(h, 'ellip', 'MatchExactly', match);

% fvtool(OUT)

if ~isempty(varargin)
    if any(strcmpi(varargin{:}),'plot')
        fvtool(OUT)
    end
end

end












