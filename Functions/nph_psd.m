

% power spectral density

function [psd,f] = nph_psd(IN,point_spacing,varargin)

switch nargin
    case 1
        point_spacing = 1;
end

F = fft(IN);
m = nph_hilbertmask(size(IN));

psd = abs(reverse(F(m == 2)) + conj(F(m == 0)));

periods = (length(IN) ./ (1:length(psd))) * point_spacing;

f = 1./periods;






















