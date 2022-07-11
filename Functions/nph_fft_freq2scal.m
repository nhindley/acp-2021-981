
% nph_fft_freq2scal.m

% For use with the 3DST! Convert real frequencies you want to analyse into
% the best fitting integer fractions of the total length of the timeseries.
% 
% scales = the integer fractions of the total length of the thing your
% trying to fft, ie 1:10 etc
%  
% freqs = the real frequencies you want to analysise
%
% L = real total length of the input timeseries, ie not just number of
% elements but the real physical distance they represent
% 


% note that due to the unique(), scales might not always be the same length
% as freqs. This is to remove duplication of frequencies in the 3DST. For
% low frequencies where integer fractions are sparse, this might be quite common.

function scales = nph_fft_freq2scal(freqs,L)

% Get fractions:
scales = round(L./freqs);

% Outliers:
% scales(scales == 0) = [];

% unique:
scales = unique(scales);







