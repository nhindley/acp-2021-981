
function s = frq2scal(freqs,wname,delta)


% FRQ2SCAL



if isempty(freqs) , s = freqs; return; end
if nargin == 2, delta = 1; end
err = (min(size(freqs))>1) | (min(freqs)<eps);
if err
    error(message('Wavelet:FunctionArgVal:Invalid_ScaVal'))
end
if delta <= 0
    error(message('Wavelet:FunctionArgVal:Invalid_DeltaVal'))
end

% Compute pseudo-frequencies
%f = centfrq(wname)./(a.*delta);     

s = centfrq(wname)./(freqs.*delta);

