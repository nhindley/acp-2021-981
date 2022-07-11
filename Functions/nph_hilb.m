
% nph_hilb ================================================================

% A newer, and much more simple, N-D approach to obtaining the analytic
% signal.
% find all the
% complex-conjugate pairs, set one to be zero and double the other. All
% coefficients not in a complex-conjugate pair are left unchanged.
% This approach uses linear indeces, so is N-D! Easy!
%


function varargout = nph_hilb(IN)

F = fftn(IN);

% mask template of NaNs
m = nan(size(F));

% first, put in the zero freqs:
m(imag(F) == 0) = 1;

% now, find pairs:
[a,ib] = sort(abs(imag(F(:))));

% get rid of the zeros before reshaping:
ib = ib(a ~= 0);
a = a(a ~= 0);

% now reshape:
% this should always work - after the zeros are taken out there should
% always be an even number of complex conjugate pairs remaining.
% note: you need the transpose ' here due to the way reshape re-lists things.
ar = reshape(a',2,length(a)/2);
ibr = reshape(ib',2,length(ib)/2);

% now assign 2s and 0s (doesn't matter which order):
m(ibr(1,:)) = 0;
m(ibr(2,:)) = 2;

% apply the hilbert mask:
switch nargout
    case 1
        varargout{1} = ifftn(F .* m);
    case 2
        varargout{1} = ifftn(F .* m);
        varargout{2} = m;
end

% and you're done!!




end
