

% sliding lomb-scargle periodogram of a 1D time series

%%%% INPUTS:
% IN = input timeseries to be analysed
% x = time coordinates of the input values
% samplepoints = points in the x coords on which to compute the LS
% windowsize = size of sliding window in the UNITS OF THE X COORDS
% % % % % stepsize = size of sliding steps in the UNITS OF THE X COORDS
% freqs = the frequencies to be analysed (units are 1/time coords)
% varargin = any extra inputs are passed to plomb.

%%%% OUTPUTS:
% OUT = 2D matrix of power spectral density, time against frequency.
% t = vector of time coordinates the same size
% f = vector of the frequencies analysed




function varargout = nph_sliding_lombscargle(IN,x,samplepoints,windowsize,freqs,varargin)

% arrange outputs
OUT = nan(length(freqs),length(samplepoints));
f = freqs;
t = samplepoints;

% nsteps = floor(minmax(x) ./ stepsize);

% begin loopy loop
for i = 1:length(samplepoints)
    
    disp(i)
    
    win = floor(windowsize/2);
    xrange = inrange(x,samplepoints(i)-win,samplepoints(i)+win);
    
    yy = IN(xrange);
    xx = x(xrange);
    
    % give zero mean?
    yy = yy - nanmean(yy);
    
    % sort out nans
    nanlocs = isnan(yy) | isnan(xx);
    
    if ~isempty(yy) && ~isempty(xx)
        [pxx,f] = plomb(yy(~nanlocs),xx(~nanlocs),freqs,varargin{:});
        OUT(:,i) = pxx;
    else
        OUT(:,i) = NaN;
    end
    
    
end



% OUTPUTS
switch nargout
    case {0,1}
        varargout{1} = OUT;
    case 2
        varargout{1} = OUT;
        varargout{2} = f;
    case 3
        varargout{1} = OUT;
        varargout{2} = f;
        varargout{3} = t;
end











