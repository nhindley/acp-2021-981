

%%%% N-D

%%%% My version of our beloved bin2mat function, only this can be up to an
%%%% order of magnitude faster. Uses histcounts to do the binning, which is
%%%% rapid, but unfortunately I can't think of a way to include a function
%%%% handle for the binning (mean, sum, etc).

%%%% EDIT: Note that this function has been adjusted to use bin CENTRES
%%%% rather than bin EDGES. Histcounts is great, but it only specifies bin
%%%% edges as inputs. This is annoying because if the bins edges are 1 2 3
%%%% 4, then 1.9 would be binned into element 1, even though it's much
%%%% closer to 2. So I'm gonna make the inputs bin centres.

%%%% EDIT: NOW SUPPORTING 1D, 2D and 3D binning! Nice!


function OUT = nph_bin2matn(varargin)


% INPUTS
% x = x coords of input data
% y = y coords of input data
% z = input data to be binned (must be same size as x and y)
% xi = x bin CENTRES of the output grid
% yi = y bin CENTRES of the output grid

% Split varargin into inputs and options:
inputs = {}; options = {};
for i = 1:length(varargin)
    switch isnumeric(varargin{i})
        case 1
            inputs{end+1} = varargin{i};
        otherwise
            options{end+1} = varargin{i};
    end
end

% now go through the inputs and determine whether we're using 1D, 2D, or 3D
% cases:
switch length(inputs)
    case 3 % 1D: x, data, xi
        type = 1;
        x = inputs{1}; data = inputs{2}; xi = inputs{3};
        X = {x}; XI = {xi};
    case 5 % 2D: x, y, data, xi, yi
        type = 2;
        x = inputs{1}; y = inputs{2}; data = inputs{3}; xi = inputs{4}; yi = inputs{5};
        X = {x,y}; XI = {xi,yi};
    case 7 % 3D: x, y, z, data, xi, yi, zi
        type = 3;
        x = inputs{1}; y = inputs{2}; z = inputs{3}; data = inputs{4}; xi = inputs{5}; yi = inputs{6}; zi = inputs{7};
        X = {x,y,z}; XI = {xi,yi,zi};
    otherwise
        error('Expected either 3 (1D), 5 (2D) or 7 (3D) input data arguments.')
end
        
% detemine method:
if ~isempty(options)
    switch options{1}
        case {'mean','nanmean'}
            method = 'nanmean';
        case {'sum','nansum'}
            method = 'nansum';
        otherwise
            method = 'nanmean';
    end
else
    method = 'nanmean';
end

% check for monotonics
for i = 1:type
    if ~ismonotonic(XI{i})
        error('Sorry, output grids must be monotonically increasing.')
    end
end
% if ~ismonotonic(xi) || ~ismonotonic(yi)
%     error('Sorry, output grids must be monotonically increasing.')
% end

% linearise everything:
x = x(:);
y = y(:);
data = data(:);
xi = xi(:);
yi = yi(:);

%%%% Convert bin CENTRES to bin EDGES for use with histcounts:
dxi = diff(xi); dyi = diff(yi);
xi_edges = [xi - dxi([1 1:end])./2 ; xi(end) + dxi(end)/2];
yi_edges = [yi - dyi([1 1:end])./2 ; yi(end) + dyi(end)/2];

% I *think*, using histcounts, there will never be anything allocated to
% the final bin edge index unless you specify ...,'includedege','right')?
% This means we can trim it off in my bin centres method?

% linearise everything:
data = data(:);
[~,~,~,xbin,ybin] = histcounts2(x(:),y(:),xi_edges(:),yi_edges(:));

%%%% Important point: although the bin locations here correspond to the
%%%% xi_edges/yi_edges data, they should actually correspond perfectly to
%%%% the input bin CENTRES data, because nothing with ever be put in the
%%%% final "edge" bin of xi_edges/yi_edges. I think.
%%%% This means our output should always be the same size as the input, and
%%%% we shouldn't get any cases of allocating outside the data limits. In
%%%% theory, anyway.

% remove nans and zeros (zeros are where a data point doesn't fall in any bin):
goodinds = xbin ~= 0 & ybin ~= 0 & ~isnan(data);
data = data(goodinds);
xbin = xbin(goodinds);
ybin = ybin(goodinds);

% preallocate output
osz = [length(xi) length(yi)];
OUT         = zeros(osz);
emptycount  = zeros(osz);

% % try a loop, can't think of a better way to do it.
switch method
    
    case 'nanmean'
        count = ones(osz);
        for i = 1:length(data)
            OUT(xbin(i),ybin(i)) = OUT(xbin(i),ybin(i)) + data(i);
            count(xbin(i),ybin(i)) = count(xbin(i),ybin(i)) + 1;
            emptycount(xbin(i),ybin(i)) = 1;
        end
        OUT = OUT ./ count;
        OUT(emptycount == 0) = NaN;
        
    case 'nansum'
        for i = 1:length(data)
            OUT(xbin(i),ybin(i)) = OUT(xbin(i),ybin(i)) + data(i);
            emptycount(xbin(i),ybin(i)) = 1;
        end
        OUT(emptycount == 0) = NaN;
        
end
% turns out loop is pretty fast because no functions are called in it, I
% think...










