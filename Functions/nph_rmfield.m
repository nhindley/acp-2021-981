
% an rmfield function that doesn't bloody error when you field it a field
% that doesn't exist. infuriating.

function OUT = nph_rmfield(IN,varargin)

for i = 1:length(varargin{:})
    
    fields = varargin{:};
   
    if isfield(IN,fields{i})
        
        IN = rmfield(IN,fields{i});
    
    end
end

OUT = IN;
    

