
% Designed to work exactly as save() does, only you input the actual
% variable to be saved AND the name you want it to be called.

function parsave(savename,IN,varname)

eval([varname ' = IN;' ]);

save(savename,varname)

% eval(['save(' savename ',' inputname(2) ');']);

% varname = inputname(2);
% eval([varname ' = IN;'])
% 
% eval(['save(savename,''' varname ''',varargin{:})']);

end



% %% input name demo
% 
% 
% function myfun(a,b)
% s = inputname(1);
% disp(['First calling variable is ''' s '''.'])
% end




