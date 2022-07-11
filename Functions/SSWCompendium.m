

%%%% List of SSW dates from ERA-Interim from Amy Butler's SSW Compondium at
%%%% https://csl.noaa.gov/groups/csl8/sswcompendium/majorevents.html


function SSWC = SSWCompendium

% Using ERA dates:

datestrings = {
    '22/02/1979',
    '29/02/1980',
    '04/03/1981',
    '04/12/1981',
    '24/02/1984',
    '01/01/1985',
    '23/01/1987',
    '08/12/1987',
    '14/03/1988',
    '21/02/1989',
    '15/12/1998',
    '26/02/1999',
    '20/03/2000',
    '11/02/2001',
    '30/12/2001',
    '17/02/2002',
    '18/01/2003',
    '05/01/2004',
    '21/01/2006',
    '24/02/2007',
    '22/02/2008',
    '24/01/2009',
    '09/02/2010',
    '24/03/2010',
    '06/01/2013',
    '12/02/2018',
    '02/01/2019'};


%%%% Split/Displacement as listed in Table 1 of Hall et al (2020) JGR.
%%%% Can't remember the method used but this can vary, so only take these
%%%% values as a guide
SSWC_Split           = [1 0 0 0 0 1 0 1 1 1 0 1 0 0 0 NaN 1 0 0 0 0 1 1 0 1 1 1]';
SSWC_Displacement    = [0 1 1 1 1 0 1 0 0 0 1 0 1 1 1 NaN 0 1 1 1 1 0 0 1 0 0 0]';
% the 2002 one is a NaN because it didn't occur in ERA, so Rich didn't
% classify it as S/D.


dates = datenum(datestrings,'dd/mm/yyyy');

SSWC = struct;
SSWC.Date           = dates;
SSWC.yyyymmdd       = cellstr(datestr(dates,'yyyymmdd'));
SSWC.Split          = SSWC_Split;
SSWC.Displacement   = SSWC_Displacement;

switch nargout
    case 0
        assignin('base','SSWC',SSWC);
end

% switch nargout
%     case 1
%         varargout{1} = dates;
%     case 2
%         varargout{1} = dates;
%         varargout{2} = SSWC_Split;
%     case 3
%         varargout{1} = dates;
%         varargout{2} = SSWC_Split;
%         varargout{3} = SSWC_Displacement;
%         
%     otherwise
%     
%     assignin('base','SSWC',dates);
%     assignin('base','SSWC_Split',SSWC_Split);
%     assignin('base','SSWC_Displacement',SSWC_Displacement);
%     
% end


end
